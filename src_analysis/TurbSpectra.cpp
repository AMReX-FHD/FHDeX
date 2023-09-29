#include "common_functions.H"
#include "TurbSpectra.H"

#include <AMReX_MultiFabUtil.H>
#include "AMReX_PlotFileUtil.H"
#include "AMReX_BoxArray.H"

#ifdef AMREX_USE_CUDA
std::string cufftError (const cufftResult& err)
{
    switch (err) {
    case CUFFT_SUCCESS:  return "CUFFT_SUCCESS";
    case CUFFT_INVALID_PLAN: return "CUFFT_INVALID_PLAN";
    case CUFFT_ALLOC_FAILED: return "CUFFT_ALLOC_FAILED";
    case CUFFT_INVALID_TYPE: return "CUFFT_INVALID_TYPE";
    case CUFFT_INVALID_VALUE: return "CUFFT_INVALID_VALUE";
    case CUFFT_INTERNAL_ERROR: return "CUFFT_INTERNAL_ERROR";
    case CUFFT_EXEC_FAILED: return "CUFFT_EXEC_FAILED";
    case CUFFT_SETUP_FAILED: return "CUFFT_SETUP_FAILED";
    case CUFFT_INVALID_SIZE: return "CUFFT_INVALID_SIZE";
    case CUFFT_UNALIGNED_DATA: return "CUFFT_UNALIGNED_DATA";
    default: return std::to_string(err) + " (unknown error code)";
    }
}
#endif

#ifdef AMREX_USE_HIP
std::string rocfftError (const rocfft_status err)
{
    if              (err == rocfft_status_success) {
        return std::string("rocfft_status_success");
    } else if       (err == rocfft_status_failure) {
        return std::string("rocfft_status_failure");
    } else if       (err == rocfft_status_invalid_arg_value) {
        return std::string("rocfft_status_invalid_arg_value");
    } else if       (err == rocfft_status_invalid_dimensions) {
        return std::string("rocfft_status_invalid_dimensions");
    } else if       (err == rocfft_status_invalid_array_type) {
        return std::string("rocfft_status_invalid_array_type");
    } else if       (err == rocfft_status_invalid_strides) {
        return std::string("rocfft_status_invalid_strides");
    } else if       (err == rocfft_status_invalid_distance) {
        return std::string("rocfft_status_invalid_distance");
    } else if       (err == rocfft_status_invalid_offset) {
        return std::string("rocfft_status_invalid_offset");
    } else {
        return std::to_string(err) + " (unknown error code)";
    }
}

void rocfft_status (std::string const& name, rocfft_status status)
{
    if (status != rocfft_status_success) {
        amrex::AllPrint() <<  name + " failed! Error: " + rocfftError(status) << "\n";;
    }
}
#endif

void TurbSpectrumScalar(const MultiFab& variables, 
                        const amrex::Geometry& geom, 
                        const int& step, 
                        const amrex::Vector<amrex::Real>& scaling,
                        const amrex::Vector< std::string >& var_names)
{
    BL_PROFILE_VAR("TurbSpectrumScalar()",TurbSpectrumScalar);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(variables.nComp() == var_names.size(), "TurbSpectrumScalar: must have same number variable names as components of input MultiFab");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(variables.nComp() == scaling.size(), "TurbSpectrumScalar: must have same number variable scaling as components of input MultiFab");
#if defined(HEFFTE)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(variables.local_size() == 1, "TurbSpectrumScalar: Must have one Box per MPI process when using heFFTe");

    int ncomp = variables.nComp();

    long npts;
    Box domain = geom.Domain();
    npts = (domain.length(0)*domain.length(1)*domain.length(2));
    Real sqrtnpts = std::sqrt(npts);
    
    // get box array and distribution map of variables
    DistributionMapping dm = variables.DistributionMap();
    BoxArray ba            = variables.boxArray();
    
    // since there is 1 MPI rank per box, each MPI rank obtains its local box and the associated boxid
    Box local_box;
    int local_boxid;
    {
        for (int i = 0; i < ba.size(); ++i) {
            Box b = ba[i];
            // each MPI rank has its own local_box Box and local_boxid ID
            if (ParallelDescriptor::MyProc() == dm[i]) {
                local_box = b;
                local_boxid = i;
            }
        }
    }

    // now each MPI rank works on its own box
    // for real->complex fft's, the fft is stored in an (nx/2+1) x ny x nz dataset

    // start by coarsening each box by 2 in the x-direction
    Box c_local_box = amrex::coarsen(local_box, IntVect(AMREX_D_DECL(2,1,1)));

    // if the coarsened box's high-x index is even, we shrink the size in 1 in x
    // this avoids overlap between coarsened boxes
    if (c_local_box.bigEnd(0) * 2 == local_box.bigEnd(0)) {
        c_local_box.setBig(0,c_local_box.bigEnd(0)-1);
    }
    // for any boxes that touch the hi-x domain we
    // increase the size of boxes by 1 in x
    // this makes the overall fft dataset have size (Nx/2+1 x Ny x Nz)
    if (local_box.bigEnd(0) == geom.Domain().bigEnd(0)) {
        c_local_box.growHi(0,1);
    }

    // each MPI rank gets storage for its piece of the fft
    BaseFab<GpuComplex<Real> > spectral_field(c_local_box, ncomp, The_Device_Arena());
    
    int r2c_direction = 0;
#ifdef AMREX_USE_CUDA
    heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif AMREX_USE_HIP
    heffte::fft3d_r2c<heffte::backend::rocfft> fft
#else
    heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
    ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
    {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
    {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
    {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
    r2c_direction, ParallelDescriptor::Communicator());

    using heffte_complex = typename heffte::fft_output<Real>::type;
    for (int comp=0; comp<ncomp; ++comp) {
        heffte_complex* spectral_data = (heffte_complex*) spectral_field.dataPtr(comp);
        fft.forward(variables[local_boxid].dataPtr(comp),spectral_data); 
        
        ParallelDescriptor::Barrier();

        // Integrate spectra over k-shells
        IntegrateKScalar(spectral_field,var_names[comp],scaling[comp],c_local_box,sqrtnpts,step);
    }
#else // not heFFTe
    int ncomp = variables.nComp();

    long npts;

    // Initialize the boxarray "ba_onegrid" from the single box "domain"
    BoxArray ba_onegrid;
    {
      Box domain = geom.Domain();
      ba_onegrid.define(domain);
      npts = (domain.length(0)*domain.length(1)*domain.length(2));
    }
    Real sqrtnpts = std::sqrt(npts);
    DistributionMapping dmap_onegrid(ba_onegrid);
    MultiFab variables_onegrid;
    variables_onegrid.define(ba_onegrid, dmap_onegrid, ncomp, 0);
    variables_onegrid.ParallelCopy(variables,0,0,ncomp);

#ifdef AMREX_USE_CUDA
    using FFTplan = cufftHandle;
    using FFTcomplex = cuDoubleComplex;
#elif AMREX_USE_HIP
    using FFTplan = rocfft_plan;
    using FFTcomplex = double2;
#else
    using FFTplan = fftw_plan;
    using FFTcomplex = fftw_complex;
#endif


    // size of box including ghost cell range
    IntVect fft_size;
    
    // contain to store FFT - note it is shrunk by "half" in x
    Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > > spectral_field;
    for (MFIter mfi(variables_onegrid); mfi.isValid(); ++mfi) {

        // grab a single box including ghost cell range
        Box realspace_bx = mfi.fabbox();
        
        // size of box including ghost cell range
        fft_size = realspace_bx.length(); // This will be different for hybrid FFT

        // this is the size of the box, except the 0th component is 'halved plus 1'
        IntVect spectral_bx_size = fft_size;
        spectral_bx_size[0] = fft_size[0]/2 + 1;

        // spectral box
        Box spectral_bx = Box(IntVect(0), spectral_bx_size - IntVect(1));

        spectral_field.emplace_back(new BaseFab<GpuComplex<Real> >(spectral_bx,ncomp,
                                                               The_Device_Arena()));
        spectral_field.back()->setVal<RunOn::Device>(0.0); // touch the memory
    }

    // for CUDA builds we only need to build the plan once; track whether we did
    Vector<FFTplan> forward_plan;
    bool built_plan = false;
    for (int comp=0; comp<ncomp; ++comp) {
        if (!built_plan) {
            for (MFIter mfi(variables_onegrid); mfi.isValid(); ++mfi) {
                FFTplan fplan;

#ifdef AMREX_USE_CUDA // CUDA
                cufftResult result = cufftPlan3d(&fplan, fft_size[2], fft_size[1], fft_size[0], CUFFT_D2Z);
                if (result != CUFFT_SUCCESS) {
                    amrex::AllPrint() << " cufftplan3d forward failed! Error: "
                                      << cufftError(result) << "\n";
                }
                built_plan = true;
#elif AMREX_USE_HIP // HIP
                const std::size_t lengths[] = {std::size_t(fft_size[0]),std::size_t(fft_size[1]),std::size_t(fft_size[2])};
                rocfft_status result = rocfft_plan_create(&fplan, rocfft_placement_notinplace, 
                                                          rocfft_transform_type_real_forward, rocfft_precision_double,
                                                          3, lengths, 1, nullptr);
                assert_rocfft_status("rocfft_plan_create", result);
                built_plan = true;
#else // host
                fplan = fftw_plan_dft_r2c_3d(fft_size[2], fft_size[1], fft_size[0],
                                              variables_onegrid[mfi].dataPtr(comp),
                                              reinterpret_cast<FFTcomplex*>
                                              (spectral_field.back()->dataPtr(comp)),
                                              FFTW_ESTIMATE);
#endif
                forward_plan.push_back(fplan);
            }
        }
        
        ParallelDescriptor::Barrier();

        // ForwardTransform
        for (MFIter mfi(variables_onegrid); mfi.isValid(); ++mfi) {
            int i = mfi.LocalIndex();
#ifdef AMREX_USE_CUDA
            cufftSetStream(forward_plan[i], amrex::Gpu::gpuStream());
            cufftResult result = cufftExecD2Z(forward_plan[i],
                                              variables_onegrid[mfi].dataPtr(comp),
                                              reinterpret_cast<FFTcomplex*>
                                                  (spectral_field[i]->dataPtr(comp)));
            if (result != CUFFT_SUCCESS) {
                amrex::AllPrint() << " forward transform using cufftExec failed! Error: "
                                  << cufftError(result) << "\n";
	    }
#elif AMREX_USE_HIP
            rocfft_execution_info execinfo = nullptr;
            rocfft_status result = rocfft_execution_info_create(&execinfo);
            assert_rocfft_status("rocfft_execution_info_create", result);

            std::size_t buffersize = 0;
            result = rocfft_plan_get_work_buffer_size(forward_plan[i], &buffersize);
            assert_rocfft_status("rocfft_plan_get_work_buffer_size", result);

            void* buffer = amrex::The_Arena()->alloc(buffersize);
            result = rocfft_execution_info_set_work_buffer(execinfo, buffer, buffersize);
            assert_rocfft_status("rocfft_execution_info_set_work_buffer", result);

            result = rocfft_execution_info_set_stream(execinfo, amrex::Gpu::gpuStream());
            assert_rocfft_status("rocfft_execution_info_set_stream", result);

	        amrex::Real* variables_onegrid_ptr = variables_onegrid[mfi].dataPtr(comp);
	        FFTcomplex* spectral_field_ptr = reinterpret_cast<FFTcomplex*>(spectral_field[i]->dataPtr(comp));
            result = rocfft_execute(forward_plan[i],
                                    (void**) &variables_onegrid_ptr, // in
                                    (void**) &spectral_field_ptr, // out
                                    execinfo);
            assert_rocfft_status("rocfft_execute", result);
            amrex::Gpu::streamSynchronize();
            amrex::The_Arena()->free(buffer);
            result = rocfft_execution_info_destroy(execinfo);
            assert_rocfft_status("rocfft_execution_info_destroy", result);
#else
            fftw_execute(forward_plan[i]);
#endif
        }

        // Integrate spectra over k-shells
        IntegrateKScalar(spectral_field,variables_onegrid,var_names[comp],scaling[comp],sqrtnpts,step);
    }
    
    // destroy fft plan
    for (int i = 0; i < forward_plan.size(); ++i) {
#ifdef AMREX_USE_CUDA
        cufftDestroy(forward_plan[i]);
#elif AMREX_USE_HIP
        rocfft_plan_destroy(forward_plan[i]);
#else
        fftw_destroy_plan(forward_plan[i]);
#endif
    }
#endif // end heFFTE
}

void TurbSpectrumVelDecomp(const MultiFab& vel,
                           MultiFab& vel_decomp,
                           const amrex::Geometry& geom,
                           const int& step,
                           const amrex::Real& scaling,
                           const amrex::Vector< std::string >& var_names)
{
    BL_PROFILE_VAR("TurbSpectrumVelDecomp()",TurbSpectrumVelDecomp);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vel.nComp() == 3, "TurbSpectrumVelDecomp: must have 3 components of input vel MultiFab");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(var_names.size() == 3, "TurbSpectrumVelDecomp: must have 3 names for output vel spectra (total, solenoidal, dilatational");
#if defined(HEFFTE)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vel.local_size() == 1, "TurbSpectrumVelDecomp: Must have one Box per MPI process when using heFFTe");

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    long npts;
    Box domain = geom.Domain();
    npts = (domain.length(0)*domain.length(1)*domain.length(2));
    Real sqrtnpts = std::sqrt(npts);
    
    // get box array and distribution map of vel
    DistributionMapping dm = vel.DistributionMap();
    BoxArray ba            = vel.boxArray();
    
    // since there is 1 MPI rank per box, each MPI rank obtains its local box and the associated boxid
    Box local_box;
    int local_boxid;
    {
        for (int i = 0; i < ba.size(); ++i) {
            Box b = ba[i];
            // each MPI rank has its own local_box Box and local_boxid ID
            if (ParallelDescriptor::MyProc() == dm[i]) {
                local_box = b;
                local_boxid = i;
            }
        }
    }

    // now each MPI rank works on its own box
    // for real->complex fft's, the fft is stored in an (nx/2+1) x ny x nz dataset

    // start by coarsening each box by 2 in the x-direction
    Box c_local_box = amrex::coarsen(local_box, IntVect(AMREX_D_DECL(2,1,1)));

    // if the coarsened box's high-x index is even, we shrink the size in 1 in x
    // this avoids overlap between coarsened boxes
    if (c_local_box.bigEnd(0) * 2 == local_box.bigEnd(0)) {
        c_local_box.setBig(0,c_local_box.bigEnd(0)-1);
    }
    // for any boxes that touch the hi-x domain we
    // increase the size of boxes by 1 in x
    // this makes the overall fft dataset have size (Nx/2+1 x Ny x Nz)
    if (local_box.bigEnd(0) == geom.Domain().bigEnd(0)) {
        c_local_box.growHi(0,1);
    }

    // each MPI rank gets storage for its piece of the fft
    BaseFab<GpuComplex<Real> > spectral_field_T(c_local_box, 3, The_Device_Arena()); // total
    BaseFab<GpuComplex<Real> > spectral_field_S(c_local_box, 3, The_Device_Arena()); // solenoidal
    BaseFab<GpuComplex<Real> > spectral_field_D(c_local_box, 3, The_Device_Arena()); // dilatational
    

    int r2c_direction = 0;
#ifdef AMREX_USE_CUDA
    heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif AMREX_USE_HIP
    heffte::fft3d_r2c<heffte::backend::rocfft> fft
#else
    heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
    ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
    {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
    {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
    {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
    r2c_direction, ParallelDescriptor::Communicator());
    
    // ForwardTransform
    using heffte_complex = typename heffte::fft_output<Real>::type;
    {
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_T.dataPtr(0);
        fft.forward(vel[local_boxid].dataPtr(0),spectral_data);
    }
    {
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_T.dataPtr(1);
        fft.forward(vel[local_boxid].dataPtr(1),spectral_data);
    }
    {
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_T.dataPtr(2);
        fft.forward(vel[local_boxid].dataPtr(2),spectral_data);
    }
    
    // Decompose velocity field into solenoidal and dilatational
    Array4< GpuComplex<Real> > spectral_tx = spectral_field_T.array(0,1);
    Array4< GpuComplex<Real> > spectral_ty = spectral_field_T.array(1,1);
    Array4< GpuComplex<Real> > spectral_tz = spectral_field_T.array(2,1);
    Array4< GpuComplex<Real> > spectral_sx = spectral_field_S.array(0,1);
    Array4< GpuComplex<Real> > spectral_sy = spectral_field_S.array(1,1);
    Array4< GpuComplex<Real> > spectral_sz = spectral_field_S.array(2,1);
    Array4< GpuComplex<Real> > spectral_dx = spectral_field_D.array(0,1);
    Array4< GpuComplex<Real> > spectral_dy = spectral_field_D.array(1,1);
    Array4< GpuComplex<Real> > spectral_dz = spectral_field_D.array(2,1);
    ParallelFor(c_local_box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {

       int nx = n_cells[0]; 
       int ny = n_cells[1]; 
       int nz = n_cells[2];

       Real GxR, GxC, GyR, GyC, GzR, GzC;
       
       if (i <= nx/2) { 
           // Gradient Operators
           GxR = (cos(2.0*M_PI*i/nx)-1.0)/dx[0];
           GxC = (sin(2.0*M_PI*i/nx)-0.0)/dx[0];
           GyR = (cos(2.0*M_PI*j/ny)-1.0)/dx[1];
           GyC = (sin(2.0*M_PI*j/ny)-0.0)/dx[1];
           GzR = (cos(2.0*M_PI*k/nz)-1.0)/dx[2];
           GzC = (sin(2.0*M_PI*k/nz)-0.0)/dx[2];
       }
       else { // conjugate
            amrex::Abort("check the code; i should not go beyond bx.length(0)/2");
       }

       // Scale Total velocity FFT components
       spectral_tx(i,j,k) *= (1.0/sqrtnpts);
       spectral_ty(i,j,k) *= (1.0/sqrtnpts);
       spectral_tz(i,j,k) *= (1.0/sqrtnpts);

       // Inverse Laplacian
       Real Lap = GxR*GxR + GxC*GxC + GyR*GyR + GyC*GyC + GzR*GzR + GzC*GzC;

       // Divergence of vel
       Real divR = spectral_tx(i,j,k).real()*GxR - spectral_tx(i,j,k).imag()*GxC +
                   spectral_ty(i,j,k).real()*GyR - spectral_ty(i,j,k).imag()*GyC +
                   spectral_tz(i,j,k).real()*GzR - spectral_tz(i,j,k).imag()*GzC ;
       Real divC = spectral_tx(i,j,k).real()*GxC + spectral_tx(i,j,k).imag()*GxR +
                   spectral_ty(i,j,k).real()*GyC + spectral_ty(i,j,k).imag()*GyR +
                   spectral_tz(i,j,k).real()*GzC + spectral_tz(i,j,k).imag()*GzR ;

       if (Lap < 1.0e-12) { // zero mode for no bulk motion
           spectral_dx(i,j,k) *= 0.0;
           spectral_dy(i,j,k) *= 0.0;
           spectral_dz(i,j,k) *= 0.0;
       }
       else {

           // Dilatational velocity 
           GpuComplex<Real> copy_dx((divR*GxR + divC*GxC) / Lap, 
                                    (divC*GxR - divR*GxC) / Lap);
           spectral_dx(i,j,k) = copy_dx;

           GpuComplex<Real> copy_dy((divR*GyR + divC*GyC) / Lap,
                                    (divC*GyR - divR*GyC) / Lap);
           spectral_dy(i,j,k) = copy_dy;

           GpuComplex<Real> copy_dz((divR*GzR + divC*GzC) / Lap,
                                    (divC*GzR - divR*GzC) / Lap);
           spectral_dz(i,j,k) = copy_dz;
           
           // Solenoidal velocity
           spectral_sx(i,j,k) = spectral_tx(i,j,k) - spectral_dx(i,j,k);
           spectral_sy(i,j,k) = spectral_ty(i,j,k) - spectral_dy(i,j,k); 
           spectral_sz(i,j,k) = spectral_tz(i,j,k) - spectral_dz(i,j,k);
       }

    });

    ParallelDescriptor::Barrier();

    // Integrate K spectrum for velocities
    IntegrateKVelocity(spectral_field_T,"turb_total",     scaling,c_local_box,step);
    IntegrateKVelocity(spectral_field_S,"turb_solenoidal",scaling,c_local_box,step);
    IntegrateKVelocity(spectral_field_D,"turb_dilational",scaling,c_local_box,step);
    
    // inverse Fourier transform solenoidal and dilatational components 
    {
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_S.dataPtr(0);
        fft.backward(spectral_data, vel_decomp[local_boxid].dataPtr(0));

    }
    {
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_S.dataPtr(1);
        fft.backward(spectral_data, vel_decomp[local_boxid].dataPtr(1));

    }
    {
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_S.dataPtr(2);
        fft.backward(spectral_data, vel_decomp[local_boxid].dataPtr(2));

    }
    {
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_D.dataPtr(0);
        fft.backward(spectral_data, vel_decomp[local_boxid].dataPtr(3));

    }
    {
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_D.dataPtr(1);
        fft.backward(spectral_data, vel_decomp[local_boxid].dataPtr(4));

    }
    {
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_D.dataPtr(2);
        fft.backward(spectral_data, vel_decomp[local_boxid].dataPtr(5));

    }

    vel_decomp.mult(1.0/sqrtnpts);
#else // not heFFTe
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    long npts;

    // Initialize the boxarray "ba_onegrid" from the single box "domain"
    BoxArray ba_onegrid;
    {
      Box domain = geom.Domain();
      ba_onegrid.define(domain);
      npts = (domain.length(0)*domain.length(1)*domain.length(2));
    }
    Real sqrtnpts = std::sqrt(npts);
    DistributionMapping dmap_onegrid(ba_onegrid);
    MultiFab vel_onegrid;
    vel_onegrid.define(ba_onegrid, dmap_onegrid, 3, 0);
    vel_onegrid.ParallelCopy(vel,0,0,3);

#ifdef AMREX_USE_CUDA
    using FFTplan = cufftHandle;
    using FFTcomplex = cuDoubleComplex;
#elif AMREX_USE_HIP
    using FFTplan = rocfft_plan;
    using FFTcomplex = double2;
#else
    using FFTplan = fftw_plan;
    using FFTcomplex = fftw_complex;
#endif

    // size of box including ghost cell range
    IntVect fft_size;
    
    // contain to store FFT - note it is shrunk by "half" in x
    Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > > spectral_field;
    Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > > spectral_field_S;
    Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > > spectral_field_D;
    for (MFIter mfi(vel_onegrid); mfi.isValid(); ++mfi) {

        // grab a single box including ghost cell range
        Box realspace_bx = mfi.fabbox();

        // size of box including ghost cell range
        fft_size = realspace_bx.length(); // This will be different for hybrid FFT

        // this is the size of the box, except the 0th component is 'halved plus 1'
        IntVect spectral_bx_size = fft_size;
        spectral_bx_size[0] = fft_size[0]/2 + 1;

        // spectral box
        Box spectral_bx = Box(IntVect(0), spectral_bx_size - IntVect(1));

        spectral_field.emplace_back(new BaseFab<GpuComplex<Real> >(spectral_bx,3,
                                                               The_Device_Arena()));
        spectral_field.back()->setVal<RunOn::Device>(0.0); // touch the memory

        spectral_field_S.emplace_back(new BaseFab<GpuComplex<Real> >(spectral_bx,3,
                                                               The_Device_Arena()));
        spectral_field_S.back()->setVal<RunOn::Device>(0.0); // touch the memory

        spectral_field_D.emplace_back(new BaseFab<GpuComplex<Real> >(spectral_bx,3,
                                                               The_Device_Arena()));
        spectral_field_D.back()->setVal<RunOn::Device>(0.0); // touch the memory
    }

    // for CUDA builds we only need to build the plan once; track whether we did
    Vector<FFTplan> forward_plan;
    bool built_plan = false;
    for (int comp=0; comp<3; comp++) {
        if (!built_plan) {
            for (MFIter mfi(vel_onegrid); mfi.isValid(); ++mfi) {
                FFTplan fplan;

#ifdef AMREX_USE_CUDA // CUDA
                cufftResult result = cufftPlan3d(&fplan, fft_size[2], fft_size[1], fft_size[0], CUFFT_D2Z);
                if (result != CUFFT_SUCCESS) {
                    amrex::AllPrint() << " cufftplan3d forward failed! Error: "
                                      << cufftError(result) << "\n";
                }
                built_plan = true;
#elif AMREX_USE_HIP // HIP
                const std::size_t lengths[] = {std::size_t(fft_size[0]),std::size_t(fft_size[1]),std::size_t(fft_size[2])};
                rocfft_status result = rocfft_plan_create(&fplan, rocfft_placement_notinplace, 
                                                          rocfft_transform_type_real_forward, rocfft_precision_double,
                                                          3, lengths, 1, nullptr);
                assert_rocfft_status("rocfft_plan_create", result);
                built_plan = true;
#else // host
                fplan = fftw_plan_dft_r2c_3d(fft_size[2], fft_size[1], fft_size[0],
                                                  vel_onegrid[mfi].dataPtr(comp),
                                                  reinterpret_cast<FFTcomplex*>
                                                  (spectral_field.back()->dataPtr(comp)),
                                                  FFTW_ESTIMATE);
#endif
               forward_plan.push_back(fplan);
            }
        }

        ParallelDescriptor::Barrier();

        // ForwardTransform
        for (MFIter mfi(vel_onegrid); mfi.isValid(); ++mfi) {
            int i = mfi.LocalIndex();
#ifdef AMREX_USE_CUDA
            cufftSetStream(forward_plan[i], amrex::Gpu::gpuStream());
            cufftResult result = cufftExecD2Z(forward_plan[i],
                                              vel_onegrid[mfi].dataPtr(comp),
                                              reinterpret_cast<FFTcomplex*>
                                                  (spectral_field[i]->dataPtr(comp)));
            if (result != CUFFT_SUCCESS) {
                amrex::AllPrint() << " forward transform using cufftExec failed! Error: "
                                  << cufftError(result) << "\n";
	        }
#elif AMREX_USE_HIP
            rocfft_execution_info execinfo = nullptr;
            rocfft_status result = rocfft_execution_info_create(&execinfo);
            assert_rocfft_status("rocfft_execution_info_create", result);

            std::size_t buffersize = 0;
            result = rocfft_plan_get_work_buffer_size(forward_plan[i], &buffersize);
            assert_rocfft_status("rocfft_plan_get_work_buffer_size", result);

            void* buffer = amrex::The_Arena()->alloc(buffersize);
            result = rocfft_execution_info_set_work_buffer(execinfo, buffer, buffersize);
            assert_rocfft_status("rocfft_execution_info_set_work_buffer", result);

            result = rocfft_execution_info_set_stream(execinfo, amrex::Gpu::gpuStream());
            assert_rocfft_status("rocfft_execution_info_set_stream", result);

	        amrex::Real* vel_onegrid_ptr = vel_onegrid[mfi].dataPtr(comp);
	        FFTcomplex* spectral_field_ptr = reinterpret_cast<FFTcomplex*>(spectral_field[i]->dataPtr(comp));
            result = rocfft_execute(forward_plan[i],
                                    (void**) &vel_onegrid_ptr, // in
                                    (void**) &spectral_field_ptr, // out
                                    execinfo);
            assert_rocfft_status("rocfft_execute", result);
            amrex::Gpu::streamSynchronize();
            amrex::The_Arena()->free(buffer);
            result = rocfft_execution_info_destroy(execinfo);
            assert_rocfft_status("rocfft_execution_info_destroy", result);
#else
            fftw_execute(forward_plan[i]);
#endif
        }
    }
    
    // destroy fft plan
    for (int i = 0; i < forward_plan.size(); ++i) {
#ifdef AMREX_USE_CUDA
        cufftDestroy(forward_plan[i]);
#elif AMREX_USE_HIP
        rocfft_plan_destroy(forward_plan[i]);
#else
        fftw_destroy_plan(forward_plan[i]);
#endif
    }

    // Decompose velocity field into solenoidal and dilatational
    for ( MFIter mfi(vel_onegrid,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();
        Array4< GpuComplex<Real> > spectral_tx = (*spectral_field[0])  .array(0,1);
        Array4< GpuComplex<Real> > spectral_ty = (*spectral_field[0])  .array(1,1);
        Array4< GpuComplex<Real> > spectral_tz = (*spectral_field[0])  .array(2,1);
        Array4< GpuComplex<Real> > spectral_sx = (*spectral_field_S[0]).array(0,1);
        Array4< GpuComplex<Real> > spectral_sy = (*spectral_field_S[0]).array(1,1);
        Array4< GpuComplex<Real> > spectral_sz = (*spectral_field_S[0]).array(2,1);
        Array4< GpuComplex<Real> > spectral_dx = (*spectral_field_D[0]).array(0,1);
        Array4< GpuComplex<Real> > spectral_dy = (*spectral_field_D[0]).array(1,1);
        Array4< GpuComplex<Real> > spectral_dz = (*spectral_field_D[0]).array(2,1);
            
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           int nx = n_cells[0]; 
           int ny = n_cells[1]; 
           int nz = n_cells[2];

           Real GxR, GxC, GyR, GyC, GzR, GzC;
           
           if (i <= nx/2) { 
               // Gradient Operators
               GxR = (cos(2.0*M_PI*i/nx)-1.0)/dx[0];
               GxC = (sin(2.0*M_PI*i/nx)-0.0)/dx[0];
               GyR = (cos(2.0*M_PI*j/ny)-1.0)/dx[1];
               GyC = (sin(2.0*M_PI*j/ny)-0.0)/dx[1];
               GzR = (cos(2.0*M_PI*k/nz)-1.0)/dx[2];
               GzC = (sin(2.0*M_PI*k/nz)-0.0)/dx[2];
           }
           else { // conjugate
                amrex::Abort("check the code; i should not go beyond bx.length(0)/2");
           }

           // Scale Total velocity FFT components
           spectral_tx(i,j,k) *= (1.0/sqrtnpts);
           spectral_ty(i,j,k) *= (1.0/sqrtnpts);
           spectral_tz(i,j,k) *= (1.0/sqrtnpts);

           // Inverse Laplacian
           Real Lap = GxR*GxR + GxC*GxC + GyR*GyR + GyC*GyC + GzR*GzR + GzC*GzC;

           // Divergence of vel
           Real divR = spectral_tx(i,j,k).real()*GxR - spectral_tx(i,j,k).imag()*GxC +
                       spectral_ty(i,j,k).real()*GyR - spectral_ty(i,j,k).imag()*GyC +
                       spectral_tz(i,j,k).real()*GzR - spectral_tz(i,j,k).imag()*GzC ;
           Real divC = spectral_tx(i,j,k).real()*GxC + spectral_tx(i,j,k).imag()*GxR +
                       spectral_ty(i,j,k).real()*GyC + spectral_ty(i,j,k).imag()*GyR +
                       spectral_tz(i,j,k).real()*GzC + spectral_tz(i,j,k).imag()*GzR ;

           if (Lap < 1.0e-12) { // zero mode for no bulk motion
               spectral_dx(i,j,k) *= 0.0;
               spectral_dy(i,j,k) *= 0.0;
               spectral_dz(i,j,k) *= 0.0;
           }
           else {
               // Dilatational velocity 
               GpuComplex<Real> copy_dx((divR*GxR + divC*GxC) / Lap, 
                                        (divC*GxR - divR*GxC) / Lap);
               spectral_dx(i,j,k) = copy_dx;

               GpuComplex<Real> copy_dy((divR*GyR + divC*GyC) / Lap,
                                        (divC*GyR - divR*GyC) / Lap);
               spectral_dy(i,j,k) = copy_dy;

               GpuComplex<Real> copy_dz((divR*GzR + divC*GzC) / Lap,
                                        (divC*GzR - divR*GzC) / Lap);
               spectral_dz(i,j,k) = copy_dz;
               
               // Solenoidal velocity
               spectral_sx(i,j,k) = spectral_tx(i,j,k) - spectral_dx(i,j,k);
               spectral_sy(i,j,k) = spectral_ty(i,j,k) - spectral_dy(i,j,k); 
               spectral_sz(i,j,k) = spectral_tz(i,j,k) - spectral_dz(i,j,k);
           }
        });
    }

    ParallelDescriptor::Barrier();

    // Integrate K spectrum for velocities
    IntegrateKVelocity(spectral_field  , vel_onegrid, "turb_total"     ,scaling,step);
    IntegrateKVelocity(spectral_field_S, vel_onegrid, "turb_solenoidal",scaling,step);
    IntegrateKVelocity(spectral_field_D, vel_onegrid, "turb_dilational",scaling,step);

    // Create one-grid array to store IFFT velocities
    MultiFab vel_decomp_onegrid;
    vel_decomp_onegrid.define(ba_onegrid, dmap_onegrid, 6, 0);
    vel_decomp_onegrid.setVal(0.0);;
    
    // Inverse FFT to get solenoidal velocity 
    Vector<FFTplan> backward_planS;
    built_plan = false;
    for (int comp=0; comp<3; comp++) {
        if (!built_plan) {
            for (MFIter mfi(vel_decomp_onegrid); mfi.isValid(); ++mfi) {
                FFTplan fplan;

#ifdef AMREX_USE_CUDA // CUDA
                cufftResult result = cufftPlan3d(&fplan, fft_size[2], fft_size[1], fft_size[0], CUFFT_Z2D);
                if (result != CUFFT_SUCCESS) {
                    amrex::AllPrint() << " cufftplan3d forward failed! Error: "
                                      << cufftError(result) << "\n";
                }
                built_plan = true;
#elif AMREX_USE_HIP // HIP
                const std::size_t lengths[] = {std::size_t(fft_size[0]),std::size_t(fft_size[1]),std::size_t(fft_size[2])};
                rocfft_status result = rocfft_plan_create(&fplan, rocfft_placement_notinplace, 
                                                          rocfft_transform_type_real_inverse, rocfft_precision_double,
                                                          3, lengths, 1, nullptr);
                assert_rocfft_status("rocfft_plan_create", result);
                built_plan = true;
#else // host
                fplan = fftw_plan_dft_c2r_3d(fft_size[2], fft_size[1], fft_size[0],
                                             reinterpret_cast<FFTcomplex*>
                                             (spectral_field_S.back()->dataPtr(comp)),
                                             vel_decomp_onegrid[mfi].dataPtr(comp),
                                             FFTW_ESTIMATE);
#endif
               backward_planS.push_back(fplan);
            }
        }

        ParallelDescriptor::Barrier();

        // Backward Transform
        for (MFIter mfi(vel_decomp_onegrid); mfi.isValid(); ++mfi) {
            int i = mfi.LocalIndex();
#ifdef AMREX_USE_CUDA
            cufftSetStream(backward_planS[i], amrex::Gpu::gpuStream());
            cufftResult result = cufftExecZ2D(backward_planS[i],
                                              reinterpret_cast<FFTcomplex*>
                                                  (spectral_field_S[i]->dataPtr(comp)),
                                              vel_decomp_onegrid[mfi].dataPtr(comp));
            if (result != CUFFT_SUCCESS) {
                amrex::AllPrint() << " forward transform using cufftExec failed! Error: "
                                  << cufftError(result) << "\n";
	        }
#elif AMREX_USE_HIP
            rocfft_execution_info execinfo = nullptr;
            rocfft_status result = rocfft_execution_info_create(&execinfo);
            assert_rocfft_status("rocfft_execution_info_create", result);

            std::size_t buffersize = 0;
            result = rocfft_plan_get_work_buffer_size(backward_planS[i], &buffersize);
            assert_rocfft_status("rocfft_plan_get_work_buffer_size", result);

            void* buffer = amrex::The_Arena()->alloc(buffersize);
            result = rocfft_execution_info_set_work_buffer(execinfo, buffer, buffersize);
            assert_rocfft_status("rocfft_execution_info_set_work_buffer", result);

            result = rocfft_execution_info_set_stream(execinfo, amrex::Gpu::gpuStream());
            assert_rocfft_status("rocfft_execution_info_set_stream", result);

	        amrex::Real* vel_onegrid_ptr = vel_decomp_onegrid[mfi].dataPtr(comp);
	        FFTcomplex* spectral_field_ptr = reinterpret_cast<FFTcomplex*>(spectral_field_S[i]->dataPtr(comp));
            result = rocfft_execute(backward_planS[i],
                                    (void**) &vel_onegrid_ptr, // in
                                    (void**) &spectral_field_ptr, // out
                                    execinfo);
            assert_rocfft_status("rocfft_execute", result);
            amrex::Gpu::streamSynchronize();
            amrex::The_Arena()->free(buffer);
            result = rocfft_execution_info_destroy(execinfo);
            assert_rocfft_status("rocfft_execution_info_destroy", result);
#else
            fftw_execute(backward_planS[i]);
#endif
        }
    }

    // destroy fft plan
    for (int i = 0; i < backward_planS.size(); ++i) {
#ifdef AMREX_USE_CUDA
        cufftDestroy(backward_planS[i]);
#elif AMREX_USE_HIP
        rocfft_plan_destroy(backward_planS[i]);
#else
        fftw_destroy_plan(backward_planS[i]);
#endif
    }

    
    // Inverse FFT to get dilatational velocity 
    Vector<FFTplan> backward_planD;
    built_plan = false;
    for (int comp=0; comp<3; comp++) {
        if (!built_plan) {
            for (MFIter mfi(vel_decomp_onegrid); mfi.isValid(); ++mfi) {
                FFTplan fplan;

#ifdef AMREX_USE_CUDA // CUDA
                cufftResult result = cufftPlan3d(&fplan, fft_size[2], fft_size[1], fft_size[0], CUFFT_Z2D);
                if (result != CUFFT_SUCCESS) {
                    amrex::AllPrint() << " cufftplan3d forward failed! Error: "
                                      << cufftError(result) << "\n";
                }
                built_plan = true;
#elif AMREX_USE_HIP // HIP
                const std::size_t lengths[] = {std::size_t(fft_size[0]),std::size_t(fft_size[1]),std::size_t(fft_size[2])};
                rocfft_status result = rocfft_plan_create(&fplan, rocfft_placement_notinplace, 
                                                          rocfft_transform_type_real_inverse, rocfft_precision_double,
                                                          3, lengths, 1, nullptr);
                assert_rocfft_status("rocfft_plan_create", result);
                built_plan = true;
#else // host
                fplan = fftw_plan_dft_c2r_3d(fft_size[2], fft_size[1], fft_size[0],
                                             reinterpret_cast<FFTcomplex*>
                                             (spectral_field_D.back()->dataPtr(comp)),
                                             vel_decomp_onegrid[mfi].dataPtr(comp+3),
                                             FFTW_ESTIMATE);
#endif
               backward_planD.push_back(fplan);
            }
        }

        ParallelDescriptor::Barrier();

        // Backward Transform
        for (MFIter mfi(vel_decomp_onegrid); mfi.isValid(); ++mfi) {
            int i = mfi.LocalIndex();
#ifdef AMREX_USE_CUDA
            cufftSetStream(backward_planD[i], amrex::Gpu::gpuStream());
            cufftResult result = cufftExecZ2D(backward_planD[i],
                                              reinterpret_cast<FFTcomplex*>
                                                  (spectral_field_D[i]->dataPtr(comp)),
                                              vel_decomp_onegrid[mfi].dataPtr(comp+3));
            if (result != CUFFT_SUCCESS) {
                amrex::AllPrint() << " forward transform using cufftExec failed! Error: "
                                  << cufftError(result) << "\n";
	        }
#elif AMREX_USE_HIP
            rocfft_execution_info execinfo = nullptr;
            rocfft_status result = rocfft_execution_info_create(&execinfo);
            assert_rocfft_status("rocfft_execution_info_create", result);

            std::size_t buffersize = 0;
            result = rocfft_plan_get_work_buffer_size(backward_planD[i], &buffersize);
            assert_rocfft_status("rocfft_plan_get_work_buffer_size", result);

            void* buffer = amrex::The_Arena()->alloc(buffersize);
            result = rocfft_execution_info_set_work_buffer(execinfo, buffer, buffersize);
            assert_rocfft_status("rocfft_execution_info_set_work_buffer", result);

            result = rocfft_execution_info_set_stream(execinfo, amrex::Gpu::gpuStream());
            assert_rocfft_status("rocfft_execution_info_set_stream", result);

	        amrex::Real* vel_onegrid_ptr = vel_decomp_onegrid[mfi].dataPtr(comp+3);
	        FFTcomplex* spectral_field_ptr = reinterpret_cast<FFTcomplex*>(spectral_field_D[i]->dataPtr(comp));
            result = rocfft_execute(backward_planD[i],
                                    (void**) &vel_onegrid_ptr, // in
                                    (void**) &spectral_field_ptr, // out
                                    execinfo);
            assert_rocfft_status("rocfft_execute", result);
            amrex::Gpu::streamSynchronize();
            amrex::The_Arena()->free(buffer);
            result = rocfft_execution_info_destroy(execinfo);
            assert_rocfft_status("rocfft_execution_info_destroy", result);
#else
            fftw_execute(backward_planD[i]);
#endif
        }
    }

    // destroy fft plan
    for (int i = 0; i < backward_planD.size(); ++i) {
#ifdef AMREX_USE_CUDA
        cufftDestroy(backward_planD[i]);
#elif AMREX_USE_HIP
        rocfft_plan_destroy(backward_planD[i]);
#else
        fftw_destroy_plan(backward_planD[i]);
#endif
    }

    // copy into external multifab
    vel_decomp.ParallelCopy(vel_decomp_onegrid,0,0,6);
    vel_decomp.mult(1.0/sqrtnpts);
#endif // end heFFTe
}

#if defined(HEFFTE)
void IntegrateKScalar(BaseFab<GpuComplex<Real> >& spectral_field,
                      const std::string& name, const Real& scaling,
                      const Box& c_local_box,
                      const Real& sqrtnpts,
                      const int& step)

{
    int npts = n_cells[0]/2;
    
    Gpu::DeviceVector<Real> phisum_device(npts);
    Gpu::DeviceVector<int>  phicnt_device(npts);
    Gpu::HostVector<Real> phisum_host(npts);
    Real* phisum_ptr = phisum_device.dataPtr();  // pointer to data
    int*  phicnt_ptr = phicnt_device.dataPtr();  // pointer to data
    
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
      phisum_ptr[d] = 0.;
      phicnt_ptr[d] = 0;
    });

    Array4< GpuComplex<Real> > spectral = spectral_field.array();
    ParallelFor(c_local_box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        if (i <= n_cells[0]/2) { // only half of kx-domain
            int ki = i;
            int kj = j;
            int kk = k;
            if (j >= n_cells[1]/2) kj = n_cells[1]-j;
            if (k >= n_cells[2]/2) kk = n_cells[2]-k;

            Real dist = (ki*ki + kj*kj + kk*kk);
            dist = std::sqrt(dist);
        
            if ( dist <= n_cells[0]/2-0.5) {
                dist = dist+0.5;
                int cell = int(dist);
                Real real = spectral(i,j,k).real();
                Real imag = spectral(i,j,k).imag();
                Real cov  = scaling*(1.0/(sqrtnpts*sqrtnpts))*(real*real + imag*imag); 
                amrex::HostDevice::Atomic::Add(&(phisum_ptr[cell]), cov);
                amrex::HostDevice::Atomic::Add(&(phicnt_ptr[cell]),1);
            }
        }
    });
    
    ParallelDescriptor::Barrier();
        
    for (int d=1; d<npts; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_device[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_device[d]);
    }
        
    Real dk = 1.;
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0) {
        phisum_ptr[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_ptr[d];
        }
    });
    
    Gpu::copy(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());
    
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb;
        std::string turbBaseName = "turb_"+name;
        std::string turbName = Concatenate(turbBaseName,step,7);
        turbName += ".txt";
        
        turb.open(turbName);
        for (int d=1; d<npts; ++d) {
            turb << d << " " << phisum_host[d] << std::endl;
        }
        turb.close();
    }
}
#else
void IntegrateKScalar(Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > >& spectral_field,
                      const MultiFab& variables_onegrid,
                      const std::string& name,
                      const Real& scaling,
                      const Real& sqrtnpts,
                      const int& step)

{
    int npts = n_cells[0]/2;
    Gpu::DeviceVector<Real> phisum_device(npts);
    Gpu::DeviceVector<int>  phicnt_device(npts);

    Gpu::HostVector<Real> phisum_host(npts);
    
    Real* phisum_ptr = phisum_device.dataPtr();  // pointer to data
    int*  phicnt_ptr = phicnt_device.dataPtr();  // pointer to data

    // Integrate spectra over k-shells
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
      phisum_ptr[d] = 0.;
      phicnt_ptr[d] = 0;
    });

    for ( MFIter mfi(variables_onegrid,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.tilebox();

        Array4< GpuComplex<Real> > spectral = (*spectral_field[0]).array();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i <= bx.length(0)/2) { // only half of kx-domain
                int ki = i;
                int kj = j;
                int kk = k;
                if (j >= bx.length(1)/2) kj = bx.length(1)-j;
                if (k >= bx.length(2)/2) kk = bx.length(2)-k;

                Real dist = (ki*ki + kj*kj + kk*kk);
                dist = std::sqrt(dist);
            
                if ( dist <= n_cells[0]/2-0.5) {
                    dist = dist+0.5;
                    int cell = int(dist);
                    Real real = spectral(i,j,k).real();
                    Real imag = spectral(i,j,k).imag();
                    Real cov  = scaling*(1.0/(sqrtnpts*sqrtnpts))*(real*real + imag*imag); 
                    amrex::HostDevice::Atomic::Add(&(phisum_ptr[cell]), cov);
                    amrex::HostDevice::Atomic::Add(&(phicnt_ptr[cell]),1);
                }
            }
        });
    }
        
    for (int d=1; d<npts; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_device[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_device[d]);
    }
        
    Real dk = 1.;
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0) {
        phisum_ptr[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_ptr[d];
        }
    });
    
    Gpu::copy(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());
    
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb;
        std::string turbBaseName = "turb_"+name;
        std::string turbName = Concatenate(turbBaseName,step,7);
        turbName += ".txt";
        
        turb.open(turbName);
        for (int d=1; d<npts; ++d) {
            turb << d << " " << phisum_host[d] << std::endl;
        }
        turb.close();
    }
}
#endif

#if defined(HEFFTE)
void IntegrateKVelocity(BaseFab<GpuComplex<Real> >& spectral_field,
                        const std::string& name, const Real& scaling,
                        const Box& c_local_box,
                        const int& step)

{
    int npts = n_cells[0]/2;
    
    Gpu::DeviceVector<Real> phisum_device(npts);
    Gpu::DeviceVector<int>  phicnt_device(npts);
    Gpu::HostVector<Real> phisum_host(npts);
    Real* phisum_ptr = phisum_device.dataPtr();  // pointer to data
    int*  phicnt_ptr = phicnt_device.dataPtr();  // pointer to data
    
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
      phisum_ptr[d] = 0.;
      phicnt_ptr[d] = 0;
    });

    Array4< GpuComplex<Real> > spectralx = spectral_field.array(0,1);
    Array4< GpuComplex<Real> > spectraly = spectral_field.array(1,1);
    Array4< GpuComplex<Real> > spectralz = spectral_field.array(2,1);
    ParallelFor(c_local_box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        if (i <= n_cells[0]/2) { // only half of kx-domain
            int ki = i;
            int kj = j;
            int kk = k;
            if (j >= n_cells[1]/2) kj = n_cells[1]-j;
            if (k >= n_cells[2]/2) kk = n_cells[2]-k;

            Real dist = (ki*ki + kj*kj + kk*kk);
            dist = std::sqrt(dist);
        
            if ( dist <= n_cells[0]/2-0.5) {
                dist = dist+0.5;
                int cell = int(dist);
                Real real, imag, cov_x, cov_y, cov_z, cov;
                real = spectralx(i,j,k).real();
                imag = spectralx(i,j,k).imag();
                cov_x  = scaling*(real*real + imag*imag); 
                real = spectraly(i,j,k).real();
                imag = spectraly(i,j,k).imag();
                cov_y  = scaling*(real*real + imag*imag); 
                real = spectralz(i,j,k).real();
                imag = spectralz(i,j,k).imag();
                cov_z  = scaling*(real*real + imag*imag); 
                cov = cov_x + cov_y + cov_z;
                amrex::HostDevice::Atomic::Add(&(phisum_ptr[cell]), cov);
                amrex::HostDevice::Atomic::Add(&(phicnt_ptr[cell]),1);
            }
        }
    });
    
    ParallelDescriptor::Barrier();
        
    for (int d=1; d<npts; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_device[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_device[d]);
    }
        
    Real dk = 1.;
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0) {
        phisum_ptr[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_ptr[d];
        }
    });
    
    Gpu::copy(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());
    
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb;
        std::string turbBaseName = "turb_"+name;
        std::string turbName = Concatenate(turbBaseName,step,7);
        turbName += ".txt";
        
        turb.open(turbName);
        for (int d=1; d<npts; ++d) {
            turb << d << " " << phisum_host[d] << std::endl;
        }
        turb.close();
    }
}
#else
void IntegrateKVelocity(Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > >& spectral_field,
                        const MultiFab& vel_onegrid,
                        const std::string& name,
                        const Real& scaling,
                        const int& step)
{
    int npts = n_cells[0]/2;
    
    Gpu::DeviceVector<Real> phisum_device(npts);
    Gpu::DeviceVector<int>  phicnt_device(npts);
    Gpu::HostVector<Real> phisum_host(npts);
    Real* phisum_ptr = phisum_device.dataPtr();  // pointer to data
    int*  phicnt_ptr = phicnt_device.dataPtr();  // pointer to data

    // Integrate spectra over k-shells
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
      phisum_ptr[d] = 0.;
      phicnt_ptr[d] = 0;
    });

    for ( MFIter mfi(vel_onegrid,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.tilebox();

        Array4< GpuComplex<Real> > spectralx = (*spectral_field[0]).array(0,1);
        Array4< GpuComplex<Real> > spectraly = (*spectral_field[0]).array(1,1);
        Array4< GpuComplex<Real> > spectralz = (*spectral_field[0]).array(2,1);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i <= bx.length(0)/2) { // only half of kx-domain
                int ki = i;
                int kj = j;
                int kk = k;
                if (j >= bx.length(1)/2) kj = bx.length(1)-j;
                if (k >= bx.length(2)/2) kk = bx.length(2)-k;

                Real dist = (ki*ki + kj*kj + kk*kk);
                dist = std::sqrt(dist);
            
                if ( dist <= n_cells[0]/2-0.5) {
                    dist = dist+0.5;
                    int cell = int(dist);
                    Real real, imag, cov_x, cov_y, cov_z, cov;
                    real = spectralx(i,j,k).real();
                    imag = spectralx(i,j,k).imag();
                    cov_x  = scaling*(real*real + imag*imag); 
                    real = spectraly(i,j,k).real();
                    imag = spectraly(i,j,k).imag();
                    cov_y  = scaling*(real*real + imag*imag); 
                    real = spectralz(i,j,k).real();
                    imag = spectralz(i,j,k).imag();
                    cov_z  = scaling*(real*real + imag*imag); 
                    cov = cov_x + cov_y + cov_z;
                    amrex::HostDevice::Atomic::Add(&(phisum_ptr[cell]), cov);
                    amrex::HostDevice::Atomic::Add(&(phicnt_ptr[cell]),1);
                }
            }
        });
    }
        
    for (int d=1; d<npts; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_device[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_device[d]);
    }
        
    Real dk = 1.;
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0) {
        phisum_ptr[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_ptr[d];
        }
    });
    
    Gpu::copy(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());
    
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb;
        std::string turbBaseName = "turb_"+name;
        std::string turbName = Concatenate(turbBaseName,step,7);
        turbName += ".txt";
        
        turb.open(turbName);
        for (int d=1; d<npts; ++d) {
            turb << d << " " << phisum_host[d] << std::endl;
        }
        turb.close();
    }
}
#endif
