#include "TurbSpectra.H"
#include "common_functions.H"

#include <AMReX_MultiFabUtil.H>
#include "AMReX_PlotFileUtil.H"
#include "AMReX_BoxArray.H"

void TurbSpectrumScalar(const MultiFab& variables, 
                              const amrex::Geometry& geom, 
                              const int& step, 
                              const amrex::Vector<amrex::Real>& scaling,
                              const amrex::Vector< std::string >& var_names)
{
    BL_PROFILE_VAR("TurbSpectrumScalar()",TurbSpectrumScalar);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(variables.nComp() == var_names.size(), 
        "TurbSpectrumScalar: must have same number variable names as components of input MultiFab");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(variables.nComp() == scaling.size(), 
        "TurbSpectrumScalar: must have same number variable scaling as components of input MultiFab");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(variables.local_size() == 1, 
        "TurbSpectrumScalar: Must have one Box per MPI process when using heFFTe");

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

    // BOX ARRAY TO STORE COVARIANCE MATRIX IN A MFAB
    // create a BoxArray containing the fft boxes
    // by construction, these boxes correlate to the associated spectral_data
    // this we can copy the spectral data into this multifab since we know they are owned by the same MPI rank
    BoxArray fft_ba;
    {
        BoxList bl;
        bl.reserve(ba.size());

        for (int i = 0; i < ba.size(); ++i) {
            Box b = ba[i];

            Box r_box = b;
            Box c_box = amrex::coarsen(r_box, IntVect(AMREX_D_DECL(2,1,1)));

            // this avoids overlap for the cases when one or more r_box's
            // have an even cell index in the hi-x cell
            if (c_box.bigEnd(0) * 2 == r_box.bigEnd(0)) {
                c_box.setBig(0,c_box.bigEnd(0)-1);
            }

            // increase the size of boxes touching the hi-x domain by 1 in x
            // this is an (Nx x Ny x Nz) -> (Nx/2+1 x Ny x Nz) real-to-complex sizing
            if (b.bigEnd(0) == geom.Domain().bigEnd(0)) {
                c_box.growHi(0,1);
            }
            bl.push_back(c_box);

        }
        fft_ba.define(std::move(bl));
    }
    MultiFab cov(fft_ba, dm, ncomp, 0);
		
    // each MPI rank gets storage for its piece of the fft
    BaseFab<GpuComplex<Real> > spectral_field(c_local_box, 1, The_Device_Arena());
    MultiFab variables_single(ba, dm, 1, 0);
    using heffte_complex = typename heffte::fft_output<Real>::type;
    
    int r2c_direction = 0;
    for (int comp=0; comp<ncomp; ++comp) {    
#if defined(HEFFTE_CUFFT)
    	  heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
    	  heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
    	  heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        heffte_complex* spectral_data = (heffte_complex*) spectral_field.dataPtr();        
        variables_single.ParallelCopy(variables,comp,0,1);
        fft.forward(variables_single[local_boxid].dataPtr(),spectral_data);        
        Gpu::streamSynchronize();

        // Fill in the covariance multifab
        int comp_gpu = comp;
        Real sqrtnpts_gpu = sqrtnpts;
        Real scaling_i_gpu = scaling[comp];
        std::string name_gpu = var_names[comp];
        for (MFIter mfi(cov); mfi.isValid(); ++mfi) {
            Array4<Real> const& data = cov.array(mfi);
            Array4<const GpuComplex<Real> > spectral = spectral_field.const_array();
            const Box& bx = mfi.fabbox();
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real re = spectral(i,j,k).real();
                Real im = spectral(i,j,k).imag();
                data(i,j,k,comp_gpu) = (re*re + im*im)/(sqrtnpts_gpu*sqrtnpts_gpu*scaling_i_gpu);
            });
        }
        
        // Integrate spectra over k-shells
        IntegrateKScalar(cov,name_gpu,step,comp_gpu);
    }
}

void TurbSpectrumVelDecomp(const MultiFab& vel,
                                 MultiFab& vel_decomp,
                                 const amrex::Geometry& geom,
                                 const int& step,
                                 const amrex::Real& scaling,
                                 const amrex::Vector< std::string >& var_names)
{
    BL_PROFILE_VAR("TurbSpectrumVelDecomp()",TurbSpectrumVelDecomp);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vel.nComp() == 3, 
        "TurbSpectrumVelDecomp: must have 3 components of input vel MultiFab");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(var_names.size() == 3, 
        "TurbSpectrumVelDecomp: must have 3 names for output vel spectra (total, solenoidal, dilatational");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vel.local_size() == 1, 
        "TurbSpectrumVelDecomp: Must have one Box per MPI process when using heFFTe");

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
    BaseFab<GpuComplex<Real> > spectral_field_Tx(c_local_box, 1, The_Device_Arena()); // totalx
    BaseFab<GpuComplex<Real> > spectral_field_Ty(c_local_box, 1, The_Device_Arena()); // totaly
    BaseFab<GpuComplex<Real> > spectral_field_Tz(c_local_box, 1, The_Device_Arena()); // totalz
    BaseFab<GpuComplex<Real> > spectral_field_Sx(c_local_box, 1, The_Device_Arena()); // solenoidalx
    BaseFab<GpuComplex<Real> > spectral_field_Sy(c_local_box, 1, The_Device_Arena()); // solenoidaly
    BaseFab<GpuComplex<Real> > spectral_field_Sz(c_local_box, 1, The_Device_Arena()); // solenoidalz
    BaseFab<GpuComplex<Real> > spectral_field_Dx(c_local_box, 1, The_Device_Arena()); // dilatationalx
    BaseFab<GpuComplex<Real> > spectral_field_Dy(c_local_box, 1, The_Device_Arena()); // dilatationaly
    BaseFab<GpuComplex<Real> > spectral_field_Dz(c_local_box, 1, The_Device_Arena()); // dilatationalz
	  MultiFab vel_single(ba, dm, 1, 0);
    	
    int r2c_direction = 0;
    
    // ForwardTransform
    // X
    using heffte_complex = typename heffte::fft_output<Real>::type;
    {
#if defined(HEFFTE_CUFFT)
      heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
      heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
      heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
      ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
      {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
      {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
      {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
      r2c_direction, ParallelDescriptor::Communicator());
	    
      vel_single.ParallelCopy(vel, 0, 0, 1);
	    heffte_complex* spectral_data = (heffte_complex*) spectral_field_Tx.dataPtr();
      fft.forward(vel_single[local_boxid].dataPtr(),spectral_data);
    }
    // Y
    {
#if defined(HEFFTE_CUFFT)
      heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
      heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
      heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
      ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
      {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
      {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
      {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
      r2c_direction, ParallelDescriptor::Communicator());
	    
      vel_single.ParallelCopy(vel, 1, 0, 1);
      heffte_complex* spectral_data = (heffte_complex*) spectral_field_Ty.dataPtr();
      fft.forward(vel_single[local_boxid].dataPtr(),spectral_data);
    }
    // Z
    {
#if defined(HEFFTE_CUFFT)
      heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
      heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
      heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
      ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
      {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
      {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
      {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
      r2c_direction, ParallelDescriptor::Communicator());
	    
      vel_single.ParallelCopy(vel, 2, 0, 1);
      heffte_complex* spectral_data = (heffte_complex*) spectral_field_Tz.dataPtr();
      fft.forward(vel_single[local_boxid].dataPtr(),spectral_data);
    }
    
    Gpu::streamSynchronize();
    
    // Decompose velocity field into solenoidal and dilatational
    Array4< GpuComplex<Real> > spectral_tx = spectral_field_Tx.array();
    Array4< GpuComplex<Real> > spectral_ty = spectral_field_Ty.array();
    Array4< GpuComplex<Real> > spectral_tz = spectral_field_Tz.array();
    Array4< GpuComplex<Real> > spectral_sx = spectral_field_Sx.array();
    Array4< GpuComplex<Real> > spectral_sy = spectral_field_Sy.array();
    Array4< GpuComplex<Real> > spectral_sz = spectral_field_Sz.array();
    Array4< GpuComplex<Real> > spectral_dx = spectral_field_Dx.array();
    Array4< GpuComplex<Real> > spectral_dy = spectral_field_Dy.array();
    Array4< GpuComplex<Real> > spectral_dz = spectral_field_Dz.array();
    ParallelFor(c_local_box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {

       int nx = n_cells[0]; 
       int ny = n_cells[1]; 
       int nz = n_cells[2];

       Real GxR = 0.0, GxC = 0.0, GyR = 0.0, GyC = 0.0, GzR = 0.0, GzC = 0.0;
       
       if (i <= nx/2) { 
           
           // Get the wavevector
           int ki = i;
           int kj = j;
           if (j >= ny/2) kj = ny - j;
           int kk = k;
           if (k >= nz/2) kk = nz - k;

           // Gradient Operators
           GxR = (cos(2.0*M_PI*ki/nx)-1.0)/dx[0];
           GxC = (sin(2.0*M_PI*ki/nx)-0.0)/dx[0];
           GyR = (cos(2.0*M_PI*kj/ny)-1.0)/dx[1];
           GyC = (sin(2.0*M_PI*kj/ny)-0.0)/dx[1];
           GzR = (cos(2.0*M_PI*kk/nz)-1.0)/dx[2];
           GzC = (sin(2.0*M_PI*kk/nz)-0.0)/dx[2];
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
       }
           
       // Solenoidal velocity
       spectral_sx(i,j,k) = spectral_tx(i,j,k) - spectral_dx(i,j,k);
       spectral_sy(i,j,k) = spectral_ty(i,j,k) - spectral_dy(i,j,k); 
       spectral_sz(i,j,k) = spectral_tz(i,j,k) - spectral_dz(i,j,k);

    });

    Gpu::streamSynchronize();
    
    // BOX ARRAY TO STORE COVARIANCE MATRIX IN A MFAB
    // create a BoxArray containing the fft boxes
    // by construction, these boxes correlate to the associated spectral_data
    // this we can copy the spectral data into this multifab since we know they are owned by the same MPI rank
    BoxArray fft_ba;
    {
        BoxList bl;
        bl.reserve(ba.size());

        for (int i = 0; i < ba.size(); ++i) {
            Box b = ba[i];

            Box r_box = b;
            Box c_box = amrex::coarsen(r_box, IntVect(AMREX_D_DECL(2,1,1)));

            // this avoids overlap for the cases when one or more r_box's
            // have an even cell index in the hi-x cell
            if (c_box.bigEnd(0) * 2 == r_box.bigEnd(0)) {
                c_box.setBig(0,c_box.bigEnd(0)-1);
            }

            // increase the size of boxes touching the hi-x domain by 1 in x
            // this is an (Nx x Ny x Nz) -> (Nx/2+1 x Ny x Nz) real-to-complex sizing
            if (b.bigEnd(0) == geom.Domain().bigEnd(0)) {
                c_box.growHi(0,1);
            }
            bl.push_back(c_box);

        }
        fft_ba.define(std::move(bl));
    }
    MultiFab cov(fft_ba, dm, 3, 0); // total, solenoidal, dilatational
    
    // Fill in the covariance multifab
    Real sqrtnpts_gpu = sqrtnpts;
    Real scaling_gpu = scaling;
    for (MFIter mfi(cov); mfi.isValid(); ++mfi) {
        Array4<Real> const& data = cov.array(mfi);
        Array4<const GpuComplex<Real> > spec_tx = spectral_field_Tx.const_array();
        Array4<const GpuComplex<Real> > spec_ty = spectral_field_Ty.const_array();
        Array4<const GpuComplex<Real> > spec_tz = spectral_field_Tz.const_array();
        Array4<const GpuComplex<Real> > spec_sx = spectral_field_Sx.const_array();
        Array4<const GpuComplex<Real> > spec_sy = spectral_field_Sy.const_array();
        Array4<const GpuComplex<Real> > spec_sz = spectral_field_Sz.const_array();
        Array4<const GpuComplex<Real> > spec_dx = spectral_field_Dx.const_array();
        Array4<const GpuComplex<Real> > spec_dy = spectral_field_Dy.const_array();
        Array4<const GpuComplex<Real> > spec_dz = spectral_field_Dz.const_array();
        const Box& bx = mfi.fabbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real re_x, re_y, re_z, im_x, im_y, im_z;
            
            re_x = spec_tx(i,j,k).real();
            im_x = spec_tx(i,j,k).imag();
            re_y = spec_ty(i,j,k).real();
            im_y = spec_ty(i,j,k).imag();
            re_z = spec_tz(i,j,k).real();
            im_z = spec_tz(i,j,k).imag();
            data(i,j,k,0) = (re_x*re_x + im_x*im_x + 
                             re_y*re_y + im_y*im_y +
                             re_z*re_z + im_z*im_z)/(scaling_gpu);
            re_x = spec_sx(i,j,k).real();
            im_x = spec_sx(i,j,k).imag();
            re_y = spec_sy(i,j,k).real();
            im_y = spec_sy(i,j,k).imag();
            re_z = spec_sz(i,j,k).real();
            im_z = spec_sz(i,j,k).imag();
            data(i,j,k,1) = (re_x*re_x + im_x*im_x + 
                             re_y*re_y + im_y*im_y +
                             re_z*re_z + im_z*im_z)/(scaling_gpu);
            re_x = spec_dx(i,j,k).real();
            im_x = spec_dx(i,j,k).imag();
            re_y = spec_dy(i,j,k).real();
            im_y = spec_dy(i,j,k).imag();
            re_z = spec_dz(i,j,k).real();
            im_z = spec_dz(i,j,k).imag();
            data(i,j,k,2) = (re_x*re_x + im_x*im_x + 
                             re_y*re_y + im_y*im_y +
                             re_z*re_z + im_z*im_z)/(scaling_gpu);
        });
    }

    // Integrate K spectrum for velocities
    IntegrateKVelocity(cov,"vel_total"     ,step,0);
    IntegrateKVelocity(cov,"vel_solenoidal",step,1);
    IntegrateKVelocity(cov,"vel_dilational",step,2);
    
	  MultiFab vel_decomp_single(ba, dm, 1, 0);
    // inverse Fourier transform solenoidal and dilatational components 
    {
#if defined(HEFFTE_CUFFT)
      heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
      heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
      heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
      ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
      {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
      {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
      {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
      r2c_direction, ParallelDescriptor::Communicator());
	    
      heffte_complex* spectral_data = (heffte_complex*) spectral_field_Sx.dataPtr();
      fft.backward(spectral_data, vel_decomp_single[local_boxid].dataPtr());
    
      Gpu::streamSynchronize();
      vel_decomp.ParallelCopy(vel_decomp_single, 0, 0, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
      heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
      heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
      heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
      ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
      {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
      {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
      {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
      r2c_direction, ParallelDescriptor::Communicator());
	    
      heffte_complex* spectral_data = (heffte_complex*) spectral_field_Sy.dataPtr();
      fft.backward(spectral_data, vel_decomp_single[local_boxid].dataPtr());
    
      Gpu::streamSynchronize();
      vel_decomp.ParallelCopy(vel_decomp_single, 0, 1, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
      heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
      heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
      heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
      ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
      {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
      {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
      {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
      r2c_direction, ParallelDescriptor::Communicator());
	    
      heffte_complex* spectral_data = (heffte_complex*) spectral_field_Sz.dataPtr();
      fft.backward(spectral_data, vel_decomp_single[local_boxid].dataPtr());
    
      Gpu::streamSynchronize();
      vel_decomp.ParallelCopy(vel_decomp_single, 0, 2, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
      heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
      heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
      heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
      ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
      {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
      {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
      {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
      r2c_direction, ParallelDescriptor::Communicator());
	    
      heffte_complex* spectral_data = (heffte_complex*) spectral_field_Dx.dataPtr();
      fft.backward(spectral_data, vel_decomp_single[local_boxid].dataPtr());
    
      Gpu::streamSynchronize();
      vel_decomp.ParallelCopy(vel_decomp_single, 0, 3, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
      heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
      heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
      heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
      ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
      {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
      {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
      {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
      r2c_direction, ParallelDescriptor::Communicator());
	    
      heffte_complex* spectral_data = (heffte_complex*) spectral_field_Dy.dataPtr();
      fft.backward(spectral_data, vel_decomp_single[local_boxid].dataPtr());
    
      Gpu::streamSynchronize();
      vel_decomp.ParallelCopy(vel_decomp_single, 0, 4, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
      heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
      heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
      heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
      ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
      {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
      {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
      {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
      r2c_direction, ParallelDescriptor::Communicator());
	    
      heffte_complex* spectral_data = (heffte_complex*) spectral_field_Dz.dataPtr();
      fft.backward(spectral_data, vel_decomp_single[local_boxid].dataPtr());
    
      Gpu::streamSynchronize();
      vel_decomp.ParallelCopy(vel_decomp_single, 0, 5, 1);
    }

    
    vel_decomp.mult(1.0/sqrtnpts);

}

void IntegrateKScalar(const MultiFab& cov_mag,
                      	    const std::string& name,
                            const int& step,
                            const int& comp)

{
    int npts = n_cells[0]/2;
    
    Gpu::DeviceVector<Real> phisum_device(npts);
    Gpu::DeviceVector<int>  phicnt_device(npts);
//    Gpu::HostVector<Real> phisum_host(npts);
//    Gpu::HostVector<int>  phicnt_host(npts);
    
    Gpu::HostVector<Real> phisum_host(npts);
    
    Real* phisum_ptr = phisum_device.dataPtr();  // pointer to data
    int*  phicnt_ptr = phicnt_device.dataPtr();  // pointer to data
    
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
      phisum_ptr[d] = 0.;
      phicnt_ptr[d] = 0;
    });
//    for (int d=0; d<npts; ++d) {
//	phisum_host[d] = 0.;
//	phicnt_host[d] = 0;
//    }

    int comp_gpu = comp;
    int nx = n_cells[0]; 
    int ny = n_cells[1]; 
    int nz = n_cells[2];
    for ( MFIter mfi(cov_mag,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.tilebox();

        const Array4<const Real> & cov = cov_mag.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ki = i; 
            int kj = j;
            if (j >= ny/2) kj = ny - j;
            int kk = k;
            if (k >= nz/2) kk = nz - k;

            Real dist = (ki*ki + kj*kj + kk*kk);
            dist = std::sqrt(dist);
            
            if ( dist <=  n_cells[0]/2-0.5) {
	              dist = dist+0.5;
                int cell = int(dist);
		            amrex::Gpu::Atomic::Add(&(phisum_ptr[cell]), cov(i,j,k,comp_gpu));
		            amrex::Gpu::Atomic::Add(&(phicnt_ptr[cell]),1);
            }
        });
    }
    
    Gpu::streamSynchronize();
        
    ParallelDescriptor::ReduceRealSum(phisum_device.dataPtr(),npts);
    ParallelDescriptor::ReduceIntSum(phicnt_device.dataPtr(),npts);
        
    Real dk = 1.;
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0) {
        phisum_ptr[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_ptr[d];
        }
    });
    
    Gpu::copyAsync(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());
    Gpu::streamSynchronize();
    
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

void IntegrateKVelocity(const MultiFab& cov_mag,
                              const std::string& name,
                              const int& step,
                              const int& comp)

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
    
    int comp_gpu = comp;
    int nx = n_cells[0]; 
    int ny = n_cells[1]; 
    int nz = n_cells[2];
    for ( MFIter mfi(cov_mag,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.tilebox();

        const Array4<const Real> & cov = cov_mag.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ki = i; 
            int kj = j;
            if (j >= ny/2) kj = ny - j;
            int kk = k;
            if (k >= nz/2) kk = nz - k;

            Real dist = (ki*ki + kj*kj + kk*kk);
            dist = std::sqrt(dist);
            
            if ( dist <=  n_cells[0]/2-0.5) {
	              dist = dist+0.5;
                int cell = int(dist);
		            amrex::Gpu::Atomic::Add(&(phisum_ptr[cell]), cov(i,j,k,comp_gpu));
		            amrex::Gpu::Atomic::Add(&(phicnt_ptr[cell]),1);
            }
        });
    }
    
    Gpu::streamSynchronize();

    ParallelDescriptor::ReduceRealSum(phisum_device.dataPtr(),npts);
    ParallelDescriptor::ReduceIntSum(phicnt_device.dataPtr(),npts);
        
    Real dk = 1.;
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0) {
        phisum_ptr[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_ptr[d];
        }
    });
    
    Gpu::copyAsync(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());
    Gpu::streamSynchronize();
    
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

