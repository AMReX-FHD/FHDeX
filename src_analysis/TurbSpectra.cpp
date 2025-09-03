#include <AMReX_FFT.H>

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
        "TurbSpectrumScalar: Must have one Box per MPI process");

    int ncomp = variables.nComp();

    Box domain = geom.Domain();
    auto npts = domain.numPts();
    Real sqrtnpts = std::sqrt(npts);

    amrex::FFT::R2C<Real,FFT::Direction::forward> r2c(geom.Domain());

    auto const& [cba, cdm] = r2c.getSpectralDataLayout();

    MultiFab cov(cba, cdm, ncomp, 0);

    for (int comp=0; comp<ncomp; ++comp) {

        MultiFab mf(variables, amrex::make_alias, comp, 1);
        cMultiFab cmf(cba, cdm, 1, 0);

        r2c.forward(mf,cmf);

        // Fill in the covariance multifab
        int comp_gpu = comp;
        Real sqrtnpts_gpu = sqrtnpts;
        Real scaling_i_gpu = scaling[comp];
        std::string name_gpu = var_names[comp];
        for (MFIter mfi(cov); mfi.isValid(); ++mfi) {
            Array4<Real> const& data = cov.array(mfi);
            Array4<const GpuComplex<Real> > spectral = cmf.const_array(mfi);
            const Box& bx = mfi.validbox();
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
        "TurbSpectrumVelDecomp: Must have one Box per MPI process");

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    Box domain = geom.Domain();
    auto npts = domain.numPts();
    Real sqrtnpts = std::sqrt(npts);

    // get box array and distribution map of vel
    DistributionMapping dm = vel.DistributionMap();
    BoxArray ba            = vel.boxArray();

    amrex::FFT::R2C<Real,FFT::Direction::both> r2c(geom.Domain());

    // box array and dmap for FFT
    auto const& [cba, cdm] = r2c.getSpectralDataLayout();

    // each MPI rank gets storage for its piece of the fft
    cMultiFab spectral_field_Tx(cba,cdm,1,0); // totalx
    cMultiFab spectral_field_Ty(cba,cdm,1,0); // totaly
    cMultiFab spectral_field_Tz(cba,cdm,1,0); // totalz
    cMultiFab spectral_field_Sx(cba,cdm,1,0); // solenoidalx
    cMultiFab spectral_field_Sy(cba,cdm,1,0); // solenoidaly
    cMultiFab spectral_field_Sz(cba,cdm,1,0); // solenoidalz
    cMultiFab spectral_field_Dx(cba,cdm,1,0); // dilatationalx
    cMultiFab spectral_field_Dy(cba,cdm,1,0); // dilatationaly
    cMultiFab spectral_field_Dz(cba,cdm,1,0); // dilatationalz

    // ForwardTransform
    // X
    {
        MultiFab vel_single(vel, amrex::make_alias, 0, 1);
        r2c.forward(vel_single,spectral_field_Tx);
    }
    // Y
    {
        MultiFab vel_single(vel, amrex::make_alias, 1, 1);
        r2c.forward(vel_single,spectral_field_Ty);
    }
    // Z
    {
        MultiFab vel_single(vel, amrex::make_alias, 2, 1);
        r2c.forward(vel_single,spectral_field_Tz);
    }

    // Decompose velocity field into solenoidal and dilatational
    for (MFIter mfi(spectral_field_Tx); mfi.isValid(); ++mfi) {
        Array4< GpuComplex<Real> > spectral_tx = spectral_field_Tx.array(mfi);
        Array4< GpuComplex<Real> > spectral_ty = spectral_field_Ty.array(mfi);
        Array4< GpuComplex<Real> > spectral_tz = spectral_field_Tz.array(mfi);
        Array4< GpuComplex<Real> > spectral_sx = spectral_field_Sx.array(mfi);
        Array4< GpuComplex<Real> > spectral_sy = spectral_field_Sy.array(mfi);
        Array4< GpuComplex<Real> > spectral_sz = spectral_field_Sz.array(mfi);
        Array4< GpuComplex<Real> > spectral_dx = spectral_field_Dx.array(mfi);
        Array4< GpuComplex<Real> > spectral_dy = spectral_field_Dy.array(mfi);
        Array4< GpuComplex<Real> > spectral_dz = spectral_field_Dz.array(mfi);
        const Box& bx = mfi.validbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
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
           }
           else { // conjugate
                amrex::Abort("check the code; i should not go beyond bx.length(0)/2");
           }

        });
    }

    MultiFab cov(cba, cdm, 3, 0); // total, solenoidal, dilatational

    // Fill in the covariance multifab
    Real scaling_gpu = scaling;
    for (MFIter mfi(cov); mfi.isValid(); ++mfi) {
        Array4<Real> const& data = cov.array(mfi);
        Array4<const GpuComplex<Real> > spec_tx = spectral_field_Tx.const_array(mfi);
        Array4<const GpuComplex<Real> > spec_ty = spectral_field_Ty.const_array(mfi);
        Array4<const GpuComplex<Real> > spec_tz = spectral_field_Tz.const_array(mfi);
        Array4<const GpuComplex<Real> > spec_sx = spectral_field_Sx.const_array(mfi);
        Array4<const GpuComplex<Real> > spec_sy = spectral_field_Sy.const_array(mfi);
        Array4<const GpuComplex<Real> > spec_sz = spectral_field_Sz.const_array(mfi);
        Array4<const GpuComplex<Real> > spec_dx = spectral_field_Dx.const_array(mfi);
        Array4<const GpuComplex<Real> > spec_dy = spectral_field_Dy.const_array(mfi);
        Array4<const GpuComplex<Real> > spec_dz = spectral_field_Dz.const_array(mfi);
        const Box& bx = mfi.validbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i <= n_cells[0]/2) {
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
            }
            else {
                amrex::Abort("check the code; i should not go beyond n_cells[0]/2");
            }
        });
    }

    // Integrate K spectrum for velocities
    IntegrateKVelocity(cov,"vel_total"     ,step,0);
    IntegrateKVelocity(cov,"vel_solenoidal",step,1);
    IntegrateKVelocity(cov,"vel_dilational",step,2);

    // inverse Fourier transform solenoidal and dilatational components
    {
        MultiFab vel_decomp_single(vel_decomp, amrex::make_alias, 0, 1);
        r2c.backward(spectral_field_Sx,vel_decomp_single);
    }
    {
        MultiFab vel_decomp_single(vel_decomp, amrex::make_alias, 1, 1);
        r2c.backward(spectral_field_Sy,vel_decomp_single);
    }
    {
        MultiFab vel_decomp_single(vel_decomp, amrex::make_alias, 2, 1);
        r2c.backward(spectral_field_Sz,vel_decomp_single);
    }
    {
        MultiFab vel_decomp_single(vel_decomp, amrex::make_alias, 3, 1);
        r2c.backward(spectral_field_Dx,vel_decomp_single);
    }
    {
        MultiFab vel_decomp_single(vel_decomp, amrex::make_alias, 4, 1);
        r2c.backward(spectral_field_Dy,vel_decomp_single);
    }
    {
        MultiFab vel_decomp_single(vel_decomp, amrex::make_alias, 5, 1);
        r2c.backward(spectral_field_Dz,vel_decomp_single);
    }

    vel_decomp.mult(1.0/sqrtnpts);

}

void IntegrateKScalar(const MultiFab& cov_mag,
                            const std::string& name,
                            const int& step,
                            const int& comp)

{
    int npts = n_cells[0]/2;

    Gpu::DeviceVector<Real> phisum_device(npts, 0);
    Gpu::DeviceVector<int>  phicnt_device(npts, 0);
    Real* phisum_ptr = phisum_device.dataPtr();  // pointer to data
    int*  phicnt_ptr = phicnt_device.dataPtr();  // pointer to data

    int comp_gpu = comp;
    int nx = n_cells[0];
    int ny = n_cells[1];
    int nz = n_cells[2];
    for ( MFIter mfi(cov_mag,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        const Array4<const Real> & cov = cov_mag.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i <= n_cells[0]/2) {
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
            }
            else {
                amrex::Abort("check the code; i should not go beyond n_cells[0]/2");
            }
        });
    }

    Gpu::HostVector<Real> phisum_host(npts);
    Gpu::HostVector<int>  phicnt_host(npts);
    Gpu::copyAsync(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());
    Gpu::copyAsync(Gpu::deviceToHost, phicnt_device.begin(), phicnt_device.end(), phicnt_host.begin());
    Gpu::streamSynchronize();

    ParallelDescriptor::ReduceRealSum(phisum_host.dataPtr(),npts);
    ParallelDescriptor::ReduceIntSum(phicnt_host.dataPtr(),npts);

    Real dk = 1.;
    for (int d = 1; d < npts; ++d) {
        phisum_host[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_host[d];
    }

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

    Gpu::DeviceVector<Real> phisum_device(npts, 0);
    Gpu::DeviceVector<int>  phicnt_device(npts, 0);
    Real* phisum_ptr = phisum_device.dataPtr();  // pointer to data
    int*  phicnt_ptr = phicnt_device.dataPtr();  // pointer to data

    int comp_gpu = comp;
    int nx = n_cells[0];
    int ny = n_cells[1];
    int nz = n_cells[2];
    for ( MFIter mfi(cov_mag,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        const Array4<const Real> & cov = cov_mag.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i <= n_cells[0]/2) {
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
            }
            else {
                amrex::Abort("check the code; i should not go beyond n_cells[0]/2");
            }
        });
    }

    Gpu::HostVector<Real> phisum_host(npts);
    Gpu::HostVector<int>  phicnt_host(npts);
    Gpu::copyAsync(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());
    Gpu::copyAsync(Gpu::deviceToHost, phicnt_device.begin(), phicnt_device.end(), phicnt_host.begin());
    Gpu::streamSynchronize();

    ParallelDescriptor::ReduceRealSum(phisum_host.dataPtr(),npts);
    ParallelDescriptor::ReduceIntSum(phicnt_host.dataPtr(),npts);

    Real dk = 1.;
    for (int d = 1; d < npts; ++d) {
        phisum_host[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_host[d];
    }

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
