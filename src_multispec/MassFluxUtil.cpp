#include "multispec_functions.H"

void ComputeMolconcMolmtot(const MultiFab& rho_in,
    const MultiFab& rhotot_in,
    MultiFab& molarconc_in,
    MultiFab& molmtot_in)
{

    BL_PROFILE_VAR("ComputeMolconcMolmtot()",ComputeMolconcMolmtot);

    int ng = molarconc_in.nGrow();

    // Loop over boxes
    for (MFIter mfi(rho_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        const Array4<const Real>& rho = rho_in.array(mfi);
        const Array4<const Real>& rhotot = rhotot_in.array(mfi);
        const Array4<      Real>& molarconc = molarconc_in.array(mfi);
        const Array4<      Real>& molmtot = molmtot_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real, MAX_SPECIES> RhoN;
            GpuArray<Real, MAX_SPECIES> MolarConcN;

            for (int n=0; n<nspecies; ++n ){
                RhoN[n] = rho(i,j,k,n);
            }

            ComputeMolconcMolmtotLocal(nspecies,
                            molmass,
                            RhoN,
                            rhotot(i,j,k),
                            MolarConcN,
                            molmtot(i,j,k));

            for (int n=0; n<nspecies; ++n ){
                molarconc(i,j,k,n) = MolarConcN[n] ;
            }

        });
    }
}

void ComputeMassfrac(const MultiFab& rho_in,
    const MultiFab& rhotot_in,
    MultiFab& massfrac_in)
{

    BL_PROFILE_VAR("ComputeMassfrac()",ComputeMassfrac);

    int ng = massfrac_in.nGrow();

    // Loop over boxes
    for (MFIter mfi(rho_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        const Array4<const Real>& rho = rho_in.array(mfi);
        const Array4<const Real>& rhotot = rhotot_in.array(mfi);
        const Array4<      Real>& massfrac = massfrac_in.array(mfi);

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {

            massfrac(i,j,k,n) = rho(i,j,k,n)/rhotot(i,j,k);

        });
    }
}

void ComputeGamma(const MultiFab& molarconc_in,
    const MultiFab& Hessian_in,
    MultiFab& Gamma_in)
{

    BL_PROFILE_VAR("ComputeGamma()",ComputeGamma);

    int ng = Gamma_in.nGrow();

    // Loop over boxes
    for (MFIter mfi(Gamma_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        const Array4<const Real>& molarconc = molarconc_in.array(mfi);
        const Array4<const Real>& Hessian = Hessian_in.array(mfi);
        const Array4<      Real>& Gamma = Gamma_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real, MAX_SPECIES> MolarConcN;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> GammaN;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> HessianN;

            // Read MultiFab data into arrays
            for (int n=0; n<nspecies; ++n ){
                MolarConcN[n] = molarconc(i,j,k,n);

                for (int m=0; m<nspecies; ++m){
                    GammaN(m+1,n+1) = Gamma(i,j,k,n*nspecies+m);
                    HessianN(m+1,n+1) = Hessian(i,j,k,n*nspecies+m);
                }
            }

            ComputeGammaLocal(MolarConcN, HessianN, GammaN);

            // Write back to MultiFab
            for (int n=0; n<nspecies; ++n ){
                for (int m=0; m<nspecies; ++m){
                    Gamma(i,j,k,n*nspecies+m) = GammaN(m+1,n+1);
                }
            }
        });
    }
}

void ComputeFHGamma(const MultiFab& massfrac_in,
    MultiFab& Gamma_in)
{

    BL_PROFILE_VAR("ComputeFHGamma()",ComputeFHGamma);

    int ng = Gamma_in.nGrow();

    // Loop over boxes
    for (MFIter mfi(Gamma_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        const Array4<const Real>& massfrac = massfrac_in.array(mfi);
        const Array4<      Real>& Gamma = Gamma_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real, MAX_SPECIES> MassfracN;
            Array2D<Real, 0, MAX_SPECIES-1, 0, MAX_SPECIES-1> GammaN;

            // Read MultiFab data into arrays
            for (int n=0; n<nspecies; ++n ){
                MassfracN[n] = massfrac(i,j,k,n);

                for (int m=0; m<nspecies; ++m){
                    GammaN(n,m) = Gamma(i,j,k,n*nspecies+m);
                }
            }

            ComputeFHGammaLocal(MassfracN, GammaN);

            // Write back to MultiFab
            for (int n=0; n<nspecies; ++n ){
                for (int m=0; m<nspecies; ++m){
                    // Gamma(i,j,k,n*nspecies+m) = GammaN(n,m);
                    Gamma(i,j,k,m*nspecies+n) = GammaN(n,m);
                }
            }
        });
    }
}


void ComputeRhoWChi(const MultiFab& rho_in,
    const MultiFab& rhotot_in,
    const MultiFab& molarconc_in,
    MultiFab& rhoWchi_in,
    const MultiFab& D_bar_in)
{
    BL_PROFILE_VAR("ComputeRhoWChi()",ComputeRhoWChi);

    int ng = rhoWchi_in.nGrow();

    // Loop over boxes
    for (MFIter mfi(rhoWchi_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        const Array4<const Real>& rho = rho_in.array(mfi);
        const Array4<const Real>& rhotot = rhotot_in.array(mfi);
        const Array4<const Real>& molarconc = molarconc_in.array(mfi);
        const Array4<      Real>& rhoWchi = rhoWchi_in.array(mfi);
        const Array4<const Real>& D_bar = D_bar_in.array(mfi);


        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real, MAX_SPECIES> rhoN;
            GpuArray<Real, MAX_SPECIES> MolarConcN;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> rhoWchiN;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> D_barN;

            // Read MultiFab data into arrays
            for (int n=0; n<nspecies; ++n){

                rhoN[n] = rho(i,j,k,n);
                MolarConcN[n] = molarconc(i,j,k,n);
                for (int m=0; m<nspecies; ++m){
                    rhoWchiN(m+1,n+1) = rhoWchi(i,j,k,n*nspecies+m);
                    D_barN(m+1,n+1) = D_bar(i,j,k,n*nspecies+m);
                }
            }

            ComputeRhoWChiLocal(rhoN,
                            rhotot(i,j,k),
                            MolarConcN,
                            rhoWchiN,
                            D_barN,
                            molmass);

            // Write back to MultiFab
            for (int n=0; n<nspecies; ++n ){
                for (int m=0; m<nspecies; ++m){
                    rhoWchi(i,j,k,n*nspecies+m) = rhoWchiN(m+1,n+1);
                }
            }
        });

    }
}

void ComputeZetaByTemp(const MultiFab& molarconc_in,
    const MultiFab& D_bar_in,
    const MultiFab& Temp_in,
    MultiFab& zeta_by_Temp_in,
    const MultiFab& D_therm_in)
{
    BL_PROFILE_VAR("ComputeZetaByTemp()",ComputeZetaByTemp);

    int ng = zeta_by_Temp_in.nGrow();

    // Loop over boxes
    for (MFIter mfi(zeta_by_Temp_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        const Array4<const Real>& molarconc = molarconc_in.array(mfi);
        const Array4<const Real>& D_bar = D_bar_in.array(mfi);
        const Array4<const Real>& Temp = Temp_in.array(mfi);
        const Array4<      Real>& zeta_by_Temp = zeta_by_Temp_in.array(mfi);
        const Array4<const Real>& D_therm = D_therm_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real, MAX_SPECIES> MolarConcN;
            GpuArray<Real, MAX_SPECIES> ZetaByTemp;
            GpuArray<Real, MAX_SPECIES> DTherm;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> DBarN;

            // Read MultiFab data into arrays
            for (int n=0; n<nspecies; ++n){
                MolarConcN[n] = molarconc(i,j,k,n);
                ZetaByTemp[n] = zeta_by_Temp(i,j,k,n);
                DTherm[n] = D_therm(i,j,k,n);
                for (int m=0; m<nspecies; ++m){
                    DBarN(m+1,n+1) = D_bar(i,j,k,n*nspecies+m);
                }
            }

            ComputeZetaByTempLocal( MolarConcN,
                                    DBarN,
                                    Temp(i,j,k),
                                    ZetaByTemp,
                                    DTherm);

            //write data back to MultiFabs
            for (int n=0; n<nspecies; ++n ){
                zeta_by_Temp(i,j,k,n) = ZetaByTemp[n] ;
            }

        });
    }
}

void ComputeSqrtLonsagerFC(const MultiFab& rho_in,
    const MultiFab& rhotot_in,
    std::array< MultiFab, AMREX_SPACEDIM >& sqrtLonsager_fc,
    const Geometry& geom)
{
    BL_PROFILE_VAR("ComputeSqrtLonsagerFC()",ComputeSqrtLonsagerFC);

    const Real* dx_old = geom.CellSize();

    // for GPU later
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // Loop over boxes
    for (MFIter mfi(rho_in); mfi.isValid(); ++mfi) {
        // note: tiling or GPU-ing requires nodal tileboxes and changes to
        // loop indices in fortran

        // Create cell-centered box
        const Box& validBox = mfi.validbox();

        const Array4<const Real>& rho = rho_in.array(mfi);
        const Array4<const Real>& rhotot = rhotot_in.array(mfi);

        AMREX_D_TERM(const Array4<      Real>& sqrtLOnsager_X = sqrtLonsager_fc[0].array(mfi);,
                     const Array4<      Real>& sqrtLOnsager_Y = sqrtLonsager_fc[1].array(mfi);,
                     const Array4<      Real>& sqrtLOnsager_Z = sqrtLonsager_fc[2].array(mfi););

        AMREX_D_TERM(const Box& box_x = mfi.nodaltilebox(0);,
                     const Box& box_y = mfi.nodaltilebox(1);,
                     const Box& box_z = mfi.nodaltilebox(2););

        amrex::ParallelFor(box_x, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real, MAX_SPECIES> RhoN;
            GpuArray<Real, MAX_SPECIES> RhoAv;
            GpuArray<Real, MAX_SPECIES> RhoNXShift;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> sqrtLOnsager_XN;

            for (int n=0; n<nspecies; ++n ){
                RhoN[n] = rho(i,j,k,n);
                RhoNXShift[n] = rho(i-1,j,k,n);
                for (int m=0; m<nspecies; ++m){
                    sqrtLOnsager_XN(m+1,n+1) = sqrtLOnsager_X(i,j,k,n*nspecies+m);
                }
            }

            ComputeNonnegativeRhoAv(RhoNXShift, RhoN, dx, molmass, RhoAv);

            //update RhoAv for SqrtLOnsager
            Real RhoAvSum = 0.0;
            for (int n=0; n<nspecies; ++n ){
                RhoAvSum += RhoAv[n];
            }

            ComputeSqrtLOnsagerLocal(molmass, RhoAv, RhoAvSum, sqrtLOnsager_XN);

            //copy data back
            for (int n=0; n<nspecies; ++n ){
                for (int m=0; m<nspecies; ++m){
                    sqrtLOnsager_X(i,j,k,n*nspecies+m) = sqrtLOnsager_XN(m+1,n+1);
                }
            }
        }); //end box_x

        amrex::ParallelFor(box_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            GpuArray<Real, MAX_SPECIES> RhoN;
            GpuArray<Real, MAX_SPECIES> RhoAv;
            GpuArray<Real, MAX_SPECIES> RhoNYShift;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> sqrtLOnsager_YN;

            for (int n=0; n<nspecies; ++n ){
                RhoN[n] = rho(i,j,k,n);
                RhoNYShift[n] = rho(i,j-1,k,n);

                for (int m=0; m<nspecies; ++m){
                    sqrtLOnsager_YN(m+1,n+1) = sqrtLOnsager_Y(i,j,k,n*nspecies+m);
                }
            }

            ComputeNonnegativeRhoAv(RhoNYShift, RhoN, dx, molmass, RhoAv);

            Real RhoAvSum = 0.0;
            for (int n=0; n<nspecies; ++n ){
                RhoAvSum += RhoAv[n];
            }

            ComputeSqrtLOnsagerLocal(molmass, RhoAv, RhoAvSum, sqrtLOnsager_YN);

            //copy data back
            for (int n=0; n<nspecies; ++n ){
                for (int m=0; m<nspecies; ++m){
                    sqrtLOnsager_Y(i,j,k,n*nspecies+m) = sqrtLOnsager_YN(m+1,n+1);
                }
            }
        }); //end box_y
#if (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(box_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            GpuArray<Real, MAX_SPECIES> RhoN;
            GpuArray<Real, MAX_SPECIES> RhoAv;
            GpuArray<Real, MAX_SPECIES> RhoNZShift;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> sqrtLOnsager_ZN;

            for (int n=0; n<nspecies; ++n ){
                RhoN[n] = rho(i,j,k,n);
                RhoNZShift[n] = rho(i,j,k-1,n);

                for (int m=0; m<nspecies; ++m){
                    sqrtLOnsager_ZN(m+1,n+1) = sqrtLOnsager_Z(i,j,k,n*nspecies+m);
                }
            }

            ComputeNonnegativeRhoAv(RhoNZShift, RhoN, dx, molmass, RhoAv);

            //update RhoAv for SqrtLOnsager
            Real RhoAvSum = 0.0;
            for (int n=0; n<nspecies; ++n ){
                RhoAvSum += RhoAv[n];
            }

            ComputeSqrtLOnsagerLocal(molmass, RhoAv, RhoAvSum, sqrtLOnsager_ZN);

            //copy data back
            for (int n=0; n<nspecies; ++n ){
                for (int m=0; m<nspecies; ++m){
                    sqrtLOnsager_Z(i,j,k,n*nspecies+m) = sqrtLOnsager_ZN(m+1,n+1);
                }
            }

        }); //end box_z
#endif
    }

}