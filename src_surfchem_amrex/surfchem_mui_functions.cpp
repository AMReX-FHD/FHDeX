#include "surfchem_mui_functions.H"
#include "common_namespace.H"

#include <AMReX_iMultiFab.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int surfchem_mui::nspec_mui;

void InitializeSurfChemMUINamespace()
{
    // extract inputs parameters
    ParmParse pp;

    // number of species involved in mui (via adsorption/desorption)
    pp.get("nspec_mui",surfchem_mui::nspec_mui);

    return;
}

void amrex_fetch_Ntot(MultiFab& Ntot, MPMD::Copier const& copier)
{
    iMultiFab one(Ntot.boxArray(), Ntot.DistributionMap(), 1, 0);
    copier.recv(one, 0, 1);

    for (MFIter mfi(Ntot); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        if (bx.smallEnd(2) == 0) {
            const Box& b2d = amrex::makeSlab(bx, 2, 0);
            auto const& Ntot_arr = Ntot.array(mfi);
            auto const& one_arr = one.const_array(mfi);
            amrex::ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Ntot_arr(i,j,k) = one_arr(i,j,k);
            });
        }
    }
}

void amrex_fetch_surfcov(MultiFab const& Ntot, MultiFab& surfcov,
                         MPMD::Copier const& copier)
{
    if (surfchem_mui::nspec_mui == 1) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(surfchem_mui::nspec_mui == 1,
                                         "Assuming nspec_mui is 1 for now");
        iMultiFab occ1(Ntot.boxArray(), Ntot.DistributionMap(), 1, 0);
        copier.recv(occ1, 0, 1);

        for (MFIter mfi(surfcov); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (bx.smallEnd(2) == 0) {
                const Box& b2d = amrex::makeSlab(bx, 2, 0);
                auto const& Ntot_arr = Ntot.const_array(mfi);
                auto const& occ1_arr = occ1.const_array(mfi);
                auto const& surfcov_arr = surfcov.array(mfi);
                amrex::ParallelFor(b2d, surfchem_mui::nspec_mui,
                                   [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    surfcov_arr(i,j,k,n) = occ1_arr(i,j,k,n) / Ntot_arr(i,j,k);
                });
            }
        }
    }
    else if (surfchem_mui::nspec_mui == 2) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(surfchem_mui::nspec_mui == 2,
                                         "Assuming nspec_mui is 2 for now");
        iMultiFab occ(Ntot.boxArray(), Ntot.DistributionMap(), surfchem_mui::nspec_mui, 0);
        copier.recv(occ, 0, 1);
        copier.recv(occ, 1, 1);

        for (MFIter mfi(surfcov); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (bx.smallEnd(2) == 0) {
                const Box& b2d = amrex::makeSlab(bx, 2, 0);
                auto const& Ntot_arr = Ntot.const_array(mfi);
                auto const& occ_arr = occ.const_array(mfi);
                auto const& surfcov_arr = surfcov.array(mfi);
                amrex::ParallelFor(b2d, surfchem_mui::nspec_mui,
                                   [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    surfcov_arr(i,j,k,n) = occ_arr(i,j,k,n) / Ntot_arr(i,j,k);
                });
            }
        }
    }
    else if (surfchem_mui::nspec_mui == 3) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(surfchem_mui::nspec_mui == 3,
                                         "Assuming nspec_mui is 3 for now");
        iMultiFab occ(Ntot.boxArray(), Ntot.DistributionMap(), surfchem_mui::nspec_mui, 0);
        copier.recv(occ, 0, 1);
        copier.recv(occ, 1, 1);
        copier.recv(occ, 2, 1);

        for (MFIter mfi(surfcov); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (bx.smallEnd(2) == 0) {
                const Box& b2d = amrex::makeSlab(bx, 2, 0);
                auto const& Ntot_arr = Ntot.const_array(mfi);
                auto const& occ_arr = occ.const_array(mfi);
                auto const& surfcov_arr = surfcov.array(mfi);
                amrex::ParallelFor(b2d, surfchem_mui::nspec_mui,
                                   [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    surfcov_arr(i,j,k,n) = occ_arr(i,j,k,n) / Ntot_arr(i,j,k);
                });
            }
        }
    }
}

void amrex_push(MultiFab const& cu, MultiFab const& prim, MPMD::Copier const& copier)
// this routine pushes the following information to KMC
// - species number densities and temperature of FHD cells contacting the interface
// it assumes that the interface is perpendicular to the z-axis
// and includes cells with the smallest value of z (i.e. k=0)
{
    if (surfchem_mui::nspec_mui == 1) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(surfchem_mui::nspec_mui == 1,
                                         "Assuming nspec_mui is 1 for now");

        constexpr Real AVONUM = 6.02214076e23;

        MultiFab dens(cu.boxArray(), cu.DistributionMap(), surfchem_mui::nspec_mui, 0);
        for (MFIter mfi(dens); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (bx.smallEnd(2) == 0) {
                const Box& b2d = amrex::makeSlab(bx, 2, 0);
                auto const& dens_arr = dens.array(mfi);
                auto const& cu_arr = cu.const_array(mfi);
                amrex::ParallelFor(b2d, surfchem_mui::nspec_mui,
                                   [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    dens_arr(i,j,k,n) = cu_arr(i,j,k,5+n) * AVONUM / common::molmass[n];
                });
            }
        }
        copier.send(dens, 0, 1);
        copier.send(prim, 4, 1);
    }
    else if (surfchem_mui::nspec_mui == 2) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(surfchem_mui::nspec_mui == 2,
                                         "Assuming nspec_mui is 2 for now");

        constexpr Real AVONUM = 6.02214076e23;
        MultiFab dens(cu.boxArray(), cu.DistributionMap(), surfchem_mui::nspec_mui, 0);
        for (MFIter mfi(dens); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (bx.smallEnd(2) == 0) {
                const Box& b2d = amrex::makeSlab(bx, 2, 0);
                auto const& dens_arr = dens.array(mfi);
                auto const& cu_arr = cu.const_array(mfi);
                amrex::ParallelFor(b2d, surfchem_mui::nspec_mui,
                                   [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    dens_arr(i,j,k,n) = cu_arr(i,j,k,5+n) * AVONUM / common::molmass[n];
                });
            }
        }
        copier.send(dens, 0, 1);
        copier.send(dens, 1, 1);
        copier.send(prim, 4, 1);
    }
    else if (surfchem_mui::nspec_mui == 3) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(surfchem_mui::nspec_mui == 3,
                                         "Assuming nspec_mui is 3 for now");

        constexpr Real AVONUM = 6.02214076e23;
        MultiFab dens(cu.boxArray(), cu.DistributionMap(), surfchem_mui::nspec_mui, 0);
        for (MFIter mfi(dens); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (bx.smallEnd(2) == 0) {
                const Box& b2d = amrex::makeSlab(bx, 2, 0);
                auto const& dens_arr = dens.array(mfi);
                auto const& cu_arr = cu.const_array(mfi);
                amrex::ParallelFor(b2d, surfchem_mui::nspec_mui,
                                   [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    dens_arr(i,j,k,n) = cu_arr(i,j,k,5+n) * AVONUM / common::molmass[n];
                });
            }
        }
        copier.send(dens, 0, 1);
        copier.send(dens, 1, 1);
        copier.send(dens, 2, 1);
        copier.send(prim, 4, 1);
    }
}

void amrex_fetch(MultiFab& cu, MultiFab const& prim, GpuArray<Real,3> const& dx,
                 MPMD::Copier const& copier)
// this routine fetches the following information from KMC:
// - adsoprtion and desoprtion counts of each species between time points
// it assumes that the interface is perpendicular to the z-axis
// and includes cells with the smallest value of z (i.e. k=0)
{
    if (surfchem_mui::nspec_mui == 1) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(surfchem_mui::nspec_mui == 1,
                                         "Assuming nspec_mui is 1 for now");

        constexpr Real AVONUM = 6.02214076e23;
        constexpr Real BETA = 0.5;

        iMultiFab acdc(cu.boxArray(), cu.DistributionMap(), 2, 0);
        copier.recv(acdc,0,1);
        copier.recv(acdc,1,1);

        for (MFIter mfi(cu); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            if (bx.smallEnd(2) == 0){
                const Box& b2d = amrex::makeSlab(bx, 2, 0);
                auto const& acdc_arr = acdc.const_array(mfi);
                auto const& prim_arr = prim.const_array(mfi);
                auto const& cu_arr = cu.array(mfi);
                amrex::ParallelFor(b2d, surfchem_mui::nspec_mui,
                                   [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {

                    int ac = acdc_arr(i,j,k,0);
                    int dc = acdc_arr(i,j,k,1);
                    amrex::Real dN = static_cast<Real>(ac-dc);
                    amrex::Real T_inst = prim_arr(i,j,k,4);
                    amrex::Real factor1 = common::molmass[n]/AVONUM/(dx[0]*dx[1]*dx[2]);
                    amrex::Real factor2 = (BETA*common::k_B*T_inst+(common::e0[n]+common::hcv[n]*T_inst)*common::molmass[n]/AVONUM)/(dx[0]*dx[1]*dx[2]);

                    cu_arr(i,j,k,0) -= factor1*dN;
                    cu_arr(i,j,k,5+n) -= factor1*dN;
                    cu_arr(i,j,k,4) -= factor2*dN;
                });
            }
        }
    }
    else if (surfchem_mui::nspec_mui == 2) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(surfchem_mui::nspec_mui == 2,
                                         "Assuming nspec_mui is 2 for now");

        constexpr Real AVONUM = 6.02214076e23;
        constexpr Real BETA = 0.5;

        iMultiFab acdc(cu.boxArray(), cu.DistributionMap(), 4, 0);
        copier.recv(acdc,0,1);
        copier.recv(acdc,1,1);
        copier.recv(acdc,2,1);
        copier.recv(acdc,3,1);

        for (MFIter mfi(cu); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            if (bx.smallEnd(2) == 0) {
                const Box& b2d = amrex::makeSlab(bx, 2, 0);
                auto const& acdc_arr = acdc.const_array(mfi);
                auto const& prim_arr = prim.const_array(mfi);
                auto const& cu_arr = cu.array(mfi);
                amrex::ParallelFor(b2d, surfchem_mui::nspec_mui,
                                   [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {

                    int ac1 = acdc_arr(i,j,k,0);
                    int dc1 = acdc_arr(i,j,k,1);
                    int ac2 = acdc_arr(i,j,k,2);
                    int dc2 = acdc_arr(i,j,k,3);
                    amrex::Real dN1 = static_cast<Real>(ac1-dc1);
                    amrex::Real dN2 = static_cast<Real>(ac2-dc2);
                    double dN[3] = {dN1, dN2, 0};
                    amrex::Real T_inst = prim_arr(i,j,k,4);
                    amrex::Real factor1 = common::molmass[n]/AVONUM/(dx[0]*dx[1]*dx[2]);
                    amrex::Real factor2 = (BETA*common::k_B*T_inst+(common::e0[n]+common::hcv[n]*T_inst)*common::molmass[n]/AVONUM)/(dx[0]*dx[1]*dx[2]);

                    cu_arr(i,j,k,0) -= factor1*dN[n];
                    cu_arr(i,j,k,5+n) -= factor1*dN[n];
                    cu_arr(i,j,k,4) -= factor2*dN[n];
                });
            }
        }
    }
    else if (surfchem_mui::nspec_mui == 3) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(surfchem_mui::nspec_mui == 3,
                                         "Assuming nspec_mui is 3 for now");

        constexpr Real AVONUM = 6.02214076e23;
        constexpr Real BETA = 0.5;

        iMultiFab acdc(cu.boxArray(), cu.DistributionMap(), 5, 0);
        copier.recv(acdc,0,1);
        copier.recv(acdc,1,1);
        copier.recv(acdc,2,1);
        copier.recv(acdc,3,1);
        copier.recv(acdc,4,1);

        for (MFIter mfi(cu); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            if (bx.smallEnd(2) == 0) {
                const Box& b2d = amrex::makeSlab(bx, 2, 0);
                auto const& acdc_arr = acdc.const_array(mfi);
                auto const& prim_arr = prim.const_array(mfi);
                auto const& cu_arr = cu.array(mfi);
                amrex::ParallelFor(b2d, surfchem_mui::nspec_mui,
                                   [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {

                    int ac1 = acdc_arr(i,j,k,0);
                    int dc1 = acdc_arr(i,j,k,1);
                    int ac2 = acdc_arr(i,j,k,2);
                    int dc2 = acdc_arr(i,j,k,3);
                    int adc3 = acdc_arr(i,j,k,4);
                    amrex::Real dN1 = static_cast<Real>(ac1-dc1);
                    amrex::Real dN2 = static_cast<Real>(ac2-dc2);
                    amrex::Real dN3 = static_cast<Real>(-1*adc3);
                    double dN[4] = {dN1, dN2, dN3, 0};
                    amrex::Real T_inst = prim_arr(i,j,k,4);
                    amrex::Real factor1 = common::molmass[n]/AVONUM/(dx[0]*dx[1]*dx[2]);
                    amrex::Real factor2 = (BETA*common::k_B*T_inst+(common::e0[n]+common::hcv[n]*T_inst)*common::molmass[n]/AVONUM)/(dx[0]*dx[1]*dx[2]);

                    cu_arr(i,j,k,0) -= factor1*dN[n];
                    cu_arr(i,j,k,5+n) -= factor1*dN[n];
                    cu_arr(i,j,k,4) -= factor2*dN[n];
                });
            }
        }
    }
}