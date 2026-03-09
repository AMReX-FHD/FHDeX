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
    iMultiFab occ(Ntot.boxArray(), Ntot.DistributionMap(), surfchem_mui::nspec_mui, 0);
    for (int m=0;m<surfchem_mui::nspec_mui;m++) copier.recv(occ, m, 1);

    for (MFIter mfi(surfcov); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        if (bx.smallEnd(2) == 0) {
            const Box& b2d = amrex::makeSlab(bx, 2, 0);
            auto const& Ntot_arr = Ntot.const_array(mfi);
            auto const& occ_arr = occ.const_array(mfi);
            auto const& surfcov_arr = surfcov.array(mfi);
            amrex::ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                for (int m=0;m<surfchem_mui::nspec_mui;m++) surfcov_arr(i,j,k,m) = occ_arr(i,j,k,m) / Ntot_arr(i,j,k);
            });
        }
    }
}

void amrex_push(MultiFab const& cu, MultiFab const& prim, MPMD::Copier const& copier)
// this routine pushes the following information to KMC
// - species number densities and temperature of FHD cells contacting the interface
// it assumes that the interface is perpendicular to the z-axis
// and includes cells with the smallest value of z (i.e. k=0)
{
    MultiFab pres(cu.boxArray(), cu.DistributionMap(), surfchem_mui::nspec_mui, 0);
    for (MFIter mfi(pres); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        if (bx.smallEnd(2) == 0) {
            const Box& b2d = amrex::makeSlab(bx, 2, 0);
            auto const& pres_arr = pres.array(mfi);
            auto const& prim_arr = prim.array(mfi);
            amrex::ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                for (int m=0;m<surfchem_mui::nspec_mui;m++) pres_arr(i,j,k,m) = prim_arr(i,j,k,5)*prim_arr(i,j,k,6+common::nspecies+m);
            });
        }
    }
    for (int m=0;m<surfchem_mui::nspec_mui;m++) copier.send(pres, m, 1);
    copier.send(prim, 4, 1);
}

void amrex_fetch(MultiFab& cu, MultiFab const& prim, GpuArray<Real,3> const& dx,
                 MPMD::Copier const& copier)
// this routine fetches the following information from KMC:
// - adsoprtion and desoprtion counts of each species between time points
// it assumes that the interface is perpendicular to the z-axis
// and includes cells with the smallest value of z (i.e. k=0)
{
    constexpr Real AVONUM = 6.02214076e23;
    constexpr Real BETA = 0.5;

    iMultiFab ac(cu.boxArray(), cu.DistributionMap(), surfchem_mui::nspec_mui, 0);
    iMultiFab dc(cu.boxArray(), cu.DistributionMap(), surfchem_mui::nspec_mui, 0);
    for (int m=0;m<surfchem_mui::nspec_mui;m++) {
        copier.recv(ac, m, 1);
        copier.recv(dc, m, 1);
    }

    for (MFIter mfi(cu); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        if (bx.smallEnd(2) == 0){
            const Box& b2d = amrex::makeSlab(bx, 2, 0);
            auto const& ac_arr = ac.const_array(mfi);
            auto const& dc_arr = dc.const_array(mfi);
            auto const& prim_arr = prim.const_array(mfi);
            auto const& cu_arr = cu.array(mfi);
            amrex::ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                for (int m=0;m<surfchem_mui::nspec_mui;m++) {
                    int ac = ac_arr(i,j,k,m);
                    int dc = dc_arr(i,j,k,m);
                    amrex::Real dN = static_cast<Real>(ac-dc);
                    amrex::Real T_inst = prim_arr(i,j,k,4);
                    amrex::Real factor1 = common::molmass[m]/AVONUM/(dx[0]*dx[1]*dx[2]);
                    amrex::Real factor2 = (BETA*common::k_B*T_inst+(common::e0[m]+common::hcv[m]*T_inst)*common::molmass[m]/AVONUM)/(dx[0]*dx[1]*dx[2]);

                    cu_arr(i,j,k,0) -= factor1*dN;
                    cu_arr(i,j,k,5+m) -= factor1*dN;
                    cu_arr(i,j,k,4) -= factor2*dN;
                }
            });
        }
    }
}
