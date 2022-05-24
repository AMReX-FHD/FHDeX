#include "MFsurfchem_functions.H"
#include "AMReX_ParmParse.H"

AMREX_GPU_MANAGED int MFsurfchem::ads_spec;
AMREX_GPU_MANAGED amrex::Real MFsurfchem::surfcov0;
AMREX_GPU_MANAGED amrex::Real MFsurfchem::surf_site_num_dens;

AMREX_GPU_MANAGED amrex::Real MFsurfchem::ads_rate_const;
AMREX_GPU_MANAGED amrex::Real MFsurfchem::des_rate;

void InitializeMFSurfchemNamespace()
{
    // extract inputs parameters
    ParmParse pp;

    ads_spec = -1;
    // get the species number that undergoes adsoprtion/desorption 
    pp.query("ads_spec",ads_spec);

    // if ads_spec is set to -1 or not defined in the inputs file, quit the routine
    if (ads_spec==-1) return;

    surfcov0 = 0.;
    // get initial surface coverage
    pp.query("surfcov0",surfcov0);

    surf_site_num_dens = 0.;
    // get number density of adsorption sites on the lattice
    pp.query("surf_site_num_dens",surf_site_num_dens);

    ads_rate_const = 0.;
    // get adsorption rate const 
    pp.query("ads_rate_const",ads_rate_const);

    des_rate = 0.;
    // get desoprtion rate
    pp.query("des_rate",des_rate);

    return;
}

void init_surfcov(MultiFab& surfcov, const amrex::Real* dx)
{
    for (MFIter mfi(surfcov,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        const Array4<Real> & surfcov_arr = surfcov.array(mfi);

        amrex::Real Ntot = surf_site_num_dens*dx[0]*dx[1];  // total number of reactive sites

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k==0) {
                surfcov_arr(i,j,k,0) = RandomPoisson(Ntot*surfcov0)/Ntot;
            } else {
                surfcov_arr(i,j,k,0) = 0.;
            }
        });
    }

    return;
}

void sample_MFsurfchem(MultiFab& cu, MultiFab& prim, MultiFab& surfcov, MultiFab& dNadsdes,
                       const amrex::Real* dx, const amrex::Real dt)
{
    for (MFIter mfi(cu,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        const Array4<Real> & cu_arr = cu.array(mfi);
        const Array4<Real> & prim_arr = prim.array(mfi);
        const Array4<Real> & surfcov_arr = surfcov.array(mfi);
        const Array4<Real> & dNadsdes_arr = dNadsdes.array(mfi);

        amrex::Real Ntot = surf_site_num_dens*dx[0]*dx[1];  // total number of reactive sites

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k==0) {
                amrex::Real dens = cu_arr(i,j,k,5+ads_spec);   // mass density
                dens *= AVONUM/molmass[ads_spec];         // number density

                amrex::Real theta = surfcov_arr(i,j,k,0);

                amrex::Real meanNads = ads_rate_const*dens*(1-theta)*Ntot*dt;
                //amrex::Real Nads = meanNads + sqrt(meanNads)*RandomNormal(0.,1.);
                amrex::Real Nads = RandomPoisson(meanNads);

                amrex::Real meanNdes = des_rate*theta*Ntot*dt;
                //amrex::Real Ndes = meanNdes + sqrt(meanNdes)*RandomNormal(0.,1.);
                amrex::Real Ndes = RandomPoisson(meanNdes);

                dNadsdes_arr(i,j,k,0) = Nads-Ndes;
            }
        });
    }

    return;
}

void update_MFsurfchem(MultiFab& cu, MultiFab& surfcov, MultiFab& dNadsdes,
                       const amrex::Real* dx, const amrex::Real dt)
{
    for (MFIter mfi(cu,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        const Array4<Real> & cu_arr = cu.array(mfi);
        const Array4<Real> & surfcov_arr = surfcov.array(mfi);
        const Array4<Real> & dNadsdes_arr = dNadsdes.array(mfi);

        amrex::Real Ntot = surf_site_num_dens*dx[0]*dx[1];  // total number of reactive sites
        amrex::Real factor = molmass[ads_spec]/AVONUM/(dx[0]*dx[1]*dx[2]);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k==0) {
                amrex::Real dN = dNadsdes_arr(i,j,k,0);

                surfcov_arr(i,j,k,0) += dN/Ntot;
                cu_arr(i,j,k,0) -= factor*dN;
                cu_arr(i,j,k,5+ads_spec) -= factor*dN;
            }
        });
    }

    return;
}
