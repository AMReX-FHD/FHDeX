#include "MFsurfchem_functions.H"
#include "AMReX_ParmParse.H"

AMREX_GPU_MANAGED int MFsurfchem::n_ads_spec;
AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_SPECIES> MFsurfchem::surfcov0;

AMREX_GPU_MANAGED amrex::Real MFsurfchem::surf_site_num_dens;

AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_SPECIES> MFsurfchem::ads_rate_const;
AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_SPECIES> MFsurfchem::des_rate;

AMREX_GPU_MANAGED int MFsurfchem::stoch_surfcov0;
AMREX_GPU_MANAGED int MFsurfchem::stoch_MFsurfchem;

// temperature correction exponent
#define BETA    0.5

void InitializeMFSurfchemNamespace()
{
    // extract inputs parameters
    ParmParse pp;

    n_ads_spec = 0;
    // get the number of species that undergoes adsoprtion/desorption
    pp.query("n_ads_spec",n_ads_spec);

    // if n_ads_spec is set to 0 or not defined in the inputs file, quit the routine
    if (n_ads_spec==0) return;

    // load default values to surfcov0 array
    for (int m=0;m<n_ads_spec;m++) surfcov0[m] = 0.;
    
    std::vector<amrex::Real> surfcov0_tmp(MAX_SPECIES);
    // get initial surface coverage
    if (pp.queryarr("surfcov0",surfcov0_tmp,0,n_ads_spec)){
        for (int m=0;m<n_ads_spec;m++) surfcov0[m] = surfcov0_tmp[m];
    }
    
    surf_site_num_dens = 0.;
    // get number density of adsorption sites on the lattice
    pp.query("surf_site_num_dens",surf_site_num_dens);

    // load default values to ads_rate_const array
    for (int m=0;m<n_ads_spec;m++) ads_rate_const[m] = 0.;

    std::vector<amrex::Real> ads_rate_const_tmp(MAX_SPECIES);
    // get adsorption rate const 
    if (pp.queryarr("ads_rate_const",ads_rate_const_tmp,0,n_ads_spec)){
        for (int m=0;m<n_ads_spec;m++) ads_rate_const[m] = ads_rate_const_tmp[m];
    }
    
    // load default values to des_rate array
    for (int m=0;m<n_ads_spec;m++) des_rate[m] = 0.;

    std::vector<amrex::Real> des_rate_tmp(MAX_SPECIES);
    // get desorption rate
    if (pp.queryarr("des_rate",des_rate_tmp,0,n_ads_spec)){
        for (int m=0;m<n_ads_spec;m++) des_rate[m] = des_rate_tmp[m];
    }

    stoch_surfcov0 = 1; // default value
    pp.query("stoch_surfcov0",stoch_surfcov0);

    stoch_MFsurfchem = 1; // default value
    pp.query("stoch_MFsurfchem",stoch_MFsurfchem);
    return;
}

void init_surfcov(MultiFab& surfcov, const amrex::Geometry& geom)
{
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    GpuArray<Real,MAX_SPECIES> sum_surfcov0;
    if (stoch_surfcov0==1) {
        sum_surfcov0[0] = surfcov0[0];
        for (int m=1;m<n_ads_spec;m++)
            sum_surfcov0[m] = sum_surfcov0[m-1] + surfcov0[m];
    }

    for (MFIter mfi(surfcov,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<Real> & surfcov_arr = surfcov.array(mfi);

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
        {
            if (k==0) {
                if (stoch_surfcov0==1) {
                    amrex::Real Ntot = rint(surf_site_num_dens*dx[0]*dx[1]);  // total number of reactive sites
                    GpuArray<int,MAX_SPECIES> Nocc;
               
                    for (int m=0;m<n_ads_spec;m++) Nocc[m] = 0;

                    for (int n=0;n<Ntot;n++) {
                        amrex::Real u = amrex::Random(engine);
                        for (int m=0;m<n_ads_spec;m++) {
                            if (u<sum_surfcov0[m]) {
                                Nocc[m]++;
                                break;
                            }
                        }
                    }

                    for (int m=0;m<n_ads_spec;m++) surfcov_arr(i,j,k,m) = Nocc[m]/Ntot;
                } else {
                    for (int m=0;m<n_ads_spec;m++) surfcov_arr(i,j,k,m) = surfcov0[m];
                }
            } else {
                for (int m=0;m<n_ads_spec;m++) surfcov_arr(i,j,k,m) = 0.;
            }
        });
    }

    return;
}

void sample_MFsurfchem(MultiFab& cu, MultiFab& prim, MultiFab& surfcov, MultiFab& dNadsdes,
                       const amrex::Geometry& geom, const amrex::Real dt)
{
    
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    for (MFIter mfi(cu,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<Real> & cu_arr = cu.array(mfi);
        const Array4<Real> & prim_arr = prim.array(mfi);
        const Array4<Real> & surfcov_arr = surfcov.array(mfi);
        const Array4<Real> & dNadsdes_arr = dNadsdes.array(mfi);

        amrex::Real Ntot = surf_site_num_dens*dx[0]*dx[1];  // total number of reactive sites

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
        {
            if (k==0) {
                amrex::Real sumtheta = 0.;
                for (int m=0;m<n_ads_spec;m++) {
                    sumtheta += surfcov_arr(i,j,k,m);
                }

                amrex:: Real tempratio = prim_arr(i,j,k,4)/T_init[0];

                for (int m=0;m<n_ads_spec;m++) {
                    amrex::Real dens = cu_arr(i,j,k,5+m);   // mass density
                    dens *= AVONUM/molmass[m];              // number density

                    amrex::Real theta = surfcov_arr(i,j,k,m);

                    amrex::Real meanNads = ads_rate_const[m]*dens*(1-sumtheta)*Ntot*dt*pow(tempratio,BETA);
                    amrex::Real meanNdes = des_rate[m]*theta*Ntot*dt;

                    amrex::Real Nads;
                    amrex::Real Ndes;

                    if (stoch_MFsurfchem==0) {
                        Nads = meanNads;
                        Ndes = meanNdes;
                    }
                    else {
                        Nads = RandomPoisson(meanNads,engine);
                        Ndes = RandomPoisson(meanNdes,engine);
                    }

                    dNadsdes_arr(i,j,k,m) = Nads-Ndes;
                }
            }
        });
    }

    return;
}

void update_MFsurfchem(MultiFab& cu, MultiFab& prim, MultiFab& surfcov, MultiFab& dNadsdes,
                       const amrex::Geometry& geom)
{
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for (MFIter mfi(cu,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<Real> & cu_arr = cu.array(mfi);
        const Array4<Real> & prim_arr = prim.array(mfi);
        const Array4<Real> & surfcov_arr = surfcov.array(mfi);
        const Array4<Real> & dNadsdes_arr = dNadsdes.array(mfi);

        amrex::Real Ntot = surf_site_num_dens*dx[0]*dx[1];  // total number of reactive sites
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k==0) {

                amrex::Real T_inst = prim_arr(i,j,k,4);

                for (int m=0;m<n_ads_spec;m++) {
                    amrex::Real dN = dNadsdes_arr(i,j,k,m);
                    amrex::Real factor1 = molmass[m]/AVONUM/(dx[0]*dx[1]*dx[2]);
                    amrex::Real factor2 = (BETA*k_B*T_inst+(e0[m]+hcv[m]*T_inst)*molmass[m]/AVONUM)/(dx[0]*dx[1]*dx[2]);

                    surfcov_arr(i,j,k,m) += dN/Ntot;
                    cu_arr(i,j,k,0) -= factor1*dN;
                    cu_arr(i,j,k,5+m) -= factor1*dN;
                    cu_arr(i,j,k,4) -= factor2*dN;
                }
            }
        });
    }

    return;
}
