#include "MFsurfchem_functions.H"
#include "AMReX_ParmParse.H"

AMREX_GPU_MANAGED int MFsurfchem::n_ads_spec;
AMREX_GPU_MANAGED int MFsurfchem::ads_wall_dir;
AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_SPECIES> MFsurfchem::surfcov0;

AMREX_GPU_MANAGED amrex::Real MFsurfchem::surf_site_num_dens;

AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_SPECIES> MFsurfchem::ads_rate_const;
AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_SPECIES> MFsurfchem::des_rate;

AMREX_GPU_MANAGED int MFsurfchem::stoch_surfcov0;
AMREX_GPU_MANAGED int MFsurfchem::stoch_MFsurfchem;

AMREX_GPU_MANAGED amrex::Real MFsurfchem::k_beta;
AMREX_GPU_MANAGED amrex::Real MFsurfchem::e_beta;

AMREX_GPU_MANAGED int MFsurfchem::splitting_MFsurfchem;
AMREX_GPU_MANAGED int MFsurfchem::conversion_MFsurfchem;

AMREX_GPU_MANAGED int MFsurfchem::mean_MFsurfchem;
AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_SPECIES> MFsurfchem::mean_pressure;

void InitializeMFSurfchemNamespace()
{
    // extract inputs parameters
    ParmParse pp;

    n_ads_spec = 0;
    // get the number of species that undergoes adsoprtion/desorption
    pp.query("n_ads_spec",n_ads_spec);

    // if n_ads_spec is set to 0 or not defined in the inputs file, quit the routine
    if (n_ads_spec==0) return;

    ads_wall_dir = 2;
    pp.query("ads_wall_dir",ads_wall_dir);

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

    k_beta = -0.5; // default value
    pp.query("k_beta",k_beta);

    e_beta = 0.5; // default value
    pp.query("e_beta",e_beta);

    // get splitting type: first order (0), strang (1)
    splitting_MFsurfchem = 0; // default value
    pp.query("splitting_MFsurfchem",splitting_MFsurfchem);

    // ads/des of species 1 -> ads of species 1 + desorption of species (1+n)
    conversion_MFsurfchem = 0; // default value
    pp.query("conversion_MFsurfchem",conversion_MFsurfchem);
    if ( (n_ads_spec + conversion_MFsurfchem) > nspecies) {
        Abort("ERROR: desorption species is not included in nspecies");
    }

    // Use equilibrium (mean) values of pressure and temperature to calculate adsorption rate
    mean_MFsurfchem = 0; // default value
    pp.query("mean_MFsurfchem",mean_MFsurfchem);
    if (mean_MFsurfchem > 0) {
        std::vector<amrex::Real> mean_pressure_tmp(MAX_SPECIES);
        pp.queryarr("mean_pressure",mean_pressure_tmp,0,n_ads_spec); // mean partial pressure of adsorption species
        for (int m=0;m<n_ads_spec;m++) {
            mean_pressure[m] = mean_pressure_tmp[m];
        }
    }
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
            if ( (ads_wall_dir == 0 && i == 0) || (ads_wall_dir == 1 && j == 0) || (ads_wall_dir == 2 && k == 0) ) {
                if (stoch_surfcov0==1) {
                    amrex::Real Ntot = rint(surf_site_num_dens*dx[(ads_wall_dir+1)%3]*dx[(ads_wall_dir+2)%3]);  // total number of reactive sites
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

void sample_MFsurfchem(MultiFab& cu, MultiFab& prim, MultiFab& surfcov, MultiFab& dNadsdes, MultiFab& dNads, MultiFab& dNdes,
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
        const Array4<Real> & dNads_arr = dNads.array(mfi);
        const Array4<Real> & dNdes_arr = dNdes.array(mfi);

        amrex::Real Ntot = surf_site_num_dens*dx[(ads_wall_dir+1)%3]*dx[(ads_wall_dir+2)%3];  // total number of reactive sites

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
        {
            if ( (ads_wall_dir == 0 && i == 0) || (ads_wall_dir == 1 && j == 0) || (ads_wall_dir == 2 && k == 0) ) {
                amrex::Real sumtheta = 0.;
                for (int m=0;m<n_ads_spec;m++) {
                    sumtheta += surfcov_arr(i,j,k,m);
                }

                amrex:: Real tempratio = prim_arr(i,j,k,4)/T_init[0];

                for (int m=0;m<n_ads_spec;m++) {
                    amrex::Real pres = prim_arr(i,j,k,5);   // total pressure
                    pres *= prim_arr(i,j,k,6+nspecies+m);   // partial pressure

                    amrex::Real theta = surfcov_arr(i,j,k,m);

                    amrex::Real meanNads;
                    amrex::Real meanNdes;

                    if (mean_MFsurfchem==0) {
                        meanNads = ads_rate_const[m]*pres*(1-sumtheta)*Ntot*dt*pow(tempratio,k_beta);
                        meanNdes = des_rate[m]*theta*Ntot*dt;
                    }
                    else {
                        meanNads = ads_rate_const[m]*mean_pressure[m]*(1-sumtheta)*Ntot*dt; // tempratio = 1
                        meanNdes = des_rate[m]*theta*Ntot*dt;
                    }

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

                    if (conversion_MFsurfchem > 0) {
                        dNads_arr(i,j,k,m) = Nads;
                        dNdes_arr(i,j,k,m) = Ndes;
                    }
                    else {
                        dNadsdes_arr(i,j,k,m) = Nads-Ndes;
                    }
                }
            }
        });
    }

    return;
}

void update_MFsurfchem(MultiFab& cu, MultiFab& prim, MultiFab& surfcov, MultiFab& dNadsdes, MultiFab& dNads, MultiFab& dNdes,
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
        const Array4<Real> & dNads_arr = dNads.array(mfi);
        const Array4<Real> & dNdes_arr = dNdes.array(mfi);

        amrex::Real Ntot = surf_site_num_dens*dx[(ads_wall_dir+1)%3]*dx[(ads_wall_dir+2)%3];  // total number of reactive sites
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( (ads_wall_dir == 0 && i == 0) || (ads_wall_dir == 1 && j == 0) || (ads_wall_dir == 2 && k == 0) ) {

                amrex::Real T_inst = prim_arr(i,j,k,4);

                for (int m=0;m<n_ads_spec;m++) {
                    if (conversion_MFsurfchem > 0) {
                        amrex::Real mconv = m + conversion_MFsurfchem;
                        amrex::Real dNads = dNads_arr(i,j,k,m);
                        amrex::Real dNdes = dNdes_arr(i,j,k,m);
                        amrex::Real factor1 = molmass[m]/AVONUM/(dx[0]*dx[1]*dx[2]);
                        amrex::Real factor1conv = molmass[mconv]/AVONUM/(dx[0]*dx[1]*dx[2]);
                        amrex::Real factor2 = (e_beta*k_B*T_inst+(e0[m]+hcv[m]*T_inst)*molmass[m]/AVONUM)/(dx[0]*dx[1]*dx[2]);
                        amrex::Real factor2conv = (e_beta*k_B*T_inst+(e0[mconv]+hcv[mconv]*T_inst)*molmass[mconv]/AVONUM)/(dx[0]*dx[1]*dx[2]);
                        surfcov_arr(i,j,k,m) += (dNads-dNdes)/Ntot;
                        cu_arr(i,j,k,0) -= factor1*dNads - factor1conv*dNdes;
                        cu_arr(i,j,k,5+m) -= factor1*dNads;
                        cu_arr(i,j,k,5+mconv) += factor1conv*dNdes;
                        cu_arr(i,j,k,4) -= factor2*dNads - factor2conv*dNdes;
                    }
                    else {
                        amrex::Real dN = dNadsdes_arr(i,j,k,m);
                        amrex::Real factor1 = molmass[m]/AVONUM/(dx[0]*dx[1]*dx[2]);
                        amrex::Real factor2 = (e_beta*k_B*T_inst+(e0[m]+hcv[m]*T_inst)*molmass[m]/AVONUM)/(dx[0]*dx[1]*dx[2]);
                        surfcov_arr(i,j,k,m) += dN/Ntot;
                        cu_arr(i,j,k,0) -= factor1*dN;
                        cu_arr(i,j,k,5+m) -= factor1*dN;
                        cu_arr(i,j,k,4) -= factor2*dN;
                    }
                }
            }
        });
    }

    return;
}
