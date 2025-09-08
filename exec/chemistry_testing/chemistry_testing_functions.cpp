#include "common_functions.H"
#include "chemistry_functions.H"
#include "rng_functions.H"
#include "chemistry_testing_functions.H"

void EMstep_chem_only(MultiFab& rho_old, MultiFab& rho_new,
                       const amrex::Geometry geom, const amrex::Real dt)
{
    if (reaction_type!=1) amrex::Abort("EMstep_chem_only assumes reaction_type=1");

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    BoxArray ba = rho_old.boxArray();
    DistributionMapping dm = rho_old.DistributionMap();

    MultiFab source(ba,dm,nspecies,0);;

    MultiFab ranchem(ba,dm,nreaction,0);

    // initialize white noise field
    for (int m=0;m<nreaction;m++) {
        MultiFabFillRandom(ranchem,m,1.,geom);
    }

    compute_chemistry_source_CLE_2(dt,dx[0]*dx[1]*dx[2],rho_old,0,source,0,ranchem);

    for ( MFIter mfi(rho_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & rho_old_fab = rho_old.array(mfi);
        const Array4<Real> & rho_new_fab = rho_new.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            rho_new_fab(i,j,k,n) = rho_old_fab(i,j,k,n) + dt*source_fab(i,j,k,n);
        });
    }

}


void RK3step_chem_only(MultiFab& rho_old, MultiFab& rho_new,
                       const amrex::Geometry geom, const amrex::Real dt)
{
    if (reaction_type!=1) amrex::Abort("RK3step_chem_only assumes reaction_type=1");

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    BoxArray ba = rho_old.boxArray();
    DistributionMapping dm = rho_old.DistributionMap();

    MultiFab rhop(ba,dm,nspecies,0);
    MultiFab rhop2(ba,dm,nspecies,0);
    MultiFab rhop3(ba,dm,nspecies,0);;

    MultiFab source(ba,dm,nspecies,0);;

    MultiFab ranchem(ba,dm,nreaction,0);
    MultiFab ranchem_A(ba,dm,nreaction,0);
    MultiFab ranchem_B(ba,dm,nreaction,0);

    // weights for stochastic fluxes; swgt2 changes each stage
    amrex::Real swgt1, swgt2;
    swgt1 = 1.;

    // initialize white noise fields
    for (int m=0;m<nreaction;m++) {
        MultiFabFillRandom(ranchem_A,m,1.,geom);
        MultiFabFillRandom(ranchem_B,m,1.,geom);
    }

    // stage1
    swgt2 = (2.*std::sqrt(2.)+std::sqrt(3.))/5.;

    MultiFab::LinComb(ranchem,
                swgt1, ranchem_A, 0,
                swgt2, ranchem_B, 0,
                0, nreaction, 0);

    compute_chemistry_source_CLE_2(dt,dx[0]*dx[1]*dx[2],rho_old,0,source,0,ranchem);

    for ( MFIter mfi(rho_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & rho_old_fab = rho_old.array(mfi);
        const Array4<Real> & rhop_fab = rhop.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            rhop_fab(i,j,k,n) = rho_old_fab(i,j,k,n) + dt*source_fab(i,j,k,n);
        });
    }

    // stage2
    swgt2 = (-4.*std::sqrt(2.)+3.*std::sqrt(3.))/5.;

    MultiFab::LinComb(ranchem,
                swgt1, ranchem_A, 0,
                swgt2, ranchem_B, 0,
                0, nreaction, 0);

    compute_chemistry_source_CLE_2(dt,dx[0]*dx[1]*dx[2],rhop,0,source,0,ranchem);

    for ( MFIter mfi(rho_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & rho_old_fab = rho_old.array(mfi);
        const Array4<Real> & rhop_fab = rhop.array(mfi);
        const Array4<Real> & rhop2_fab = rhop2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            rhop2_fab(i,j,k,n) = 0.25*( 3.*rho_old_fab(i,j,k,n) + rhop_fab(i,j,k,n) + dt*source_fab(i,j,k,n) );
        });
    }

    // stage3
    swgt2 = (std::sqrt(2.)-2.*std::sqrt(3.))/10.;

    MultiFab::LinComb(ranchem,
                swgt1, ranchem_A, 0,
                swgt2, ranchem_B, 0,
                0, nreaction, 0);

    compute_chemistry_source_CLE_2(dt,dx[0]*dx[1]*dx[2],rhop2,0,source,0,ranchem);

    for ( MFIter mfi(rho_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & rho_old_fab = rho_old.array(mfi);
        const Array4<Real> & rho_new_fab = rho_new.array(mfi);
        const Array4<Real> & rhop2_fab = rhop2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            rho_new_fab(i,j,k,n) = (2./3.)*( 0.5*rho_old_fab(i,j,k,n) + rhop2_fab(i,j,k,n) + dt*source_fab(i,j,k,n) );
        });
    }
}



void compute_chemistry_source_CLE_1(amrex::Real dt, amrex::Real dV,
                              MultiFab& mf_in, int startComp_in,
                              MultiFab& source, int startComp_out)
// mf_in: input MultiFab containing mass densitities rho1, rho2, ..., rho_nspecies
// startComp_in: position of rho1 in mf_in
// source: output MultiFab containing source terms corresponding to rho1, rho2, ..., rho_nspecies
// startComp_out: position of the first source term corresponding to rho1 in MultiFab source
{
    if (reaction_type<0 || reaction_type>2) amrex::Abort("ERROR: invalid reaction_type");

    if (reaction_type==2) amrex::Abort("ERROR: reaction_type=2 not implemented yet");

    GpuArray<amrex::Real,MAX_SPECIES> m_s;
    for (int n=0; n<nspecies; n++) m_s[n] = molmass[n]/(avogadro);

    for (MFIter mfi(mf_in); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& rho_arr = mf_in.array(mfi);
        const Array4<Real>& source_arr = source.array(mfi);

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, RandomEngine const& engine) noexcept
        {
            GpuArray<amrex::Real,MAX_SPECIES> n_dens;
            for (int n=0; n<nspecies; n++) n_dens[n] = rho_arr(i,j,k,n+startComp_in)/m_s[n];

            GpuArray<amrex::Real,MAX_REACTION> avg_react_rate;
            compute_reaction_rates(n_dens,avg_react_rate);

            GpuArray<amrex::Real,MAX_SPECIES> sourceArr;
            for (int n=0; n<nspecies; n++) sourceArr[n] = 0.;

            for (int m=0; m<nreaction; m++)
            {
                avg_react_rate[m] = std::max(0.,avg_react_rate[m]);
                for (int n=0; n<nspecies; n++)
                    sourceArr[n] += m_s[n]*stoich_coeffs_PR(m,n)*avg_react_rate[m];
            }

            if (reaction_type==1)
            {
                for (int m=0; m<nreaction; m++)
                {
                    amrex::Real W = RandomNormal(0.,1.,engine)/sqrt(dt*dV);
                    for (int n=0; n<nspecies; n++)
                        sourceArr[n] += m_s[n]*stoich_coeffs_PR(m,n)*sqrt(avg_react_rate[m])*W;
                }
            }

            for (int n=0; n<nspecies; n++) source_arr(i,j,k,n+startComp_out) = sourceArr[n];
        });
    }
}

void compute_chemistry_source_CLE_2(amrex::Real dt, amrex::Real dV,
                              MultiFab& mf_in, int startComp_in,
                              MultiFab& source, int startComp_out, MultiFab& ranchem)
// mf_in: input MultiFab containing mass densitities rho1, rho2, ..., rho_nspecies
// startComp_in: position of rho1 in mf_in
// source: output MultiFab containing source terms corresponding to rho1, rho2, ..., rho_nspecies
// startComp_out: position of the first source term corresponding to rho1 in MultiFab source
{
    if (reaction_type!=1) amrex::Abort("ERROR: compute_chemistry_source_CLE assumes reaction_type=1");

    GpuArray<amrex::Real,MAX_SPECIES> m_s;
    for (int n=0; n<nspecies; n++) m_s[n] = molmass[n]/(avogadro);

    for (MFIter mfi(mf_in); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& rho_arr = mf_in.array(mfi);
        const Array4<Real>& source_arr = source.array(mfi);
        const Array4<Real>& ranchem_arr = ranchem.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            GpuArray<amrex::Real,MAX_SPECIES> n_dens;
            for (int n=0; n<nspecies; n++) n_dens[n] = rho_arr(i,j,k,n+startComp_in)/m_s[n];

            GpuArray<amrex::Real,MAX_REACTION> avg_react_rate;
            compute_reaction_rates(n_dens,avg_react_rate);

            GpuArray<amrex::Real,MAX_SPECIES> sourceArr;
            for (int n=0; n<nspecies; n++) sourceArr[n] = 0.;

            for (int m=0; m<nreaction; m++)
            {
                avg_react_rate[m] = std::max(0.,avg_react_rate[m]);

                amrex::Real W = ranchem_arr(i,j,k,m)/sqrt(dt*dV);

                for (int n=0; n<nspecies; n++)
                {
                    sourceArr[n] += m_s[n]*stoich_coeffs_PR(m,n)*avg_react_rate[m];
                    sourceArr[n] += m_s[n]*stoich_coeffs_PR(m,n)*sqrt(avg_react_rate[m])*W;
                }
            }

            for (int n=0; n<nspecies; n++) source_arr(i,j,k,n+startComp_out) = sourceArr[n];
        });
    }
}

AMREX_GPU_HOST_DEVICE void compute_reaction_rates(GpuArray<Real,MAX_SPECIES>& n_dens,
                                                  GpuArray<Real,MAX_REACTION>& a_r)
{
    for (int m=0; m<nreaction; m++)
    {
        a_r[m] = rate_const[m];

        for (int n=0; n<nspecies; n++)
            a_r[m] *= pow(n_dens[n],stoich_coeffs_R(m,n));
    }

    return;
}

AMREX_GPU_HOST_DEVICE void advance_reaction_det_cell(GpuArray<amrex::Real,MAX_SPECIES>& n_old,
                                                     GpuArray<amrex::Real,MAX_SPECIES>& n_new,amrex::Real dt)
{
    for (int n=0; n<nspecies; n++) n_new[n] = n_old[n];

    GpuArray<amrex::Real,MAX_REACTION> avg_react_rate;
    compute_reaction_rates(n_new,avg_react_rate);

    for (int m=0; m<nreaction; m++)
    {
        avg_react_rate[m] = std::max(0.,avg_react_rate[m]);
        for (int n=0; n<nspecies; n++)
            n_new[n] += dt*stoich_coeffs_PR(m,n)*avg_react_rate[m];
    }

    return;
}

AMREX_GPU_HOST_DEVICE void advance_reaction_CLE_cell(GpuArray<amrex::Real,MAX_SPECIES>& n_old,
                                                     GpuArray<amrex::Real,MAX_SPECIES>& n_new,
                                                     amrex::Real dt, amrex::Real dV,
                                                     RandomEngine const& engine)
{
    for (int n=0; n<nspecies; n++) n_new[n] = n_old[n];

    GpuArray<amrex::Real,MAX_REACTION> avg_react_rate;
    compute_reaction_rates(n_new,avg_react_rate);

    for (int m=0; m<nreaction; m++)
    {
        avg_react_rate[m] = std::max(0.,avg_react_rate[m]);

        amrex::Real W = sqrt(dt/dV)*RandomNormal(0.,1.,engine);

        for (int n=0; n<nspecies; n++)
        {
            n_new[n] += dt*stoich_coeffs_PR(m,n)*avg_react_rate[m];
            n_new[n] += stoich_coeffs_PR(m,n)*sqrt(avg_react_rate[m])*W;
        }
    }

    return;
}

AMREX_GPU_HOST_DEVICE void advance_reaction_SSA_cell(GpuArray<amrex::Real,MAX_SPECIES>& n_old,
                                                     GpuArray<amrex::Real,MAX_SPECIES>& n_new,
                                                     amrex::Real dt, amrex::Real dV,
                                                     RandomEngine const& engine)
{
    amrex::Real t_local = 0.;

    for (int n=0; n<nspecies; n++) n_new[n] = n_old[n];

    while(true)
    {
        GpuArray<amrex::Real,MAX_REACTION> avg_react_rate;
        compute_reaction_rates(n_new,avg_react_rate);

        amrex::Real rTotal = 0.;
        for (int m=0; m<nreaction; m++)
        {
            // convert reation rates to propensities
            avg_react_rate[m] = std::max(0.,avg_react_rate[m]*dV);
            rTotal += avg_react_rate[m];
        }

        if (rTotal==0.) break;

        amrex::Real u1 = amrex::Random(engine);
        amrex::Real tau = -log(1-u1)/rTotal;
        t_local += tau; // update t_local

        if (t_local > dt) break;

        amrex::Real u2 = amrex::Random(engine);
        u2 *= rTotal;

        // find which reaction has occured
        int which_reaction=0;
        amrex::Real rSum = 0.;
        for (int m=0; m<nreaction; m++)
        {
            rSum = rSum + avg_react_rate[m];
            which_reaction = m;
            if (rSum >= u2) break;
        }

        // update number densities for the reaction that has occured
        for (int n=0; n<nspecies; n++)
            n_new[n] += stoich_coeffs_PR(which_reaction,n)/dV;
    }

    return;
}
