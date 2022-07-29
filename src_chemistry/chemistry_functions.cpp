#include "chemistry_functions.H"
#include "AMReX_ParmParse.H"

AMREX_GPU_MANAGED int chemistry::nreaction;

AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> chemistry::rate_const;
AMREX_GPU_MANAGED amrex::Real chemistry::T0_chem;

AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> chemistry::stoich_coeffs_R; 
AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> chemistry::stoich_coeffs_P;
AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> chemistry::stoich_coeffs_PR;

AMREX_GPU_MANAGED int chemistry::reaction_type;

void InitializeChemistryNamespace()
{
    // extract inputs parameters
    ParmParse pp;

    nreaction = 0;
    // get number of reactions
    pp.query("nreaction",nreaction);

    // if nreaction is set to zero or not defined in the inputs file, quit the routine
    if (nreaction==0) return;

    // get rate constants
    std::vector<amrex::Real> k_tmp(MAX_REACTION);
    pp.getarr("rate_const",k_tmp,0,nreaction);
    for (int m=0; m<nreaction; m++) rate_const[m] = k_tmp[m];

    T0_chem = 0.;
    // get temperature T0 for rate constants
    pp.query("T0_chem",T0_chem);

    // get stoich coeffs for reactants
    for (int m=0; m<nreaction; m++)
    {
        // keyword to extract data from input files
        char keyword[128];
        sprintf(keyword,"stoich_%dR",m+1);
        // temporary vector to extract data from input files
        std::vector<int> s_tmp(MAX_SPECIES);
        // get stoich coeffs from input files
        pp.getarr(keyword,s_tmp,0,nspecies);
        // assign them to stoich_coeffs_R
        for (int n=0; n<nspecies; n++) stoich_coeffs_R(m,n) = s_tmp[n];
    }

    // get stoich coeffs for products
    for (int m=0; m<nreaction; m++)
    {
        // keyword to extract data from input files
        char keyword[128];
        sprintf(keyword,"stoich_%dP",m+1);
        // temporary vector to extract data from input files
        std::vector<int> s_tmp(MAX_SPECIES);
        // get stoich coeffs from input files
        pp.getarr(keyword,s_tmp,0,nspecies);
        // assign them to stoich_coeffs_P
        for (int n=0; n<nspecies; n++) stoich_coeffs_P(m,n) = s_tmp[n];
    }

    // stoich coeffs to update number density
    for (int m=0; m<nreaction; m++)
        for (int n=0; n<nspecies; n++)
            stoich_coeffs_PR(m,n) = stoich_coeffs_P(m,n)-stoich_coeffs_R(m,n);

    // get reaction type: Deterministic, CLE or SSA
    pp.get("reaction_type",reaction_type);
    
    return;
}

void compute_chemistry_source_CLE(amrex::Real dt, amrex::Real dV,
                                  MultiFab& prim, MultiFab& source, MultiFab& ranchem)
{
    if (reaction_type!=1) amrex::Abort("ERROR: compute_chemistry_source_CLE assumes reaction_type=1");

    if (T0_chem<=0.) amrex::Abort("ERROR: T0_chem>0 expected");

    GpuArray<amrex::Real,MAX_SPECIES> m_s;
    for (int n=0; n<nspecies; n++) m_s[n] = molmass[n]/(Runiv/k_B);

    for (MFIter mfi(prim); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& prim_arr = prim.array(mfi);
        const Array4<Real>& source_arr = source.array(mfi);
        const Array4<Real>& ranchem_arr = ranchem.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            amrex::Real pres = prim_arr(i,j,k,5);
            amrex::Real pres0 = 1.013250e6;

            amrex::Real T = prim_arr(i,j,k,4);
            amrex::Real T0 = T0_chem;

            GpuArray<amrex::Real,MAX_SPECIES> Xk;
            for (int n=0; n<nspecies; n++) Xk[n] = prim_arr(i,j,k,6+nspecies+n);

            GpuArray<amrex::Real,MAX_REACTION> avg_react_rate;
            for (int m=0; m<nreaction; m++)
            {
                avg_react_rate[m] = rate_const[m];

                for (int n=0; n<nspecies; n++)
                {
                    // corrections for fluctuating temperature
                    avg_react_rate[m] *= exp(stoich_coeffs_R(m,n)*m_s[n]*e0[n]/k_B*(1/T-1/T0));
                    avg_react_rate[m] *= pow(T/T0,-stoich_coeffs_R(m,n)*m_s[n]*hcp[n]/k_B);
                    // rate in terms of pressure (more precisely activity)
                    avg_react_rate[m] *= pow(pres/pres0*Xk[n],stoich_coeffs_R(m,n));
                }
            }

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

            for (int n=0; n<nspecies; n++) source_arr(i,j,k,5+n) = sourceArr[n];
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
    for (int n=0; n<nspecies; n++) m_s[n] = molmass[n]/(Runiv/k_B);

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
    for (int n=0; n<nspecies; n++) m_s[n] = molmass[n]/(Runiv/k_B);

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
