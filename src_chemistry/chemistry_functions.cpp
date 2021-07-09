#include "chemistry_functions.H"
#include <iostream>
#include "AMReX_ParmParse.H"

void InitializeChemistryNamespace()
{
    // extract inputs parameters
    ParmParse pp;

    // get number of reactions
    pp.get("nreaction",nreaction);
    // get rate constants from input files and assign them to rate_const
    std::vector<amrex::Real> k_tmp(MAX_REACTION);               // temporary vector to extract data from input files
    pp.getarr("rate_const",k_tmp,0,nreaction);                  // get reaction constants from input files
    for (int m=0;m<nreaction;m++) rate_const[m] = k_tmp[m];     // assign them to reac_k

    // get stoich coeffs for reactants
    for (int m=0;m<nreaction;m++)
    {
        // keyword to extract data from input files
        char keyword[128];
        sprintf(keyword,"stoich_%dR",m+1);
        // temporary vector to extract data from input files
        std::vector<int> s_tmp(MAX_SPECIES);
        // get stoich coeffs from input files
        pp.getarr(keyword,s_tmp,0,nspecies);
        // assign them to stoich_coeffs_R
        for (int n=0;n<nspecies;n++) stoich_coeffs_R[m][n] = s_tmp[n];
    }

    // get stoich coeffs for products
    for (int m=0;m<nreaction;m++)
    {
        // keyword to extract data from input files
        char keyword[128];
        sprintf(keyword,"stoich_%dP",m+1);
        // temporary vector to extract data from input files
        std::vector<int> s_tmp(MAX_SPECIES);
        // get stoich coeffs from input files
        pp.getarr(keyword,s_tmp,0,nspecies);
        // assign them to stoich_coeffs_P
        for (int n=0;n<nspecies;n++) stoich_coeffs_P[m][n] = s_tmp[n];
    }

    // stoich coeffs to update number density
    for (int m=0;m<nreaction;m++)
    {
        for (int n=0;n<nspecies;n++)
        {
            stoich_coeffs_PR[m][n] = stoich_coeffs_P[m][n]-stoich_coeffs_R[m][n];
        }
    }


    // get reaction type: Deterministic, CLE or SSA
    pp.get("reaction_type",reaction_type);
    
    return;
}

void compute_reaction_rates(amrex::Real n_dens[MAX_SPECIES], amrex::Real a_r[MAX_REACTION])
{
    for (int m=0; m<nreaction; m++)
    {
        a_r[m] = rate_const[m];
        for (int n=0; n<nspecies; n++) a_r[m] *= pow(n_dens[n],stoich_coeffs_R[m][n]);
    }

    return;
}

void advance_reaction_SSA_cell(amrex::Real n_old[MAX_SPECIES],amrex::Real n_new[MAX_SPECIES],amrex::Real dt,amrex::Real dV,RandomEngine const& engine)
{
    amrex::Real t_local = 0.;

    for (int n=0; n<nspecies; n++) n_new[n] = n_old[n];

    while(true)
    {
        amrex::Real avg_react_rate[MAX_REACTION];
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
        int which_reaction;
        amrex::Real rSum = 0.;
        for (int m=0; m<nreaction; m++)
        {
            rSum = rSum + avg_react_rate[m];
            which_reaction = m;
            if (rSum >= u2) break;
        }

        // update number densities for the reaction that has occured
        for (int n=0; n<nspecies; n++)
        {
            n_new[n] += stoich_coeffs_PR[which_reaction][n]/dV;
        }
    }

    return;
}

void advance_reaction_det_cell(amrex::Real n_old[MAX_SPECIES],amrex::Real n_new[MAX_SPECIES],amrex::Real dt)
{
    for (int n=0; n<nspecies; n++) n_new[n] = n_old[n];

    amrex::Real avg_react_rate[MAX_REACTION];
    compute_reaction_rates(n_new,avg_react_rate);
    
    for (int m=0; m<nreaction; m++)
    {
        avg_react_rate[m] = std::max(0.,avg_react_rate[m]);

        for (int n=0;n<nspecies;n++)
        {
            n_new[n] += dt*stoich_coeffs_PR[m][n]*avg_react_rate[m];
        }
    }

    return;
}

void advance_reaction_CLE_cell(amrex::Real n_old[MAX_SPECIES],amrex::Real n_new[MAX_SPECIES],amrex::Real dt,amrex::Real dV,RandomEngine const& engine)
{
    for (int n=0; n<nspecies; n++) n_new[n] = n_old[n];

    amrex::Real avg_react_rate[MAX_REACTION];
    compute_reaction_rates(n_new,avg_react_rate);
    
    for (int m=0; m<nreaction; m++)
    {
        avg_react_rate[m] = std::max(0.,avg_react_rate[m]);

        amrex::Real W = sqrt(dt/dV)*RandomNormal(0.,1.,engine);

        for (int n=0;n<nspecies;n++)
        {
            n_new[n] += dt*stoich_coeffs_PR[m][n]*avg_react_rate[m];
            n_new[n] += stoich_coeffs_PR[m][n]*sqrt(avg_react_rate[m])*W;
        }
    }

    return;
}
