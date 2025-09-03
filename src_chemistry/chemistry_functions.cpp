#include "chemistry_functions.H"
#include "AMReX_ParmParse.H"

AMREX_GPU_MANAGED int chemistry::nreaction;

// from the fortran code, stoich_coeffs_R = stoichiometric_factors(spec,1,reac)
// from the fortran code, stoich_coeffs_P = stoichiometric_factors(spec,2,reac)
// stoich_coeffs_PR = stoich_coeffs_P - stoich_coeffs_R
AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> chemistry::stoich_coeffs_R;
AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> chemistry::stoich_coeffs_P;
AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> chemistry::stoich_coeffs_PR;

// reaction rate constant for each reaction (assuming Law of Mass Action holds)
// using rate_multiplier, reaction rates can be changed by the same factor
// if include_discrete_LMA_correction, n^2 and n^3 in rate expressions become
// n*(n-1/dv) and n*(n-1/dv)*(n-2/dv).
AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> chemistry::rate_const;
AMREX_GPU_MANAGED amrex::Real chemistry::rate_multiplier;
AMREX_GPU_MANAGED int chemistry::include_discrete_LMA_correction;

// if n is positive, exclude species n (=solvent) when computing reaction rates
// in this case, the concentration of the solvent is assumed to be constant,
// which should be reflected on rate constants.
// if 0, no species is excluded
// e.g. U + S -> 2U, if exclude_solvent_comput_rates=0, rate=k*n_U*n_S
//                   if exclude_solvent_comput_rates=2, rate=k_new*n_U where k_new=k*n_S
AMREX_GPU_MANAGED int chemistry::exclude_solvent_comput_rates;

// from the fortran code this was use_Poisson_rng (0=CLE; 1=tau leaping; -1=deterministic; 2=SSA)
// here it's being used as reaction_type (0=deterministic; 1=CLE; 2=SSA; 3=tau leap)
AMREX_GPU_MANAGED int chemistry::reaction_type;

// use mole fraction based LMA
AMREX_GPU_MANAGED int chemistry::use_mole_frac_LMA;

// specific to compressible codes
AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> chemistry::alpha_param;
AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> chemistry::beta_param;
AMREX_GPU_MANAGED amrex::Real chemistry::T0_chem;


void InitializeChemistryNamespace()
{
    // extract inputs parameters
    ParmParse pp;

    nreaction = 0;
    // get number of reactions
    pp.query("nreaction",nreaction);
    if (nreaction > MAX_REACTION) {
        Abort("nreaction > MAX_REACTION; recompile with a new MAX_REAC in the GNUmakefile");
    }

    // if nreaction is set to zero or not defined in the inputs file, quit the routine
    if (nreaction==0) return;

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

    // get rate constants
    std::vector<amrex::Real> k_tmp(MAX_REACTION);
    pp.getarr("rate_const",k_tmp,0,nreaction);
    for (int m=0; m<nreaction; m++) rate_const[m] = k_tmp[m];

    rate_multiplier = 1.;
    pp.query("rate_multiplier",rate_multiplier);

    include_discrete_LMA_correction = 0;
    pp.query("include_discrete_LMA_correction",include_discrete_LMA_correction);

    exclude_solvent_comput_rates = -1;
    pp.query("exclude_solvent_comput_rates",exclude_solvent_comput_rates);

    // get reaction type (0=deterministic; 1=CLE; 2=SSA; 3=tau leap)
    pp.get("reaction_type",reaction_type);

    use_mole_frac_LMA = 0;
    pp.query("use_mole_frac_LMA",use_mole_frac_LMA);

    // get alpha parameter for compressible code
    std::vector<amrex::Real> alpha_tmp(MAX_REACTION);
    pp.queryarr("alpha_param",alpha_tmp,0,nreaction);
    for (int m=0; m<nreaction; m++) alpha_param[m] = alpha_tmp[m];

    // get beta parameter for compressible code
    std::vector<amrex::Real> beta_tmp(MAX_REACTION);
    pp.queryarr("beta_param",beta_tmp,0,nreaction);
    for (int m=0; m<nreaction; m++) beta_param[m] = beta_tmp[m];

    T0_chem = 0.;
    // get temperature T0 for rate constants for compressible code
    pp.query("T0_chem",T0_chem);

    return;
}

// used in compressible code only
void compute_compressible_chemistry_source_CLE(amrex::Real dt, amrex::Real dV,
                                               MultiFab& prim, MultiFab& source, MultiFab& ranchem)
{
    if (reaction_type!=0 && reaction_type!=1) {
        amrex::Abort("ERROR: compute_compressible_chemistry_source_CLE only works for reaction_type = 0 or 1");
    }

    if (T0_chem<=0.) amrex::Abort("ERROR: T0_chem>0 expected");

    for (MFIter mfi(prim); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& prim_arr = prim.array(mfi);
        const Array4<Real>& source_arr = source.array(mfi);
        const Array4<Real>& ranchem_arr = ranchem.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            amrex::Real pres = prim_arr(i,j,k,5);
            amrex::Real T = prim_arr(i,j,k,4);
            amrex::Real T0 = T0_chem;
            amrex::Real ctot = pres/Runiv/T;
            amrex::Real Navo = avogadro;

            GpuArray<amrex::Real,MAX_SPECIES> ck;   // molar concentrations
            for (int n=0; n<nspecies; n++) ck[n] = prim_arr(i,j,k,6+nspecies+n)*ctot;

            GpuArray<amrex::Real,MAX_REACTION> avg_react_rate;
            for (int m=0; m<nreaction; m++)
            {
                // rate constants
                avg_react_rate[m] = rate_const[m];
                avg_react_rate[m] *= exp(-alpha_param[m]/Runiv*(1/T-1/T0));
                avg_react_rate[m] *= pow(T/T0,beta_param[m]);

                for (int n=0; n<nspecies; n++)
                {
                    // rate in terms of molar concentrations
                    avg_react_rate[m] *= pow(ck[n],stoich_coeffs_R(m,n));
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
                    // mean
                    sourceArr[n] += molmass[n]*stoich_coeffs_PR(m,n)*avg_react_rate[m];

                    // fluctuation
                    if (reaction_type==1) sourceArr[n] += molmass[n]*stoich_coeffs_PR(m,n)*sqrt(avg_react_rate[m]/Navo)*W;
                }
            }

            for (int n=0; n<nspecies; n++) source_arr(i,j,k,5+n) = sourceArr[n];
        });
    }
}


void ChemicalRates(const MultiFab& n_cc, MultiFab& chem_rate, const amrex::Geometry& geom, const amrex::Real& dt,
                   const MultiFab& n_interm, Vector<Real>& lin_comb_coef_in, Real volume_factor_in)
{
    if (nreaction == 1) {
        chem_rate.setVal(0.);
        return;
    }

    int lin_comb_avg_react_rate = 1;
    if (lin_comb_coef_in[0] == 1. && lin_comb_coef_in[1] == 0.) {
        lin_comb_avg_react_rate = 0;
    }

    GpuArray<Real,2> lin_comb_coef;
    lin_comb_coef[0] = lin_comb_coef_in[0];
    lin_comb_coef[1] = lin_comb_coef_in[1];

    const Real* dx = geom.CellSize();

    Real dv = (AMREX_SPACEDIM == 3) ? dx[0]*dx[1]*dx[2]*cell_depth : dx[0]*dx[1]*cell_depth;
    dv *= volume_factor_in;

    for (MFIter mfi(n_cc); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<const Real>& n_arr = n_cc.array(mfi);
        const Array4<const Real>& n_int = n_interm.array(mfi);

        const Array4<Real>& rate = chem_rate.array(mfi);

        if (reaction_type == 2) { // SSA

            amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, amrex::RandomEngine const& engine) noexcept
            {
                GpuArray<Real,MAX_SPECIES> n_old;
                GpuArray<Real,MAX_SPECIES> n_new;
                GpuArray<Real,MAX_REACTION> avg_reaction_rate;

                Real t_local = 0.;

                for (int n=0; n<nspecies; ++n) {
                    n_old[n] = n_arr(i,j,k,n);
                    n_new[n] = n_arr(i,j,k,n);
                }

                while(true)
                {
                    compute_reaction_rates(n_new,avg_reaction_rate,dv);

                    Real rTotal = 0.;
                    for (int m=0; m<nreaction; m++)
                    {
                        // convert reation rates to propensities
                        avg_reaction_rate[m] = std::max(0.,avg_reaction_rate[m]*dv);
                        rTotal += avg_reaction_rate[m];
                    }

                    if (rTotal==0.) break;

                    Real u1 = amrex::Random(engine);
                    Real tau = -log(1-u1)/rTotal;
                    t_local += tau; // update t_local

                    if (t_local > dt) break;

                    Real u2 = amrex::Random(engine);
                    u2 *= rTotal;

                    // find which reaction has occured
                    int which_reaction=0;
                    Real rSum = 0.;
                    for (int m=0; m<nreaction; m++)
                    {
                        rSum = rSum + avg_reaction_rate[m];
                        which_reaction = m;
                        if (rSum >= u2) break;
                    }

                    // update number densities for the reaction that has occured
                    for (int n=0; n<nspecies; n++) {
                        n_new[n] += stoich_coeffs_PR(which_reaction,n)/dv;
                    }
                }

                for (int n=0; n<nspecies; ++n) {
                    rate(i,j,k,n) = (n_new[n] - n_old[n] ) / dt;
                }
            });

        } else {
            amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, amrex::RandomEngine const& engine) noexcept
            {
                GpuArray<Real,MAX_SPECIES> n_in;
                GpuArray<Real,MAX_SPECIES> n_int_in;
                GpuArray<Real,MAX_REACTION> avg_reaction_rate;
                GpuArray<Real,MAX_REACTION> avg_reaction_rate_interm;
                GpuArray<Real,MAX_REACTION> avg_num_reactions;
                GpuArray<Real,MAX_REACTION> num_reactions;

                for (int n=0; n<nspecies; ++n) {
                    rate(i,j,k,n) = 0.;
                    n_in[n]     = n_arr(i,j,k,n);
                    n_int_in[n] = n_int(i,j,k,n);
                }

                if (lin_comb_avg_react_rate == 1) {
                    compute_reaction_rates(n_in    , avg_reaction_rate       , dv);
                    compute_reaction_rates(n_int_in, avg_reaction_rate_interm, dv);
                    for (int r=0; r<nreaction; ++r) {
                        avg_reaction_rate[r] = lin_comb_coef[0]*avg_reaction_rate[r] + lin_comb_coef[1]*avg_reaction_rate_interm[r];
                    }
                } else {
                    compute_reaction_rates(n_in, avg_reaction_rate, dv);
                }

                for (int r=0; r<nreaction; ++r) {
                    avg_num_reactions[r] = std::max(0.,avg_reaction_rate[r]*dv*dt);
                }
                sample_num_reactions(n_in,num_reactions,avg_num_reactions,engine);
                for (int r=0; r<nreaction; ++r) {
                    for (int n=0; n<nspecies; ++n) {
                        rate(i,j,k,n) += num_reactions[r]/dv/dt * stoich_coeffs_PR(r,n);
                    }
                }
            });
        }
    }
}

AMREX_GPU_HOST_DEVICE void compute_reaction_rates(GpuArray<Real,MAX_SPECIES>& n_in,
                                                  GpuArray<Real,MAX_REACTION>& reaction_rates,
                                                  const amrex::Real& dv)
{
    GpuArray<Real,MAX_SPECIES> n_nonneg;

    Real n_sum = 0.;
    for (int n=0; n<nspecies; ++n) {
        n_nonneg[n] = std::max(0.,n_in[n]);
        n_sum += n_nonneg[n];
    }
    if (n_sum < 0.) {
        n_sum = 1./dv;
        Abort("compute_reaction_rates() - n_sum < 0, is this right?");
    }

    if (use_mole_frac_LMA && include_discrete_LMA_correction) {

        Abort("compute_reaction_rates() - use_mole_frac_LMA && include_discrete_LMA_correction not supported yet");

/*
      ! Use mole-fraction based LMA (general ideal mixtures) with integer corrections

      do reaction=1, nreactions
        reaction_rates(reaction) = rate_multiplier*rate_const(reaction)
        do species=1, nspecies
          ! Donev: Replaced case statement by if here
          ! Donev: Made sure n_sum is never zero for empty cells to avoid division by zero

          if(stoichiometric_factors(species,1,reaction)>=1) then
            ! rate ~ N/N_sum
            if(n_nonneg(species)>0.0d0) then ! This species is present in this cell
               reaction_rates(reaction) = reaction_rates(reaction) * n_nonneg(species)/n_sum
            else
               reaction_rates(reaction) = 0.0d0
            end if
          end if
          if(stoichiometric_factors(species,1,reaction)>=2) then
            ! rate ~ (N/N_sum)*((N-1)/(N_sum-1))
            ! Donev: Avoid division by zero or negative rates
            if(n_nonneg(species)>1.0d0/dv) then ! There is at least one molecule of this species in this cell
               reaction_rates(reaction) = reaction_rates(reaction) * (n_nonneg(species)-1.0d0/dv)/(n_sum-1.0d0/dv)
            else
               reaction_rates(reaction) = 0.0d0
            end if
          end if
          if(stoichiometric_factors(species,1,reaction)>=3) then ! Donev added ternary reactions here
            ! rate ~ (N/N_sum)*((N-1)/(N_sum-1))*((N-2)/(N_sum-2))
            if(n_nonneg(species)>2.0d0/dv) then ! There is at least two molecules of this species in this cell
              reaction_rates(reaction) = reaction_rates(reaction) * (n_nonneg(species)-2.0d0/dv)/(n_sum-2.0d0/dv)
            else
               reaction_rates(reaction) = 0.0d0
            end if
          end if
          if(stoichiometric_factors(species,1,reaction)>=4) then
            ! This is essentially impossible in practice and won't happen
            call bl_error("Stochiometric coefficients larger then 3 not supported")
          end if
        end do
      end do
*/

    } else if (include_discrete_LMA_correction == 0 && exclude_solvent_comput_rates == -1) {

        if (use_mole_frac_LMA) {
            for (int n=0; n<nspecies; ++n) {
                n_nonneg[n] /= n_sum;
            }
        }

        for (int r=0; r<nreaction; ++r) {
            reaction_rates[r] = rate_multiplier*rate_const[r];
            for (int n=0; n<nspecies; ++n) {
                reaction_rates[r] *= std::pow(n_nonneg[n],stoich_coeffs_R(r,n));
            }
        }

    } else { // General case of number-density based LMA is handled by slower code that includes species by species

        for (int r=0; r<nreaction; ++r) {
            reaction_rates[r] = rate_multiplier*rate_const[r];

            for (int n=0; n<nspecies; ++n) {
                if (n == exclude_solvent_comput_rates) {
                    continue;
                }
                if (include_discrete_LMA_correction) {

                    int coef = stoich_coeffs_R(r,n);
                    if (coef == 0) {
                        // Species doe not participate in reaction
                    } else if (coef == 1) {
                        reaction_rates[r] *= n_nonneg[n];
                    } else if (coef == 2) {
                        reaction_rates[r] *= n_nonneg[n]*std::max(0.,n_nonneg[n]-1./dv);
                    } else if (coef == 3) {
                        reaction_rates[r] *= n_nonneg[n]*std::max(0.,n_nonneg[n]-1./dv)*std::max(0.,n_nonneg[n]-2./dv);
                    } else {
                        // This is essentially impossible in practice and won't happen
                        Abort("Stochiometric coefficients larger then 3 not supported");
                    }

                } else {
                    reaction_rates[r] *= std::pow(n_nonneg[n],stoich_coeffs_R(r,n));
                }
            } // end loop over species
        } // end loop over reaction
    }

}

AMREX_GPU_HOST_DEVICE void sample_num_reactions(GpuArray<Real,MAX_SPECIES>& n_in,
                                                GpuArray<Real,MAX_REACTION>& num_reactions,
                                                GpuArray<Real,MAX_REACTION>& avg_num_reactions,
                                                const amrex::RandomEngine& engine)
{
    if (reaction_type == 0) { // deterministic
        for (int n=0; n<nreaction; ++n) {
            num_reactions[n] = avg_num_reactions[n];
        }
    } else if (reaction_type == 1) { // CLE
        for (int n=0; n<nreaction; ++n) {
            Real rand = RandomNormal(0.,1.,engine);
            num_reactions[n] = avg_num_reactions[n] + std::sqrt(avg_num_reactions[n])*rand;
        }
    } else if (reaction_type == 3) { // tau leaping
        for (int n=0; n<nreaction; ++n) {
            num_reactions[n] = RandomPoisson(avg_num_reactions[n], engine);
        }
    } else {
        Abort("sample_num_reactions() - reaction_type not supported");
    }
}

