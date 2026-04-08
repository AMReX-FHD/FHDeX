#include "multispec_functions.H"

#include "AMReX_ParmParse.H"

int                                                         multispec::inverse_type;
int                                                         multispec::temp_type;
int                                                         multispec::chi_iterations;
amrex::Real                                                 multispec::start_time;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_ELEMENT> multispec::Dbar;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> multispec::Dtherm;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_ELEMENT> multispec::H_offdiag;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> multispec::H_diag;
amrex::Real                                                 multispec::fraction_tolerance;
int                                                         multispec::correct_flux;
int                                                         multispec::print_error_norms;
AMREX_GPU_MANAGED int                                       multispec::is_nonisothermal;
AMREX_GPU_MANAGED int                                       multispec::is_ideal_mixture;
int                                                         multispec::use_lapack;
AMREX_GPU_MANAGED int                                       multispec::use_multiphase;
AMREX_GPU_MANAGED int                                       multispec::use_flory_huggins;
AMREX_GPU_MANAGED amrex::Real                               multispec::kc_tension;
AMREX_GPU_MANAGED amrex::Real                               multispec::alpha_gex;
AMREX_GPU_MANAGED amrex::Array2D<Real,0,MAX_SPECIES-1,0,MAX_SPECIES-1> multispec::fh_kappa;
AMREX_GPU_MANAGED amrex::Array2D<Real,0,MAX_SPECIES-1,0,MAX_SPECIES-1> multispec::fh_chi;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> multispec::fh_monomers;
AMREX_GPU_MANAGED amrex::Real                               multispec::monomer_mass;
AMREX_GPU_MANAGED int                                       multispec::n_gex;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> multispec::c_init_1;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> multispec::c_init_2;

int                                                         multispec::midpoint_stoch_mass_flux_type;
AMREX_GPU_MANAGED int                                       multispec::avg_type;
int                                                         multispec::mixture_type;

// charged fluid
int                                                         multispec::use_charged_fluid;
int                                                         multispec::print_debye_len;
amrex::Real                                                 multispec::dielectric_const;
int                                                         multispec::dielectric_type;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> multispec::charge_per_mass;
amrex::Real                                                 multispec::theta_pot;
int                                                         multispec::num_pot_iters;
amrex::Real                                                 multispec::dpdt_factor;
amrex::Real                                                 multispec::relxn_param_charge;
int                                                         multispec::E_ext_type;
amrex::Vector<amrex::Real>                                  multispec::E_ext_value;
int                                                         multispec::electroneutral;
int                                                         multispec::induced_charge_eo;
int                                                         multispec::zero_eps_on_wall_type;
int                                                         multispec::zero_charge_on_wall_type;
amrex::Real                                                 multispec::zero_eps_on_wall_left_end;
amrex::Real                                                 multispec::zero_eps_on_wall_right_start;
int                                                         multispec::bc_function_type;
amrex::Real                                                 multispec::L_pos;
amrex::Real                                                 multispec::L_trans;
amrex::Real                                                 multispec::L_zero;

void InitializeMultispecNamespace() {

    if (MAX_ELEMENT != MAX_SPECIES*(MAX_SPECIES-1)/2) {
        Print() << "MAX_ELEMENT " << MAX_ELEMENT << std::endl;
        Print() << "MAX_SPECIES " << MAX_SPECIES << std::endl;
        Abort("MAX_ELEMENT != (MAX_SPECIES*MAX_SPECIES-1)/2)");
    }

    E_ext_value.resize(AMREX_SPACEDIM);

    // specify default values first, then read in values from inputs file

    // Physical properties:
    //----------------------
    fraction_tolerance = 1.e-14; // For roundoff errors in mass and mole fractions
                                 // must be larger than machine eps of else the W=(1,0) case fails)
    start_time = 0.;
    inverse_type = 1;       // Only for LAPACK:  1=inverse, 2=pseudo inverse
    correct_flux = 1;       // Manually ensure mass is conserved to roundoff
    print_error_norms = 1;
    is_ideal_mixture = 1;   // If T assume Gamma=I (H=0) and simplify
    is_nonisothermal = 0;   // If T Soret effect will be included
    use_lapack = 0;         // Use LAPACK or iterative method for diffusion matrix (recommend False)
    use_multiphase = 0;     // for RTIL
    use_flory_huggins = 0;   // for flory huggins
    kc_tension = 0;         // for RTIL
    alpha_gex = 0;          // for RTIL
    n_gex = 1;              // for RTIL
    chi_iterations = 10;    // number of iterations used in Dbar2chi_iterative

    // Initial and boundary conditions
    //----------------------

    temp_type = 0;  // for initializing temperature
    for (int i=0; i<MAX_SPECIES; ++i) {
        c_init_1[i] = 1.;   // initial values for c
        c_init_2[i] = 1.;
    }

    // Thermodynamic and transport properties:
    //----------------------

    // These are lower-triangules of symmetric matrices represented as vectors
    // Number of elements is (nspecies*(nspecies-1)/2)
    // The values are read row by row starting from top going down (this allows easy addition/deletion of new species/rows)
    // So D_12; D_13, D_23; D_14, D_24, D_34; ...
    for (int i=0; i<MAX_SPECIES; ++i) {
        Dtherm[i] = 0.; // thermo-diffusion coefficients, only differences among elements matter
        H_diag[i] = 0.; // Diagonal of H=d^2F/dx^2, these are vectors of length nspecies
    }
    for (int i=0; i<MAX_ELEMENT; ++i) {
        Dbar[i] = 1.;      // Maxwell-Stefan diffusion constant
        H_offdiag[i] = 0.;
    }

    // Algorithm control
    //----------------------
    midpoint_stoch_mass_flux_type = 1; // 1 = Strato
                                       // 2 = Ito

    avg_type = 1;  // how to compute stochastc_mass_fluxdiv
                   // 1=arithmetic (with C0-Heaviside), 2=geometric, 3=harmonic
                   // 10=arithmetic average with discontinuous Heaviside function
                   // 11=arithmetic average with C1-smoothed Heaviside function
                   // 12=arithmetic average with C2-smoothed Heaviside function

    mixture_type = 0; // Model for how transport and thermodynamic coefficients depend on composition
                      // See compute_mixture_properties.f90 for values supported at present
                      // The default mixture_type=0 means no dependence on composition

    // for charged fluid
    use_charged_fluid = 0;
    print_debye_len = 0;
    dielectric_const = 1.;
    dielectric_type = 0;   // 0 = assumes constant epsilon
                           // 1 = (1+c1)*dielectric_const
                           // see fluid_charge.f90:compute_permittivity()
    for (int i=0; i<MAX_SPECIES; ++i) {
        charge_per_mass[i] = 0.;
    }
    bc_function_type = 0;  // 0 = constant
                           // 1 = cubic, see description below

    L_pos = 0.;            // length of part of boundary where there is positive charge flux, if cubic function is imposed
    L_trans = 0.;          // length of transition part of boundary, where the value varies like a cubic
    L_zero = 0.;           // length of part of boundary where there is zero charge flux, if cubic function is imposed

    theta_pot = 0.5;       // for implicit algorithm_type=3, controls
                           // temporal discretization for potential term
    num_pot_iters = 2;
    dpdt_factor = 0.;
    E_ext_type = 0;         // if 1, sets an external E field to E_ext_value
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        E_ext_value[i] = 0.;        // spacedim-vector specifying external E field
    }
    electroneutral = 0;              // use electroneutral diffusion fluxes
    induced_charge_eo = 0;           // are we simulating ICEO?
    relxn_param_charge = 1.;         // Used to prevent slow buildup of charge for electroneutral, keep at 1.0
    zero_eps_on_wall_type = 0;       // set eps=0 on certain Dirichlet walls
    // if we want homogeneous Neumann bc's on
    // phi for part of a Dirichlet wall
    zero_charge_on_wall_type = 0;    // set sigma=0 on certain Neumann walls
    // if we want homogeneous Neumann bc's on
    // phi for part of a Dirichlet wall
    zero_eps_on_wall_left_end = 0.25;    // eg if set to 0.25, then eps will be set to 0 on the wall from 0*Lx --> 0.25*Lx
    zero_eps_on_wall_right_start = 0.75; // eg if set to 0.75, then eps will be set to 0 on the wall from 0.75*Lx --> 1.*Lx

    ParmParse pp;

    int temp_max = std::max(3,MAX_SPECIES);
    temp_max = std::max(temp_max,MAX_ELEMENT);

    amrex::Vector<amrex::Real> temp(temp_max,0.);

    // pp.query searches for optional parameters
    // pp.get aborts if the parameter is not found
    // pp.getarr and queryarr("string",inputs,start_indx,count); can be used for arrays

    pp.query("fraction_tolerance",fraction_tolerance);
    pp.query("start_time",start_time);
    pp.query("inverse_type",inverse_type);
    pp.query("correct_flux",correct_flux);
    pp.query("print_error_norms",print_error_norms);
    pp.query("is_ideal_mixture",is_ideal_mixture);
    pp.query("is_nonisothermal",is_nonisothermal);
    pp.query("use_lapack",use_lapack);
    pp.query("use_multiphase",use_multiphase);
    pp.query("use_flory_huggins",use_flory_huggins);
    pp.query("monomer_mass",monomer_mass);
    pp.query("kc_tension",kc_tension);
    pp.query("alpha_gex",alpha_gex);
    pp.query("n_gex",n_gex);
    pp.query("chi_iterations",chi_iterations);
    pp.query("temp_type",temp_type);
    if(pp.queryarr("fh_kappa",temp)) {
        for (int i=0; i<nspecies; ++i) {
        for (int j=0; j<nspecies; ++j) {
            fh_kappa(i,j) = temp[i*nspecies+j];
        }
        }
    }
    if(pp.queryarr("fh_chi",temp)) {
        for (int i=0; i<nspecies; ++i) {
        for (int j=0; j<nspecies; ++j) {
            fh_chi(i,j) = temp[i*nspecies+j];
        }
        }
    }
    if(pp.queryarr("fh_monomers",temp)) {
        for (int i=0; i<nspecies; ++i) {
            fh_monomers[i] = temp[i];
        }
    }
    if(use_flory_huggins == 1) {
        for (int i=0; i<nspecies; ++i) {
            molmass[i] = fh_monomers[i]*monomer_mass;
        }
    }
    if(pp.queryarr("c_init_1",temp)) {
        for (int i=0; i<nspecies; ++i) {
            c_init_1[i] = temp[i];
        }
    }
    if(pp.queryarr("c_init_2",temp)) {
        for (int i=0; i<nspecies; ++i) {
            c_init_2[i] = temp[i];
        }
    }
    if(pp.queryarr("Dtherm",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            Dtherm[i] = temp[i];
        }
    }
    if(pp.queryarr("H_diag",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            H_diag[i] = temp[i];
        }
    }
    if(pp.queryarr("Dbar",temp,0,nspecies*(nspecies-1)/2)) {
        for (int i=0; i<nspecies*(nspecies-1)/2; ++i) {
            Dbar[i] = temp[i];
        }
    }
    if(pp.queryarr("H_offdiag",temp,0,nspecies*(nspecies-1)/2)) {
        for (int i=0; i<nspecies*(nspecies-1)/2; ++i) {
            H_offdiag[i] = temp[i];
        }
    }
    pp.query("midpoint_stoch_mass_flux_type",midpoint_stoch_mass_flux_type);
    pp.query("avg_type",avg_type);
    pp.query("mixture_type",mixture_type);
    pp.query("use_charged_fluid",use_charged_fluid);
    pp.query("print_debye_len",print_debye_len);
    pp.query("dielectric_const",dielectric_const);
    pp.query("dielectric_type",dielectric_type);
    if(pp.queryarr("charge_per_mass",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            charge_per_mass[i] = temp[i];
        }
    }
    pp.query("bc_function_type",bc_function_type);
    pp.query("L_pos",L_pos);
    pp.query("L_trans",L_trans);
    pp.query("L_zero",L_zero);
    pp.query("theta_pot",theta_pot);
    pp.query("num_pot_iters",num_pot_iters);
    pp.query("dpdt_factor",dpdt_factor);
    pp.query("E_ext_type",E_ext_type);
    if(pp.queryarr("E_ext_value",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            E_ext_value[i] = temp[i];
        }
    }
    pp.query("electroneutral",electroneutral);
    pp.query("induced_charge_eo",induced_charge_eo);
    pp.query("relxn_param_charge",relxn_param_charge);
    pp.query("zero_eps_on_wall_type",zero_eps_on_wall_type);
    pp.query("zero_charge_on_wall_type",zero_charge_on_wall_type);
    pp.query("zero_eps_on_wall_left_end",zero_eps_on_wall_left_end);
    pp.query("zero_eps_on_wall_right_start",zero_eps_on_wall_right_start);
}
