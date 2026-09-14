#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include "chemistry_functions.H"
#include "common_functions.H"

#include "chemistry_testing_functions.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    if (argc==1) amrex::Abort("ERROR: inputs file expected");

    main_main(argv[1]);

    amrex::Finalize();
    return 0;
}

void main_main(const char* argv)
{
    // **********************************
    // variables defined in src_common

    std::string inputs_file = argv;


    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();

    // print variable values read from inputs file
    amrex::Print() << "\n";

    amrex::Print() << "(src_common param) n_cells: " << n_cells[0] << " " << n_cells[1];
#if AMREX_SPACEDIM==3
    amrex::Print() << " " << n_cells[2] << "\n";
#else
    amrex::Print() << "\n";
#endif

    amrex::Print() << "(src_common param) prob_lo: " << prob_lo[0] << " " << prob_lo[1];
#if AMREX_SPACEDIM==3
    amrex::Print() << " " << prob_lo[2] << "\n";
#else
    amrex::Print() << "\n";
#endif

    amrex::Print() << "(src_common param) prob_hi: " << prob_hi[0] << " " << prob_hi[1];
#if AMREX_SPACEDIM==3
    amrex::Print() << " " << prob_hi[2] << "\n";
#else
    amrex::Print() << "\n";
#endif

    amrex::Print() << "(src_common param) max_step = "  << max_step   << "\n";
    amrex::Print() << "(src_common param) fixed_dt = "  << fixed_dt   << "\n";
    amrex::Print() << "(src_common param) plot_int = "  << plot_int   << "\n";
    amrex::Print() << "(src_common param) prob_type = " << prob_type  << "\n";
    amrex::Print() << "(src_common param) nspecies = "  << nspecies   << "\n";

    amrex::Print() << "(src_common param) molmass: ";
    for (int n=0; n<nspecies; n++)
        amrex::Print() << molmass[n] << " ";
    amrex::Print() << "\n";

    // **********************************
    // variables defined in src_chemistry

    InitializeChemistryNamespace();

    // print variable values read from inputs file
    amrex::Print() << "\n";

    amrex::Print() << "(src_chemistry param) nreaction = " << nreaction << "\n";

    amrex::Print() << "(src_chemistry param) rate_const: ";
    for (int m=0; m<nreaction; m++)
        amrex::Print() << rate_const[m] << " ";
    amrex::Print() << "\n";

    for (int m=0; m<nreaction; m++)
    {
        amrex::Print() << "(src_chemistry param) stoich_coeffs_R_" << m+1 << ": ";
        for (int n=0; n<nspecies; n++)
            amrex::Print() << stoich_coeffs_R(m,n) << " ";
        amrex::Print() << "\n";
    }

    for (int m=0;m<nreaction;m++)
    {
        amrex::Print() << "(src_chemistry param) stoich_coeffs_P_" << m+1 << ": ";
        for (int n=0; n<nspecies; n++)
            amrex::Print() << stoich_coeffs_P(m,n) << " ";
        amrex::Print() << "\n";
    }

    amrex::Print() << "(src_chemistry param) reaction_type = " << reaction_type << "\n";

    amrex::Print() << "\n";

    // **********************************
    // SIMULATION SETUP

    // make BoxArray and Geometry
    // ba will contain a list of boxes that cover the domain
    // geom contains information such as the physical domain size,
    //               number of points in the domain, and periodicity
    BoxArray ba;
    Geometry geom;

    // AMREX_D_DECL means "do the first X of these, where X is the dimensionality of the simulation"
    IntVect dom_lo(AMREX_D_DECL(0, 0, 0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(IntVect(max_grid_size));

    // This defines the physical box, [0,1] in each direction.
    RealBox real_box({AMREX_D_DECL( prob_lo[0], prob_lo[1], prob_lo[2])},
                     {AMREX_D_DECL( prob_hi[0], prob_hi[1], prob_hi[2])});

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // Nghost = number of ghost cells for each array
    int Nghost = 1;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two rho multifabs; one will store the old state, the other the new.
    MultiFab rho_old(ba, dm, nspecies, Nghost);
    MultiFab rho_new(ba, dm, nspecies, Nghost);

    // allocate source MultiFab
    MultiFab source(ba, dm, nspecies, Nghost);

    // time = starting time in the simulation
    amrex::Real time = 0.0;
    amrex::Real dt = fixed_dt;

    amrex::Real dV ;
#if AMREX_SPACEDIM==3
    dV = dx[0]*dx[1]*dx[2];
    amrex::Print() << "dx: " << dx[0] << " " << dx[1] << " " << dx[2] << "\n";
#else
    dV = dx[0]*dx[1];
    amrex::Print() << "dx: " << dx[0] << " " << dx[1] << "\n";
#endif
    amrex::Print() << "dV = " << dV << "\n";
    amrex::Print() << "\n";

    // **********************************
    // INITIALIZE DATA

    // loop over boxes
    for (MFIter mfi(rho_old); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& rhoOld = rho_old.array(mfi);

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
        {
            rhoOld(i,j,k,n) = rho0*rhobar[n];
        });
    }

    // vector to store the name of the components
    // NOTE: its size must be equal to the number of components
    Vector<std::string> var_names(nspecies);
    for (int n=0; n<nspecies; n++) var_names[n] = "spec" + std::to_string(n+1);

    // Write a plotfile of the initial data if plot_int > 0
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,5);
        WriteSingleLevelPlotfile(pltfile, rho_old, var_names, geom, time, 0);
    }

    // mean and variance computed numerically at initial time

    amrex::Print() << "Stats ";
    amrex::Print() << 0 << " ";

    for (int n=0; n<nspecies; n++)
    {
        //amrex::Print() << ComputeSpatialMean(rho_old,n)*(avogadro)/molmass[n] << " ";
        amrex::Print() << ComputeSpatialMean(rho_old,n) << " ";
    }


    for (int n=0; n<nspecies; n++)
    {
        //amrex::Print() << ComputeSpatialVariance(rho_old,n)*((avogadro)/molmass[n])*((avogadro)/molmass[n]) << " ";
        amrex::Print() << ComputeSpatialVariance(rho_old,n) << " ";
    }

    amrex::Print() << "\n";

    // **********************************
    // MAIN LOOP: Time advancement

    for (int step = 1; step <= max_step; ++step)
    {
        // fill periodic ghost cells
        rho_old.FillBoundary(geom.periodicity());

        if (prob_type==1)   // cell-based routines
        {
            for ( MFIter mfi(rho_old); mfi.isValid(); ++mfi )
            {
                const Box& bx = mfi.validbox();

                const Array4<Real>& rhoOld = rho_old.array(mfi);
                const Array4<Real>& rhoNew = rho_new.array(mfi);

                amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
                {
                    GpuArray<amrex::Real,MAX_SPECIES> n_old;
                    GpuArray<amrex::Real,MAX_SPECIES> n_new;
                    for (int n=0; n<nspecies; n++) n_old[n] = rhoOld(i,j,k,n)*(avogadro)/molmass[n];
                    switch(reaction_type){
                        case 0: // deterministic case
                            advance_reaction_det_cell(n_old,n_new,dt);
                            break;
                        case 1: // CLE case
                            advance_reaction_CLE_cell(n_old,n_new,dt,dV,engine);
                            break;
                        case 2: // SSA case
                            advance_reaction_SSA_cell(n_old,n_new,dt,dV,engine);
                            break;
                        default:
                            amrex::Abort("ERROR: invalid reaction_type");
                    }
                    for (int n=0; n<nspecies; n++) rhoNew(i,j,k,n) = n_new[n]*(k_B/Runiv)*molmass[n];
                });
            }
        }
        else if (prob_type==2)  // MultiFab-based routine
        {
            // compute source
            compute_chemistry_source_CLE_1(dt,dV,rho_old,0,source,0);

            for ( MFIter mfi(rho_old); mfi.isValid(); ++mfi )
            {
                const Box& bx = mfi.validbox();

                const Array4<Real>& rhoOld = rho_old.array(mfi);
                const Array4<Real>& rhoNew = rho_new.array(mfi);

                const Array4<Real>& sourceArr = source.array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int n=0; n<nspecies; n++) rhoNew(i,j,k,n) = rhoOld(i,j,k,n) + dt*sourceArr(i,j,k,n);
                });
            }
        }
        else if (prob_type==3)
        {
            RK3step_chem_only(rho_old,rho_new,geom,dt);
        }
        else if (prob_type==4)
        {
            EMstep_chem_only(rho_old,rho_new,geom,dt);
        }
        else
        {
            amrex::Abort("ERROR: invalid prob_type");
        }

        // update time
        time = time + dt;

        // copy new solution into old solution
        MultiFab::Copy(rho_old, rho_new, 0, 0, nspecies, 0);

        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << step << "\n";

        // print out mean and variance of both species in one single line at each time step
        amrex::Print()  << "Stats ";
        amrex::Print()  << dt*step << " ";
        for (int n=0; n<nspecies; n++)
        {
            //amrex::Print()  << ComputeSpatialMean(rho_new,n)*(avogadro)/molmass[n] << " ";
            amrex::Print()  << ComputeSpatialMean(rho_new,n) << " ";
        }
        for (int n=0; n<nspecies; n++)
        {
            //amrex::Print()  << ComputeSpatialVariance(rho_new,n)*((avogadro)/molmass[n])*((avogadro)/molmass[n]) << " ";
            amrex::Print()  << ComputeSpatialVariance(rho_new,n) << " ";
        }
        amrex::Print() << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && step%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",step,5);
            WriteSingleLevelPlotfile(pltfile, rho_new, var_names, geom, time, step);
        }
    }

    return;
}
