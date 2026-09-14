
//#include "hydro_test_functions.H"
//#include "hydro_test_functions_F.H"

//#include "hydro_functions.H"
//#include "StochMomFlux.H"
//#include "StructFact.H"


#include "common_functions.H"


//#include "compressible_functions.H"

#include "exec_functions.H"
#include "GL_functions.H"
#include "GL_functions_F.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{


    amrex::Real umbrella; //spring constant
    amrex::Real phi0; // umbrella center
    amrex::Real alpha; // spring constant scaling parameter
    amrex::Real r1; //  phi0 step parameter, see overleaf notes
    amrex::Real r2; //  overlap parameter, see overleaf notes
    int N_Burn; // This corresponds to  the integer "Equil" input. This is the number of time-steps taken before data collection begins
    int L; //Integer corresponding to the "Number_of_Samples" input. This is the number of samples considered in each umbrella
    int Plot_Skip; //An integer used to limit the number of plot files saved. See input file and overleaf for more information.
    int adaptive; //  1= adaptive method used, 0= uniform phi_0 step method. Specified at input
    int Reverse; //! 1= "forward" direction (i.e phi_0) is increased, 0= "backward" direction (i.e phi_0) is decreased. Specified at input

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    // read in parameters from inputs file into F90 modules
    // we use "+1" because of amrex_string_c_to_f expects a null char termination
    read_common_namelist(inputs_file.c_str(),inputs_file.size()+1);
    read_GL_namelist(inputs_file.c_str(),inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    set_inputs(&Plot_Skip,&L,&N_Burn,&alpha,&r1,&r2,&adaptive,&Reverse); //read-in inputs added by Kevin

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1) {
            is_periodic[i] = 1;
            Print() << "Periodic: " << is_periodic[i] << "\n";
        }
        //is_periodic[i] = 0;
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
        IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        // note we are converting "Vector<int> max_grid_size" to an IntVect
        ba.maxSize(IntVect(max_grid_size));

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                         {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    Real dt = fixed_dt;

    const Real* dx = geom.CellSize();

    setdt(dx, &dt);

    Real dtinv = 1.0/dt;

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////
    const int n_rngs = 1;

    const int proc = ParallelDescriptor::MyProc();

    int fhdSeed = 0;

    fhdSeed += 10000*proc;

    //Initialise rngs
    rng_init(&fhdSeed);
    /////////////////////////////////////////

    //conserved quantaties
    MultiFab phi  (ba,dmap,nvars,ngc);
    MultiFab phin (ba,dmap,nvars,ngc);

    //possibly for later
    MultiFab rannums(ba,dmap,nvars,ngc);


    Real time = 0;

    int step, statsCount;

    step = 0;

    statsCount = 1;

     Init_Phi(phi,dx);
     Param_Output(&umbrella,&phi0);
     WritePlotFile(step, time, geom, phi,umbrella,phi0);

     if( proc == 0 ) {

        std::cout << "DT = " << dt << std::endl;
        std::cout << "Mesh spacing = " << dx[0] << " , " << dx[1] << std::endl;
        std::cout << " Domain = [ 0 , " << prob_hi[0] << " ] x [ 0 , " << prob_hi[1] <<" ]" << std::endl;

     }


    //Time stepping loop

    bool Make_PltFiles = true; //boolean used as flag for making plot files
    int umbrella_size=0; //Serves as a flag for the spring constant strength (1=has attained smallest allowable value,2=has attained largest allowable value,0=neither too large or small)
    amrex::Real Expec;// The average of samples collected in PREVIOUS accepted umbrella
    amrex::Real MAD; // The median absolute deviation ( measure of spread) of the data in PREVIOUS accepted umbrella
    amrex::Real Expec2; // The average of samples collected in CURRENT candidate umbrella
    amrex::Real MAD2;// The median absolute deviation ( measure of spread) of the data in candidate accepted umbrella
    bool sucessful_compare=false; // boolean that states whether there was sufficient overlap with PREVIOUS accepted umbrella
    int Shift_Flag=0; //-Integer that indicates what "approach" is used to compute phi_0. See fortran "inc_phi0_Adapt" subroutine comments for more information
    bool while_loop_comp=true;// boolean that is used to break the while loop under certain conditions. ( set to false to break the loop)
    bool First_Loop_Step=true;// boolean that is used to differentiate the situation when entering the while loop after a previous sucessful umbrella comparison. (true= in first step, false=NOT in first step)
    bool weak_umb=false; // boolean that is used to skip comparisons IF both the previous umbrella comparison and current comparison were sucessful with weakest possible spring constant
    int Plot_Num=0;//  An integer keeping track of the contigous plot file number count
    int umbrella_number=0; //An integer keeping track of the contigous umbrella file number count

    if(adaptive==1 and Reverse==0) // adaptive algorithm in "forward" direction (i.e increasing phi_0)
    {
        Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                N_Burn,L,Expec,MAD,Plot_Num,Plot_Skip,umbrella_number);
        inc_phi0_Adapt(&Expec,&MAD,&r1,&Shift_Flag);
        Make_PltFiles = false;

        while((Expec+r1*MAD) <1.0) // cut-off condition for adaptive method
        {
                Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                        N_Burn,L,Expec2,MAD2,Plot_Num,Plot_Skip,umbrella_number);
                Check_Overlap(Expec,MAD,Expec2,MAD2,r2,alpha,sucessful_compare,umbrella_size,Shift_Flag,while_loop_comp,First_Loop_Step,weak_umb);
                First_Loop_Step=false;
            if(sucessful_compare and  while_loop_comp and !weak_umb)// If first comparson is sucessful, the adaptive algorithm tries to optimize spring constant values
            {
                while(sucessful_compare and  while_loop_comp and !weak_umb)// This while loop continues until a spring constant leading to an unsucessful comparsion is found
                {
                    Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                        N_Burn,L,Expec2,MAD2,Plot_Num,Plot_Skip,umbrella_number);
                    Check_Overlap(Expec,MAD,Expec2,MAD2,r2,alpha,sucessful_compare,umbrella_size,Shift_Flag,while_loop_comp,First_Loop_Step,weak_umb);
                }
            }else // If first comparson is NOT sucessful, the adaptive algorithm tries to find a spring constant and phi_0 that does work. The first such pair is used.
            {
                while(!sucessful_compare)// This while loop continues until a successful combination of phi_0 and spring constant is found
                {
                    Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                        N_Burn,L,Expec2,MAD2,Plot_Num,Plot_Skip,umbrella_number);
                    Check_Overlap(Expec,MAD,Expec2,MAD2,r2,alpha,sucessful_compare,umbrella_size,Shift_Flag,while_loop_comp,First_Loop_Step,weak_umb);
                }
            }
                //Below, the sucessful phi_0 and spring constant values are used to compute the accepted umbrella data. This data is written to an umbrellaxxxxx.txt file.
                //After this, phi_0 is shifted, and all flags are reset appropriately
                Make_PltFiles = true;
                Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                        N_Burn,L,Expec,MAD,Plot_Num,Plot_Skip,umbrella_number);
                Shift_Flag=1;
                inc_phi0_Adapt(&Expec,&MAD,&r1,&Shift_Flag);
                umbrella_size=0;
                Make_PltFiles = false;
                First_Loop_Step=true;
                while_loop_comp=true;
                weak_umb=false;
                sucessful_compare=true;

        }
    }
    if(adaptive==1 and Reverse==1) // adaptive algorithm in "backward" direction (i.e decreasing phi_0)
    {
        Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                N_Burn,L,Expec,MAD,Plot_Num,Plot_Skip,umbrella_number);
        r1=-r1;
        inc_phi0_Adapt(&Expec,&MAD,&r1,&Shift_Flag);
        r1=-r1;
        Make_PltFiles = false;

        while((Expec+r1*MAD)>0)
        {
                Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                        N_Burn,L,Expec2,MAD2,Plot_Num,Plot_Skip,umbrella_number);
                Check_Overlap_Backwards(Expec,MAD,Expec2,MAD2,r2,alpha,sucessful_compare,umbrella_size,Shift_Flag,while_loop_comp,First_Loop_Step,weak_umb);
                First_Loop_Step=false;
            if(sucessful_compare and  while_loop_comp and !weak_umb)// If first comparson is sucessful, the adaptive algorithm tries to optimize spring constant values
            {
                while(sucessful_compare and  while_loop_comp and !weak_umb)// This while loop continues until a spring constant leading to an unsucessful comparsion is found
                {
                    Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                        N_Burn,L,Expec2,MAD2,Plot_Num,Plot_Skip,umbrella_number);
                    Check_Overlap_Backwards(Expec,MAD,Expec2,MAD2,r2,alpha,sucessful_compare,umbrella_size,Shift_Flag,while_loop_comp,First_Loop_Step,weak_umb);
                }
            }else// If first comparson is NOT sucessful, the adaptive algorithm tries to find a spring constant and phi_0 that does work. The first such pair is used.
            {
                while(!sucessful_compare)// This while loop continues until a successful combination of phi_0 and spring constant is found
                {
                    Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                        N_Burn,L,Expec2,MAD2,Plot_Num,Plot_Skip,umbrella_number);
                    Check_Overlap_Backwards(Expec,MAD,Expec2,MAD2,r2,alpha,sucessful_compare,umbrella_size,Shift_Flag,while_loop_comp,First_Loop_Step,weak_umb);
                }
            }
                //Below, the sucessful phi_0 and spring constant values are used to compute the accepted umbrella data. This data is written to an umbrellaxxxxx.txt file.
                //After this, phi_0 is shifted, and all flags are reset appropriately
                Make_PltFiles = true;
                Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                        N_Burn,L,Expec,MAD,Plot_Num,Plot_Skip,umbrella_number);
                Shift_Flag=1;
                r1=-r1;
                inc_phi0_Adapt(&Expec,&MAD,&r1,&Shift_Flag);
                r1=-r1;
                umbrella_size=0;
                Make_PltFiles = false;
                First_Loop_Step=true;
                while_loop_comp=true;
                weak_umb=false;
                sucessful_compare=true;
        }
    }
    if(adaptive==0)// direct computation of umbrella data with uniform phi_0 steps and fixed spring constant
    {
        Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                N_Burn,L,Expec,MAD,Plot_Num,Plot_Skip,umbrella_number);
        int inc_dir=1; // forward direction
        fixed_inc_phi0(&inc_dir);

        while((Expec+r1*MAD) <0.24)
        {
        Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                N_Burn,L,Expec,MAD,Plot_Num,Plot_Skip,umbrella_number);
        fixed_inc_phi0(&inc_dir);
        }
        inc_dir=0; // forward direction
        while((Expec+r1*MAD) >0.00286481020)
        {
        Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
                N_Burn,L,Expec,MAD,Plot_Num,Plot_Skip,umbrella_number);
        fixed_inc_phi0(&inc_dir);
        }
    }

    // ****NOTE: Another portion can be added to the total "trip" by copying and pasting the while loops as done below
    //*****NOTE: the while loop conditions should all be consistent with each other when taking as a sequence from top to bottom (i.e the exit condition for a prior while loop should allow you to enter the next loop)


    // while((Expec+r1*MAD)>0.004)
    // {
    //         Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
    //                 N_Burn,L,Expec2,MAD2,Plot_Num,Plot_Skip,umbrella_number);
    //         Check_Overlap_Backwards(Expec,MAD,Expec2,MAD2,r2,alpha,sucessful_compare,umbrella_size,Shift_Flag,while_loop_comp,First_Loop_Step,weak_umb);
    //         First_Loop_Step=false;
    //     if(sucessful_compare and  while_loop_comp and !weak_umb)
    //     {
    //         while(sucessful_compare and  while_loop_comp and !weak_umb)
    //         {
    //             Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
    //                 N_Burn,L,Expec2,MAD2,Plot_Num,Plot_Skip,umbrella_number);
    //             Check_Overlap_Backwards(Expec,MAD,Expec2,MAD2,r2,alpha,sucessful_compare,umbrella_size,Shift_Flag,while_loop_comp,First_Loop_Step,weak_umb);
    //         }
    //     }else
    //     {
    //         while(!sucessful_compare)
    //         {
    //             Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
    //                 N_Burn,L,Expec2,MAD2,Plot_Num,Plot_Skip,umbrella_number);
    //             Check_Overlap_Backwards(Expec,MAD,Expec2,MAD2,r2,alpha,sucessful_compare,umbrella_size,Shift_Flag,while_loop_comp,First_Loop_Step,weak_umb);
    //         }
    //     }

    //     Make_PltFiles = true;
    //     Run_Steps(phi,phin,rannums,geom,dx,dt, time,plot_int,Make_PltFiles,
    //             N_Burn,L,Expec,MAD,Plot_Num,Plot_Skip,umbrella_number);
    //     Shift_Flag=1;
    //     r1=-r1;
    //     inc_phi0_Adapt(&Expec,&MAD,&r1,&Shift_Flag);
    //     r1=-r1;
    //     umbrella_size=0;
    //     Make_PltFiles = false;
    //     First_Loop_Step=true;
    //     while_loop_comp=true;
    //     weak_umb=false;
    //     sucessful_compare=true;
    // }


    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

//    amrex::Finalize();

}

