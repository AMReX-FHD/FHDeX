#include "GL_functions.H"
#include "GL_functions_F.H"

#include "common_functions.H"


#include "exec_functions.H"

#include "AMReX_ArrayLim.H"

#include <algorithm>

using namespace amrex;

void RK2step(MultiFab& phi, MultiFab& phin, MultiFab& rannums,
               const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt, amrex::Real& integral, int n, amrex::Real& phi_avg,amrex::Real& energy,amrex::Real& teng,amrex::Real& H1_semi_norm)
{

    const int comp=0;
    int ncomp=1;


    for (MFIter mfi(rannums); mfi.isValid(); ++mfi) {

        const Box& validBox = mfi.validbox();

        //multifab_fill_random(ARLIM_2D(validBox.loVect()), ARLIM_2D(validBox.hiVect()),
        multifab_fill_random(BL_TO_FORTRAN_BOX(validBox),
                             BL_TO_FORTRAN_ANYD(rannums[mfi]), &ncomp, &comp);
    }

    phi.FillBoundary(geom.periodicity());

    // iterating over the data in phi multifab to compute integral term
    integral = 0.;
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        integrate_1D(BL_TO_FORTRAN_BOX(bx),
                phi[mfi].dataPtr(),
                dx,
                &integral);

    }
    ParallelDescriptor::ReduceRealSum(integral);

                //    if(ParallelDescriptor::MyProc() == 0 ){
                //           std::cout << "integral = "<< integral << std::endl;
                //    }


    // iterating over the data in phi multifab to compute phi for next time step.
    //The function called also computes the averages of phi and some energy quantitites
    energy =0.;
    teng =0.;
    phi_avg=0.0;
    H1_semi_norm = 0.;
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        rk2_stage1_1D(BL_TO_FORTRAN_BOX(bx),
                phi[mfi].dataPtr(),
                phin[mfi].dataPtr(),
                rannums[mfi].dataPtr(),
                &integral,
                &energy, &teng,&H1_semi_norm,
                dx, &dt,&phi_avg);
    }
    ParallelDescriptor::ReduceRealSum(phi_avg);
    phi_avg = phi_avg/(n_cells[0]);
    ParallelDescriptor::ReduceRealSum(energy);
    ParallelDescriptor::ReduceRealSum(teng);
    ParallelDescriptor::ReduceRealSum(H1_semi_norm);


                        //if(ParallelDescriptor::MyProc() == 0 ){
                        //       std::cout << n << " " << energy << "  energy  " << std::endl;
                        //       std::cout << n << " " << teng << "  teng  " << std::endl;
                        //      // std::cout << time << " " << energy << "  energy = "<< energy << "  dphi/dt  = " << phit << std::endl;
                    // }

                    //    phin.FillBoundary(geom.periodicity());
                    //
                    //
                    //    for ( MFIter mfi(cu); mfi.isValid(); ++mfi)
                    //    {
                    //        const Box& bx = mfi.validbox();
                    //
                    //        rk2_stage2(ARLIM_2D(bx.loVect()), ARLIM_2D(bx.hiVect()),
                    //                   phi[mfi].dataPtr(),
                    //                   phin[mfi].dataPtr(),
                    //                   rannums[mfi].dataPtr(),
                    //                     &integral,
                    //                     ZFILL(dx), &dt);
                    //    }
}

void Init_Phi(MultiFab& phi, const amrex::Real* dx )
{
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        initphi_1D(BL_TO_FORTRAN_BOX(bx),
                phi[mfi].dataPtr(),
                dx);
    }
}

void Run_Steps(MultiFab& phi, MultiFab& phin, MultiFab& rannums, const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt,
                        amrex::Real& time, int plot_int, bool Make_PltFiles,int N_Burn,int L, amrex::Real& Expec, amrex::Real& MAD, int& Plot_Num,int& Plot_Skip, int& umbrella_number)
{
    // This function computes all the data collected in a single umbrella by marching the explcit method forward.

    // After some equilibration period, the data is collected, with plot files saved on occasion. The umbrella data is saved in a
    // .txt" file labled with the current umbrella number.  After this, statistics on the data are computed and set as output.

    // INPUT:
    // phi -- current time step phi field multifab
    // phin -- current time step phi field multifab
    // rannums -- random number multifab
    // geom -- amrex Geometry object
    // dx -- an array of doubles that specifies the spatial grid width
    // dt -- a double that specifies the time step
    // time -- time quantity needed in plot file function
    // plot_int -- An integer specified during input. It is the number of steps that determines a solution "snapshot"/"sample"
    // Make_PltFiles -- A boolean that specifies when plot files are to be saved
    // N_Burn -- This corresponds to  the integer "Equil" input. This is the number of time-steps taken before data collection begins
    // L -- Integer corresponding to the "Number_of_Samples" input. This is the number of samples considered in each umbrella
    // Plot_Num -- An integer keeping track of the contigous plot file number count
    // Plot_Skip -- An integer used to limit the number of plot files saved. See input file and overleaf for more information.
    // umbrella_number -- An integer keeping track of the contigous umbrella file number count

    // OUTPUT:
    // Expec -- The average of samples collected in umbrella
    // MAD -- The median absolute deviation of the data in umbrella (i.e a measure of spread, see overleaf for more info.)




    int Steps_Counted;  // counter that tracks all steps taken
    int Total_Steps=L*plot_int; // Total steps to be taken not considering equilibration steps
    int j=0;

    // Below, quantities evaluated in the RK function are initialized
    amrex::Real integral=0.0;
    amrex::Real Expec_Temp=0.0;
    amrex::Real Phi_Avg=0.0;
    amrex::Real energy=0.0;
    amrex::Real teng=0.0;
    amrex::Real H1_semi_norm=0.0;


    // Below, quantitites used to compute statistics are intitialized
    amrex::Real median;
    double *Avg_collect= new double[L]; //initialize array that will hold spatial averages
    Avg_collect[L]={0};

    // Below, quantitites specifing the umbrella are set so that they may be written to a file
    amrex::Real umbrella;
    amrex::Real phi0;
    Param_Output(&umbrella,&phi0);


    // The block below creates a "umbrellaxxxxxxx.txt" where the x's are the current umbrella number.
    // This is used as the name for the umbrella file
    std::stringstream ss;
    ss << std::setw(8) << std::setfill('0') << umbrella_number;
    std::string s = ss.str();
    std::string umb_Num;
    umb_Num="umbrella"+ s;
    umb_Num=umb_Num+".txt";

    // The block below sets the formatting of the output for the two output files created in this function
    // The two files are the individual "umbrella" files and a file that keeps all the umbrella data in one contigous file
    std::ofstream ofs;
    std::ofstream ofs2;

    ofs.setf(std::ios::scientific);
    ofs.setf(std::ios::showpos);
    ofs.precision(13);
    ofs2.setf(std::ios::scientific);
    ofs2.setf(std::ios::showpos);
    ofs2.precision(13);

                // std::ofstream ofs3;
                // ofs3.open("Console_output_Fortran.txt", std::ofstream::out | std::ofstream::app);



    if(ParallelDescriptor::MyProc() == 0  and Make_PltFiles)
    {
        // At the top of every "umbrella" file, the first line is the spring constant, and the second line is phi_0
        ofs.open (umb_Num.c_str(), std::ofstream::out | std::ofstream::app);
        ofs << umbrella<< "\n";
        ofs << phi0<< "\n";
        ofs.close();
    }





    // The loop below runs through all steps including equilibration steps and desired sample steps
    for(int step=1;step<=Total_Steps+N_Burn;++step)
    {
        // Rk2 step computes the field at the next time step and the spatial average of the field
        H1_semi_norm = 0.;
        RK2step(phi, phin, rannums, geom, dx, dt, integral, step,Phi_Avg,energy,teng,H1_semi_norm);

        if(step>N_Burn) // statistics and data are taken only after equilibration period
        {
            Steps_Counted=step-N_Burn;
            if(step%500 == 0)
            {
                // supressed output of older version

                //amrex::Print() << "Integral " << integral << "\n";
                //amrex::Print() << "Advanced step " << step << "\n";
            }
            if (plot_int > 0 && Steps_Counted%plot_int == 0) // collect data every "plot_int" steps
            {
                Avg_collect[j]=Phi_Avg;
                Expec_Temp=Expec_Temp+Phi_Avg; // sum all averages
                j=j+1;
                amrex::Print() << "The average of phi is " << Phi_Avg<< "\n";

                if(Make_PltFiles)
                {
                    if(ParallelDescriptor::MyProc() == 0 )
                    {
                        ofs2.open("AVG_DATA.txt", std::ofstream::out | std::ofstream::app);// write a file that holds all output of every umbrella
                        ofs2 << Phi_Avg << " " << energy << " " << teng  << " " << H1_semi_norm << "\n";
                        ofs2.close();

                        ofs.open (umb_Num.c_str(), std::ofstream::out | std::ofstream::app);// write data to an "umbrella" file
                        ofs << Phi_Avg << " " << energy << " " << teng  << " " << H1_semi_norm << "\n";
                        ofs.close();
                    }
                    if((j+1)%Plot_Skip==0) // only save plot file at multiples of "Plot_Skip" parameter
                    {
                        Plot_Num+=1; //increment index used in plot file name
                        WritePlotFile(Plot_Num, time, geom, phi, umbrella,phi0);
                    }
                }

            }
        }
        time = time + dt;
    }

    Expec=Expec_Temp/L; //compute average of field averages in umbrella
    std::sort(Avg_collect, Avg_collect+L); //sort field averages in umbrella
    median=Avg_collect[L/2];  //find median of field averages in umbrella
                    // if(ParallelDescriptor::MyProc() == 0)
                    // {
                    //     ofs3 << "median of Phi in umbrella is " << median << "\n";
                    // }
    amrex::Print() << "median of Phi in umbrella is " << median << "\n";
    for(j=0;j<L;j++)
    {
        Avg_collect[j]=amrex::Math::abs(Avg_collect[j]-median) ;// compute the distance from the median for each data point
    }
    std::sort(Avg_collect, Avg_collect+L); // Sort the distances from the median--the "median absolute deviation"
    MAD=Avg_collect[L/2]; //The median of all such distances is the measure of spread used--the "median absolute deviation"
    amrex::Print() << "Average of Phi in umbrella is " << Expec << "\n";
    amrex::Print() << "MAD of Umbrella is " << MAD << "\n";
                    // if(ParallelDescriptor::MyProc() == 0)
                    // {
                    //     ofs3 << "Average of Phi in umbrella is " << Expec << "\n";
                    //     ofs3 << "MAD of Umbrella is " << MAD << "\n";
                    // }
    if(Make_PltFiles)
    {
        umbrella_number=umbrella_number+1; // when an umbrella file is saved, increment the index used to name umbrella files
    }
                    //ofs3.close();
    delete [] Avg_collect;
}


void Check_Overlap(amrex::Real& Expec,amrex::Real& MAD,amrex::Real& Expec2,amrex::Real& MAD2,amrex::Real& r2,amrex::Real& alpha, bool& sucessful_compare, int& umbrella_size, int& Shift_Flag, bool& while_loop_comp, bool& First_Loop_Step, bool& weak_umb)
{
    // This function is the implementation of the algorithm detailed in the overleaf notes. This function is meant for the scenario
    // where we are INCREASING values of phi_0. (i.e moving "forward")
    // This function checks the overlap of the current umbrella with the previous umbrella and adjusts spring constant kappa and phi_0 accordingly

    // INPUT:
    // Expec -- The average of samples collected in PREVIOUS umbrella
    // MAD -- The median absolute deviation of the data in PREVIOUS umbrella
    // Expec2 -- The average of samples collected in CURRENT umbrella
    // MAD2 -- The median absolute deviation of the data in CURRENT umbrella
    // r2 -- double that tweaks the overlap condition
    // alpha -- double that scales the spring constant
    // sucessful_compare -- boolean that states whether there was sufficient overlap with PREVIOUS umbrella
    // umbrella_size -- integer that is set to 0,1,or 2. Serves as a flag for the spring constant strength (1=has attained smallest allowable value,2=has attained largest allowable value,0=neither too large or small)
    // Shift_Flag --Integer that indicates what "approach" is used to compute phi_0. See fortran "inc_phi0_Adapt" subroutine comments for more information
    // while_loop_comp -- boolean that is used to break the while loop under certain conditions. ( set to false to break the loop)
    // First_Loop_Step -- boolean that is used to differentiate the situation when entering the while loop after a previous sucessful umbrella comparison. (true= in first step, false=NOT in first step)
    // weak_umb -- boolean that is used to skip comparisons IF both the previous umbrella comparison and current comparison were sucessful with weakest possible spring constant

    // OUTPUT:
    // sucessful_compare -- boolean that states whether there was sufficient overlap with CURRENT umbrella
    // umbrella_size -- integer that is set to 0,1,or 2. Serves as a flag for the spring constant strength (1=has attained smallest allowable value,2=has attained largest allowable value,0=neither too large or small)
    // Shift_Flag --Integer that indicates what "approach" is used to compute phi_0. See fortran "inc_phi0_Adapt" subroutine comments for more information
    // weak_umb -- boolean that is used to skip comparisons IF both the previous umbrella comparison and current comparison were sucessful with weakest possible spring constant

    int sucessful_iter; //integer version of sucessful_compare for passing to fortran (1=sucessful comparison, 0=NOT sucessful comparison)
    int sucessful_iter_prev; //integer version of sucessful_compare FOR PREVIOUS umbrella for passing to fortran (1=sucessful comparison, 0=NOT sucessful comparison)
    int Umbrella_Size_Prev; // umbrella_size parameter for previous umbrella
    amrex::Real r_temp; // temorary double for computations
    amrex::Real umbrella_reset_val; //CHANGE TO INPUT

    Umbrella_Size_Prev=umbrella_size;

    if(sucessful_compare) // convert umbrella overlap boolean for PREVIOUS umbrella comparison to integer for passing to Fortran
    {
        sucessful_iter_prev=1;
    }else
    {
        sucessful_iter_prev=0;
    }

                // std::ofstream ofs;
                // ofs.open("Console_output_Fortran.txt", std::ofstream::out | std::ofstream::app);
                // if(ParallelDescriptor::MyProc() == 0 )
                // {
                //     ofs << "E2-r2*S2 "  << Expec2-r2*MAD2 << "\n";
                //     ofs << "E2+r2*S2 "  << Expec2+r2*MAD2 << "\n";
                //     ofs << "E1+r2*S1 "  << Expec + r2*MAD <<  "\n";
                // }

    // print umbrella statistics to screen
    amrex::Print() << "E2-r2*S2 "  << Expec2-r2*MAD2 << "\n";
    amrex::Print() << "E2+r2*S2 "  << Expec2+r2*MAD2 << "\n";
    amrex::Print() << "E1+r2*S1 "  << Expec + r2*MAD <<  "\n";



    if(Expec2-r2*MAD2 < Expec + r2*MAD && Expec2+r2*MAD2 > Expec + r2*MAD)// check if there was sufficient overlap AND new umbrella does not overlap too much
    {

        sucessful_compare=true;
        sucessful_iter=1;
        Shift_Flag=0;
        Umbrella_Adjust(&sucessful_iter,&alpha,&umbrella_size,&sucessful_iter_prev);
        amrex::Print() << "Overlap occured "  << "\n";
                // if(ParallelDescriptor::MyProc() == 0 )
                // {
                //     ofs << "Overlap occured "  << "\n";
                // }
        if(umbrella_size==1 &&  Umbrella_Size_Prev==1)
        {
            weak_umb=true;
        }


    }else if (Expec2-r2*MAD2 > Expec + r2*MAD) // Consider the case where new umbrella does not overlap enough
    {
        sucessful_compare=false;
        Shift_Flag=0;
        sucessful_iter=0;
        Umbrella_Adjust (&sucessful_iter,&alpha ,&umbrella_size,&sucessful_iter_prev);
        amrex::Print() << "Overlap did not occur"  << "\n";
        amrex::Print() << umbrella_size  << "\n";
                    // if(ParallelDescriptor::MyProc() == 0 )
                    // {
                    //     ofs << "Overlap did not occur"  << "\n";
                    //     ofs << umbrella_size  << "\n";
                    // }
        if(umbrella_size==2)// shift phi_0 down and reset the spring constant to a low value if there was insufficient overlap when the spring constant upper limit is reached
        {
            r_temp=r2;
            r2=-1.25*r2;
            inc_phi0_Adapt(&Expec,&MAD,&r2,&Shift_Flag);
            r2=r_temp;
            umbrella_reset_val=65.0;
            umbrella_reset(&umbrella_reset_val);
            umbrella_size=0;
        }
   }
   else if (Expec2+r2*MAD2 < Expec + r2*MAD) // Consider the case where too much of the current umbrella overlaps
    {
        sucessful_compare=false;
        Shift_Flag=1;
        amrex::Print() << "New Phi_0 resulted in lower average"  << "\n";
        amrex::Print() << "Increasing Phi_0"  << "\n";
                    // if(ParallelDescriptor::MyProc() == 0 )
                    // {
                    //     ofs << "New Phi_0 resulted in lower average"  << "\n";
                    //     ofs << "Increasing Phi_0"  << "\n";
                    // }
        inc_phi0_Adapt(&Expec,&MAD,&r2,&Shift_Flag);
    }
    if(sucessful_iter_prev==1 and sucessful_iter==0 and !First_Loop_Step) // break the "sucess" while loop when a comparison fails(after many sucessful comparisons)
    {
        while_loop_comp=false;
    }
    if(sucessful_iter_prev==1 and sucessful_iter==1 and  umbrella_size==1) // break the "sucess" while loop after contigously sucessful comparisons with spring constants at lower limit
    {
        weak_umb=true;
    }
                    // ofs.close();
}





void Check_Overlap_Backwards(amrex::Real& Expec,amrex::Real& MAD,amrex::Real& Expec2,amrex::Real& MAD2,amrex::Real& r2,amrex::Real& alpha, bool& sucessful_compare, int& umbrella_size, int& Shift_Flag, bool& while_loop_comp, bool& First_Loop_Step, bool& weak_umb)
{
    // This function is the implementation of the algorithm detailed in the overleaf notes. This function is meant for the scenario
    // where we are DECREASING values of phi_0. (i.e moving "backward")
    // This function checks the overlap of the current umbrella with the previous umbrella and adjusts spring constant kappa and phi_0 accordingly

    // INPUT:
    // Expec -- The average of samples collected in PREVIOUS umbrella
    // MAD -- The median absolute deviation of the data in PREVIOUS umbrella
    // Expec2 -- The average of samples collected in CURRENT umbrella
    // MAD2 -- The median absolute deviation of the data in CURRENT umbrella
    // r2 -- double that tweaks the overlap condition
    // alpha -- double that scales the spring constant
    // sucessful_compare -- boolean that states whether there was sufficient overlap with PREVIOUS umbrella
    // umbrella_size -- integer that is set to 0,1,or 2. Serves as a flag for the spring constant strength (1=has attained smallest allowable value,2=has attained largest allowable value,0=neither too large or small)
    // Shift_Flag --Integer that indicates what "approach" is used to compute phi_0. See fortran "inc_phi0_Adapt" subroutine comments for more information
    // while_loop_comp -- boolean that is used to break the while loop under certain conditions. ( set to false to break the loop)
    // First_Loop_Step -- boolean that is used to differentiate the situation when entering the while loop after a previous sucessful umbrella comparison. (true= in first step, false=NOT in first step)
    // weak_umb -- boolean that is used to skip comparisons IF both the previous umbrella comparison and current comparison were sucessful with weakest possible spring constant

    // OUTPUT:
    // sucessful_compare -- boolean that states whether there was sufficient overlap with CURRENT umbrella
    // umbrella_size -- integer that is set to 0,1,or 2. Serves as a flag for the spring constant strength (1=has attained smallest allowable value,2=has attained largest allowable value,0=neither too large or small)
    // Shift_Flag --Integer that indicates what "approach" is used to compute phi_0. See fortran "inc_phi0_Adapt" subroutine comments for more information
    // weak_umb -- boolean that is used to skip comparisons IF both the previous umbrella comparison and current comparison were sucessful with weakest possible spring constant

    int sucessful_iter; //integer version of sucessful_compare for passing to fortran (1=sucessful comparison, 0=NOT sucessful comparison)
    int sucessful_iter_prev; //integer version of sucessful_compare FOR PREVIOUS umbrella for passing to fortran (1=sucessful comparison, 0=NOT sucessful comparison)
    int Umbrella_Size_Prev; // umbrella_size parameter for previous umbrella
    amrex::Real r_temp; // temorary double for computations
    amrex::Real umbrella_reset_val; //CHANGE TO INPUT

    Umbrella_Size_Prev=umbrella_size;

    if(sucessful_compare)// convert umbrella overlap boolean for PREVIOUS umbrella comparison to integer for passing to Fortran
    {
        sucessful_iter_prev=1;
    }else
    {
        sucessful_iter_prev=0;
    }


                // std::ofstream ofs;
                // ofs.open("Console_output_Fortran.txt", std::ofstream::out | std::ofstream::app);
                // if(ParallelDescriptor::MyProc() == 0 )
                // {
                //     ofs << "E2-r2*S2 "  << Expec2-r2*MAD2 << "\n";
                //     ofs << "E2+r2*S2 "  << Expec2+r2*MAD2 << "\n";
                //     ofs << "E1+r2*S1 "  << Expec + r2*MAD <<  "\n";
                // }
    amrex::Print() << "E2-r2*S2 "  << Expec2-r2*MAD2 << "\n";
    amrex::Print() << "E2+r2*S2 "  << Expec2+r2*MAD2 << "\n";
    amrex::Print() << "E1-r2*S1 "  << Expec - r2*MAD <<  "\n";

    // print umbrella statistics to screen
    if(Expec-r2*MAD < Expec2 + r2*MAD2 && Expec2-r2*MAD2 < Expec - r2*MAD)// check if there was sufficient overlap AND new umbrella does not overlap too much
    {

        sucessful_compare=true;
        sucessful_iter=1;
        Shift_Flag=0;
        Umbrella_Adjust(&sucessful_iter,&alpha,&umbrella_size,&sucessful_iter_prev);
        amrex::Print() << "Overlap occured "  << "\n";
                    // if(ParallelDescriptor::MyProc() == 0 )
                    // {
                    //     ofs << "Overlap occured "  << "\n";
                    // }
        if(umbrella_size==1 &&  Umbrella_Size_Prev==1)
        {
            weak_umb=true;
        }


    }else if (Expec-r2*MAD > Expec2 + r2*MAD2)// Consider the case where new umbrella does not overlap enough
    {
        sucessful_compare=false;
        Shift_Flag=0;
        sucessful_iter=0;
        Umbrella_Adjust (&sucessful_iter,&alpha ,&umbrella_size,&sucessful_iter_prev);
        amrex::Print() << "Overlap did not occur"  << "\n";
        amrex::Print() << umbrella_size  << "\n";
                    // if(ParallelDescriptor::MyProc() == 0 )
                    // {
                    //     ofs << "Overlap did not occur"  << "\n";
                    //     ofs << umbrella_size  << "\n";
                    // }
        if(umbrella_size==2)// shift phi_0 up and reset the spring constant to a low value if there was insufficient overlap when the spring constant upper limit is reached
        {
            r_temp=r2;
            r2=1.25*r2;
            inc_phi0_Adapt(&Expec,&MAD,&r2,&Shift_Flag);
            r2=r_temp;
            umbrella_reset_val=65.0;
            umbrella_reset(&umbrella_reset_val);
            umbrella_size=0;
        }
   }
   else if (Expec2-r2*MAD2 > Expec - r2*MAD)// Consider the case where too much of the current umbrella overlaps
    {
        sucessful_compare=false;
        Shift_Flag=1;
        amrex::Print() << "New Phi_0 resulted in higher average"  << "\n";
        amrex::Print() << "decreasing Phi_0"  << "\n";
                    // if(ParallelDescriptor::MyProc() == 0 )
                    // {
                    //     ofs << "New Phi_0 resulted in higher average"  << "\n";
                    //     ofs << "decreasing Phi_0"  << "\n";
                    // }
        r2=-r2;
        inc_phi0_Adapt(&Expec,&MAD,&r2,&Shift_Flag);
        r2=-r2;
    }
    if(sucessful_iter_prev==1 and sucessful_iter==0 and !First_Loop_Step) // break the "sucess" while loop when a comparison fails(after many sucessful comparisons)
    {
        while_loop_comp=false;
    }
    if(sucessful_iter_prev==1 and sucessful_iter==1 and  umbrella_size==1)// break the "sucess" while loop after contigously sucessful comparisons with spring constants at lower limit
    {
        weak_umb=true;
    }
                // ofs.close();
}