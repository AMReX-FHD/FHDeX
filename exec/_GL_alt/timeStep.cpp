
#include "GL_functions.H"
#include "GL_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "exec_functions.H"

#include "AMReX_ArrayLim.H"

#include <algorithm>

void RK2step(MultiFab& phi, MultiFab& phin, MultiFab& rannums, 
               const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt, amrex::Real& integral, int n, amrex::Real& phi_avg)
{

    const int comp=0;
    int ncomp=1;

    amrex::Real energy, teng;

    for (MFIter mfi(rannums); mfi.isValid(); ++mfi) {

        const Box& validBox = mfi.validbox();

        //multifab_fill_random(ARLIM_2D(validBox.loVect()), ARLIM_2D(validBox.hiVect()),
        multifab_fill_random(BL_TO_FORTRAN_BOX(validBox),
                             BL_TO_FORTRAN_ANYD(rannums[mfi]), &ncomp, &comp);
    }

    phi.FillBoundary(geom.periodicity());

    integral = 0.;

    for ( MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        integrate(BL_TO_FORTRAN_BOX(bx),
                   phi[mfi].dataPtr(),  
      	           dx,
                   &integral);   

    }

    ParallelDescriptor::ReduceRealSum(integral);

//    if(ParallelDescriptor::MyProc() == 0 ){
//           std::cout << "integral = "<< integral << std::endl;
//    }

    energy =0.;
    teng =0.;
    phi_avg=0.0;
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        rk2_stage1(BL_TO_FORTRAN_BOX(bx),
                   phi[mfi].dataPtr(),  
                   phin[mfi].dataPtr(),  
                   rannums[mfi].dataPtr(),
                   &integral,
                   &energy, &teng,
      	           dx, &dt,&phi_avg);   
    }
    ParallelDescriptor::ReduceRealSum(phi_avg);
    ParallelDescriptor::ReduceRealSum(energy);
    ParallelDescriptor::ReduceRealSum(teng);
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
//      	           ZFILL(dx), &dt);
//    }
}

void Init_Phi(MultiFab& phi, const amrex::Real* dx )
  {


    for ( MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        initphi(BL_TO_FORTRAN_BOX(bx),
                   phi[mfi].dataPtr(),
                   dx);
    }

   }

void Run_Steps(MultiFab& phi, MultiFab& phin, MultiFab& rannums, const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt, amrex::Real& integral,int step,
                        amrex::Real& time, int plot_int,int n_steps_skip, bool Make_PltFiles,int N_Burn,int L, amrex::Real& Expec, amrex::Real& MAD, int max_step, int& Plot_Num,int& Plot_Skip, int& umbrella_number)
{
    int Steps_Counted;
    int Total_Steps;
    int j=0;
    amrex::Real Expec_Temp;
    amrex::Real Phi_Avg;
    amrex::Real median;
    amrex::Real umbrella;
    amrex::Real phi0;
    double Avg_collect[L]={0};


    Total_Steps=L*plot_int;
    Expec_Temp=0.0;
    Phi_Avg=0.0;
    Param_Output(&umbrella,&phi0);
    std::stringstream ss;
    ss << std::setw(8) << std::setfill('0') << umbrella_number;
    std::string s = ss.str();
    std::string umb_Num;
    umb_Num="umbrella"+ s;
    umb_Num=umb_Num+".txt";
    std::ofstream ofs;
    std::ofstream ofs2;
 if(ParallelDescriptor::MyProc() == 0  and Make_PltFiles)                    
 {
    ofs.open (umb_Num.c_str(), std::ofstream::out | std::ofstream::app);
    ofs << umbrella<< "\n";
    ofs << phi0<< "\n";
    ofs.close();
 }






    for(step=1;step<=Total_Steps+N_Burn;++step)
    {

        RK2step(phi, phin, rannums, geom, dx, dt, integral, step,Phi_Avg);

        if(step>N_Burn)
        {
            Steps_Counted=step-N_Burn;
            if(step%500 == 0)
            {
                    //amrex::Print() << "Integral " << integral << "\n";
                    //amrex::Print() << "Advanced step " << step << "\n";
            }
            if (plot_int > 0 && Steps_Counted > n_steps_skip && Steps_Counted%plot_int == 0)
            {
                
                Avg_collect[j]=Phi_Avg;
                Expec_Temp=Expec_Temp+Phi_Avg;
                j=j+1;
                amrex::Print() << "The average of phi is " << Phi_Avg<< "\n";

                if(Make_PltFiles)
                    {
                        if(ParallelDescriptor::MyProc() == 0 )
                        {
                            ofs2.open("AVG_DATA.txt", std::ofstream::out | std::ofstream::app);
                            ofs2 << Phi_Avg<< "\n";
                            ofs2.close();

                            ofs.open (umb_Num.c_str(), std::ofstream::out | std::ofstream::app);
                            ofs << Phi_Avg<< "\n";
                            ofs.close();
                        }
                        if((j+1)%Plot_Skip==0)
                        {
                            Plot_Num+=1;
                            WritePlotFile(Plot_Num, time, geom, phi, umbrella,phi0);
                        }
                    }
                
            }
        }
        time = time + dt;
    }
    Expec=Expec_Temp/L;
    std::sort(Avg_collect, Avg_collect+L);
    median=Avg_collect[L/2];
    amrex::Print() << "median of Phi in umbrella is " << median << "\n";
    for(j=0;j<L;j++)
    {
        Avg_collect[j]=std::abs(Avg_collect[j]-median) ;
    }
    std::sort(Avg_collect, Avg_collect+L);
    MAD=Avg_collect[L/2];
    amrex::Print() << "Average of Phi in umbrella is " << Expec << "\n";
    amrex::Print() << "MAD of Umbrella is " << MAD << "\n";
    //Expec=median; // using median instead of average for comparison
    if(Make_PltFiles)
    {
        umbrella_number=umbrella_number+1;
    }

}

void Check_Overlap(amrex::Real& Expec,amrex::Real& MAD,amrex::Real& Expec2,amrex::Real& MAD2,amrex::Real& r2,amrex::Real& alpha, bool& sucessful_compare, int& umbrella_size, int& Shift_Flag, bool& while_loop_comp, bool& First_Loop_Step, bool& weak_phi)
{ 
      
    int sucessful_iter;
    int sucessful_iter_prev;
    int Umbrella_Size_Prev;
    amrex::Real r_temp;
    amrex::Real umbrella_reset_val;

    Umbrella_Size_Prev=umbrella_size;

    if(sucessful_compare)
    {
    sucessful_iter_prev=1;
    }else
    {
    sucessful_iter_prev=0;
    }   
    amrex::Print() << "E2-r2*S2 "  << Expec2-r2*MAD2 << "\n";
    amrex::Print() << "E2+r2*S2 "  << Expec2+r2*MAD2 << "\n";
    amrex::Print() << "E1+r2*S1 "  << Expec + r2*MAD <<  "\n";

    if(Expec2-r2*MAD2 < Expec + r2*MAD && Expec2+r2*MAD2 > Expec + r2*MAD)
    {    

        sucessful_compare=true;
        sucessful_iter=1;
        Shift_Flag=0;
        Umbrella_Adjust(&sucessful_iter,&alpha,&umbrella_size,&sucessful_iter_prev);
        amrex::Print() << "Overlap occured "  << "\n";
        if(umbrella_size==1 &&  Umbrella_Size_Prev==1)
        {
        //    sucessful_compare=false;
        //    sucessful_iter=0;
        //    amrex::Print() << "Next umbrella is too close with weak k, shifting phi0 up "  << "\n";
        //    Shift_Flag=1;
        //    r_temp=0.5*r2;
        //    inc_phi0_Adapt(&Expec,&MAD,&r_temp,&Shift_Flag);
            weak_phi=true;
        //    umbrella_reset_val=150.0;
        //    umbrella_reset(&umbrella_reset_val);
        //    umbrella_size=0;
        }
 

    }else if (Expec2-r2*MAD2 > Expec + r2*MAD)
    {   
        sucessful_compare=false;
        Shift_Flag=0;
        sucessful_iter=0;
        Umbrella_Adjust (&sucessful_iter,&alpha ,&umbrella_size,&sucessful_iter_prev);
        amrex::Print() << "Overlap did not occur"  << "\n";
        amrex::Print() << umbrella_size  << "\n";
        if(umbrella_size==2)
        {   
            r_temp=r2;
            r2=-0.5*r2;
            inc_phi0_Adapt(&Expec,&MAD,&r2,&Shift_Flag);
            r2=r_temp;  
            umbrella_reset_val=150.0;
            umbrella_reset(&umbrella_reset_val);
            umbrella_size=0;
        } 
   }
   else if (Expec2+r2*MAD2 < Expec + r2*MAD)
    {
        sucessful_compare=false;
        Shift_Flag=1;
        amrex::Print() << "New Phi_0 resulted in lower average"  << "\n";
        amrex::Print() << "Increasing Phi_0"  << "\n";
        inc_phi0_Adapt(&Expec,&MAD,&r2,&Shift_Flag);
    }
    if(sucessful_iter_prev==1 and sucessful_iter==0 and !First_Loop_Step) 
    {
        while_loop_comp=false;
    }
    if(sucessful_iter_prev==1 and sucessful_iter==1 and  umbrella_size==1)
    {
        weak_phi=true;
    }
}





void Check_Overlap_Backwards(amrex::Real& Expec,amrex::Real& MAD,amrex::Real& Expec2,amrex::Real& MAD2,amrex::Real& r2,amrex::Real& alpha, bool& sucessful_compare, int& umbrella_size, int& Shift_Flag, bool& while_loop_comp, bool& First_Loop_Step, bool& weak_phi)
{ 
      
    int sucessful_iter;
    int sucessful_iter_prev;
    int Umbrella_Size_Prev;
    amrex::Real r_temp;
    amrex::Real umbrella_reset_val;

    Umbrella_Size_Prev=umbrella_size;

    if(sucessful_compare)
    {
    sucessful_iter_prev=1;
    }else
    {
    sucessful_iter_prev=0;
    }   
    amrex::Print() << "E2-r2*S2 "  << Expec2-r2*MAD2 << "\n";
    amrex::Print() << "E2+r2*S2 "  << Expec2+r2*MAD2 << "\n";
    amrex::Print() << "E1-r2*S1 "  << Expec - r2*MAD <<  "\n";

    if(Expec-r2*MAD < Expec2 + r2*MAD2 && Expec2-r2*MAD2 < Expec - r2*MAD)
    {    

        sucessful_compare=true;
        sucessful_iter=1;
        Shift_Flag=0;
        Umbrella_Adjust(&sucessful_iter,&alpha,&umbrella_size,&sucessful_iter_prev);
        amrex::Print() << "Overlap occured "  << "\n";
        if(umbrella_size==1 &&  Umbrella_Size_Prev==1)
        {
        //    sucessful_compare=false;
        //    sucessful_iter=0;
        //    amrex::Print() << "Next umbrella is too close with weak k, shifting phi0 up "  << "\n";
        //    Shift_Flag=1;
        //    r_temp=0.5*r2;
        //    inc_phi0_Adapt(&Expec,&MAD,&r_temp,&Shift_Flag);
            weak_phi=true;
        //    umbrella_reset_val=150.0;
        //    umbrella_reset(&umbrella_reset_val);
        //    umbrella_size=0;
        }
 

    }else if (Expec-r2*MAD > Expec2 + r2*MAD2)
    {   
        sucessful_compare=false;
        Shift_Flag=0;
        sucessful_iter=0;
        Umbrella_Adjust (&sucessful_iter,&alpha ,&umbrella_size,&sucessful_iter_prev);
        amrex::Print() << "Overlap did not occur"  << "\n";
        amrex::Print() << umbrella_size  << "\n";
        if(umbrella_size==2)
        {   
            r_temp=r2;
            r2=0.5*r2;
            inc_phi0_Adapt(&Expec,&MAD,&r2,&Shift_Flag);
            r2=r_temp;  
            umbrella_reset_val=150.0;
            umbrella_reset(&umbrella_reset_val);
            umbrella_size=0;
        } 
   }
   else if (Expec2-r2*MAD2 > Expec - r2*MAD)
    {
        sucessful_compare=false;
        Shift_Flag=1;
        amrex::Print() << "New Phi_0 resulted in higher average"  << "\n";
        amrex::Print() << "decreasing Phi_0"  << "\n";
        r2=-r2;
        inc_phi0_Adapt(&Expec,&MAD,&r2,&Shift_Flag);
        r2=-r2;
    }
    if(sucessful_iter_prev==1 and sucessful_iter==0 and !First_Loop_Step) 
    {
        while_loop_comp=false;
    }
    if(sucessful_iter_prev==1 and sucessful_iter==1 and  umbrella_size==1)
    {
        weak_phi=true;
    }
}