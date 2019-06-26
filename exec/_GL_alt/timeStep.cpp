
#include "GL_functions.H"
#include "GL_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "exec_functions.H"

#include "AMReX_ArrayLim.H"

#include <algorithm>

void RK2step(MultiFab& phi, MultiFab& phin, MultiFab& rannums, 
               const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt, amrex::Real& integral, int n, amrex::Real& phi_avg, amrex::Real& phi_sq_avg)
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
    phi_sq_avg=0.0;
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        rk2_stage1(BL_TO_FORTRAN_BOX(bx),
                   phi[mfi].dataPtr(),  
                   phin[mfi].dataPtr(),  
                   rannums[mfi].dataPtr(),
                   &integral,
                   &energy, &teng,
      	           dx, &dt,&phi_avg,&phi_sq_avg);   
    }
    ParallelDescriptor::ReduceRealSum(phi_avg);
    ParallelDescriptor::ReduceRealSum(phi_sq_avg);
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
                        amrex::Real& time, int plot_int,int n_steps_skip, bool Make_PltFiles, amrex::Real& alpha,amrex::Real& r1, amrex::Real& r2,int N_Burn,int L, amrex::Real& Expec, amrex::Real& MAD, int max_step)
{
    int Steps_Counted;
    int Total_Steps;
    int j=0;
    amrex::Real Expec_Temp;
    amrex::Real Phi_Avg;
    amrex::Real Phi_sqrd_Avg;
    amrex::Real median;
    double Avg_collect[L]={0};


    Total_Steps=L*plot_int;
    Expec_Temp=0.0;
    Phi_Avg=0.0;
    Phi_Avg=0.0;

    for(step=1;step<=Total_Steps+N_Burn;++step)
    {

        RK2step(phi, phin, rannums, geom, dx, dt, integral, step,Phi_Avg,Phi_sqrd_Avg);

//        if(step == n_steps_skip)
//     {
//     on the fly statistics setup?
//      }

//        evaluateStats(phi, dx);

//        statsCount++;

        //inc_phi0(&step);
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
                        std::ofstream ofs;
                        ofs.open ("AVG_DATA.txt", std::ofstream::out | std::ofstream::app);
                         ofs << Phi_Avg<< "\n";
                        ofs.close();
                        //WritePlotFile(step, time, geom, phi);
                    }
            }
        }
        time = time + dt;
    }
    Expec=Expec_Temp/L;
    std::sort(Avg_collect, Avg_collect+L);
    median=Avg_collect[L/2];
    amrex::Print() << " median of Phi in umbrella  is " << median << "\n";
    for(j=0;j<L;j++)
    {
        Avg_collect[j]=std::abs(Avg_collect[j]-median) ;
    }
    std::sort(Avg_collect, Avg_collect+L);
    MAD=Avg_collect[L/2];
    amrex::Print() << " Average of Phi in umbrella  is " << Expec << "\n";
    amrex::Print() << "MAD of Umbrella is " << MAD << "\n";

}

void Check_Overlap(amrex::Real& Expec,amrex::Real& MAD,amrex::Real& Expec2,amrex::Real& MAD2,amrex::Real& r2,amrex::Real& alpha, bool& sucessful_compare, int& umbrella_size)
{   
    int sucessful_iter;
    int sucessful_iter_prev;

    if(sucessful_compare)
    {
    sucessful_iter_prev=1;
    }else
    {
    sucessful_iter_prev=0;
    }
    


    if(Expec2-r2*MAD2 < Expec + r2*MAD && Expec2+r2*MAD2 > Expec + r2*MAD)
    {    
        sucessful_compare=true;
        sucessful_iter=1;
        Umbrella_Adjust(&sucessful_iter,&alpha,&umbrella_size,&sucessful_iter_prev);
        amrex::Print() << "Overlap occured "  << "\n";

    }else{
        sucessful_compare=false;
        sucessful_iter=0;
        Umbrella_Adjust (&sucessful_iter,&alpha ,&umbrella_size,&sucessful_iter_prev);
        amrex::Print() << "Overlap did not occur"  << "\n";

   }
}