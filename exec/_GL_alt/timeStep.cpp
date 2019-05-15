#include "GL_functions.H"
#include "GL_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "AMReX_ArrayLim.H"

void RK2step(MultiFab& phi, MultiFab& phin, MultiFab& rannums, 
               const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt, amrex::Real& integral, int n)
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

    for ( MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        rk2_stage1(BL_TO_FORTRAN_BOX(bx),
                   phi[mfi].dataPtr(),  
                   phin[mfi].dataPtr(),  
                   rannums[mfi].dataPtr(),
                   &integral,
                   &energy, &teng,
      	           dx, &dt);   
    }

    ParallelDescriptor::ReduceRealSum(energy);
    ParallelDescriptor::ReduceRealSum(teng);
    if(ParallelDescriptor::MyProc() == 0 ){
           std::cout << n << " " << energy << "  energy  " << std::endl;
           std::cout << n << " " << teng << "  teng  " << std::endl;
          // std::cout << time << " " << energy << "  energy = "<< energy << "  dphi/dt  = " << phit << std::endl;
    }

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


