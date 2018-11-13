#include "common_functions.H"
#include "common_functions_F.H"
#include "rng_functions_F.H"
#include "common_namespace.H"

using namespace common;
using namespace amrex;

void MultiFABFillRandom(MultiFab& mf, const int& comp, const amrex::Real& variance, const Geometry& geom)
{

    BL_PROFILE_VAR("MultiFABFillRandom()",MultiFABFillRandom);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

        const Box& validBox = mfi.validbox();

	multifab_fill_random(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
			     BL_TO_FORTRAN_FAB(mf[mfi]), &comp);
    }
    
//----------------------------------------

//error in here?
    // Scale standard gaussian samples by standard deviation
    mf.mult(sqrt(variance), comp, 1, 0);

    // Enforce boundary conditions on nodal boundaries & ghost cells
    mf.OverrideSync(geom.periodicity());

    mf.FillBoundary(geom.periodicity());

    

//----------------------------------------
}

void MultiFABFillRandomHack(MultiFab& mf, const int& comp, const amrex::Real& variance, const Geometry& geom)
{

    BL_PROFILE_VAR("MultiFABFillRandom()",MultiFABFillRandom);

    double myTopNumber, myBottomNumber;

    int rightFriend = (ParallelDescriptor::MyProc()+1)%ParallelDescriptor::NProcs();
    int leftFriend = (ParallelDescriptor::MyProc()+ParallelDescriptor::NProcs()-1)%ParallelDescriptor::NProcs();

    if(ParallelDescriptor::MyProc()%2 == 0)
    {
        myTopNumber = get_fhd_normal_func();

        ParallelDescriptor::Send(&myTopNumber, 1, rightFriend, 1);

    }else
    {
        ParallelDescriptor::Recv(&myBottomNumber, 1, leftFriend, 1);
    }

    if(ParallelDescriptor::MyProc()%2 != 0)
    {
        myTopNumber = get_fhd_normal_func();
        ParallelDescriptor::Send(&myTopNumber, 1, rightFriend, 1);

    }else
    {
        ParallelDescriptor::Recv(&myBottomNumber, 1, leftFriend, 1);

    }

    if(ParallelDescriptor::MyProc() == 0)
    {
        myBottomNumber = get_fhd_normal_func();

    }

    if(ParallelDescriptor::MyProc() == ParallelDescriptor::NProcs()-1)
    {
        myTopNumber = get_fhd_normal_func();

    }

    //std::cout << "proc: " << ParallelDescriptor::MyProc() << ", top: " << myTopNumber << ", bottom: " << myBottomNumber << std::endl;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

        const Box& validBox = mfi.validbox();

	    multifab_fill_random_hack(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
			         BL_TO_FORTRAN_FAB(mf[mfi]), &comp, &myTopNumber, &myBottomNumber);
    }
    
//----------------------------------------

//error in here?
    // Scale standard gaussian samples by standard deviation
    mf.mult(sqrt(variance), comp, 1, 0);

    // Enforce boundary conditions on nodal boundaries & ghost cells
    //mf.OverrideSync(geom.periodicity());

    mf.FillBoundary(geom.periodicity());    

//----------------------------------------
}
