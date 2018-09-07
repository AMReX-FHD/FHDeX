#include "common_functions.H"
#include "common_functions_F.H"


void FindNodalValues(const MultiFab& inFab, MultiFab& outFab, const MultiFab& ccFab)
{

    BL_PROFILE_VAR("FindNodalValues()",FindNodalValues);

    const int* hiVectIn;
    const int* hiVectCC;

    const int* loVectIn;
    const int* loVectCC;

    int xCheck, yCheck, zCheck;
   
    for (MFIter mfi(ccFab); mfi.isValid(); ++mfi) 
    {
        const Box& validBox = mfi.validbox();

        hiVectCC = validBox.hiVect();
        loVectCC = validBox.loVect();

        hiVectIn = inFab[mfi].hiVect();
        loVectIn = inFab[mfi].loVect();

//This can be simplified with macros..

#if (AMREX_SPACEDIM == 3)
        xCheck = (hiVectIn[0] - loVectIn[0]) - (hiVectCC[0] - loVectCC[0]) - 2;  //Minus 2 is for ghost cells
        yCheck = (hiVectIn[1] - loVectIn[1]) - (hiVectCC[1] - loVectCC[1]) - 2;  
        zCheck = (hiVectIn[2] - loVectIn[2]) - (hiVectCC[2] - loVectCC[2]) - 2;
#endif

#if (AMREX_SPACEDIM == 2)
        xCheck = (hiVectIn[0] - loVectIn[0]) - (hiVectCC[0] - loVectCC[0]) - 2;
        yCheck = (hiVectIn[1] - loVectIn[1]) - (hiVectCC[1] - loVectCC[1]) - 2;  
#endif
        find_nodal_values(loVectCC, hiVectCC, 
                          BL_TO_FORTRAN_3D(inFab[mfi]),
                          BL_TO_FORTRAN_3D(outFab[mfi]),
                          &xCheck, &yCheck, &zCheck);
    }

}
