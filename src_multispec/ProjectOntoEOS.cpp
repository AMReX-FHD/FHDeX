#include "multispec_functions.H"

void ProjectOntoEOS(MultiFab& rho)
{
    if (algorithm_type == 4 || algorithm_type == 6) {

        Abort("ProjectOntoEOS algorithm_type = 4, 6 not written");
        
    } else {

        Abort("ProjectOntoEOS algorithm_type = 0, 2, 3, 5 not written");
        
    }
    
}
