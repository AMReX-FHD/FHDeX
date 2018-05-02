#include "gmres_functions.H"
#include "gmres_params_F.H"
#include "gmres_params.H"

using namespace gmres;

void initialize_gmres_namespace() {
    
    copy_gmres_params_to_c();

}
