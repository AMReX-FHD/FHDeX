#include "rng_functions.H"
#include "common_functions.H"
#include "hydro_functions.H"

#include "TurbForcing.H"

// initialize n_rngs, geom
// build MultiFabs to hold random numbers
TurbForcing::TurbForcing(BoxArray ba_in, DistributionMapping dmap_in, Geometry geom_in)
{

    BL_PROFILE_VAR("TurbForcing()",TurbForcing);

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        sines  [i].define(convert(ba_in,nodal_flag_dir[i]), dmap_in, 22, 0);
        cosines[i].define(convert(ba_in,nodal_flag_dir[i]), dmap_in, 22, 0);
    }
    
    
}
