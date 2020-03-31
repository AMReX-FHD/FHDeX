#include "common_functions.H"
#include "multispec_functions.H"

// compute rhotot from rho in VALID REGION
void ComputeRhotot(const MultiFab& rho,
		   MultiFab& rhotot)
{
  
    BL_PROFILE_VAR("ComputeRhotot()",ComputeRhotot);

    rhotot.setVal(0.0);
    for (int i=0; i<nspecies; i++) {
        MultiFab::Add(rhotot,rho,i,0,1,0);
    }

}

void ConvertRhoCToC(MultiFab& rho, const MultiFab& rhotot, MultiFab conc, int rho_to_c)
{
    Abort("ConvertRhoCToC: write me!");
    
    if (rho_to_c == 1) {
        // rho to conc - NO GHOST CELLS

    }
    else if (rho_to_c == 0) {
        // conc to rho - VALID + GHOST

    }
    else {
        Abort("ConvertRhoCToC: invalid rho_to_c");
    }
}

// void FillRhoRhototGhost


