#include "common_functions.H"
#include "multispec_functions.H"

// compute rhotot from rho in VALID REGION
void ComputeRhotot(const MultiFab& rho,
		   MultiFab& rhotot,
                   int include_ghost) // include_ghost=0 by default
{
  
    BL_PROFILE_VAR("ComputeRhotot()",ComputeRhotot);

    int ng = (include_ghost == 1) ? rhotot.nGrow() : 0;

    if (include_ghost == 1) {
        int ng_r = rho.nGrow();
        if (ng > ng_r) {
            Abort("ComputeRhotot: rho needs as many ghost cells as rhotot");
        }
    }
    
    rhotot.setVal(0.0);
    for (int i=0; i<nspecies; i++) {
        MultiFab::Add(rhotot,rho,i,0,1,ng);
    }

}

void ConvertRhoCToC(MultiFab& rho, const MultiFab& rhotot, MultiFab conc, int rho_to_c)
{
    Abort("ConvertRhoCToC: write me!");
    
    if (rho_to_c == 1) {
        // rho to conc - NO GHOST CELLS
        MultiFab::Copy(conc,rho,0,0,nspecies,0);
        for (int i=0; i<nspecies; ++i) {
            MultiFab::Divide(conc,rhotot,0,i,1,0);
        }

    }
    else if (rho_to_c == 0) {

        int ng = rho.nGrow();
        int ng_c = conc.nGrow();
        int ng_r = rhotot.nGrow();

        if (ng > ng_c || ng > ng_r) {
            Abort("ConvertRhoCToC: conc needs as many ghost cells as rho or rhotot");
        }
        
        // conc to rho - VALID + GHOST
        MultiFab::Copy(rho,conc,0,0,nspecies,ng);
        for (int i=0; i<nspecies; ++i) {
            MultiFab::Multiply(rho,rhotot,0,i,1,ng);
        }

    }
    else {
        Abort("ConvertRhoCToC: invalid rho_to_c");
    }
}

// void FillRhoRhototGhost


