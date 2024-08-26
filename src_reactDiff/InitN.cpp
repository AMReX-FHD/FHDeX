#include "reactDiff_functions.H"

void InitN(MultiFab& n_in,
           const Geometry& geom,
           const Real& time) {

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(n_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.tilebox();

        const Array4<Real> & ninit = n_in.array(mfi);

        if (prob_type == 0) {
            //============================================================
            // Thermodynamic equilibrium
            //============================================================
        
            amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                ninit(i,j,k,n) = 0.;
            });

        }
    }
    
    n_in.FillBoundary(geom.periodicity());
    MultiFabPhysBC(n_in, geom, 0, nspecies, SPEC_BC_COMP, time);
}
