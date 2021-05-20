#include "INS_functions.H"
#include "common_functions.H"

#include "gmres_functions.H"


//Explicit Euler
void eulerStep(const MultiFab& betaCC, const MultiFab& gammaCC,
                 std::array<MultiFab, NUM_EDGE>& betaEdge,
                 std::array<MultiFab, AMREX_SPACEDIM>& umacIn, 
                 std::array<MultiFab, AMREX_SPACEDIM>& umacOut,
                 std::array<MultiFab, AMREX_SPACEDIM>& umacNew,
                 std::array<MultiFab, AMREX_SPACEDIM>& alpha,
                 const Geometry geom,
                 Real* dt)
{

    StagApplyOp(geom,betaCC, gammaCC,
                betaEdge,
                umacIn, umacOut, alpha, geom.CellSize(), 1.);

    const int xOff[3] = {1,0,0};
    const int yOff[3] = {0,1,0};
    const int zOff[3] = {0,0,1};

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(betaCC,TilingIfNotGPU()); mfi.isValid(); ++mfi) 
    {
        const Box& validBox = mfi.validbox();
        
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););

        AMREX_D_TERM(const Array4<Real> & oldx = umacIn[0].array(mfi);,
                     const Array4<Real> & oldy = umacIn[1].array(mfi);,
                     const Array4<Real> & oldz = umacIn[2].array(mfi););

        AMREX_D_TERM(const Array4<Real> & stagopx = umacOut[0].array(mfi);,
                     const Array4<Real> & stagopy = umacOut[1].array(mfi);,
                     const Array4<Real> & stagopz = umacOut[2].array(mfi););

        AMREX_D_TERM(const Array4<Real> & newdatax = umacNew[0].array(mfi);,
                     const Array4<Real> & newdatay = umacNew[1].array(mfi);,
                     const Array4<Real> & newdataz = umacNew[2].array(mfi););

#if (AMREX_SPACEDIM == 2)

        amrex::ParallelFor(bx_x,bx_y,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            newdatax(i,j,k) = oldx(i,j,k)-stagopx(i,j,k)*dt;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            newdatay(i,j,k) = oldy(i,j,k)-stagopy(i,j,k)*dt;
        });

#elif (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(bx_x,bx_y,bx_z,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            newdatax(i,j,k) = oldx(i,j,k)-stagopx(i,j,k)*dt;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            newdatay(i,j,k) = oldy(i,j,k)-stagopy(i,j,k)*dt;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            newdataz(i,j,k) = oldz(i,j,k)-stagopz(i,j,k)*dt;        
        });
#endif
    }
}

