#include "DsmcParticleContainer.H"
using namespace std;
void FhdParticleContainer::TimeCorrelation(MultiFab& mfcuInst,
						MultiFab& mfprimInst,
            MultiFab& mftimeCross,
            MultiFab& mft0Cross,
						const int steps)
{
    BL_PROFILE_VAR("EvaluateStats()",EvaluateStats);
    const Real osteps = 1.0/steps;
    const Real stepsMinusOne = steps-1.;

    const int lev = 0;
    if (plot_cross)
    {
        for ( MFIter mfi(mfcuInst); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();

            const Array4<const Real> prim      = mfprimInst.array(mfi);
            const Array4<const Real> cu        = mfcuInst.array(mfi);

						const Array4<      Real> t0Cross   = mftimeCross.array(mfi);
            const Array4<      Real> timeCross = mft0Cross.array(mfi);

						if(steps==1)
						{
							amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            	{
								t0Cross(i,j,k,0)  = cu(i,j,k,0);  // rho-instant
								t0Cross(i,j,k,1)  = cu(i,j,k,4);  // energy-instant
								t0Cross(i,j,k,2)  = cu(i,j,k,1);  // jx-instant
								t0Cross(i,j,k,3)  = cu(i,j,k,2);  // jy-instant
								t0Cross(i,j,k,4)  = cu(i,j,k,3);  // jz-instant
								t0Cross(i,j,k,5) = prim(i,j,k,2); // velx-instant
								t0Cross(i,j,k,6) = prim(i,j,k,3); // vely-instant
								t0Cross(i,j,k,7) = prim(i,j,k,4); // velz-instant
								t0Cross(i,j,k,8) = prim(i,j,k,6); // T-instant
							});
						}

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
            	timeCross(i,j,k,0) = (cu(i,j,k,0)*t0Cross(i,j,k,0)+timeCross(i,j,k,0)*stepsMinusOne)*osteps;
            	timeCross(i,j,k,1) = (cu(i,j,k,4)*t0Cross(i,j,k,1)+timeCross(i,j,k,1)*stepsMinusOne)*osteps;
            	timeCross(i,j,k,2) = (cu(i,j,k,1)*t0Cross(i,j,k,2)+timeCross(i,j,k,2)*stepsMinusOne)*osteps;
            	timeCross(i,j,k,3) = (cu(i,j,k,2)*t0Cross(i,j,k,3)+timeCross(i,j,k,3)*stepsMinusOne)*osteps;
            	timeCross(i,j,k,4) = (cu(i,j,k,3)*t0Cross(i,j,k,4)+timeCross(i,j,k,4)*stepsMinusOne)*osteps;
            	timeCross(i,j,k,5) = (prim(i,j,k,2)*t0Cross(i,j,k,5)+timeCross(i,j,k,5)*stepsMinusOne)*osteps;
            	timeCross(i,j,k,6) = (prim(i,j,k,3)*t0Cross(i,j,k,6)+timeCross(i,j,k,6)*stepsMinusOne)*osteps;
            	timeCross(i,j,k,7) = (prim(i,j,k,4)*t0Cross(i,j,k,7)+timeCross(i,j,k,7)*stepsMinusOne)*osteps;
            	timeCross(i,j,k,8) = (prim(i,j,k,6)*t0Cross(i,j,k,8)+timeCross(i,j,k,8)*stepsMinusOne)*osteps;
            });
        }
    }
}
