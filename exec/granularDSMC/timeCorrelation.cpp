#include "DsmcParticleContainer.H"
using namespace std;
void FhdParticleContainer::TimeCorrelation(
	const MultiFab& mfrho_time,
	const MultiFab& mfu_time,
	const MultiFab& mfK_time,
	MultiFab& mfrhotimeCross,
	MultiFab& mfutimeCross,
	MultiFab& mfKtimeCross,
	const int nCor,
	const int steps)
{
	BL_PROFILE_VAR("EvaluateStats()",EvaluateStats);
	Print() << "step " << steps << "\n";
	const Real osteps = 1.0/steps;
	const Real stepsMinusOne = steps-1.;
	const int lev = 0;
	for ( MFIter mfi(mfrhotimeCross); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.validbox();

		const Array4<const Real> rho_time  = mfrho_time.array(mfi);
		const Array4<const Real> u_time    = mfu_time.array(mfi);
		const Array4<const Real> K_time    = mfK_time.array(mfi);

		Array4<      Real> rhoCross  = mfrhotimeCross.array(mfi);
		Array4<      Real> uCross    = mfutimeCross.array(mfi);
		Array4<      Real> KCross    = mfKtimeCross.array(mfi);

		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		{
			for(int l=0;l<nCor;l++){
				rhoCross(i,j,k,l) = (rho_time(i,j,k,0)*rho_time(i,j,k,l)
					+rhoCross(i,j,k,l)*stepsMinusOne)*osteps; 
				uCross(i,j,k,l) = (u_time(i,j,k,0)*u_time(i,j,k,l)
					+uCross(i,j,k,l)*stepsMinusOne)*osteps;
				KCross(i,j,k,l) = (K_time(i,j,k,0)*K_time(i,j,k,l)
					+KCross(i,j,k,l)*stepsMinusOne)*osteps;
//				Print() << uCross(i,j,k,0) << " " << u_time(i,j,k,0) << "\n";
			}
		});
	}
}

// Update data points
// Assumes equally spaced
void FhdParticleContainer::updateTimeData(
	MultiFab& mfcuInst,
	MultiFab& mfprimInst,
	MultiFab& mfrho_time,
	MultiFab& mfu_time,
	MultiFab& mfK_time,
	const int nCor)
{
	BL_PROFILE_VAR("updateTimeData()",EvaluateStats);
	const int lev = 0;
	for ( MFIter mfi(mfcuInst); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.validbox();

		Array4<Real> rho_time              = mfrho_time.array(mfi);
		Array4<Real> u_time                = mfu_time.array(mfi);
		Array4<Real> K_time                = mfK_time.array(mfi);
		const Array4<const Real> cuInst    = mfcuInst.array(mfi);
		const Array4<const Real> primInst  = mfprimInst.array(mfi);

		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		{
			// shift all values one to left
			for (int l=1;l<nCor;l++) {
				rho_time(i,j,k,l-1) = rho_time(i,j,k,l);
				u_time(i,j,k,l-1) = u_time(i,j,k,l);
				K_time(i,j,k,l-1) = K_time(i,j,k,l);
			}
			rho_time(i,j,k,nCor-1) = cuInst(i,j,k,0);
			u_time(i,j,k,nCor-1)   = primInst(i,j,k,2);
			K_time(i,j,k,nCor-1)   = cuInst(i,j,k,4);
		});
	}
}
