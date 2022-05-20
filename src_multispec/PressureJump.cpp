#include "multispec_functions.H"
#include "common_functions.H"

void PressureJump(const MultiFab& press,
                  const MultiFab& rho_in,
                  const MultiFab& rhotot_in,
	 	  const MultiFab& Temp,
                  MultiFab& pressure_jump_in,
		  const Geometry& geom)
{
    BL_PROFILE_VAR("PressureJump()",PressureJump);
    for(MFIter mfi(pressure_jump_in,TilingIfNotGPU()); mfi.isValid(); ++mfi){

        const Box& bx = mfi.validbox();
	const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
        const Array4<const Real>& pi = press.array(mfi);
        const Array4<const Real>& rho = rho_in.array(mfi);
        const Array4<const Real>& rhotot = rhotot_in.array(mfi);
        const Array4<const Real>& T = Temp.array(mfi);
        const Array4<      Real>& pressure_jump = pressure_jump_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
 
            Real phi_alpha, lap_phi_beta;
	    Real summation = 0.;
            for (int alpha = 0; alpha < nspecies; ++alpha){
                for (int beta = 0; beta < nspecies; ++beta){
		        phi_alpha = rho(i,j,k,alpha) / rhotot(i,j,k);
                	lap_phi_beta = ( (rho(i+1,j,k,beta) -2*rho(i,j,k,beta) + rho(i-1,j,k,beta))/(dx[0]*dx[0])
            		               + (rho(i,j+1,k,beta) -2*rho(i,j,k,beta) + rho(i,j-1,k,beta))/(dx[1]*dx[1])
				       ) / rhotot(i,j,k);
            	        summation += phi_alpha * fh_kappa(alpha,beta) * lap_phi_beta;
                }
            }
            pressure_jump(i,j,k) = pi(i,j,k) - rho0*k_B*T(i,j,k) / monomer_mass * summation;

        });
    }
} 
