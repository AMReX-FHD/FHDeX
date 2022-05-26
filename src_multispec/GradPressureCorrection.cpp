#include "multispec_functions.H"
#include "common_functions.H"

void GradPressureCorrection(const MultiFab& rho_in,
                            const MultiFab& rhotot_in,
	 	            const MultiFab& Temp,
			    MultiFab& lap_phi_beta_in,
			    std::array< MultiFab, AMREX_SPACEDIM >& grad_rhs,
		            const Geometry& geom,
			    int increment)
{
    BL_PROFILE_VAR("GradPressureCorrection()",GradPressureCorrection);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // make sure grad rhs is zero initially if not incrementing
    if (increment==0)
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir)
            grad_rhs[dir].setVal(0.);

    for(MFIter mfi(rho_in,TilingIfNotGPU()); mfi.isValid(); ++mfi){

        const Box& bx = mfi.validbox();
        const Array4<const Real>& rho = rho_in.array(mfi);
        const Array4<const Real>& rhotot = rhotot_in.array(mfi);
        const Array4<const Real>& T = Temp.array(mfi);
        const Array4<      Real>& lap_phi_beta = lap_phi_beta_in.array(mfi);

        AMREX_D_TERM(const Array4<Real> & grad_rhs_x = grad_rhs[0].array(mfi);,
                     const Array4<Real> & grad_rhs_y = grad_rhs[1].array(mfi);,
                     const Array4<Real> & grad_rhs_z = grad_rhs[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        // compute laplacian
        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n){

            lap_phi_beta(i,j,k,n) = ( (rho(i+1,j,k,n) -2*rho(i,j,k,n) + rho(i-1,j,k,n))/(dx[0]*dx[0])
                                    + (rho(i,j+1,k,n) -2*rho(i,j,k,n) + rho(i,j-1,k,n))/(dx[1]*dx[1])
	    		            ) / rhotot(i,j,k);
#if (AMREX_SPACEDIM == 3)
            lap_phi_beta(i,j,k,n) += (rho(i,j,k+1,n)-2*rho(i,j,k,n)+rho(i,j,k-1,n))
	                             /rhotot(i,j,k)/(dx[2]*dx[2]);
#endif
        });

	amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real phi_alpha;
            Real grad_phi_alpha_x;
	    Real grad_lap_phi_beta_x;
	    Real grad_summation = 0.;

	    // temporarily fill grad_rhs with summation values
            for (int alpha = 0; alpha < nspecies; ++alpha)
	    {
                for (int beta = 0; beta < nspecies; ++beta)
		{
		        phi_alpha = rho(i,j,k,alpha) / rhotot(i,j,k);
			grad_phi_alpha_x = ( rho(i,j,k,alpha) - rho(i-1,j,k,alpha) ) / dx[0] 
			                   / rhotot(i,j,k);
			grad_lap_phi_beta_x = (lap_phi_beta(i,j,k,beta)-lap_phi_beta(i-1,j,k,beta))/dx[0];
			grad_summation += fh_kappa(alpha,beta) * (
                                          phi_alpha*grad_lap_phi_beta_x +
                                          grad_phi_alpha_x*lap_phi_beta(i,j,k,beta) );
                }
            }
	    // compute correction term
	    grad_rhs_x(i,j,k) += rho0*k_B*T(i,j,k) / monomer_mass * grad_summation;
        });
	amrex::ParallelFor(bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real phi_alpha;
            Real grad_phi_alpha_y;
	    Real grad_lap_phi_beta_y;
	    Real grad_summation = 0.;

	    // temporarily fill grad_rhs with summation values
            for (int alpha = 0; alpha < nspecies; ++alpha)
	    {
                for (int beta = 0; beta < nspecies; ++beta)
		{
		        phi_alpha = rho(i,j,k,alpha) / rhotot(i,j,k);
			grad_phi_alpha_y = ( rho(i,j,k,alpha) - rho(i,j-1,k,alpha) ) / dx[1] 
			                   / rhotot(i,j,k);
			grad_lap_phi_beta_y = (lap_phi_beta(i,j,k,beta)-lap_phi_beta(i,j-1,k,beta))/dx[1];
			grad_summation += fh_kappa(alpha,beta) * (
                                          phi_alpha*grad_lap_phi_beta_y +
                                          grad_phi_alpha_y*lap_phi_beta(i,j,k,beta) );
                }
            }
	    // compute correction term
	    grad_rhs_y(i,j,k) += rho0*k_B*T(i,j,k) / monomer_mass * grad_summation;
        });
#if (AMREX_SPACEDIM == 3)
	amrex::ParallelFor(bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real phi_alpha;
            Real grad_phi_alpha;
	    Real grad_lap_phi_beta_z;
	    Real grad_summation = 0.;

	    // temporarily fill grad_rhs with summation values
            for (int alpha = 0; alpha < nspecies; ++alpha)
	    {
                for (int beta = 0; beta < nspecies; ++beta)
		{
		        phi_alpha = rho(i,j,k,alpha) / rhotot(i,j,k);
			grad_phi_alpha_z = ( rho(i,j,k,alpha) - rho(i,j,k-1,alpha) ) / dx[2] 
			                   / rhotot(i,j,k);
			grad_lap_phi_beta_z = (lap_phi_beta(i,j,k,beta)-lap_phi_beta(i,j,k-1,beta))/dx[2];
			grad_summation += fh_kappa(alpha,beta) * (
                                          phi_alpha*grad_lap_phi_beta_y +
                                          grad_phi_alpha_z*lap_phi_beta(i,j,k,beta) );
                }
            }
	    // compute correction term
	    grad_rhs_z(i,j,k) += rho0*k_B*T(i,j,k) / monomer_mass * grad_summation;
        });
#endif


    }

} 
