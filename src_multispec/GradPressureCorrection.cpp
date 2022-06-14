#include "multispec_functions.H"
#include "common_functions.H"

void GradPressureCorrection(const MultiFab& rho_in,
                            const MultiFab& rhotot_in,
	 	            const MultiFab& Temp,
			    std::array< MultiFab, AMREX_SPACEDIM >& grad_rhs,
		            const Geometry& geom,
			    int increment)
{
    BL_PROFILE_VAR("GradPressureCorrection()",GradPressureCorrection);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    int gibbs_duhem = 1-increment;

    // create temp multifab for pressure correction
    BoxArray ba = rho_in.boxArray();
    DistributionMapping dmap = rho_in.DistributionMap();
    MultiFab correction_mf(ba, dmap, 1, 1);
    MultiFab chemicalpotential_mf(ba, dmap, nspecies, 1);

    // make sure grad rhs is zero initially if not incrementing
    if (increment==0)
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir)
            grad_rhs[dir].setVal(0.);

    for(MFIter mfi(rho_in,TilingIfNotGPU()); mfi.isValid(); ++mfi){

        const Box& bx = mfi.growntilebox(1);
        const Array4<const Real>& rho = rho_in.array(mfi);
        const Array4<const Real>& rhotot = rhotot_in.array(mfi);
        const Array4<const Real>& T = Temp.array(mfi);
        const Array4<      Real>& correction = correction_mf.array(mfi);
        const Array4<      Real>& mu = chemicalpotential_mf.array(mfi);

        // compute correction term
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){

            Real phi_alpha, phi_beta, lap_phi_beta;
            Real summation = 0.;
            for (int alpha = 0; alpha < nspecies; ++alpha){
                phi_alpha = rho(i,j,k,alpha) / rhotot(i,j,k);
		
		if (gibbs_duhem)
		    mu(i,j,k,alpha) = ( log(phi_alpha) + 1 ) / fh_monomers[alpha];

                for (int beta = 0; beta < nspecies; ++beta){
                   phi_beta = rho(i,j,k,beta) / rhotot(i,j,k);
                   lap_phi_beta = ( (rho(i+1,j,k,beta) -2*rho(i,j,k,beta) + rho(i-1,j,k,beta))/(dx[0]*dx[0])
                                  + (rho(i,j+1,k,beta) -2*rho(i,j,k,beta) + rho(i,j-1,k,beta))/(dx[1]*dx[1])
                                  ) / rhotot(i,j,k);
#if (AMREX_SPACEDIM == 3)
                   lap_phi_beta += (rho(i,j,k+1,beta)-2*rho(i,j,k,beta)+rho(i,j,k-1,beta))
                                   /rhotot(i,j,k)/(dx[2]*dx[2]);
#endif

                   if (gibbs_duhem)
		       mu(i,j,k,alpha) += fh_chi(alpha,beta)*phi_beta + fh_kappa(alpha,beta)*lap_phi_beta;

                   else
                       summation += phi_alpha * fh_kappa(alpha,beta) * lap_phi_beta;
                }
		if (gibbs_duhem)
		    mu(i,j,k,alpha) *= k_B*T(i,j,k) / monomer_mass;
            }
            if (increment)
		correction(i,j,k) = rho0*k_B*T(i,j,k) / monomer_mass * summation;

        });

    }


    for(MFIter mfi(rho_in,TilingIfNotGPU()); mfi.isValid(); ++mfi){
        AMREX_D_TERM(const Array4<Real> & grad_rhs_x = grad_rhs[0].array(mfi);,
                     const Array4<Real> & grad_rhs_y = grad_rhs[1].array(mfi);,
                     const Array4<Real> & grad_rhs_z = grad_rhs[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        const Array4<Real>& correction = correction_mf.array(mfi);
        const Array4<const Real>& rho = rho_in.array(mfi);
        const Array4<      Real>& mu = chemicalpotential_mf.array(mfi);

	// take gradient
	amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
	    if (gibbs_duhem){
	        for (int n = 0; n < nspecies; n++)
		    grad_rhs_x(i,j,k) += (rho(i,j,k,n)+rho(i-1,j,k,n))/2.
		                         * (mu(i,j,k,n)-mu(i-1,j,k,n))/dx[0];
            }
	    else{
	        grad_rhs_x(i,j,k) += (correction(i,j,k) - correction(i-1,j,k))/dx[0];
	    }
        });
	amrex::ParallelFor(bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
	    if (gibbs_duhem){
	        for (int n = 0; n < nspecies; n++)
		    grad_rhs_y(i,j,k) += (rho(i,j,k,n)+rho(i,j-1,k,n))/2.
		                         * (mu(i,j,k,n)-mu(i,j-1,k,n))/dx[1];
            }
	    else{
	        grad_rhs_y(i,j,k) += (correction(i,j,k) - correction(i,j-1,k))/dx[1];
	    }
        });
#if (AMREX_SPACEDIM == 3)
	amrex::ParallelFor(bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
	    if (gibbs_duhem){
	        for (int n = 0; n < nspecies; n++)
		    grad_rhs_z(i,j,k) += (rho(i,j,k,n)+rho(i,j,k-1,n))/2.
		                         * (mu(i,j,k,n)-mu(i,j,k-1,n))/dx[2];
            }
	    else{
	        grad_rhs_z(i,j,k) += (correction(i,j,k) - correction(i,j,k-1))/dx[2];
	    }
        });
#endif


    }

} 
