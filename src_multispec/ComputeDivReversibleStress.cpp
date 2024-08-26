#include "multispec_functions.H"

void ComputeDivFHReversibleStress(std::array<MultiFab,AMREX_SPACEDIM>& div_reversible_stress,
                                const MultiFab& rhotot_in,
//                                Array2D<Real, 0, MAX_SPECIES-1, 0, MAX_SPECIES-1>& fh_kappa,
                                MultiFab& rho_in,
                                const Geometry& geom)
{
    BL_PROFILE_VAR("ComputeDivFHReversibleStress()",ComputeDivReversibleStress);
    
    BoxArray ba = rho_in.boxArray();
    DistributionMapping dmap = rho_in.DistributionMap();
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    MultiFab node_grad_x_mf(convert(ba,nodal_flag), dmap, nspecies, 1);
    MultiFab node_grad_y_mf(convert(ba,nodal_flag), dmap, nspecies, 1);
    MultiFab node_grad_z_mf(convert(ba,nodal_flag), dmap, nspecies, 1);

    MultiFab conc(ba, dmap, nspecies, 2);

    // rho to conc - VALID REGION ONLY
    ConvertRhoCToC(rho_in,rhotot_in,conc,1);

    Real scale_factor = rhobar[0]*k_B*T_init[0]/monomer_mass;

    // fill conc ghost cells
    conc.FillBoundary(geom.periodicity());
    MultiFabPhysBCFH(conc,geom,0,nspecies,scale_factor);    

//  JBB  check what's in this
//    Real scale_factor = rhobar[0]*k_B*T_init[0]/molmass[0];
    
    for ( MFIter mfi(node_grad_x_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.growntilebox(1);

        const Array4<Real>& node_grad_x = node_grad_x_mf.array(mfi);
        const Array4<Real>& node_grad_y = node_grad_y_mf.array(mfi);
#if (AMREX_SPACEDIM == 3)            
        const Array4<Real>& node_grad_z = node_grad_z_mf.array(mfi);
#endif
        const Array4<Real const>& c = conc.array(mfi);
        
        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
#if (AMREX_SPACEDIM == 2)            

            node_grad_x(i,j,k,n) = (c(i,j,k,n)-c(i-1,j,k,n)+c(i,j-1,k,n)-c(i-1,j-1,k,n))/(2*dx[0]);
            node_grad_y(i,j,k,n) = (c(i,j,k,n)-c(i,j-1,k,n)+c(i-1,j,k,n)-c(i-1,j-1,k,n))/(2*dx[1]);

#elif (AMREX_SPACEDIM == 3)

            node_grad_x(i,j,k,n) = (c(i,j,k  ,n)-c(i-1,j,k  ,n)+c(i,j-1,k  ,n)-c(i-1,j-1,k  ,n)
                                   +c(i,j,k-1,n)-c(i-1,j,k-1,n)+c(i,j-1,k-1,n)-c(i-1,j-1,k-1,n))/(4*dx[0]);
            node_grad_y(i,j,k,n) = (c(i,j,k  ,n)-c(i,j-1,k  ,n)+c(i-1,j,k  ,n)-c(i-1,j-1,k  ,n)
                                   +c(i,j,k-1,n)-c(i,j-1,k-1,n)+c(i-1,j,k-1,n)-c(i-1,j-1,k-1,n))/(4*dx[1]);
            node_grad_z(i,j,k,n) = (c(i,j  ,k,n)-c(i,j  ,k-1,n)+c(i-1,j  ,k,n)-c(i-1,j  ,k-1,n)
                                   +c(i,j-1,k,n)-c(i,j-1,k-1,n)+c(i-1,j-1,k,n)-c(i-1,j-1,k-1,n))/(4*dx[2]);
#endif
        });
    }

// JBB is this needed?
    for (int n = 0 ; n < AMREX_SPACEDIM; n++)
        div_reversible_stress[n].setVal(0.);

    for ( MFIter mfi(rhotot_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        AMREX_D_TERM(const Array4<Real>& node_grad_x = node_grad_x_mf.array(mfi);,
                     const Array4<Real>& node_grad_y = node_grad_y_mf.array(mfi);,
                     const Array4<Real>& node_grad_z = node_grad_z_mf.array(mfi););

        AMREX_D_TERM(const Array4<Real> & forcex = div_reversible_stress[0].array(mfi);,
                     const Array4<Real> & forcey = div_reversible_stress[1].array(mfi);,
                     const Array4<Real> & forcez = div_reversible_stress[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

#if (AMREX_SPACEDIM == 2)
        
        amrex::ParallelFor(bx_x, bx_y,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for( int n=0 ; n< nspecies; n++){
               for( int m=0 ; m< nspecies; m++){
                  Real cax_local_plus  = 0.25*(node_grad_x(i,j,k,n)+node_grad_x(i,j+1,k,n)+node_grad_x(i+1,j+1,k,n)+node_grad_x(i+1,j,k,n));
                  Real cay_local_plus  = 0.25*(node_grad_y(i,j,k,n)+node_grad_y(i,j+1,k,n)+node_grad_y(i+1,j+1,k,n)+node_grad_y(i+1,j,k,n));
                  Real cax_local_minus = 0.25*(node_grad_x(i,j,k,n)+node_grad_x(i,j+1,k,n)+node_grad_x(i-1,j+1,k,n)+node_grad_x(i-1,j,k,n));
                  Real cay_local_minus = 0.25*(node_grad_y(i,j,k,n)+node_grad_y(i,j+1,k,n)+node_grad_y(i-1,j+1,k,n)+node_grad_y(i-1,j,k,n));
                  Real cbx_local_plus  = 0.25*(node_grad_x(i,j,k,m)+node_grad_x(i,j+1,k,m)+node_grad_x(i+1,j+1,k,m)+node_grad_x(i+1,j,k,m));
                  Real cby_local_plus  = 0.25*(node_grad_y(i,j,k,m)+node_grad_y(i,j+1,k,m)+node_grad_y(i+1,j+1,k,m)+node_grad_y(i+1,j,k,m));
                  Real cbx_local_minus = 0.25*(node_grad_x(i,j,k,m)+node_grad_x(i,j+1,k,m)+node_grad_x(i-1,j+1,k,m)+node_grad_x(i-1,j,k,m));
                  Real cby_local_minus = 0.25*(node_grad_y(i,j,k,m)+node_grad_y(i,j+1,k,m)+node_grad_y(i-1,j+1,k,m)+node_grad_y(i-1,j,k,m));

                  forcex(i,j,k) += scale_factor * fh_kappa(n,m)*( (node_grad_x(i,j+1,k,n)*node_grad_y(i,j+1,k,m) - node_grad_x(i,j,k,n)*node_grad_y(i,j,k,m))/dx[1] 
                      +0.5*(cax_local_plus*cbx_local_plus-cay_local_plus*cby_local_plus -cax_local_minus*cbx_local_minus+cay_local_minus*cby_local_minus)/dx[0]);
               }
            }
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for( int n=0 ; n< nspecies; n++){
               for( int m=0 ; m< nspecies; m++){
                  Real cax_local_plus  = 0.25*(node_grad_x(i,j,k,n)+node_grad_x(i,j+1,k,n)+node_grad_x(i+1,j+1,k,n)+node_grad_x(i+1,j,k,n));
                  Real cay_local_plus  = 0.25*(node_grad_y(i,j,k,n)+node_grad_y(i,j+1,k,n)+node_grad_y(i+1,j+1,k,n)+node_grad_y(i+1,j,k,n));
                  Real cax_local_minus = 0.25*(node_grad_x(i,j,k,n)+node_grad_x(i,j-1,k,n)+node_grad_x(i+1,j-1,k,n)+node_grad_x(i+1,j,k,n));
                  Real cay_local_minus = 0.25*(node_grad_y(i,j,k,n)+node_grad_y(i,j-1,k,n)+node_grad_y(i+1,j-1,k,n)+node_grad_y(i+1,j,k,n));
                  Real cbx_local_plus  = 0.25*(node_grad_x(i,j,k,m)+node_grad_x(i,j+1,k,m)+node_grad_x(i+1,j+1,k,m)+node_grad_x(i+1,j,k,m));
                  Real cby_local_plus  = 0.25*(node_grad_y(i,j,k,m)+node_grad_y(i,j+1,k,m)+node_grad_y(i+1,j+1,k,m)+node_grad_y(i+1,j,k,m));
                  Real cbx_local_minus = 0.25*(node_grad_x(i,j,k,m)+node_grad_x(i,j-1,k,m)+node_grad_x(i+1,j-1,k,m)+node_grad_x(i+1,j,k,m));
                  Real cby_local_minus = 0.25*(node_grad_y(i,j,k,m)+node_grad_y(i,j-1,k,m)+node_grad_y(i+1,j-1,k,m)+node_grad_y(i+1,j,k,m));

                  forcey(i,j,k) += scale_factor*fh_kappa(n,m)*((node_grad_x(i+1,j,k,n)*node_grad_y(i+1,j,k,m) - node_grad_x(i,j,k,n)*node_grad_y(i,j,k,m))/dx[0]
                      +0.5*(cay_local_plus*cby_local_plus- cay_local_minus*cby_local_minus-cax_local_plus*cbx_local_plus+cax_local_minus*cbx_local_minus)/dx[1]);
               }
            }
        });
        
#elif (AMREX_SPACEDIM == 3)
        
        amrex::ParallelFor(bx_x, bx_y, bx_z,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for( int n=0 ; n< nspecies; n++){
               for( int m=0 ; m< nspecies; m++){

                  Real cax_local_plus  = .125*(node_grad_x(i  ,j  ,k  ,n)+node_grad_x(i  ,j+1,k  ,n)+node_grad_x(i+1,j+1,k  ,n)+node_grad_x(i+1,j  ,k  ,n)
                                              +node_grad_x(i  ,j  ,k+1,n)+node_grad_x(i  ,j+1,k+1,n)+node_grad_x(i+1,j+1,k+1,n)+node_grad_x(i+1,j  ,k+1,n));

                  Real cax_local_minus = .125*(node_grad_x(i  ,j  ,k  ,n)+node_grad_x(i  ,j+1,k  ,n)+node_grad_x(i-1,j+1,k  ,n)+node_grad_x(i-1,j  ,k  ,n)
                                              +node_grad_x(i  ,j  ,k+1,n)+node_grad_x(i  ,j+1,k+1,n)+node_grad_x(i-1,j+1,k+1,n)+node_grad_x(i-1,j  ,k+1,n));

                  Real cay_local_plus  = .125*(node_grad_y(i  ,j  ,k  ,n)+node_grad_y(i  ,j+1,k  ,n)+node_grad_y(i+1,j+1,k  ,n)+node_grad_y(i+1,j  ,k  ,n)
                                              +node_grad_y(i  ,j  ,k+1,n)+node_grad_y(i  ,j+1,k+1,n)+node_grad_y(i+1,j+1,k+1,n)+node_grad_y(i+1,j  ,k+1,n));

                  Real cay_local_minus = .125*(node_grad_y(i  ,j  ,k  ,n)+node_grad_y(i  ,j+1,k  ,n)+node_grad_y(i-1,j+1,k  ,n)+node_grad_y(i-1,j  ,k  ,n)
                                              +node_grad_y(i  ,j  ,k+1,n)+node_grad_y(i  ,j+1,k+1,n)+node_grad_y(i-1,j+1,k+1,n)+node_grad_y(i-1,j  ,k+1,n));

                  Real caz_local_plus  = .125*(node_grad_z(i  ,j  ,k  ,n)+node_grad_z(i  ,j+1,k  ,n)+node_grad_z(i+1,j+1,k  ,n)+node_grad_z(i+1,j  ,k  ,n)
                                              +node_grad_z(i  ,j  ,k+1,n)+node_grad_z(i  ,j+1,k+1,n)+node_grad_z(i+1,j+1,k+1,n)+node_grad_z(i+1,j  ,k+1,n));

                  Real caz_local_minus = .125*(node_grad_z(i  ,j  ,k  ,n)+node_grad_z(i  ,j+1,k  ,n)+node_grad_z(i-1,j+1,k  ,n)+node_grad_z(i-1,j  ,k  ,n)
                                              +node_grad_z(i  ,j  ,k+1,n)+node_grad_z(i  ,j+1,k+1,n)+node_grad_z(i-1,j+1,k+1,n)+node_grad_z(i-1,j  ,k+1,n));

                  Real cbx_local_plus  = .125*(node_grad_x(i  ,j  ,k  ,m)+node_grad_x(i  ,j+1,k  ,m)+node_grad_x(i+1,j+1,k  ,m)+node_grad_x(i+1,j  ,k  ,m)
                                              +node_grad_x(i  ,j  ,k+1,m)+node_grad_x(i  ,j+1,k+1,m)+node_grad_x(i+1,j+1,k+1,m)+node_grad_x(i+1,j  ,k+1,m));

                  Real cbx_local_minus = .125*(node_grad_x(i  ,j  ,k  ,m)+node_grad_x(i  ,j+1,k  ,m)+node_grad_x(i-1,j+1,k  ,m)+node_grad_x(i-1,j  ,k  ,m)
                                              +node_grad_x(i  ,j  ,k+1,m)+node_grad_x(i  ,j+1,k+1,m)+node_grad_x(i-1,j+1,k+1,m)+node_grad_x(i-1,j  ,k+1,m));

                  Real cby_local_plus  = .125*(node_grad_y(i  ,j  ,k  ,m)+node_grad_y(i  ,j+1,k  ,m)+node_grad_y(i+1,j+1,k  ,m)+node_grad_y(i+1,j  ,k  ,m)
                                              +node_grad_y(i  ,j  ,k+1,m)+node_grad_y(i  ,j+1,k+1,m)+node_grad_y(i+1,j+1,k+1,m)+node_grad_y(i+1,j  ,k+1,m));

                  Real cby_local_minus = .125*(node_grad_y(i  ,j  ,k  ,m)+node_grad_y(i  ,j+1,k  ,m)+node_grad_y(i-1,j+1,k  ,m)+node_grad_y(i-1,j  ,k  ,m)
                                              +node_grad_y(i  ,j  ,k+1,m)+node_grad_y(i  ,j+1,k+1,m)+node_grad_y(i-1,j+1,k+1,m)+node_grad_y(i-1,j  ,k+1,m));

                  Real cbz_local_plus  = .125*(node_grad_z(i  ,j  ,k  ,m)+node_grad_z(i  ,j+1,k  ,m)+node_grad_z(i+1,j+1,k  ,m)+node_grad_z(i+1,j  ,k  ,m)
                                              +node_grad_z(i  ,j  ,k+1,m)+node_grad_z(i  ,j+1,k+1,m)+node_grad_z(i+1,j+1,k+1,m)+node_grad_z(i+1,j  ,k+1,m));

                  Real cbz_local_minus = .125*(node_grad_z(i  ,j  ,k  ,m)+node_grad_z(i  ,j+1,k  ,m)+node_grad_z(i-1,j+1,k  ,m)+node_grad_z(i-1,j  ,k  ,m)
                                              +node_grad_z(i  ,j  ,k+1,m)+node_grad_z(i  ,j+1,k+1,m)+node_grad_z(i-1,j+1,k+1,m)+node_grad_z(i-1,j  ,k+1,m));

                  forcex(i,j,k) += scale_factor * fh_kappa(n,m)*0.5*
                      ( ((node_grad_x(i,j+1,k,n)*node_grad_y(i,j+1,k,m) + node_grad_x(i,j+1,k+1,n)*node_grad_y(i,j+1,k+1,m)) 
                       - (node_grad_x(i,j  ,k,n)*node_grad_y(i,j  ,k,m) + node_grad_x(i,j  ,k+1,n)*node_grad_y(i,j  ,k+1,m)))/dx[1] 
                       +((cax_local_plus *cbx_local_plus - caz_local_plus *cbz_local_plus -cay_local_plus *cby_local_plus) 
                        -(cax_local_minus*cbx_local_minus- caz_local_minus*cbz_local_minus-cay_local_minus*cby_local_minus))/(dx[0])  
                      + ((node_grad_x(i,j,k+1,n)*node_grad_z(i,j,k+1,m) + node_grad_x(i,j+1,k+1,n)*node_grad_z(i,j+1,k+1,m)) -
                         (node_grad_x(i,j,k  ,n)*node_grad_z(i,j,k  ,m) + node_grad_x(i,j+1,k  ,n)*node_grad_z(i,j+1,k  ,m)))/dx[2]);
                  
               }
            }
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for( int n=0 ; n< nspecies; n++){
               for( int m=0 ; m< nspecies; m++){

                  Real cax_local_plus  = .125*(node_grad_x(i  ,j  ,k  ,n)+node_grad_x(i  ,j+1,k  ,n)+node_grad_x(i+1,j+1,k  ,n)+node_grad_x(i+1,j  ,k  ,n)
                                              +node_grad_x(i  ,j  ,k+1,n)+node_grad_x(i  ,j+1,k+1,n)+node_grad_x(i+1,j+1,k+1,n)+node_grad_x(i+1,j  ,k+1,n));

                  Real cax_local_minus = .125*(node_grad_x(i  ,j  ,k  ,n)+node_grad_x(i  ,j-1,k  ,n)+node_grad_x(i+1,j-1,k  ,n)+node_grad_x(i+1,j  ,k  ,n)
                                              +node_grad_x(i  ,j  ,k+1,n)+node_grad_x(i  ,j-1,k+1,n)+node_grad_x(i+1,j-1,k+1,n)+node_grad_x(i+1,j  ,k+1,n));

                  Real cay_local_plus  = .125*(node_grad_y(i  ,j  ,k  ,n)+node_grad_y(i  ,j+1,k  ,n)+node_grad_y(i+1,j+1,k  ,n)+node_grad_y(i+1,j  ,k  ,n)
                                              +node_grad_y(i  ,j  ,k+1,n)+node_grad_y(i  ,j+1,k+1,n)+node_grad_y(i+1,j+1,k+1,n)+node_grad_y(i+1,j  ,k+1,n));

                  Real cay_local_minus = .125*(node_grad_y(i  ,j  ,k  ,n)+node_grad_y(i  ,j-1,k  ,n)+node_grad_y(i+1,j-1,k  ,n)+node_grad_y(i+1,j  ,k  ,n)
                                              +node_grad_y(i  ,j  ,k+1,n)+node_grad_y(i  ,j-1,k+1,n)+node_grad_y(i+1,j-1,k+1,n)+node_grad_y(i+1,j  ,k+1,n));

                  Real caz_local_plus  = .125*(node_grad_z(i  ,j  ,k  ,n)+node_grad_z(i  ,j+1,k  ,n)+node_grad_z(i+1,j+1,k  ,n)+node_grad_z(i+1,j  ,k  ,n)
                                              +node_grad_z(i  ,j  ,k+1,n)+node_grad_z(i  ,j+1,k+1,n)+node_grad_z(i+1,j+1,k+1,n)+node_grad_z(i+1,j  ,k+1,n));

                  Real caz_local_minus = .125*(node_grad_z(i  ,j  ,k  ,n)+node_grad_z(i  ,j-1,k  ,n)+node_grad_z(i+1,j-1,k  ,n)+node_grad_z(i+1,j  ,k  ,n)
                                              +node_grad_z(i  ,j  ,k+1,n)+node_grad_z(i  ,j-1,k+1,n)+node_grad_z(i+1,j-1,k+1,n)+node_grad_z(i+1,j  ,k+1,n));


                  Real cbx_local_plus  = .125*(node_grad_x(i  ,j  ,k  ,m)+node_grad_x(i  ,j+1,k  ,m)+node_grad_x(i+1,j+1,k  ,m)+node_grad_x(i+1,j  ,k  ,m)
                                              +node_grad_x(i  ,j  ,k+1,m)+node_grad_x(i  ,j+1,k+1,m)+node_grad_x(i+1,j+1,k+1,m)+node_grad_x(i+1,j  ,k+1,m));

                  Real cbx_local_minus = .125*(node_grad_x(i  ,j  ,k  ,m)+node_grad_x(i  ,j-1,k  ,m)+node_grad_x(i+1,j-1,k  ,m)+node_grad_x(i+1,j  ,k  ,m)
                                              +node_grad_x(i  ,j  ,k+1,m)+node_grad_x(i  ,j-1,k+1,m)+node_grad_x(i+1,j-1,k+1,m)+node_grad_x(i+1,j  ,k+1,m));

                  Real cby_local_plus  = .125*(node_grad_y(i  ,j  ,k  ,m)+node_grad_y(i  ,j+1,k  ,m)+node_grad_y(i+1,j+1,k  ,m)+node_grad_y(i+1,j  ,k  ,m)
                                              +node_grad_y(i  ,j  ,k+1,m)+node_grad_y(i  ,j+1,k+1,m)+node_grad_y(i+1,j+1,k+1,m)+node_grad_y(i+1,j  ,k+1,m));

                  Real cby_local_minus = .125*(node_grad_y(i  ,j  ,k  ,m)+node_grad_y(i  ,j-1,k  ,m)+node_grad_y(i+1,j-1,k  ,m)+node_grad_y(i+1,j  ,k  ,m)
                                              +node_grad_y(i  ,j  ,k+1,m)+node_grad_y(i  ,j-1,k+1,m)+node_grad_y(i+1,j-1,k+1,m)+node_grad_y(i+1,j  ,k+1,m));

                  Real cbz_local_plus  = .125*(node_grad_z(i  ,j  ,k  ,m)+node_grad_z(i  ,j+1,k  ,m)+node_grad_z(i+1,j+1,k  ,m)+node_grad_z(i+1,j  ,k  ,m)
                                              +node_grad_z(i  ,j  ,k+1,m)+node_grad_z(i  ,j+1,k+1,m)+node_grad_z(i+1,j+1,k+1,m)+node_grad_z(i+1,j  ,k+1,m));

                  Real cbz_local_minus = .125*(node_grad_z(i  ,j  ,k  ,m)+node_grad_z(i  ,j-1,k  ,m)+node_grad_z(i+1,j-1,k  ,m)+node_grad_z(i+1,j  ,k  ,m)
                                              +node_grad_z(i  ,j  ,k+1,m)+node_grad_z(i  ,j-1,k+1,m)+node_grad_z(i+1,j-1,k+1,m)+node_grad_z(i+1,j  ,k+1,m));

                  forcey(i,j,k) += scale_factor * fh_kappa(n,m)*0.5* (
                       ((node_grad_y(i+1,j,k,n)*node_grad_x(i+1,j,k,m) + node_grad_y(i+1,j,k+1,n)*node_grad_x(i+1,j,k+1,m)) 
                       -(node_grad_y(i  ,j,k,n)*node_grad_x(i  ,j,k,m) + node_grad_y(i  ,j,k+1,n)*node_grad_x(i  ,j,k+1,m)))/dx[0] 
                       +((cay_local_plus *cby_local_plus -caz_local_plus *cbz_local_plus -cax_local_plus *cbx_local_plus) 
                        -(cay_local_minus*cby_local_minus-caz_local_minus*cbz_local_minus-cax_local_minus*cbx_local_minus))/dx[1]  
                      +((node_grad_y(i,j,k+1,n)*node_grad_z(i,j,k+1,m) + node_grad_y(i+1,j,k+1,n)*node_grad_z(i+1,j,k+1,m)) -
                        (node_grad_y(i,j,k  ,n)*node_grad_z(i,j,k  ,m) + node_grad_y(i+1,j,k  ,n)*node_grad_z(i+1,j,k  ,m)))/dx[2]);
                  
               }
            }
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for( int n=0 ; n< nspecies; n++){
               for( int m=0 ; m< nspecies; m++){

                  Real cax_local_plus  = .125*(node_grad_x(i  ,j  ,k  ,n)+node_grad_x(i  ,j+1,k  ,n)+node_grad_x(i+1,j+1,k  ,n)+node_grad_x(i+1,j  ,k  ,n)
                                              +node_grad_x(i  ,j  ,k+1,n)+node_grad_x(i  ,j+1,k+1,n)+node_grad_x(i+1,j+1,k+1,n)+node_grad_x(i+1,j  ,k+1,n));

                  Real cax_local_minus = .125*(node_grad_x(i  ,j  ,k  ,n)+node_grad_x(i  ,j+1,k  ,n)+node_grad_x(i+1,j+1,k  ,n)+node_grad_x(i+1,j  ,k  ,n)
                                              +node_grad_x(i  ,j  ,k-1,n)+node_grad_x(i  ,j+1,k-1,n)+node_grad_x(i+1,j+1,k-1,n)+node_grad_x(i+1,j  ,k-1,n));

                  Real cay_local_plus  = .125*(node_grad_y(i  ,j  ,k  ,n)+node_grad_y(i  ,j+1,k  ,n)+node_grad_y(i+1,j+1,k  ,n)+node_grad_y(i+1,j  ,k  ,n)
                                              +node_grad_y(i  ,j  ,k+1,n)+node_grad_y(i  ,j+1,k+1,n)+node_grad_y(i+1,j+1,k+1,n)+node_grad_y(i+1,j  ,k+1,n));

                  Real cay_local_minus = .125*(node_grad_y(i  ,j  ,k  ,n)+node_grad_y(i  ,j+1,k  ,n)+node_grad_y(i+1,j+1,k  ,n)+node_grad_y(i+1,j  ,k  ,n)
                                              +node_grad_y(i  ,j  ,k-1,n)+node_grad_y(i  ,j+1,k-1,n)+node_grad_y(i+1,j+1,k-1,n)+node_grad_y(i+1,j  ,k-1,n));

                  Real caz_local_plus  = .125*(node_grad_z(i  ,j  ,k  ,n)+node_grad_z(i  ,j+1,k  ,n)+node_grad_z(i+1,j+1,k  ,n)+node_grad_z(i+1,j  ,k  ,n)
                                              +node_grad_z(i  ,j  ,k+1,n)+node_grad_z(i  ,j+1,k+1,n)+node_grad_z(i+1,j+1,k+1,n)+node_grad_z(i+1,j  ,k+1,n));

                  Real caz_local_minus = .125*(node_grad_z(i  ,j  ,k  ,n)+node_grad_z(i  ,j+1,k  ,n)+node_grad_z(i+1,j+1,k  ,n)+node_grad_z(i+1,j  ,k  ,n)
                                              +node_grad_z(i  ,j  ,k-1,n)+node_grad_z(i  ,j+1,k-1,n)+node_grad_z(i+1,j+1,k-1,n)+node_grad_z(i+1,j  ,k-1,n));

                  Real cbx_local_plus  = .125*(node_grad_x(i  ,j  ,k  ,m)+node_grad_x(i  ,j+1,k  ,m)+node_grad_x(i+1,j+1,k  ,m)+node_grad_x(i+1,j  ,k  ,m)
                                              +node_grad_x(i  ,j  ,k+1,m)+node_grad_x(i  ,j+1,k+1,m)+node_grad_x(i+1,j+1,k+1,m)+node_grad_x(i+1,j  ,k+1,m));

                  Real cbx_local_minus = .125*(node_grad_x(i  ,j  ,k  ,m)+node_grad_x(i  ,j+1,k  ,m)+node_grad_x(i+1,j+1,k  ,m)+node_grad_x(i+1,j  ,k  ,m)
                                              +node_grad_x(i  ,j  ,k-1,m)+node_grad_x(i  ,j+1,k-1,m)+node_grad_x(i+1,j+1,k-1,m)+node_grad_x(i+1,j  ,k-1,m));

                  Real cby_local_plus  = .125*(node_grad_y(i  ,j  ,k  ,m)+node_grad_y(i  ,j+1,k  ,m)+node_grad_y(i+1,j+1,k  ,m)+node_grad_y(i+1,j  ,k  ,m)
                                              +node_grad_y(i  ,j  ,k+1,m)+node_grad_y(i  ,j+1,k+1,m)+node_grad_y(i+1,j+1,k+1,m)+node_grad_y(i+1,j  ,k+1,m));

                  Real cby_local_minus = .125*(node_grad_y(i  ,j  ,k  ,m)+node_grad_y(i  ,j+1,k  ,m)+node_grad_y(i+1,j+1,k  ,m)+node_grad_y(i+1,j  ,k  ,m)
                                              +node_grad_y(i  ,j  ,k-1,m)+node_grad_y(i  ,j+1,k-1,m)+node_grad_y(i+1,j+1,k-1,m)+node_grad_y(i+1,j  ,k-1,m));

                  Real cbz_local_plus  = .125*(node_grad_z(i  ,j  ,k  ,m)+node_grad_z(i  ,j+1,k  ,m)+node_grad_z(i+1,j+1,k  ,m)+node_grad_z(i+1,j  ,k  ,m)
                                              +node_grad_z(i  ,j  ,k+1,m)+node_grad_z(i  ,j+1,k+1,m)+node_grad_z(i+1,j+1,k+1,m)+node_grad_z(i+1,j  ,k+1,m));

                  Real cbz_local_minus = .125*(node_grad_z(i  ,j  ,k  ,m)+node_grad_z(i  ,j+1,k  ,m)+node_grad_z(i+1,j+1,k  ,m)+node_grad_z(i+1,j  ,k  ,m)
                                              +node_grad_z(i  ,j  ,k-1,m)+node_grad_z(i  ,j+1,k-1,m)+node_grad_z(i+1,j+1,k-1,m)+node_grad_z(i+1,j  ,k-1,m));



                  forcez(i,j,k) += scale_factor * fh_kappa(n,m) * 0.5 * (
                       ((node_grad_z(i+1,j,k,n)*node_grad_x(i+1,j,k,m) + node_grad_z(i+1,j+1,k,n)*node_grad_x(i+1,j+1,k,m)) 
                       -(node_grad_z(i  ,j,k,n)*node_grad_x(i  ,j,k,m) + node_grad_z(i  ,j+1,k,n)*node_grad_x(i  ,j+1,k,m)))/dx[0] 
                      +((caz_local_plus *cbz_local_plus -cay_local_plus *cby_local_plus -cax_local_plus *cbx_local_plus) 
                       -(caz_local_minus*cbz_local_minus-cay_local_minus*cby_local_minus-cax_local_minus*cbx_local_minus))/dx[2]  
                      +((node_grad_z(i,j+1,k,n)*node_grad_y(i,j+1,k,m) + node_grad_z(i+1,j+1,k,n)*node_grad_y(i+1,j+1,k,m)) -
                        (node_grad_z(i,j  ,k,n)*node_grad_y(i,j  ,k,m) + node_grad_z(i+1,j  ,k,n)*node_grad_y(i+1,j  ,k,m)))/dx[1]);
                  
               }
            }
        });
                           
#endif
    }

    // set force on walls to be zero since normal velocity is zero
    ZeroEdgevalWalls(div_reversible_stress, geom, 0, 1);
    
}

void ComputeDivReversibleStress(std::array<MultiFab,AMREX_SPACEDIM>& div_reversible_stress,
                                const MultiFab& rhotot_in,
                                MultiFab& rho_in,
                                const Geometry& geom)
{
    BL_PROFILE_VAR("ComputeDivReversibleStress()",ComputeDivReversibleStress);
    
    BoxArray ba = rho_in.boxArray();
    DistributionMapping dmap = rho_in.DistributionMap();
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    MultiFab node_grad_c_mf(convert(ba,nodal_flag), dmap, AMREX_SPACEDIM, 1);

    MultiFab conc(ba, dmap, nspecies, 2);

    // rho to conc - VALID REGION ONLY
    ConvertRhoCToC(rho_in,rhotot_in,conc,1);

    // fill conc ghost cells
    conc.FillBoundary(geom.periodicity());
    MultiFabPhysBC(conc,geom,0,nspecies,SPEC_BC_COMP);    

    Real scale_factor = rhobar[0]*k_B*T_init[0]/molmass[0];
    
    for ( MFIter mfi(node_grad_c_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.growntilebox(1);

        const Array4<Real>& node_grad_c = node_grad_c_mf.array(mfi);
        const Array4<Real const>& c = conc.array(mfi);
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if (AMREX_SPACEDIM == 2)            

            node_grad_c(i,j,k,0) = (c(i,j,k)-c(i-1,j,k)+c(i,j-1,k)-c(i-1,j-1,k))/(2*dx[0]);
            node_grad_c(i,j,k,1) = (c(i,j,k)-c(i,j-1,k)+c(i-1,j,k)-c(i-1,j-1,k))/(2*dx[1]);

#elif (AMREX_SPACEDIM == 3)

            node_grad_c(i,j,k,0) = (c(i,j,k)-c(i-1,j,k)+c(i,j-1,k)-c(i-1,j-1,k)
                                    +c(i,j,k-1)-c(i-1,j,k-1)+c(i,j-1,k-1)-c(i-1,j-1,k-1))/(4*dx[0]);
            node_grad_c(i,j,k,1) = (c(i,j,k)-c(i,j-1,k)+c(i-1,j,k)-c(i-1,j-1,k)
                                    +c(i,j,k-1)-c(i,j-1,k-1)+c(i-1,j,k-1)-c(i-1,j-1,k-1))/(4*dx[1]);
            node_grad_c(i,j,k,2) = (c(i,j,k)-c(i,j,k-1)+c(i-1,j,k)-c(i-1,j,k-1)
                                    +c(i,j-1,k)-c(i,j-1,k-1)+c(i-1,j-1,k)-c(i-1,j-1,k-1))/(4*dx[2]);
#endif
        });
    }

    for ( MFIter mfi(rhotot_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Array4<Real>& node_grad_c = node_grad_c_mf.array(mfi);

        AMREX_D_TERM(const Array4<Real> & forcex = div_reversible_stress[0].array(mfi);,
                     const Array4<Real> & forcey = div_reversible_stress[1].array(mfi);,
                     const Array4<Real> & forcez = div_reversible_stress[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

#if (AMREX_SPACEDIM == 2)
        
        amrex::ParallelFor(bx_x, bx_y,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real cx_local_plus  = 0.25*(node_grad_c(i,j,k,0)+node_grad_c(i,j+1,k,0)+node_grad_c(i+1,j+1,k,0)+node_grad_c(i+1,j,k,0));
            Real cy_local_plus  = 0.25*(node_grad_c(i,j,k,1)+node_grad_c(i,j+1,k,1)+node_grad_c(i+1,j+1,k,1)+node_grad_c(i+1,j,k,1));
            Real cx_local_minus = 0.25*(node_grad_c(i,j,k,0)+node_grad_c(i,j+1,k,0)+node_grad_c(i-1,j+1,k,0)+node_grad_c(i-1,j,k,0));
            Real cy_local_minus = 0.25*(node_grad_c(i,j,k,1)+node_grad_c(i,j+1,k,1)+node_grad_c(i-1,j+1,k,1)+node_grad_c(i-1,j,k,1));

            forcex(i,j,k) = -(node_grad_c(i,j+1,k,0)*node_grad_c(i,j+1,k,1) - node_grad_c(i,j,k,0)*node_grad_c(i,j,k,1))/dx[1] 
                +(0.5*(cy_local_plus*cy_local_plus-cx_local_plus*cx_local_plus)
                  -0.5*(cy_local_minus*cy_local_minus-cx_local_minus*cx_local_minus))/(dx[0]);
            forcex(i,j,k) = scale_factor * kc_tension*forcex(i,j,k);
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real cx_local_plus  = 0.25*(node_grad_c(i,j,k,0)+node_grad_c(i,j+1,k,0)+node_grad_c(i+1,j+1,k,0)+node_grad_c(i+1,j,k,0));
            Real cy_local_plus  = 0.25*(node_grad_c(i,j,k,1)+node_grad_c(i,j+1,k,1)+node_grad_c(i+1,j+1,k,1)+node_grad_c(i+1,j,k,1));
            Real cx_local_minus = 0.25*(node_grad_c(i,j,k,0)+node_grad_c(i,j-1,k,0)+node_grad_c(i+1,j-1,k,0)+node_grad_c(i+1,j,k,0));
            Real cy_local_minus = 0.25*(node_grad_c(i,j,k,1)+node_grad_c(i,j-1,k,1)+node_grad_c(i+1,j-1,k,1)+node_grad_c(i+1,j,k,1));

            forcey(i,j,k) = -(node_grad_c(i+1,j,k,0)*node_grad_c(i+1,j,k,1) - node_grad_c(i,j,k,0)*node_grad_c(i,j,k,1))/dx[0]
                +(0.5*(cx_local_plus*cx_local_plus-cy_local_plus*cy_local_plus)
                  -0.5*(cx_local_minus*cx_local_minus-cy_local_minus*cy_local_minus))/(dx[1]);
            forcey(i,j,k) = scale_factor * kc_tension*forcey(i,j,k);            
        });
        
#elif (AMREX_SPACEDIM == 3)
        
        amrex::ParallelFor(bx_x, bx_y, bx_z,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real cx_local_plus = (.125)*(node_grad_c(i,j,k,0)+node_grad_c(i,j+1,k,0)+node_grad_c(i+1,j+1,k,0) 
                                         +node_grad_c(i+1,j,k,0)+node_grad_c(i+1,j+1,k+1,0)+node_grad_c(i,j,k+1,0) 
                                         +node_grad_c(i,j+1,k+1,0)+node_grad_c(i+1,j,k+1,0));
            Real cy_local_plus = (.125)*(node_grad_c(i,j,k,1)+node_grad_c(i,j+1,k,1)+node_grad_c(i+1,j+1,k,1) 
                                         +node_grad_c(i+1,j,k,1)+node_grad_c(i+1,j+1,k+1,1)+node_grad_c(i,j,k+1,1) 
                                         +node_grad_c(i,j+1,k+1,1)+node_grad_c(i+1,j,k+1,1));
            Real cz_local_plus = (.125)*(node_grad_c(i,j,k,2)+node_grad_c(i,j+1,k,2)+node_grad_c(i+1,j+1,k,2) 
                                         +node_grad_c(i+1,j,k,2)+node_grad_c(i+1,j+1,k+1,2)+node_grad_c(i,j,k+1,2)+ 
                                         node_grad_c(i,j+1,k+1,2)+node_grad_c(i+1,j,k+1,2));
            Real cx_local_minus = (.125)*(node_grad_c(i,j,k,0)+node_grad_c(i,j+1,k,0)+node_grad_c(i-1,j+1,k,0)+
                                          node_grad_c(i-1,j,k,0)+node_grad_c(i,j,k+1,0)+node_grad_c(i,j+1,k+1,0)+
                                          node_grad_c(i-1,j+1,k+1,0)+node_grad_c(i-1,j,k+1,0));
            Real cy_local_minus = (.125)*(node_grad_c(i,j,k,1)+node_grad_c(i,j+1,k,1)+node_grad_c(i-1,j+1,k,1)+
                                          node_grad_c(i-1,j,k,1)+node_grad_c(i,j,k+1,1)+node_grad_c(i,j+1,k+1,1)+
                                          node_grad_c(i-1,j+1,k+1,1)+node_grad_c(i-1,j,k+1,1));
            Real cz_local_minus = (.125)*(node_grad_c(i,j,k,2)+node_grad_c(i,j+1,k,2)+node_grad_c(i-1,j+1,k,2)
                                          +node_grad_c(i-1,j,k,2)+node_grad_c(i,j,k+1,2)+node_grad_c(i,j+1,k+1,2)+
                                          node_grad_c(i-1,j+1,k+1,2)+node_grad_c(i-1,j,k+1,2));

            forcex(i,j,k) = -(.5*(node_grad_c(i,j+1,k,0)*node_grad_c(i,j+1,k,1) 
                                  +node_grad_c(i,j+1,k+1,0)*node_grad_c(i,j+1,k+1,1)) -
                              .5*(node_grad_c(i,j,k,0)*node_grad_c(i,j,k,1)+
                                  node_grad_c(i,j,k+1,0)*node_grad_c(i,j,k+1,1)))/dx[1] 
                +(0.5*(cy_local_plus*cy_local_plus+cz_local_plus*cz_local_plus-cx_local_plus*cx_local_plus) 
                  -0.5*(cy_local_minus*cy_local_minus+cz_local_minus*cz_local_minus-cx_local_minus*cx_local_minus))/(dx[0])  
                -(.5*(node_grad_c(i,j,k+1,0)*node_grad_c(i,j,k+1,2) 
                      +node_grad_c(i,j+1,k+1,0)*node_grad_c(i,j+1,k+1,2)) -
                  .5*(node_grad_c(i,j,k,0)*node_grad_c(i,j,k,2)+node_grad_c(i,j+1,k,0)*node_grad_c(i,j+1,k,2)))/dx[2];
            
            forcex(i,j,k) = scale_factor * kc_tension*forcex(i,j,k);
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real cx_local_plus = (.125)*(node_grad_c(i,j,k,0)+node_grad_c(i,j+1,k,0)+
                                         node_grad_c(i+1,j+1,k,0)+node_grad_c(i+1,j,k,0)+node_grad_c(i+1,j+1,k+1,0)+
                                         node_grad_c(i,j,k+1,0)+node_grad_c(i,j+1,k+1,0)+node_grad_c(i+1,j,k+1,0));
            Real cy_local_plus = (.125)*(node_grad_c(i,j,k,1)+node_grad_c(i,j+1,k,1)+node_grad_c(i+1,j+1,k,1)+
                                         node_grad_c(i+1,j,k,1)+node_grad_c(i+1,j+1,k+1,1)+
                                         node_grad_c(i,j,k+1,1)+node_grad_c(i,j+1,k+1,1)+node_grad_c(i+1,j,k+1,1));
            Real cz_local_plus = (.125)*(node_grad_c(i,j,k,2)+node_grad_c(i,j+1,k,2)+
                                         node_grad_c(i+1,j+1,k,2)+node_grad_c(i+1,j,k,2)+node_grad_c(i+1,j+1,k+1,2)+
                                         node_grad_c(i,j,k+1,2)+node_grad_c(i,j+1,k+1,2)+node_grad_c(i+1,j,k+1,2));

            Real cx_local_minus = (.125)*(node_grad_c(i,j,k,0)+node_grad_c(i,j-1,k,0)+
                                          node_grad_c(i+1,j-1,k,0)+node_grad_c(i+1,j,k,0)+node_grad_c(i,j,k+1,0)+
                                          node_grad_c(i,j-1,k+1,0)+node_grad_c(i+1,j-1,k+1,0)+node_grad_c(i+1,j,k+1,0));
            Real cy_local_minus = (.125)*(node_grad_c(i,j,k,1)+node_grad_c(i,j-1,k,1)+
                                          node_grad_c(i+1,j-1,k,1)+node_grad_c(i+1,j,k,1)+node_grad_c(i,j,k+1,1)+
                                          node_grad_c(i,j-1,k+1,1)+node_grad_c(i+1,j-1,k+1,1)+node_grad_c(i+1,j,k+1,1));
            Real cz_local_minus = (.125)*(node_grad_c(i,j,k,2)+node_grad_c(i,j-1,k,2)+
                                          node_grad_c(i+1,j-1,k,2)+node_grad_c(i+1,j,k,2)+node_grad_c(i,j,k+1,2)+
                                          node_grad_c(i,j-1,k+1,2)+node_grad_c(i+1,j-1,k+1,2)+node_grad_c(i+1,j,k+1,2));

            forcey(i,j,k) = -(.5*(node_grad_c(i+1,j,k,0)*node_grad_c(i+1,j,k,1) 
                                  +node_grad_c(i+1,j,k+1,0)*node_grad_c(i+1,j,k+1,1)) 
                              - .5*(node_grad_c(i,j,k,0)*node_grad_c(i,j,k,1) 
                                    +node_grad_c(i,j,k+1,0)*node_grad_c(i,j,k+1,1)))/dx[0] 
                +(0.5*(cx_local_plus*cx_local_plus+cz_local_plus*cz_local_plus-cy_local_plus*cy_local_plus) 
                  -0.5*(cx_local_minus*cx_local_minus+cz_local_minus*cz_local_minus-cy_local_minus*cy_local_minus))/(dx[1])  
                -(.5*(node_grad_c(i,j,k+1,1)*node_grad_c(i,j,k+1,2) 
                      +node_grad_c(i+1,j,k+1,1)*node_grad_c(i+1,j,k+1,2)) 
                  - .5*(node_grad_c(i,j,k,1)*node_grad_c(i,j,k,2)+node_grad_c(i+1,j,k,1)*node_grad_c(i+1,j,k,2)))/dx[2];
            
            forcey(i,j,k) = scale_factor * kc_tension*forcey(i,j,k);
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real cx_local_plus = (.125)*(node_grad_c(i,j,k,0)+node_grad_c(i,j+1,k,0)+
                                         node_grad_c(i+1,j+1,k,0)+node_grad_c(i+1,j,k,0)+node_grad_c(i+1,j+1,k+1,0)+
                                         node_grad_c(i,j,k+1,0)+node_grad_c(i,j+1,k+1,0)+node_grad_c(i+1,j,k+1,0));
            Real cy_local_plus = (.125)*(node_grad_c(i,j,k,1)+node_grad_c(i,j+1,k,1)+node_grad_c(i+1,j+1,k,1)+
                                         node_grad_c(i+1,j,k,1)+node_grad_c(i+1,j+1,k+1,1)+node_grad_c(i,j,k+1,1)+
                                         node_grad_c(i,j+1,k+1,1)+node_grad_c(i+1,j,k+1,1));
            Real cz_local_plus = (.125)*(node_grad_c(i,j,k,2)+node_grad_c(i,j+1,k,2)+node_grad_c(i+1,j+1,k,2)+
                                         node_grad_c(i+1,j,k,2)+node_grad_c(i+1,j+1,k+1,2)+node_grad_c(i,j,k+1,2)+
                                         node_grad_c(i,j+1,k+1,2)+node_grad_c(i+1,j,k+1,2));

            Real cx_local_minus = (.125)*(node_grad_c(i,j,k,0)+node_grad_c(i,j+1,k,0)+node_grad_c(i+1,j+1,k,0)+ 
                                          node_grad_c(i+1,j,k,0)+node_grad_c(i,j,k-1,0)+node_grad_c(i,j+1,k-1,0) 
                                          +node_grad_c(i+1,j+1,k-1,0)+node_grad_c(i+1,j,k-1,0));
            Real cy_local_minus = (.125)*(node_grad_c(i,j,k,1)+node_grad_c(i,j+1,k,1)+node_grad_c(i+1,j+1,k,1) 
                                          +node_grad_c(i+1,j,k,1)+node_grad_c(i,j,k-1,1)+node_grad_c(i,j+1,k-1,1) 
                                          +node_grad_c(i+1,j+1,k-1,1)+node_grad_c(i+1,j,k-1,1));
            Real cz_local_minus = (.125)*(node_grad_c(i,j,k,2)+node_grad_c(i,j+1,k,2)+node_grad_c(i+1,j+1,k,2) 
                                          +node_grad_c(i+1,j,k,2)+node_grad_c(i,j,k-1,2)+node_grad_c(i,j+1,k-1,2) 
                                          +node_grad_c(i+1,j+1,k-1,2)+node_grad_c(i+1,j,k-1,2));

            forcez(i,j,k) = -(.5*(node_grad_c(i+1,j,k,0)*node_grad_c(i+1,j,k,2) 
                                  +(node_grad_c(i+1,j+1,k,0)*node_grad_c(i+1,j+1,k,2)))  
                              - .5*(node_grad_c(i,j,k,0)*node_grad_c(i,j,k,2) 
                                    +node_grad_c(i,j+1,k,0)*node_grad_c(i,j+1,k,2)))/dx[0] 
                +(0.5*(cx_local_plus*cx_local_plus+cy_local_plus*cy_local_plus-cz_local_plus*cz_local_plus) 
                  -0.5*(cx_local_minus*cx_local_minus+cy_local_minus*cy_local_minus-cz_local_minus*cz_local_minus))/(dx[2])  
                -(.5*(node_grad_c(i,j+1,k,1)*node_grad_c(i,j+1,k,2) 
                      +node_grad_c(i+1,j+1,k,1)*node_grad_c(i+1,j+1,k,2)) 
                  - .5*(node_grad_c(i,j,k,1)*node_grad_c(i,j,k,2)+node_grad_c(i+1,j,k,1)*node_grad_c(i+1,j,k,2)))/dx[1];

            forcez(i,j,k) = scale_factor * kc_tension*forcez(i,j,k);
        });
                           
#endif
    }

    // set force on walls to be zero since normal velocity is zero
    ZeroEdgevalWalls(div_reversible_stress, geom, 0, 1);
    
}

void ComputeDisjoiningPressure(std::array<MultiFab,AMREX_SPACEDIM>& disjoining_pressure,
                                const MultiFab& rhotot_in,
                                MultiFab& rho_in,
                                const Geometry& geom)
{
    BL_PROFILE_VAR("ComputeDivFHReversibleStress()",ComputeDivReversibleStress);
    
    BoxArray ba = rho_in.boxArray();
    DistributionMapping dmap = rho_in.DistributionMap();
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    

    MultiFab conc(ba, dmap, nspecies, 2);

    MultiFab disjp(ba, dmap, nspecies, 2);
    disjp.setVal(0.);

    // rho to conc - VALID REGION ONLY
    ConvertRhoCToC(rho_in,rhotot_in,conc,1);

    // fill conc ghost cells
    conc.FillBoundary(geom.periodicity());
    MultiFabPhysBC(conc,geom,0,nspecies,SPEC_BC_COMP);    
    amrex::Real A_coeff = hamaker_A;
    //A_coeff = -3.5e-14;
    //A_coeff = 4.*1.85e-13;
    //A_coeff = 8.e-15;
    amrex::Real pi = 3.141592653589793238;

    
    for ( MFIter mfi(conc,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.growntilebox(1);
//        const Array4<const Real>& conc_arr = conc.array(mfi);
        const Array4<      Real>& djp = disjp.array(mfi);

	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

             amrex::Real y = (j+0.5)*dx[1];
 //            amrex::Real alpha = (conc_arr(i,j,k,0)-c_init_1[0])/(c_init_1[1]-c_init_1[0]);

//	     alpha = std::max(0.,std::min(1.,alpha));

	     //djp(i,j,k) = (y<dx[1]) ?  0. : A_coeff/(6*pi*y*y*y);
	     djp(i,j,k) = A_coeff/(6*pi*y*y*y);
	     // djp(i,j,k) = alpha*A_coeff/(6*pi*y*y*y);
	     /*
	     if(i == 255){
	        amrex::Print() << "j,y,djp(i,j,k)  " << j << " " << y << " "  << " " << djp(i,j,k) << std::endl; 
	     }
	     */

	});
     }

     amrex::Real F0 = 4.e12;
        for (int n = 0 ; n < AMREX_SPACEDIM; n++)
            disjoining_pressure[n].setVal(0.);

        for ( MFIter mfi(rhotot_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            const Array4<      Real>& djp = disjp.array(mfi);

            AMREX_D_TERM(const Array4<Real> & djpx = disjoining_pressure[0].array(mfi);,
                     const Array4<Real> & djpy = disjoining_pressure[1].array(mfi);,
                     const Array4<Real> & djpz = disjoining_pressure[2].array(mfi););

            const Array4<const Real>& conc_arr = conc.array(mfi);
            const Array4<const Real>& rho_arr = rho_in.array(mfi);

            AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

#if (AMREX_SPACEDIM == 2)
        
            amrex::ParallelFor(bx_x, bx_y,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                amrex::Real alphap = (conc_arr(i,j,k,0)-fh_ce)/(1.-2.*fh_ce);
                amrex::Real alpham = (conc_arr(i-1,j,k,0)-fh_ce)/(1.-2.*fh_ce);

	        alphap = std::max(0.,std::min(1.,alphap));
                alpham = std::max(0.,std::min(1.,alpham));
  //               djpx(i,j,k) = 0.;
                djpx(i,j,k) = 0.5*(alpham+alphap)*F0;
  //               djpx(i,j,k) = 0.5*(alpham+alphap)* (djp(i-1,j,k)-djp(i,j,k))/dx[0];
                 //djpx(i,j,k) = (djp(i,j,k)-djp(i-1,j,k))/dx[0];
            },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
             amrex::Real alphap = (conc_arr(i,j,k,0)-fh_ce)/(1.-2.*fh_ce);
             amrex::Real alpham = (conc_arr(i,j-1,k,0)-fh_ce)/(1.-2.*fh_ce);

	     alphap = std::max(0.,std::min(1.,alphap));
             alpham = std::max(0.,std::min(1.,alpham));
	     /*
	         amrex::Real alphap, alpham;
		 if(conc_arr(i,j,k,0) < .5){
                     alphap =  std::max(0.,2.*(conc_arr(i,j,k,0)-c_init_1[0])/(c_init_1[1]-c_init_1[0]));
		 } else {
                     alphap =  std::max(0.,(2.-2.*conc_arr(i,j,k,0)-c_init_1[0])/(c_init_1[1]-c_init_1[0]));
		 }
		 if(conc_arr(i,j-1,k,0) < .5){
                     alpham =  std::max(0.,2.*(conc_arr(i,j-1,k,0)-c_init_1[0])/(c_init_1[1]-c_init_1[0]));
		 } else {
                     alpham =  std::max(0.,(2.-2.*conc_arr(i,j-1,k,0)-c_init_1[0])/(c_init_1[1]-c_init_1[0]));
		 }
	     */

                 //djpy(i,j,k) = 0.5*(alpham+alphap)* (djp(i,j,k)-djp(i,j-1,k))/dx[1];
                 djpy(i,j,k) = 0.;
		 /*
  		 if( i == 255){
  		    amrex::Print() << " y grad of disjoining pressure " << j << " " << djpy(i,j,k) << " " << djp(i,j,k) << " " << djp(i,j-1,k) << " " << alphap << " " << alpham << std::endl;
		    amrex::Print() << j << " " << conc_arr(i,j,k,0) << " " << conc_arr(i,j-1,k,0) << std::endl;
		    amrex::Print() << j << " " << rho_arr(i,j,k,0) << " " << rho_arr(i,j-1,k,0) << std::endl;
  		 }
		 */
            });
        
#else
        
            amrex::ParallelFor(bx_x, bx_y, bx_z,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept

            {
                 //djpx(i,j,k) = (djp(i,j,k)-djp(i-1,j,k))/dx[0];
                //amrex::Real alphap = conc_arr(i,  j,k,1);
               // amrex::Real alpham = conc_arr(i-1,j,k,1);
	        amrex::Real threshold = fh_ce;
                amrex::Real alphap = (conc_arr(i,j,k,1)-threshold)/(1.-2.*threshold);
                amrex::Real alpham = (conc_arr(i-1,j,k,1)-threshold)/(1.-2.*threshold);

                alphap = std::max(0.,std::min(1.,alphap));
                alpham = std::max(0.,std::min(1.,alpham));

                amrex::Real avg = 0.5*(alphap+alpham);
//              avg = (avg > 0.03) ? avg : 0.;
                djpx(i,j,k) = avg*F0;
                //if(j == 63 && k == 63){
                //    amrex::Print() << "force " << i << " " << djpx(i,j,k) << std::endl;
                //}

            },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                 // djpy(i,j,k) = (djp(i,j,k)-djp(i,j-1,k))/dx[1];
                 djpy(i,j,k) = 0.;
            },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                 // djpz(i,j,k) = (djp(i,j,k)-djp(i,j,k-1))/dx[2];
                 djpz(i,j,k) = 0.;
            });

/*        
            {
                 djpx(i,j,k) = (djp(i,j,k)-djp(i-1,j,k))/dx[0];
            },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                 djpy(i,j,k) = (djp(i,j,k)-djp(i,j-1,k))/dx[1];
            },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                 djpz(i,j,k) = (djp(i,j,k)-djp(i,j,k-1))/dx[2];
            });
        
*/

#endif

    }

}
