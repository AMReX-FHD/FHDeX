#include "multispec_functions.H"

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
            node_grad_c(i,j,k,2) = (c(i,j,k)-c(i,j-1,k)+c(i-1,j,k)-c(i-1,j-1,k)
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
