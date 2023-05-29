#include "common_functions.H"

#include "gmres_functions.H"
using namespace amrex;
void ApplyMatrix(std::array<MultiFab, AMREX_SPACEDIM> & b_u,
                 MultiFab                             & b_p,
                 std::array<MultiFab, AMREX_SPACEDIM> & x_u,
                 MultiFab                             & x_p,
                 std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                 MultiFab                             & beta,
                 std::array<MultiFab, NUM_EDGE>       & beta_ed,
                 MultiFab                             & gamma,
                 Real                                 & theta_alpha,
                 Geometry                             & geom,
                 int                                  is_inhomogeneous)


{

    std::array<MultiFab, AMREX_SPACEDIM>* b_up = &b_u;
    std::array<MultiFab, AMREX_SPACEDIM>* x_up = &x_u;
    std::array<MultiFab, NUM_EDGE>* beta_edp = &beta_ed;    
    std::array<MultiFab, AMREX_SPACEDIM>* alpha_fcp = &alpha_fc;
         
    MultiFab* b_pp = &b_p;
    MultiFab* x_pp = &x_p;
    MultiFab* gammap = &gamma;   
    MultiFab* betap = &beta;             
    Geometry* geomp = &geom;     

    ApplyMatrix(b_up,b_pp,x_up,x_pp,alpha_fcp,betap,beta_edp,gammap,theta_alpha,geomp,1,is_inhomogeneous);

}


void ApplyMatrix(Vector<std::array<MultiFab, AMREX_SPACEDIM>> & b_u,
                 Vector<MultiFab>                             & b_p,
                 std::array<MultiFab, AMREX_SPACEDIM>* & x_u,
                 MultiFab*                             & x_p,
                 std::array<MultiFab, AMREX_SPACEDIM>* & alpha_fc,
                 MultiFab*                       & beta,
                 std::array<MultiFab, NUM_EDGE>* & beta_ed,
                 MultiFab*                       & gamma,
                 Real                            & theta_alpha,
                 Geometry*                       & geom,
                 int                                    nlevels,
                 int                                    is_inhomogeneous)
{

    std::array<MultiFab, AMREX_SPACEDIM>* b_up = &b_u[0];
    MultiFab* b_pp = &b_p[0];


    ApplyMatrix(b_up,b_pp,x_u,x_p,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom,nlevels,is_inhomogeneous);

}

void ApplyMatrix(Vector<std::array<MultiFab, AMREX_SPACEDIM>> & b_u,
                 Vector<MultiFab>                             & b_p,
                 Vector<std::array<MultiFab, AMREX_SPACEDIM>> & x_u,
                 Vector<MultiFab>                             & x_p,
                 std::array<MultiFab, AMREX_SPACEDIM>* & alpha_fc,
                 MultiFab*                       & beta,
                 std::array<MultiFab, NUM_EDGE>* & beta_ed,
                 MultiFab*                       & gamma,
                 Real                            & theta_alpha,
                 Geometry*                       & geom,
                 int                                    nlevels,
                 int                                    is_inhomogeneous)
{

    std::array<MultiFab, AMREX_SPACEDIM>* b_up = &b_u[0];
    MultiFab* b_pp = &b_p[0];
    
    std::array<MultiFab, AMREX_SPACEDIM>* x_up = &x_u[0];
    MultiFab* x_pp = &x_p[0];


    ApplyMatrix(b_up,b_pp,x_up,x_pp,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom,nlevels,is_inhomogeneous);

}            
                     
// This computes A x = b explicitly
void ApplyMatrix(std::array<MultiFab, AMREX_SPACEDIM>* & b_u,
                 MultiFab*                             & b_p,
                 std::array<MultiFab, AMREX_SPACEDIM>* & x_u,
                 MultiFab*                             & x_p,
                 std::array<MultiFab, AMREX_SPACEDIM>* & alpha_fc,
                 MultiFab*                       & beta,
                 std::array<MultiFab, NUM_EDGE>* & beta_ed,
                 MultiFab*                       & gamma,
                 Real                            & theta_alpha,
                 Geometry*                       & geom,
                 int                                    nlevels,
                 int                                    is_inhomogeneous){

    BL_PROFILE_VAR("ApplyMatrix()", ApplyMatrix);

    // check to make sure x_u and x_p have enough ghost cells
    if (gmres_spatial_order == 2) {
        if (x_u[0][0].nGrow() < 1) {
            Abort("apply_matrix.f90: x_u needs at least 1 ghost cell");
        }
        if (x_p[0].nGrow() < 1) {
            Abort("apply_matrix.f90: x_p needs at least 1 ghost cell");
        } else if (gmres_spatial_order == 4) {
            if (x_u[0][0].nGrow() < 2) {
                Abort("apply_matrix.f90: x_u needs at least 2 ghost cells");
            }
            if (x_p[0].nGrow() < 2) {
                Abort("apply_matrix.f90: x_p needs at least 2 ghost cells");
            }
        }
    }

    // fill ghost cells for x_u and x_p
    // TODO: this should not be done here
    
    for(int lev=0;lev<nlevels;++lev)
    {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            x_u[lev][i].FillBoundary(geom[lev].periodicity());
        }

        x_p[lev].FillBoundary(geom[lev].periodicity());        
    }
    
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        MultiFabPhysBCDomainVel(x_u[0][i], geom[0],i);
        MultiFabPhysBCMacVel(x_u[0][i], geom[0], i, is_inhomogeneous);
    }
    MultiFabPhysBC(x_p[0], geom[0], 0, 1, PRES_BC_COMP);

    // compute b_u = A x_u
    //b_u - gmres_rhs_u
    //x_u - umac
    //
    if(nlevels >0)
    {
        //0->div free, 1->conservative
        //FaceFillCoarse(x_u,0);
        //FaceFillGhost(x_u,0);
    } 
    if (gmres_spatial_order == 2) {
        for(int lev=0;lev<nlevels;++lev)
        {
            StagApplyOp(geom[lev], beta[lev], gamma[lev], beta_ed[lev], x_u[lev], b_u[lev], alpha_fc[lev], geom[lev].CellSize(), theta_alpha);
        }
    }
    else if (gmres_spatial_order == 4) {
        Abort("ApplyMatrix.cpp: gmres_spatial_order=4 not supported yet");
    }

    // compute G x_p and add to b_u
    if (gmres_spatial_order == 2) {
        for(int lev=0;lev<nlevels;++lev)
        {
            int bcc = (lev==0)? PRES_BC_COMP : 100;
            ComputeGrad(x_p[lev], b_u[lev], 0, 0, 1, bcc, geom[lev], 1);
        }
    }
    else if (gmres_spatial_order == 4) {
        Abort("ApplyMatrix.cpp: gmres_spatial_order=4 not supported yet");
    }

    // set b_p = -D x_u
    for(int lev=0;lev<nlevels;++lev)
    {
        ComputeDiv(b_p[lev], x_u[lev], 0, 0, 1, geom[lev], 0);
        b_p[lev].mult(-1., 0, 1, 0);
    }
}





