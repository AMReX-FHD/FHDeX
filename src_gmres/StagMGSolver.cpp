#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

#include "common_functions.H"
#include "common_namespace.H"

using namespace amrex;
using namespace gmres;
using namespace common;

// solve "(theta*alpha*I - L) phi = rhs" using multigrid with Jacobi relaxation
// if abs(visc_type) = 1, L = div beta grad
// if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
// if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
// if visc_type > 1 we assume constant coefficients
// if visc_type < 1 we assume variable coefficients
// beta_cc, and gamma_cc are cell-centered
// alpha_fc, phi_fc, and rhs_fc are face-centered
// beta_ed is nodal (2d) or edge-centered (3d)
// phi_fc must come in initialized to some value, preferably a reasonable guess
void StagMGSolver(const std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
                  const MultiFab& beta_cc,
#if (AMREX_SPACEDIM == 2)
                  const std::array< MultiFab, 1 >& beta_ed,
#elif (AMREX_SPACEDIM == 3)
                  const std::array< MultiFab, 3 >& beta_ed,
#endif
                  const MultiFab& gamma_cc,
                  const MultiFab& phi_fc,
                  const std::array< MultiFab, AMREX_SPACEDIM >& rhs_fc,
                  const Real& theta,
                  const Geometry& geom)
{

    // get the problem domain and boxarray at level 0
    Box pd_base = geom.Domain();
    BoxArray ba_base = beta_cc.boxArray();

    // compute the number of multigrid levels assuming stag_mg_minwidth is the length of the
    // smallest dimension of the smallest grid at the coarsest multigrid level
    int nlevs_mg = ComputeNlevsMG(ba_base);

    // allocate multifabs used in multigrid coarsening

    // cell-centered
    Vector<MultiFab>  beta_cc_mg(nlevs_mg);
    Vector<MultiFab> gamma_cc_mg(nlevs_mg);

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > alpha_fc_mg(nlevs_mg);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >   rhs_fc_mg(nlevs_mg);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >   phi_fc_mg(nlevs_mg);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >  Lphi_fc_mg(nlevs_mg);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > resid_fc_mg(nlevs_mg);

    const int num_beta_fc = ( AMREX_SPACEDIM == 2 ) ? 1 : 3;
    Vector<std::array< MultiFab, num_beta_fc > >  beta_ed_mg(nlevs_mg); // nodal

    const Real* dx = geom.CellSize();

    Vector<std::array< Real, AMREX_SPACEDIM > > dx_mg;

    DistributionMapping dmap = beta_cc.DistributionMap();

    for (int n=0; n<nlevs_mg; ++n) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            // compute dx at this level of multigrid
            dx_mg[n][d] = dx[d] * pow(2,n);
        }

        // create the problem domain for this multigrid level
        Box pd = pd_base.coarsen(pow(2,n));

        // create the boxarray for this multigrid level
        BoxArray ba(ba_base); 
        ba.coarsen(pow(2,n));
     
        if ( n == 0 && !(ba == ba_base) ) {
            Abort("Finest multigrid level boxarray and coarsest problem boxarrays do not match");
        }

        // build multifabs used in multigrid coarsening
         beta_cc_mg[n].define(ba,dmap,1,1);
        gamma_cc_mg[n].define(ba,dmap,1,1);

        AMREX_D_TERM(alpha_fc_mg[n][0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                     alpha_fc_mg[n][1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                     alpha_fc_mg[n][2].define(convert(ba,nodal_flag_z), dmap, 1, 0););
        AMREX_D_TERM(  rhs_fc_mg[n][0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                       rhs_fc_mg[n][1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                       rhs_fc_mg[n][2].define(convert(ba,nodal_flag_z), dmap, 1, 0););
        AMREX_D_TERM(  phi_fc_mg[n][0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                       phi_fc_mg[n][1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                       phi_fc_mg[n][2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
        AMREX_D_TERM( Lphi_fc_mg[n][0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                      Lphi_fc_mg[n][1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                      Lphi_fc_mg[n][2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
        AMREX_D_TERM(resid_fc_mg[n][0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                     resid_fc_mg[n][1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                     resid_fc_mg[n][2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

        // build beta_ed_mg
        if (AMREX_SPACEDIM == 2) {
            beta_ed_mg[n][0].define(convert(ba,nodal_flag), dmap, 1, 0);
        }
        else if (AMREX_SPACEDIM == 3) {
            beta_ed_mg[n][0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
            beta_ed_mg[n][1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
            beta_ed_mg[n][2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
        }
    } // end loop over multigrid levels

    // copy level 1 coefficients into mg array of coefficients
    MultiFab::Copy(beta_cc_mg[0],beta_cc,0,0,1,1);
    MultiFab::Copy(gamma_cc_mg[0],gamma_cc,0,0,1,1);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(alpha_fc_mg[0][d],alpha_fc[d],0,0,1,0);
        // multiply alpha_fc_mg by theta
        alpha_fc_mg[0][d].mult(theta,0,1,0);
    }
    MultiFab::Copy(beta_ed_mg[0][0],beta_ed[0],0,0,1,0);
    if (AMREX_SPACEDIM == 3) {
        MultiFab::Copy(beta_ed_mg[0][1],beta_ed[1],0,0,1,0);
        MultiFab::Copy(beta_ed_mg[0][2],beta_ed[2],0,0,1,0);
    }
    
    // coarsen coefficients
    for (int n=1; n<nlevs_mg; ++n) {
        // need ghost cells set to zero to prevent intermediate NaN states
        // that cause some compilers to fail
        beta_cc_mg[n].setVal(0.);
        gamma_cc_mg[n].setVal(0.);

        // cc_restriction on beta_cc_mg

        // cc_restriction on gamma_cc_mg

        // stag_restriction on alpha_fc_mg

        if (AMREX_SPACEDIM == 2) {
            // nodal_restriction on beta_ed_mg

        }
        else {
            // edge_restriction on beta_ed_mg

        }
    }

    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Now we wolve the homogeneous problem
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

}

// compute the number of multigrid levels assuming minwidth is the length of the
// smallest dimension of the smallest grid at the coarsest multigrid level
int ComputeNlevsMG(const BoxArray& ba) {

    int nlevs_mg = -1;

    for (int i=0; i<ba.size(); ++i) {
        Box bx = ba.get(i);
        IntVect iv = bx.bigEnd() - bx.smallEnd() + IntVect(1);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            int temp = iv[d];
            int rdir = 1;
            while (temp%2 == 0 && temp/stag_mg_minwidth != 1) {
                temp /= 2;
                ++rdir;
            }

            if (nlevs_mg == -1) {
                nlevs_mg = rdir;
            }
            else {
                nlevs_mg = std::min(rdir,nlevs_mg);
            }
        }
    }

    return nlevs_mg;
}

void CCRestriction(MultiFab& phi_c, const MultiFab& phi_f)
{
    // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for ( MFIter mfi(phi_c); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        const Box& validBox = mfi.validbox();

        cc_restriction(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                       BL_TO_FORTRAN_3D(phi_c[mfi]),
                       BL_TO_FORTRAN_3D(phi_f[mfi]));
    }
}

void FaceRestriction(std::array< MultiFab, AMREX_SPACEDIM >& phi_c, 
                     const std::array< MultiFab, AMREX_SPACEDIM >& phi_f,
                     int simple_stencil)
{

    // loop over boxes (note we are not passing in a cell-centered MultiFab)
    for ( MFIter mfi(phi_c[0]); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        // there are no cell-centered MultiFabs so use this to get
        // a cell-centered box
        const Box& validBox = amrex::enclosedCells(mfi.validbox());

        face_restriction(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                         BL_TO_FORTRAN_3D(phi_c[0][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[0][mfi]),
                         BL_TO_FORTRAN_3D(phi_c[1][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(phi_c[2][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[2][mfi]),
#endif
                         &simple_stencil);
    }
}

void EdgeRestriction(std::array< MultiFab, 3 >& phi_c, 
                      const std::array< MultiFab, 3 >& phi_f)
{
    if (AMREX_SPACEDIM != 3) {
        Abort("Edge restriction can only be called for 3D!");
    }

    // loop over boxes (note we are not passing in a cell-centered MultiFab)
    for ( MFIter mfi(phi_c[0]); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        // there are no cell-centered MultiFabs so use this to get
        // a cell-centered box
        const Box& validBox = amrex::enclosedCells(mfi.validbox());

        edge_restriction(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                         BL_TO_FORTRAN_3D(phi_c[0][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[0][mfi]),
                         BL_TO_FORTRAN_3D(phi_c[1][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[1][mfi]),
                         BL_TO_FORTRAN_3D(phi_c[2][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[2][mfi]));
    }
}

void NodalRestriction(MultiFab& phi_c, const MultiFab& phi_f)
{
    // loop over boxes (note we are not passing in a cell-centered MultiFab)
    for ( MFIter mfi(phi_c); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        // there are no cell-centered MultiFabs so use this to get
        // a cell-centered box
        const Box& validBox = amrex::enclosedCells(mfi.validbox());

        nodal_restriction(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                          BL_TO_FORTRAN_3D(phi_c[mfi]),
                          BL_TO_FORTRAN_3D(phi_f[mfi]));

    }
}
