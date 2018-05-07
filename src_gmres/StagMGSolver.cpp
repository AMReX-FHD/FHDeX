#include "gmres_functions.H"
#include "gmres_namespace.H"

using namespace amrex;
using namespace gmres;

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
void StagMGSolver(const std::array<MultiFab, AMREX_SPACEDIM>& alpha_fc,
                  const MultiFab& beta_cc,
#if (AMREX_SPACEDIM == 2)
                  const std::array<MultiFab, 1>& beta_ed,
#elif (AMREX_SPACEDIM == 3)
                  const std::array<MultiFab, 3>& beta_ed,
#endif
                  const MultiFab& gamma_cc,
                  const Real& theta,
                  const Geometry& geom)
{

    // get the problem domain and boxarray at level 0
    Box pd_base = geom.Domain();
    BoxArray ba_base = beta_cc.boxArray();

    int nlevs_mg = ComputeNlevsMG(ba_base);
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
