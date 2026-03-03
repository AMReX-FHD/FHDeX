#include "FPU.H"

#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void ComputePhiFromState(MultiFab& phi,
                         const Real& r_eq,
                         const Real& p_eq,
                         const Real& e_eq,
                         const Real& R_00,
                         const Real& R_01,
                         const Real& R_02,
                         const Real& R_10,
                         const Real& R_11,
                         const Real& R_12,
                         const Real& R_20,
                         const Real& R_21,
                         const Real& R_22) {

    // convert state to state-state_eq
    phi.plus(-r_eq,0,1,0);
    phi.plus(-p_eq,1,1,0);
    phi.plus(-e_eq,2,1,0);

    // phi = R(state-state_eq)
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real>& phi_fab = phi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real R0 = R_00*phi_fab(i,j,k,0) + R_01*phi_fab(i,j,k,1) + R_02*phi_fab(i,j,k,2);
            Real R1 = R_10*phi_fab(i,j,k,0) + R_11*phi_fab(i,j,k,1) + R_12*phi_fab(i,j,k,2);
            Real R2 = R_20*phi_fab(i,j,k,0) + R_21*phi_fab(i,j,k,1) + R_22*phi_fab(i,j,k,2);

            phi_fab(i,j,k,0) = R0;
            phi_fab(i,j,k,1) = R1;
            phi_fab(i,j,k,2) = R2;

        });
    }

}

void ComputeCalphaalpha(MultiFab& C_alphaalpha,
                        const MultiFab& phi,
                        const MultiFab& phi0) {

    for (MFIter mfi(C_alphaalpha); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<      Real>& C_alphaalpha_fab = C_alphaalpha.array(mfi);
        const Array4<const Real>& phi_fab          = phi.array(mfi);
        const Array4<const Real>& phi0_fab         = phi0.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            C_alphaalpha_fab(i,j,k,0) = phi_fab(i,j,k,0) * phi0_fab(0,j,k,0);
            C_alphaalpha_fab(i,j,k,1) = phi_fab(i,j,k,1) * phi0_fab(0,j,k,1);
            C_alphaalpha_fab(i,j,k,2) = phi_fab(i,j,k,2) * phi0_fab(0,j,k,2);
        });
    }

}


