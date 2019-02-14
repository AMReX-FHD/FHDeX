#include "hydro_functions.H"
#include "hydro_functions_F.H"

void setBC(amrex::MultiFab & u_mac, amrex::MultiFab & v_mac, amrex::MultiFab & w_mac) {

    // Loop over boxes

    for ( amrex::MFIter mfi(u_mac); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();

        set_bc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
               u_mac[mfi].dataPtr(),
               v_mac[mfi].dataPtr(),
               w_mac[mfi].dataPtr());

    }

}
