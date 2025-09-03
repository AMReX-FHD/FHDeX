#ifndef _AmrCoreAdv_F_H_
#define _AmrCoreAdv_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    // these are interfaces to fortran subroutines

    void initdata(const int* level, const amrex_real* time,
        const int* lo, const int* hi,
        BL_FORT_FAB_ARG_3D(state),
        const amrex_real* dx, const amrex_real* problo);


    void state_error(int* tag, const int* tag_lo, const int* tag_hi,
        const BL_FORT_FAB_ARG_3D(state),
        const int* tagval, const int* clearval,
        const int* lo, const int* hi,
        const amrex_real* dx, const amrex_real* problo,
        const amrex_real* time, const amrex_real* conerr);

    void advect_2d(const amrex_real* time, const int* lo, const int* hi,
        const BL_FORT_FAB_ARG_3D(statein),
        BL_FORT_FAB_ARG_3D(stateout),
        BL_FORT_FAB_ARG_3D(ptS),
        int* iface, const int* if_lo, const int* if_hi,
        AMREX_D_DECL(const BL_FORT_FAB_ARG_3D(xvel),
           const BL_FORT_FAB_ARG_3D(yvel),
           const BL_FORT_FAB_ARG_3D(zvel)),
        AMREX_D_DECL(BL_FORT_FAB_ARG_3D(fx),
           BL_FORT_FAB_ARG_3D(fy),
           BL_FORT_FAB_ARG_3D(fz)),
        const amrex_real* dx, const amrex_real* dt, const amrex_real* nu);

    void advect_3d(const amrex_real* time, const int* lo, const int* hi,
        const BL_FORT_FAB_ARG_3D(statein_p),
        const BL_FORT_FAB_ARG_3D(statein_f),
        BL_FORT_FAB_ARG_3D(stateout),
        BL_FORT_FAB_ARG_3D(ptSp),
        BL_FORT_FAB_ARG_3D(ptSf),
        int* ifacep, const int* ifp_lo, const int* ifp_hi,
        int* ifacef, const int* iff_lo, const int* iff_hi,
        AMREX_D_DECL(const BL_FORT_FAB_ARG_3D(xvel_p),
               const BL_FORT_FAB_ARG_3D(yvel_p),
               const BL_FORT_FAB_ARG_3D(zvel_p)),
        AMREX_D_DECL(const BL_FORT_FAB_ARG_3D(xvel_f),
               const BL_FORT_FAB_ARG_3D(yvel_f),
               const BL_FORT_FAB_ARG_3D(zvel_f)),
        AMREX_D_DECL(BL_FORT_FAB_ARG_3D(fx),
               BL_FORT_FAB_ARG_3D(fy),
               BL_FORT_FAB_ARG_3D(fz)),
        AMREX_D_DECL(BL_FORT_FAB_ARG_3D(fx1),
               BL_FORT_FAB_ARG_3D(fy1),
               BL_FORT_FAB_ARG_3D(fz1)),
        AMREX_D_DECL(BL_FORT_FAB_ARG_3D(fx2),
               BL_FORT_FAB_ARG_3D(fy2),
               BL_FORT_FAB_ARG_3D(fz2)),
        const amrex_real* dx, const amrex_real* dt,
        const amrex_real* diffcoeff, const int* correct);

    void get_congrad_3d( const int* lo, const int* hi,
        const BL_FORT_FAB_ARG_3D(con),
        BL_FORT_FAB_ARG_3D(con_x),
        BL_FORT_FAB_ARG_3D(con_y),
        BL_FORT_FAB_ARG_3D(con_z),
        BL_FORT_FAB_ARG_3D(xface),
        BL_FORT_FAB_ARG_3D(yface),
        BL_FORT_FAB_ARG_3D(zface),
        int* iface, const int* if_lo, const int* if_hi,
        const amrex_real* ib_cen_x, const amrex_real* ib_cen_y,
        const amrex_real* ib_cen_z,
        const amrex_real* dx, const amrex_real* prob_lo);

    void get_surfgrad_3d( const int* lo, const int* hi,
        const BL_FORT_FAB_ARG_3D(con),
        BL_FORT_FAB_ARG_3D(cons_x),
        BL_FORT_FAB_ARG_3D(cons_y),
        BL_FORT_FAB_ARG_3D(cons_z),
        BL_FORT_FAB_ARG_3D(MagDCon),
        const BL_FORT_FAB_ARG_3D(ls),
        int* iface, const int* if_lo, const int* if_hi,
        const amrex_real* dx, const amrex_real* prob_lo);


   void get_ptsource_2d( const int* lo, const int* hi,
       int* ctag, const int* c_lo, const int* c_hi,
        BL_FORT_FAB_ARG_3D(ptS),
        const amrex_real* strength,const amrex_real* dx,
        const amrex_real* prob_lo, int* Num_Source);

    void get_ptsource_3d( const int* lo, const int* hi,
        int * ctag,const int* c_lo, const int* c_hi,
        BL_FORT_FAB_ARG_3D(ptS),
        const amrex_real* strength, const amrex_real* dx,
        const amrex_real* prob_lo, int* Num_source);


     void confill (BL_FORT_FAB_ARG_3D(state),
        const int* dlo, const int* dhi,
        const amrex_real* dx, const amrex_real* glo,
        const amrex_real* time, const int* bc);


#ifdef __cplusplus
}
#endif

#endif