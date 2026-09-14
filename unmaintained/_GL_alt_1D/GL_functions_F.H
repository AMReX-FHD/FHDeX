#ifndef _GL_functions_F_H_
#define GL_functions_F_H_

#include <AMReX_MultiFab.H>

#ifdef __cplusplus
extern "C"
{
#endif



  void set_inputs(const int* Plot_Skip_out,const int* Number_of_Samples_out, const int* Equil,amrex::Real* alpha,amrex::Real* r1,amrex::Real* r2, const int* adaptive, const int* reverse);
//  void rk2_stage1(const int* lo, const int* hi,
//                amrex::Real* phi,
//                amrex::Real* phin,
//                const amrex::Real* rannums,
//                amrex::Real* integral,
//                amrex::Real* energy,
//                amrex::Real* teng,
//                amrex::Real* H1_semi_norm,
//                const amrex::Real* dx,
//                const amrex::Real* dt,
//                amrex::Real* phi_avg);




  void rk2_stage2(const int* lo, const int* hi,
                amrex::Real* phi,
                amrex::Real* phin,
                const amrex::Real* rannums,
                const amrex::Real* integral,
                const amrex::Real* dx,
                const amrex::Real* dt);

  void read_GL_namelist( const char* inputs_file, const int length);

  void integrate(const int* lo, const int* hi,
                amrex::Real* phi,
                const amrex::Real* dx,
                amrex::Real* integral);

  void  setdt ( const amrex::Real* dx, amrex::Real* dt );

  void  inc_phi0 ( const int* step );

  void  fixed_inc_phi0 ( const int* forward_input );

  void  inc_phi0_Adapt ( amrex::Real* Expec,amrex::Real* MAD,amrex::Real* r1, int* Shift_Flag);


  void  initphi (const int* lo, const int* hi,
                amrex::Real* phi,
                const amrex::Real* dx);

  void Stat_Quant(const int* lo, const int* hi,
                amrex::Real* phi,
                amrex::Real* phi_avg,
                amrex::Real* phi_sq_avg);


  void Umbrella_Adjust ( int* sucessful_iter, amrex::Real* alpha, int* umbrella_size , int* sucessful_iter_prev);

  void Param_Output (amrex::Real*  umbrella_output, amrex::Real*  phi0_output);
  void umbrella_reset( amrex::Real*  umbrella_input );



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Stat_Quant_1D(const int* lo, const int* hi,
                amrex::Real* phi,
                amrex::Real* phi_avg,
                amrex::Real* phi_sq_avg);

  void integrate_1D(const int* lo, const int* hi,
                amrex::Real* phi,
                const amrex::Real* dx,
                amrex::Real* integral);

  void Comp_H1_semi_norm_1D(const int* lo, const int* hi,
                amrex::Real* phi,
                const amrex::Real* dx,
                amrex::Real* H1_semi);

  void  initphi_1D (const int* lo, const int* hi,
                amrex::Real* phi,
                const amrex::Real* dx);

  void rk2_stage1_1D(const int* lo, const int* hi,
                amrex::Real* phi,
                amrex::Real* phin,
                const amrex::Real* rannums,
                amrex::Real* integral,
                amrex::Real* energy,
                amrex::Real* teng,
                amrex::Real* H1_semi_norm,
                const amrex::Real* dx,
                const amrex::Real* dt,
                amrex::Real* phi_avg);

  void  setdt_1D ( const amrex::Real* dx, amrex::Real* dt );


#ifdef __cplusplus
}
#endif

#endif