#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

using namespace amrex;


void RK2step(MultiFab& phi, MultiFab& phin, MultiFab& rannums, const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt, amrex::Real& integral, int n, amrex::Real& phi_avg); 

void Init_Phi(MultiFab& phi, const amrex::Real* dx);

void Run_Steps(MultiFab& phi, MultiFab& phin, MultiFab& rannums, const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt, amrex::Real& integral,int step,
                        amrex::Real& time, int plot_int,int n_steps_skip, bool Make_PltFiles,int N_Burn,int L, amrex::Real& Expec, amrex::Real& MAD, int max_step, int& Plot_Num,int& Plot_Skip, int& umbrella_number);
	
void Check_Overlap(amrex::Real& Expec,amrex::Real& MAD,amrex::Real& Expec2,amrex::Real& MAD2,amrex::Real& r2,amrex::Real& alpha, bool& sucessful_compare,int& umbrella_size, int& Shift_Flag, bool& while_loop_comp, bool& First_Loop_Step, bool& weak_phi);

void Check_Overlap_Backwards(amrex::Real& Expec,amrex::Real& MAD,amrex::Real& Expec2,amrex::Real& MAD2,amrex::Real& r2,amrex::Real& alpha, bool& sucessful_compare,int& umbrella_size, int& Shift_Flag, bool& while_loop_comp, bool& First_Loop_Step, bool& weak_phi);
