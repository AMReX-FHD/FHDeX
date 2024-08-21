#ifndef MYFUNC_H_
#define MYFUNC_H_

void main_main(const char* argv);

void EMstep_chem_only(amrex::MultiFab& rho_old, amrex::MultiFab& rho_new,
                      const amrex::Geometry geom, const amrex::Real dt);

void RK3step_chem_only(amrex::MultiFab& rho_old, amrex::MultiFab& rho_new,
                       const amrex::Geometry geom, const amrex::Real dt);

void compute_chemistry_source_CLE_1(amrex::Real dt, amrex::Real dV, MultiFab& mf_in, int startComp_in,
                                    MultiFab& source, int startComp_out);

void compute_chemistry_source_CLE_2(amrex::Real dt, amrex::Real dV, MultiFab& mf_in, int startComp_in,
                                    MultiFab& source, int startComp_out, MultiFab& ranchem);
#endif
