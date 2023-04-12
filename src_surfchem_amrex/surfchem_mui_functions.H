#ifndef surfchem_mui_functions_H_
#define surfchem_mui_functions_H_

#include <AMReX_MPMD.H>
#include <AMReX_MultiFab.H>

namespace surfchem_mui {

    extern int nspec_mui;

}

void InitializeSurfChemMUINamespace();

void amrex_fetch_Ntot(amrex::MultiFab& Ntot, amrex::MPMD::Copier const& copier);

void amrex_fetch_surfcov(amrex::MultiFab const& Ntot, amrex::MultiFab& surfcov,
                         amrex::MPMD::Copier const& copier);

void amrex_push(amrex::MultiFab const& cu, amrex::MultiFab const& prim,
                amrex::MPMD::Copier const& copier);

void amrex_fetch(amrex::MultiFab& cu, amrex::MultiFab const& prim,
                 amrex::GpuArray<amrex::Real,3> const& dx,
                 amrex::MPMD::Copier const& copier);

#endif
