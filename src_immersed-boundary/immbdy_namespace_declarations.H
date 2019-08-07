#ifndef _immbdy_namespace_declarations_H_
#define _immbdy_namespace_declarations_H_


int  immbdy::n_immbdy;
bool immbdy::contains_flagellum;

// Flagellum
amrex::Vector<int>             ib_flagellum::n_marker;
amrex::Vector<amrex::RealVect> ib_flagellum::offset_0;
amrex::Vector<amrex::Real>     ib_flagellum::amplitude;
amrex::Vector<amrex::Real>     ib_flagellum::frequency;
amrex::Vector<amrex::Real>     ib_flagellum::length;
amrex::Vector<amrex::Real>     ib_flagellum::wavelength;


#endif
