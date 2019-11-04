#ifndef _immbdy_namespace_declarations_H_
#define _immbdy_namespace_declarations_H_


int  immbdy::n_immbdy;
bool immbdy::contains_flagellum;
bool immbdy::contains_fourier;
bool immbdy::contains_colloid;


// Flagellum
amrex::Vector<int>             ib_flagellum::n_marker;
amrex::Vector<amrex::RealVect> ib_flagellum::offset_0;
amrex::Vector<amrex::Real>     ib_flagellum::amplitude;
amrex::Vector<amrex::Real>     ib_flagellum::frequency;
amrex::Vector<amrex::Real>     ib_flagellum::length;
amrex::Vector<amrex::Real>     ib_flagellum::wavelength;
amrex::Vector<amrex::Real>     ib_flagellum::k_spring;
amrex::Vector<amrex::Real>     ib_flagellum::k_driving;


// Chlamy
int         ib_flagellum::fourier_coef_len;
std::string ib_flagellum::chlamy_flagellum_datafile;

amrex::Vector<amrex::Vector<int>>                        chlamy_flagellum::N;
amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> chlamy_flagellum::a_coef;
amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> chlamy_flagellum::b_coef;


// Colloid
amrex::Vector<int>             ib_colloid::n_marker;
amrex::Vector<amrex::RealVect> ib_colloid::center;
amrex::Vector<amrex::Real>     ib_colloid::radius;
amrex::Vector<amrex::Real>     ib_colloid::rho;
amrex::Vector<amrex::Real>     ib_colloid::k_spring;

#endif
