#ifndef _immbdy_namespace_H_
#define _immbdy_namespace_H_


namespace immbdy {

    extern int n_immbdy;

    extern bool contains_flagellum;
    extern bool contains_fourier;
}



namespace ib_flagellum {

    extern amrex::Vector<int>             n_marker;
    extern amrex::Vector<amrex::RealVect> offset_0;
    extern amrex::Vector<amrex::Real>     amplitude;
    extern amrex::Vector<amrex::Real>     frequency;
    extern amrex::Vector<amrex::Real>     length;
    extern amrex::Vector<amrex::Real>     wavelength;
    extern amrex::Vector<amrex::Real>     k_spring;
    extern amrex::Vector<amrex::Real>     k_driving;

    extern int fourier_coef_len;
    extern std::string chlamy_flagellum_datafile;
}


namespace chlamy_flagellum {
    extern amrex::Vector<amrex::Vector<int>> N;
    extern amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> a_coef;
    extern amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> b_coef;
}


#endif
