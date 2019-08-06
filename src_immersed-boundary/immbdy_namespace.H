#ifndef _immbdy_namespace_H_
#define _immbdy_namespace_H_


namespace immbdy {

    extern int n_immbdy;

    extern bool contains_flagellum;
}



namespace ib_flagellum {

    extern amrex::Vector<int>             n_marker;
    extern amrex::Vector<amrex::RealVect> offset_0;
    extern amrex::Vector<amrex::Real>     amplitude;
    extern amrex::Vector<amrex::Real>     frequency;
    extern amrex::Vector<amrex::Real>     length;
    extern amrex::Vector<amrex::Real>     wavelength;

}


#endif
