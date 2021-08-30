#include <ib_functions.H>
#include <ib_functions_F.H>
#include <immbdy_namespace.H>

using namespace immbdy;
using namespace ib_flagellum;

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


void InitializeImmbdyNamespace() {

    n_immbdy = 0;
    contains_flagellum = false;
    contains_colloid   = false;

    int cf = 0, cfourier = 0, ccolloid = 0;
    initialize_immbdy_namespace(& n_immbdy, & cf, & cfourier, & ccolloid);
    if (cf == 1) contains_flagellum = true;
    if (cfourier == 1) contains_fourier = true;
    if (ccolloid == 1) contains_colloid = true;
}


void InitializeIBFlagellumNamespace() {

    ib_flagellum::n_marker.resize(n_immbdy);
    ib_flagellum::offset_0.resize(n_immbdy);
    ib_flagellum::amplitude.resize(n_immbdy);
    ib_flagellum::frequency.resize(n_immbdy);
    ib_flagellum::length.resize(n_immbdy);
    ib_flagellum::wavelength.resize(n_immbdy);
    ib_flagellum::k_spring.resize(n_immbdy);
    ib_flagellum::k_driving.resize(n_immbdy);

    int len_buffer = 0;
    chlamy_flagellum_datafile_len(& len_buffer);
    char chlamy_flagellum_datafile[len_buffer+1];

    initialize_ib_flagellum_namespace(n_immbdy,
                                      ib_flagellum::n_marker.dataPtr(),
                                      ib_flagellum::offset_0.dataPtr(),
                                      ib_flagellum::amplitude.dataPtr(),
                                      ib_flagellum::frequency.dataPtr(),
                                      ib_flagellum::length.dataPtr(),
                                      ib_flagellum::wavelength.dataPtr(),
                                      ib_flagellum::k_spring.dataPtr(),
                                      ib_flagellum::k_driving.dataPtr(),
                                      & ib_flagellum::fourier_coef_len,
                                      chlamy_flagellum_datafile);

    ib_flagellum::chlamy_flagellum_datafile = chlamy_flagellum_datafile;

    if (contains_fourier) {
        chlamy_flagellum::N.resize(n_immbdy);
        chlamy_flagellum::a_coef.resize(n_immbdy);
        chlamy_flagellum::b_coef.resize(n_immbdy);

        int max_n_markers = 0;
        flagellum_max_markers(& max_n_markers);

        for (int i=0; i<n_immbdy; ++i) {
            // NOTE: "last" marker contains no fourier mode
            chlamy_flagellum::N[i].resize(max_n_markers-1);
            chlamy_flagellum::a_coef[i].resize(max_n_markers-1);
            chlamy_flagellum::b_coef[i].resize(max_n_markers-1);

            for(int j=0; j<max_n_markers-1; ++j) {
                chlamy_flagellum::a_coef[i][j].resize(ib_flagellum::fourier_coef_len);
                chlamy_flagellum::b_coef[i][j].resize(ib_flagellum::fourier_coef_len);

                copy_ib_fourier_data(i+1, j+1, // Don't forget: fortran array indices start at 1
                                     & chlamy_flagellum::N[i][j],
                                     chlamy_flagellum::a_coef[i][j].dataPtr(),
                                     chlamy_flagellum::b_coef[i][j].dataPtr());

            }
        }


        Print() << "Loaded Fourier Data:" << std::endl;

        for (int i=0; i<n_immbdy; ++i) {
            Print() << "N " << i << " :";
            for (int elt : chlamy_flagellum::N[i])
                Print() << " " << elt ;
            Print() << std::endl;
        }

        Print() << "a_coef:" << std::endl;
        for (int i=0; i<n_immbdy; ++i) {
            for (int j=0; j<max_n_markers-1; ++j) {
                Print() << i << ", " << j << ":";
                for (Real elt : chlamy_flagellum::a_coef[i][j])
                    Print() << " " << elt ;
                Print() << std::endl;
            }
        }

        Print() << "b_coef:" << std::endl;
        for (int i=0; i<n_immbdy; ++i) {
            for (int j=0; j<max_n_markers-1; ++j) {
                Print() << i << ", " << j << ":";
                for (Real elt : chlamy_flagellum::b_coef[i][j])
                    Print() << " " << elt ;
                Print() << std::endl;
            }
        }


    }
}



void InitializeIBColloidNamespace() {

    ib_colloid::n_marker.resize(n_immbdy);
    ib_colloid::center.resize(n_immbdy);
    ib_colloid::radius.resize(n_immbdy);
    ib_colloid::rho.resize(n_immbdy);
    ib_colloid::k_spring.resize(n_immbdy);

    initialize_ib_colloid_namespace(n_immbdy,
                                    ib_colloid::n_marker.dataPtr(),
                                    ib_colloid::center.dataPtr(),
                                    ib_colloid::radius.dataPtr(),
                                    ib_colloid::rho.dataPtr(),
                                    ib_colloid::k_spring.dataPtr());

}
