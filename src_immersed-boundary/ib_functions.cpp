#include <ib_functions.H>
#include <ib_functions_F.H>
#include <immbdy_namespace.H>
#include <immbdy_namespace_declarations.H>


using namespace immbdy;
using namespace ib_flagellum;


void InitializeImmbdyNamespace() {

    n_immbdy = 0;
    contains_flagellum = false;

    int cf = 0, cfourier = 0;
    initialize_immbdy_namespace(& n_immbdy, & cf, & cfourier);
    if (cf == 1) contains_flagellum = true;
    if (cfourier == 1) contains_fourier = true;

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
        chlamy_flagellum::a_coef.resize(n_immbdy);
        chlamy_flagellum::b_coef.resize(n_immbdy);

        int max_n_markers = 0;
        flagellum_max_markers(& max_n_markers);

        for (int i=0; i<n_immbdy; ++i) {
            // NOTE: "last" marker contains no fourier mode
            chlamy_flagellum::a_coef[i].resize(max_n_markers-1);
            chlamy_flagellum::b_coef[i].resize(max_n_markers-1);

            for(int j=0; j<max_n_markers-1; ++j) {
                chlamy_flagellum::a_coef[i][j].resize(ib_flagellum::fourier_coef_len);
                chlamy_flagellum::b_coef[i][j].resize(ib_flagellum::fourier_coef_len);

                copy_ib_fourier_data(i+1, j+1, // Don't forget: fortran array indices start at 1
                                     chlamy_flagellum::a_coef[i][j].dataPtr(),
                                     chlamy_flagellum::b_coef[i][j].dataPtr());

            }
        }


        Print() << "Loaded Fourier Data:" << std::endl;
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
