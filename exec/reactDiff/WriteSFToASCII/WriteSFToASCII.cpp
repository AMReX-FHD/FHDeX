#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace std;
using namespace amrex;

int main (int argc, char* argv[]) {

    amrex::Initialize(argc,argv);

    MultiFab mf;

    // Parse command-line arguments
    std::string plt_file; 
    int n_bins; 

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            plt_file= argv[++i]; 
        } 

        if (arg == "-b" && i + 1 < argc) {
            n_bins = std::stoi(argv[++i]);
        } 
    }

    // Append Level_0/Cell for mf read
    plt_file += "/Level_0/Cell";
    std::cout << plt_file << std::endl;
     
    // read in plotfile to MultiFab
    VisMF::Read(mf, plt_file);
    //VisMF::Read(mf, "plt_SF_mag000200000/Level_0/Cell");

    std::vector<std::pair<double, double>> sf_flat;
        
    const int center_x = 32;
    const int center_y = 32;

    for ( MFIter mfi(mf,false); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.validbox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        const Array4<Real>& mfdata = mf.array(mfi);

        // assume everything is 2d for now
        for (auto n=0; n<1; ++n) {
            for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                    for (auto i = lo.x; i <= hi.x; ++i) {
                        auto dk_x = i - center_x;
                        auto dk_y = j - center_y;
                        auto dk = sqrt(dk_x*dk_x + dk_y*dk_y);

                        sf_flat.push_back(std::make_pair(dk, mfdata(i,j,k,n)));
                        //std::cout << dk << ", " << mfdata(i,j,0,n) << "\n";
                    } 
                } 
            } 
        }

    } // end MFIter

    // sort and bin data
    std::sort(sf_flat.begin(), sf_flat.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    std::cout << "... and we are done!\n\n";

}


