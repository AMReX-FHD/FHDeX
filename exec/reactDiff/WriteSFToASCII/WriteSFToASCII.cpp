#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace std;
using namespace amrex;


/*
 *   Inputs:
 *      -i  input file (this needs to be given)
 *      -b  number of bins (will default if not given)
 */  
int main (int argc, char* argv[]) {

    amrex::Initialize(argc,argv);

    MultiFab mf;

    // Parse command-line arguments
    std::string plt_file; 
    std::string out_name;
    int n_bins; 

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            plt_file = argv[++i]; 
        } 

        if (arg == "-b" && i + 1 < argc) {
            n_bins = std::stoi(argv[++i]);
        }
        
        if (arg == "-o" && i + 1 < argc) {
            out_name = argv[++i];
        }
    }

    if (plt_file == "")
    {
        std::cerr << "No input file given!\n";
        return 1;
    }

    // if no output file name is given, use the plt file name and append with _sf_raw.csv
    if (out_name == "") {
        out_name += plt_file + "_sf_raw.csv";
    } 

    // Append Level_0/Cell for mf read
    plt_file += "/Level_0/Cell";
    std::cout << plt_file << std::endl;
     
    // read in plotfile to MultiFab
    VisMF::Read(mf, plt_file);

    std::vector<std::pair<double, double>> sf_flat;
       
    // TODO: need to figure out how to know the center from the level_0 data 
    const int center_x = 32;
    const int center_y = 32;

    for ( MFIter mfi(mf,false); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.validbox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        const Array4<Real>& mfdata = mf.array(mfi);

        for (auto n=0; n<1; ++n) {                          // only look at the first component 
            for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                    for (auto i = lo.x; i <= hi.x; ++i) {
                        auto dk_x = i - center_x;
                        auto dk_y = j - center_y;
                        auto dk = sqrt(dk_x*dk_x + dk_y*dk_y);

                        sf_flat.push_back(std::make_pair(dk, mfdata(i,j,k,n)));
                    } 
                } 
            } 
        }

    } // end MFIter

    // sort 
    std::sort(sf_flat.begin(), sf_flat.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    // if no bins specified, print all the data out
    if (n_bins == 0)
    {
        std::ofstream out_file(out_name);

        if (!out_file.is_open()) {
            std::cerr << "Error: Could not open file " << out_name << " for writing." << std::endl;
            return 1; 
        }

        out_file << "dk, sf\n";

        for (const auto& data : sf_flat)
        {
            out_file << data.first << ", " << data.second << "\n";
        }
    }
    else 
    {
        // simple way to make sure there are no empty bins
        if (n_bins > 2*center_x || n_bins > 2*center_y) {
            n_bins = 2*center_x;
            std::cout << "The number of bins exceeds the number of cells in one direction!\n";
            std::cout << "Setting n_bins = " << n_bins << "\n";
        }

        // bin ranges
        auto max_dk = std::max_element(sf_flat.begin(), sf_flat.end(),
                [](const auto& a, const auto& b) {
                return a.first < b.first;
                });

        auto bin_r = max_dk->first / double(n_bins);

        // bins<bin_index<sum_values, n_values>>
        std::map<int, std::pair<double, int>> bins;
        
        for (const auto& data : sf_flat)
        {
            size_t bin_i = static_cast<size_t>(std::floor(data.first / bin_r));
            bins[bin_i].first += data.second;
            bins[bin_i].second += 1; 
        }

        std::ofstream out_file(out_name);

        if (!out_file.is_open()) {
            std::cerr << "Error: Could not open file " << out_name << " for writing." << std::endl;
            return 1; 
        }

        out_file << "dk, sf\n";

        for (const auto& [bin_i, bin_d] : bins)
        {
            auto bin_center = bin_i * bin_r + 0.5*bin_r;
            auto bin_average = double(bin_d.first / bin_d.second);

            out_file << bin_center << ", " << bin_average << "\n";
        }

    }
    
    
    std::cout << "... and we are done!\n\n";

}


