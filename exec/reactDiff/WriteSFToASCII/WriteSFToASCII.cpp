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

int main (int argc, char* argv[]) {

    amrex::Initialize(argc,argv);

    // read inputs
    ParmParse pp;

    // input plt file
    std::string input_plt_file;
    pp.query("input_plt_file", input_plt_file);
    if (input_plt_file.empty()) {
        amrex::Abort("You must specify an input plotfile.");
    }

    // output file name
    std::string output_file;
    pp.query("output_ASCII_file", output_file);
    if (output_file.empty()) {
        output_file = input_plt_file;
    }

    // component to print
    std::string component_to_print;
    pp.query("component_to_print", component_to_print);

    // number of bins to use
    int n_bins = 1;
    pp.query("n_bins", n_bins);

    // read information from header
    std::string Header = input_plt_file;
    Header += "/Header";

    std::cout << Header << std::endl;

    // open header
    ifstream x;
    x.open(Header.c_str(), ios::in);

    // read in first line of header
    string str;
    x >> str;

    // read in number of components from header
    int ncomp;
    x >> ncomp;

    // read in variable names from header
    std::vector<std::string> comp_names;

    for (int n=0; n<ncomp; ++n) {
        x >> str;
        comp_names.push_back(str);
    }

    // Read data from file
    MultiFab mf;

    // Append Level_0/Cell for mf read
    std::string mf_file = input_plt_file;

    mf_file += "/Level_0/Cell";
    std::cout << "Reading mf data from " << mf_file << std::endl;

    // read in plotfile to MultiFab
    VisMF::Read(mf, mf_file);

    // figure out the center cell
    auto domain_box = mf.boxArray().minimalBox();

    auto center_x = ceil(domain_box.length(0)/2.);
    auto center_y = ceil(domain_box.length(1)/2.);
#if (AMREX_SPACEDIM == 3)
    auto center_z = ceil(domain_box.length(2)/2.);
#endif

    // container for all component data
    std::vector<std::vector<std::pair<double, double>>> comp_SF_data(ncomp);

    for ( MFIter mfi(mf,false); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.validbox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        const Array4<Real>& mfdata = mf.array(mfi);

        for (auto n=0; n<ncomp; ++n) { 
#if (AMREX_SPACEDIM == 3)
            for (auto k = lo.z; k <= hi.z; ++k) {
#endif
                for (auto j = lo.y; j <= hi.y; ++j) {
                    for (auto i = lo.x; i <= hi.x; ++i) {
                        auto dk_x = i - center_x;
                        auto dk_y = j - center_y;
#if (AMREX_SPACEDIM == 2)
                        auto dk = sqrt(dk_x*dk_x + dk_y*dk_y);
                        comp_SF_data[n].push_back(std::make_pair(dk, mfdata(i,j,0,n)));
#elif (AMREX_SPACEDIM == 3)
                        auto dk_z = k - center_z;
                        auto dk = sqrt(dk_x*dk_x + dk_y*dk_y + dk_z*dk_z);
                        comp_SF_data[n].push_back(std::make_pair(dk, mfdata(i,j,k,n)));
#endif
                    } 
                } 
#if (AMREX_SPACEDIM == 3)
            } 
#endif
        } // end ncomp
    } // end MFIter

    // for each component, print the data to file
    for (int n = 0; n < ncomp; ++n)
    {
        // sort the data
        std::sort(comp_SF_data[n].begin(), comp_SF_data[n].end(), [](const auto& a, const auto& b) {
                return a.first < b.first;
                });

        // bin the data
        auto max_dk = std::max_element(comp_SF_data[n].begin(), comp_SF_data[n].end(),
                [](const auto& a, const auto& b) {
                return a.first < b.first;
                });

        auto bin_r = max_dk->first / double(n_bins);

        std::map<int, std::pair<double, int>> bins;

        for (const auto& data : comp_SF_data[n])
        {
            size_t bin_i = static_cast<size_t>(std::floor(data.first / bin_r));
            bins[bin_i].first += data.second;
            bins[bin_i].second += 1; 
        }

        // print the data
        std::cout << "Printing data to \n";

        std::string output_name = output_file;
        output_name += "_";
        output_name += comp_names[n];
        output_name += "_out.csv";

        std::cout << "  " << output_name << "\n";

        std::ofstream out_file(output_name);

        if (!out_file.is_open()) {
            std::cerr << "Error: Could not open file " << output_name << " for writing." << std::endl;
            return 1; 
        }

        out_file << "dk, sf\n";

        for (const auto& [bin_i, bin_d] : bins)
        {
            if (bin_d.second == 0) {
                // TODO: We could just skip this bin instead of aborting.
                amrex::Abort("The specified 'n_bins' creates empty bins.");
            } 

            auto bin_center = bin_i * bin_r + 0.5*bin_r;
            auto bin_average = double(bin_d.first / bin_d.second);
            
            out_file << bin_center << ", " << bin_average << "\n";
        }

        out_file.close();
    }

}


