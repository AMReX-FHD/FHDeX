#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include <typeinfo>

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace std;
using namespace amrex;

static void PrintUsage (const char* progName)
{
    Print() << std::endl
            << "This utility is intended to be used to convert structure factor data, via a plotfile," << std::endl
            << "into ASCII data representing a (x,y) scatter plot." << std::endl
            << "  inputs:" << std::endl
            << "    input_plt_file=inputFileName" << std::endl
            << "    output_ASCII_file=outputFileName" << std::endl
            << "    n_bins=numberOfBins (if n_bins = 0, the raw data is printed)" << std::endl; 

    Print() << "Usage:" << '\n';
    Print() << progName << " input_file=inputFileName output_ASCII_file=outputFileName n_bins=numberOfBins" << '\n' << '\n';

    exit(1);
}

int main (int argc, char* argv[]) {

    amrex::Initialize(argc,argv);

    // print usage
    if (argc == 1) {
        PrintUsage(argv[0]);
    }

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
    // TODO: Should we only print 1 component if empty?
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

    // read in dimensionality from header
    int dim;
    x >> dim;

    if (dim != AMREX_SPACEDIM) {
        Print() << "\nError: you are using a " << AMREX_SPACEDIM << "D build to open a "
                << dim << "D plotfile\n\n";
        Abort();
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
        std::map<int, std::pair<double, int>> bins;
        double bin_r = 1.e40;
                
        if (n_bins > 0) {
            auto max_dk = std::max_element(comp_SF_data[n].begin(), comp_SF_data[n].end(),
                    [](const auto& a, const auto& b) {
                    return a.first < b.first;
                    });

            bin_r = max_dk->first / double(n_bins);


            for (const auto& data : comp_SF_data[n])
            {
                size_t bin_i = static_cast<size_t>(std::floor(data.first / bin_r));
                bins[bin_i].first += data.second;
                bins[bin_i].second += 1; 
            }
        }
         
        // print the data with comp name appended to output_file
        std::string output_name = output_file;
        output_name += "_";
        output_name += comp_names[n];
        output_name += "_out.csv";


        std::ofstream out_file(output_name);
        if (!out_file.is_open()) {
            std::cerr << "Error: Could not open file " << output_name << " for writing." << std::endl;
            return 1; 
        }

        std::cout << "..Printing data to\n    " << output_name << "\n";
        
        out_file << "dk, sf\n";

        if (n_bins > 0) {
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
        } else {
            for (size_t i = 0; i < comp_SF_data[n].size(); ++i)
            {
                out_file << comp_SF_data[n][i].first << ", " << comp_SF_data[n][i].second << "\n";
            }
        }

        out_file.close();
    }
}


