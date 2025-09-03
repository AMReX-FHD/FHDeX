#include <fstream>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace std;
using namespace amrex;

static
void
PrintUsage (const char* progName)
{
    Print() << std::endl
            << "This utility computes a PDF of the differentiated velocity field," << std::endl
            << "Lap^n(u)" << std::endl << std::endl;

    Print() << "Usage:" << '\n';
    Print() << progName << " <inputs>" << std::endl
            << "OR" << std::endl
            << progName << std::endl
            << " infile=<plotFileName>" << std::endl
            << " outfile=<base_output_filename> " << std::endl
            << " nbins=<number of bins> " << std::endl
            << " range=<lo/hi end of range> " << std::endl
            << std::endl;

    exit(1);
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
{
    if (argc == 1) {
        PrintUsage(argv[0]);
    }

    ParmParse pp;

    // plotfile name
    std::string iFile;
    pp.get("infile", iFile);

    // plotfile name
    std::string oFile;
    pp.get("outfile", oFile);

    std::string oFile_save = oFile;


    // how many bins
    int nbins;
    pp.get("nbins", nbins);

    Real range;
    pp.get("range",range);

    // for the Header
    std::string iFile2 = iFile;
    iFile2 += "/Header";

    // open header
    ifstream x;
    x.open(iFile2.c_str(), ios::in);

    // read in first line of header (typically "HyperCLaw-V1.1" or similar)
    string str;
    x >> str;

    // read in number of components from header
    int ncomp;
    x >> ncomp;

    // read in variable names from header
    for (int n=0; n<ncomp; ++n) {
        x >> str;
    }

    // read in dimensionality from header
    int dim;
    x >> dim;

    if (dim != AMREX_SPACEDIM) {
        Print() << "\nError: you are using a " << AMREX_SPACEDIM
                << "D build to open a " << dim << "D plotfile\n\n";
        Abort();
    }

    // read in time
    Real time;
    x >> time;

    // read in finest level
    int finest_level;
    x >> finest_level;

    // read in prob_lo and prob_hi
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo, prob_hi;
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        x >> prob_lo[i];
    }
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        x >> prob_hi[i];
    }

    // now read in the plotfile data
    // check to see whether the user pointed to the plotfile base directory
    // or the data itself
    if (amrex::FileExists(iFile+"/Level_0/Cell_H")) {
       iFile += "/Level_0/Cell";
    }
    if (amrex::FileExists(iFile+"/Level_00/Cell_H")) {
       iFile += "/Level_00/Cell";
    }

    // storage for the input coarse and fine MultiFabs
    MultiFab mf;

    // read in plotfile mf to MultiFab
    VisMF::Read(mf, iFile);

    // get BoxArray and DistributionMapping
    BoxArray ba = mf.boxArray();
    DistributionMapping dmap = mf.DistributionMap();

    // physical dimensions of problem
    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    // single box with the enire domain
    Box domain = ba.minimalBox().enclosedCells();

    // set to 1 (periodic)
    Vector<int> is_periodic(AMREX_SPACEDIM,1);

    Geometry geom(domain,&real_box,CoordSys::cartesian,is_periodic.data());

    const Real* dx = geom.CellSize();

    MultiFab mf_grown(ba,dmap,AMREX_SPACEDIM,1);
    MultiFab laplacian(ba,dmap,AMREX_SPACEDIM,1);

    // copy shifted velocity components from mf into mf_grown
    Copy(mf_grown,mf,AMREX_SPACEDIM,0,AMREX_SPACEDIM,0);

    for (int nderivs = 0 ; nderivs <5; nderivs++){

        if(nderivs == 0){
            Copy(laplacian,mf,AMREX_SPACEDIM,0,AMREX_SPACEDIM,0);
        } else {

            // fill ghost cells of mf_grown
            mf_grown.FillBoundary(geom.periodicity());

            //    for (int m=0; m<nderivs; ++m) {

            for ( MFIter mfi(mf_grown,false); mfi.isValid(); ++mfi ) {

                const Box& bx = mfi.validbox();
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                const Array4<Real>& vel = mf_grown.array(mfi);
                const Array4<Real>& lap = laplacian.array(mfi);

                for (auto n=0; n<AMREX_SPACEDIM; ++n) {
                for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {

                    lap(i,j,k,n) = -(vel(i+1,j,k,n) - 2.*vel(i,j,k,n) + vel(i-1,j,k,n)) / (dx[0]*dx[0])
                                   -(vel(i,j+1,k,n) - 2.*vel(i,j,k,n) + vel(i,j-1,k,n)) / (dx[1]*dx[1])
#if (AMREX_SPACEDIM == 3)
                                   -(vel(i,j,k+1,n) - 2.*vel(i,j,k,n) + vel(i,j,k+1,n)) / (dx[2]*dx[2])
#endif
                        ;
                }
                }
                }
                }

            } // end MFIter


            // copy lap into mf_grown
            Copy(mf_grown,laplacian,0,0,AMREX_SPACEDIM,0);

            // fill ghost cells of mf_grown
            mf_grown.FillBoundary(geom.periodicity());

        } // end loop over nderivs

        Vector<Real> L2(AMREX_SPACEDIM,0.);
        for (int i=0; i<AMREX_SPACEDIM; i++)
            L2[i]=0.;

        for ( MFIter mfi(laplacian,false); mfi.isValid(); ++mfi ) {

            const Box& bx = mfi.validbox();
            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);

            const Array4<Real>& lap = laplacian.array(mfi);

            for (auto n=0; n<AMREX_SPACEDIM; ++n) {
            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
            for (auto i = lo.x; i <= hi.x; ++i) {

                L2[n] += lap(i,j,k,n)*lap(i,j,k,n);

            }
            }
            }
            }

        } // end MFIter

        ParallelDescriptor::ReduceRealSum(L2.dataPtr(),AMREX_SPACEDIM);
        amrex::Long totpts =  domain.numPts();
        L2[0] = sqrt(L2[0]/totpts);
        L2[1] = sqrt(L2[1]/totpts);
        L2[2] = sqrt(L2[2]/totpts);
        Print() << "L2 norm of Laplacian to power " << nderivs << " is " << L2[0] << " "  << L2[1] << " "  << L2[2] << " " << std::endl;

        for ( MFIter mfi(laplacian,false); mfi.isValid(); ++mfi ) {

            const Box& bx = mfi.validbox();
            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);

            const Array4<Real>& lap = laplacian.array(mfi);

            for (auto n=0; n<AMREX_SPACEDIM; ++n) {
            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
            for (auto i = lo.x; i <= hi.x; ++i) {

                lap(i,j,k,n) = lap(i,j,k,n)/L2[n];

            }
            }
            }
            }

        } // end MFIter

        Vector<Real> bins(nbins+1,0.);

        int halfbin = nbins/2;
        Real hbinwidth = range/nbins;
        Real binwidth = 2.*range/nbins;
        amrex::Long count=0;
        amrex::Long totbin=0;
        for (int ind=0 ; ind < nbins+1; ind++)
            bins[ind]=0;

        for ( MFIter mfi(laplacian,false); mfi.isValid(); ++mfi ) {

            const Box& bx = mfi.validbox();
            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);

            const Array4<Real>& lap = laplacian.array(mfi);

            for (auto n=0; n<AMREX_SPACEDIM; ++n) {
            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
            for (auto i = lo.x; i <= hi.x; ++i) {

                int index = floor((lap(i,j,k,n) + hbinwidth)/binwidth);
                index += halfbin;

                if( index >=0 && index <= nbins) {
                    bins[index] += 1;
                    totbin++;
                }

                count++;

            }
            }
            }
            }

        } // end MFIter

        ParallelDescriptor::ReduceRealSum(bins.dataPtr(),nbins+1);
        ParallelDescriptor::ReduceLongSum(count);
        ParallelDescriptor::ReduceLongSum(totbin);
        Print() << "Points outside of range "<< count - totbin << " " << (double)(count-totbin)/count << std::endl;

        // print out contents of bins to the screen
        for (int i=0; i<nbins+1; ++i) {
            Print() << "For  m="<< nderivs<< " " <<  (i-halfbin)*binwidth << " " << bins[i]/(count*binwidth) << std::endl;
        }
        if (ParallelDescriptor::IOProcessor()) {
            std::ofstream outfile;
            oFile = oFile_save;
            oFile +="_";
            oFile += std::to_string(nderivs);
            oFile += ".dat";

            outfile.open(oFile);
            for (int i=0; i<nbins+1; ++i) {
                outfile << (i-halfbin)*binwidth << " " << bins[i]/(count*binwidth) << std::endl;
            }
        }

    }

}
    amrex::Finalize();


}