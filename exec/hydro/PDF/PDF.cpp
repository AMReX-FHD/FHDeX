
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
            << "This utility computes... " << std::endl << std::endl;

    Print() << "Usage:" << '\n';
    Print() << progName << " <inputs>" << std::endl
            << "OR" << std::endl
            << progName << " infile=plotFileName ..." << '\n' << '\n';

    exit(1);
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    if (argc == 1) {
        PrintUsage(argv[0]);
    }

    // plotfile name
    std::string iFile;

    // how many 2nd-derivatives to take
    int nderivs;

    // read in parameters from inputs file or command line
    ParmParse pp;

    // MultiFab
    pp.get("infile", iFile);

    // how many 2nd derivatives
    pp.get("nderivs", nderivs);

    // single-level only

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
        Print() << "\nError: you are using a " << AMREX_SPACEDIM << "D build to open a "
                << dim << "D plotfile\n\n";
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

    // read in plotfiles, 'coarse' and 'fine' to MultiFabs
    // note: fine could be the same resolution as coarse
    VisMF::Read(mf, iFile);

    // get boxArray
    BoxArray ba = mf.boxArray();
    DistributionMapping dmap = mf.DistributionMap();

    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});
            
    Box domain = ba.minimalBox().enclosedCells();

    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // set to 1 (periodic)
    
    Geometry geom(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    
    const Real* dx = geom.CellSize();
  
    MultiFab mf_grown(ba,dmap,AMREX_SPACEDIM,1);
    MultiFab laplacian(ba,dmap,AMREX_SPACEDIM,1);

    // copy shifted velocity components from mf into mf_grown
    Copy(mf_grown,mf,AMREX_SPACEDIM,0,AMREX_SPACEDIM,0);

    // fill ghost cells of mf_grown
    mf_grown.FillBoundary(geom.periodicity());

    for (int i=0; i<nderivs; ++i) {
    
    
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

                lap(i,j,k,n) =  (vel(i+1,j,k,n) - 2.*vel(i,j,k,n) + vel(i-1,j,k,n)) / (dx[0]*dx[0])
                               +(vel(i,j+1,k,n) - 2.*vel(i,j,k,n) + vel(i,j-1,k,n)) / (dx[1]*dx[1])
#if (AMREX_SPACEDIM == 3)
                               +(vel(i,j,k+1,n) - 2.*vel(i,j,k,n) + vel(i,j,k+1,n)) / (dx[2]*dx[2])
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
    
    
    }

}
