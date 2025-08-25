#include <fstream>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;
using namespace std;

static
void
PrintUsage (const char* progName)
{
    Print() << std::endl
            << "This utility computes PDF of scalars, and various powers of Laplacian of velocity field," << std::endl;

    Print() << "Usage:" << '\n';
    Print() << progName << " <inputs>" << std::endl
            << "OR" << std::endl
            << progName << std::endl
            << " step=<step number of plotfile to be read>" << std::endl
            << " nbins=<number of bins> " << std::endl
            << " range=<lo/hi end of range> " << std::endl
            << std::endl;

    exit(1);
}


int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {

        if (argc == 1) {
            PrintUsage(argv[0]);
        }

        ParmParse pp;

        int step;
        pp.query("step",step);

        std::string iFile         = amrex::Concatenate("plt",step,9);

        Vector<std::string> scalar_out(3);
        scalar_out[0] = amrex::Concatenate("rho_pdf",step,9);
        scalar_out[1] = amrex::Concatenate("press_pdf",step,9);
        scalar_out[2] = amrex::Concatenate("temp_pdf",step,9);
        Vector<std::string> Lap_out(5);
        Lap_out[0] = amrex::Concatenate("L0_pdf",step,9);
        Lap_out[1] = amrex::Concatenate("L1_pdf",step,9);
        Lap_out[2] = amrex::Concatenate("L2_pdf",step,9);
        Lap_out[3] = amrex::Concatenate("L3_pdf",step,9);
        Lap_out[4] = amrex::Concatenate("L4_pdf",step,9);

        int nbins;
        pp.get("nbins", nbins);

        Real range;
        pp.get("range",range);

        amrex::Print() << "Reading from plotfile " << iFile << "\n";

        // for the Header
        std::string iFile2 = iFile;
        iFile2 += "/Header";

        // open header
        ifstream x;
        x.open(iFile2.c_str(), ios::in);

        // read in first line of header (typically "HyperCLaw-V1.1" or similar)
        std::string str;
        x >> str;

        // read in number of components from header
        int ncomp;
        x >> ncomp;

        // read in variable names from header
        int flag = 0;
        int rho_ind, press_ind, temp_ind, velx_ind;
        for (int n=0; n<ncomp; ++n) {
            x >> str;
            if (str == "rhoInstant") rho_ind = flag;
            if (str == "pInstant") press_ind = flag;
            if (str == "tInstant") temp_ind = flag;
            if (str == "uxInstantFACE") velx_ind = flag;
            flag ++;
        }

        // read in dimensionality from header
        int dim;
        x >> dim;

        // read in time
        Real time;
        x >> time;

        // read in finest level
        int finest_level;
        x >> finest_level;

        // read in prob_lo and prob_hi
        amrex::GpuArray<amrex::Real, 3> prob_lo, prob_hi;
        for (int i=0; i<3; ++i) {
            x >> prob_lo[i];
        }
        for (int i=0; i<3; ++i) {
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

        Real ncells = (double) domain.numPts();

        // set to 1 (periodic)
        Vector<int> is_periodic(3,1);

        Geometry geom(domain,&real_box,CoordSys::cartesian,is_periodic.data());

        const Real* dx = geom.CellSize();

        ////////////////////////////////////////////////////////////////////////
        ////////////// velocity Laplacian PDFs /////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        MultiFab vel_grown(ba,dmap,3,1);
        MultiFab laplacian(ba,dmap,3,1);

        // copy shifted velocity components from mf into vel_grown
        Copy(vel_grown,mf,velx_ind,0,3,0);
        Copy(laplacian,mf,velx_ind,0,3,0);

        // fill ghost cells of vel_grown
        vel_grown.FillBoundary(geom.periodicity());
        laplacian.FillBoundary(geom.periodicity());

        for (int m=0; m<5; ++m) {

            Vector<Real> L2(3,0.);
            for (int i=0; i<3; i++)
                L2[i]=0.;

            for ( MFIter mfi(laplacian,false); mfi.isValid(); ++mfi ) {

                const Box& bx = mfi.validbox();
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                const Array4<Real>& lap = laplacian.array(mfi);

                for (auto n=0; n<3; ++n) {
                for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {

                    L2[n] += lap(i,j,k,n)*lap(i,j,k,n);

                }
                }
                }
                }

            } // end MFIter

            ParallelDescriptor::ReduceRealSum(L2.dataPtr(),3);
            amrex::Long totpts =  domain.numPts();
            L2[0] = sqrt(L2[0]/totpts);
            L2[1] = sqrt(L2[1]/totpts);
            L2[2] = sqrt(L2[2]/totpts);
            Print() << "L2 norm of Laplacian to power " << m << " is " << L2[0]
                    << " "  << L2[1] << " "  << L2[2] << " " << std::endl;

            for ( MFIter mfi(laplacian,false); mfi.isValid(); ++mfi ) {

                const Box& bx = mfi.validbox();
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                const Array4<Real>& lap = laplacian.array(mfi);

                for (auto n=0; n<3; ++n) {
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
            for (int ind=0 ; ind < nbins+1; ind++) bins[ind]=0;

            for ( MFIter mfi(laplacian,false); mfi.isValid(); ++mfi ) {

                const Box& bx = mfi.validbox();
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                const Array4<Real>& lap = laplacian.array(mfi);

                for (auto n=0; n<3; ++n) {
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
            Print() << "Points outside of range "<< count - totbin << " " <<
                       (double)(count-totbin)/count << std::endl;

            // print out contents of bins to the screen
            for (int i=0; i<nbins+1; ++i) {
                Print() << "For  m= "<< m << " " <<  (i-halfbin)*binwidth << " "
                        << bins[i]/(count*binwidth) << std::endl;
            }
            if (ParallelDescriptor::IOProcessor()) {
                std::ofstream outfile;
                outfile.open(Lap_out[m]);
                for (int i=0; i<nbins+1; ++i) {
                    outfile << (i-halfbin)*binwidth << " " << bins[i]/(count*binwidth) << std::endl;
                }
                outfile.close();
            }

            for ( MFIter mfi(vel_grown,false); mfi.isValid(); ++mfi ) {

                const Box& bx = mfi.validbox();
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                const Array4<Real>& vel = vel_grown.array(mfi);
                const Array4<Real>& lap = laplacian.array(mfi);

                for (auto n=0; n<3; ++n) {
                for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {

                    lap(i,j,k,n) = -(vel(i+1,j,k,n) - 2.*vel(i,j,k,n) + vel(i-1,j,k,n)) / (dx[0]*dx[0])
                                   -(vel(i,j+1,k,n) - 2.*vel(i,j,k,n) + vel(i,j-1,k,n)) / (dx[1]*dx[1])
                                   -(vel(i,j,k+1,n) - 2.*vel(i,j,k,n) + vel(i,j,k+1,n)) / (dx[2]*dx[2]);
                }
                }
                }
                }

            } // end MFIter

            // copy lap into vel_grown
            Copy(vel_grown,laplacian,0,0,3,0);

            // fill ghost cells of vel_grown
            vel_grown.FillBoundary(geom.periodicity());

        } // end loop
        ////////////////////////////////////////////////////////////////////////
        ////////////// velocity Laplacian PDFs /////////////////////////////////
        ////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////
        ///////////////////////// scalar  PDFs /////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        MultiFab scalar(ba,dmap,3,0);
        Copy(scalar,mf,rho_ind,0,1,0);
        Copy(scalar,mf,press_ind,1,1,0);
        Copy(scalar,mf,temp_ind,2,1,0);

        // compute spatial mean
        Real mean_rho   = scalar.sum(0) / (ncells);
        Real mean_press = scalar.sum(1) / (ncells);
        Real mean_temp  = scalar.sum(2) / (ncells);

        // get fluctuations
        scalar.plus(-1.0*mean_rho,   0, 1);
        scalar.plus(-1.0*mean_press, 1, 1);
        scalar.plus(-1.0*mean_temp,  2, 1);

        // get rms
        Real rms_rho   = scalar.norm2(0) / sqrt(ncells);
        Real rms_press = scalar.norm2(1) / sqrt(ncells);
        Real rms_temp  = scalar.norm2(2) / sqrt(ncells);

        // scale by rms
        scalar.mult(1.0/rms_rho,   0, 1);
        scalar.mult(1.0/rms_press, 1, 1);
        scalar.mult(1.0/rms_temp,  2, 1);

        // now compute pdfs
        for (int m = 0; m < 3; ++m) {

            Vector<Real> bins(nbins+1,0.);

            int halfbin = nbins/2;
            Real hbinwidth = range/nbins;
            Real binwidth = 2.*range/nbins;
            amrex::Long count=0;
            amrex::Long totbin=0;
            for (int ind=0 ; ind < nbins+1; ind++) bins[ind]=0;

            for ( MFIter mfi(scalar,false); mfi.isValid(); ++mfi ) {

                const Box& bx = mfi.validbox();
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                const Array4<Real>& sca = scalar.array(mfi);

                for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {

                    int index = floor((sca(i,j,k,m) + hbinwidth)/binwidth);
                    index += halfbin;

                    if( index >=0 && index <= nbins) {
                        bins[index] += 1;
                        totbin++;
                    }

                    count++;

                }
                }
                }

            } // end MFIter

            ParallelDescriptor::ReduceRealSum(bins.dataPtr(),nbins+1);
            ParallelDescriptor::ReduceLongSum(count);
            ParallelDescriptor::ReduceLongSum(totbin);
            Print() << "Points outside of range "<< count - totbin << " " <<
                       (double)(count-totbin)/count << std::endl;

            // print out contents of bins to the screen
            for (int i=0; i<nbins+1; ++i) {
                Print() << "For scalar m = "<< m << " " <<  (i-halfbin)*binwidth << " "
                        << bins[i]/(count*binwidth) << std::endl;
            }
            if (ParallelDescriptor::IOProcessor()) {
                std::ofstream outfile;
                outfile.open(scalar_out[m]);
                for (int i=0; i<nbins+1; ++i) {
                    outfile << (i-halfbin)*binwidth << " " << bins[i]/(count*binwidth) << std::endl;
                }
                outfile.close();
            }
        }
    }

    amrex::Finalize();

}

