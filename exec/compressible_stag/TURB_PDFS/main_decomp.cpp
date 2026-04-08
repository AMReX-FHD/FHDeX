#include <fstream>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

using namespace amrex;
using namespace std;

static
void
PrintUsage (const char* progName)
{
    Print() << std::endl
            << "This utility computes PDF of vorticity and divergence, and various powers of Laplacian of solenoidal and dilatational velocity field," << std::endl;

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

        std::string iFile         = amrex::Concatenate("vel_grad_decomp",step,9);

        Vector<std::string> scalar_out(5);
        scalar_out[0] = amrex::Concatenate("div_pdf",step,9);
        scalar_out[1] = amrex::Concatenate("vortx_pdf",step,9);
        scalar_out[2] = amrex::Concatenate("vorty_pdf",step,9);
        scalar_out[3] = amrex::Concatenate("vortz_pdf",step,9);
        scalar_out[4] = amrex::Concatenate("vort_pdf",step,9);
        Vector<std::string> Lap_out_sol(5);
        Lap_out_sol[0] = amrex::Concatenate("L0_pdf_sol",step,9);
        Lap_out_sol[1] = amrex::Concatenate("L1_pdf_sol",step,9);
        Lap_out_sol[2] = amrex::Concatenate("L2_pdf_sol",step,9);
        Lap_out_sol[3] = amrex::Concatenate("L3_pdf_sol",step,9);
        Lap_out_sol[4] = amrex::Concatenate("L4_pdf_sol",step,9);
        Vector<std::string> Lap_out_dil(5);
        Lap_out_dil[0] = amrex::Concatenate("L0_pdf_dil",step,9);
        Lap_out_dil[1] = amrex::Concatenate("L1_pdf_dil",step,9);
        Lap_out_dil[2] = amrex::Concatenate("L2_pdf_dil",step,9);
        Lap_out_dil[3] = amrex::Concatenate("L3_pdf_dil",step,9);
        Lap_out_dil[4] = amrex::Concatenate("L4_pdf_dil",step,9);

        int nbins;
        pp.get("nbins", nbins);

        Real range;
        pp.get("range",range);

        amrex::Print() << "Reading from vel_grad_decomp plotfile " << iFile << "\n";

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
        int vort_ind, div_ind, velx_sol_ind, vely_sol_ind, velz_sol_ind, velx_dil_ind, vely_dil_ind, velz_dil_ind;
        for (int n=0; n<ncomp; ++n) {
            x >> str;
            if (str == "vort") vort_ind = flag;
            if (str == "div")  div_ind = flag;
            if (str == "ux_s") velx_sol_ind = flag;
            if (str == "uy_s") vely_sol_ind = flag;
            if (str == "uz_s") velz_sol_ind = flag;
            if (str == "ux_d") velx_dil_ind = flag;
            if (str == "uy_d") vely_dil_ind = flag;
            if (str == "uz_d") velz_dil_ind = flag;
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
        ////////////// velocity Laplacian PDFs///////////// ////////////////////
        ////////////////////////////////////////////////////////////////////////
        MultiFab vel_grown(ba,dmap,6,1);
        MultiFab vel_sol  (ba,dmap,3,1);
        MultiFab laplacian(ba,dmap,6,1);

        // copy shifted velocity components from mf into vel_grown
        Copy(vel_grown,mf,velx_sol_ind,0,1,0); // sol
        Copy(vel_grown,mf,vely_sol_ind,1,1,0); // sol
        Copy(vel_grown,mf,velz_sol_ind,2,1,0); // sol

        Copy(laplacian,mf,velx_sol_ind,0,1,0); // sol
        Copy(laplacian,mf,vely_sol_ind,1,1,0); // sol
        Copy(laplacian,mf,velz_sol_ind,2,1,0); // sol

        Copy(vel_grown,mf,velx_dil_ind,3,1,0); // dil
        Copy(vel_grown,mf,vely_dil_ind,4,1,0); // dil
        Copy(vel_grown,mf,velz_dil_ind,5,1,0); // dil

        Copy(laplacian,mf,velx_dil_ind,3,1,0); // dil
        Copy(laplacian,mf,vely_dil_ind,4,1,0); // dil
        Copy(laplacian,mf,velz_dil_ind,5,1,0); // dil

        Copy(vel_sol,mf,velx_sol_ind,0,1,0); // sol
        Copy(vel_sol,mf,vely_sol_ind,1,1,0); // sol
        Copy(vel_sol,mf,velz_sol_ind,2,1,0); // sol

        // fill ghost cells of vel_grown
        vel_grown.FillBoundary(geom.periodicity());
        laplacian.FillBoundary(geom.periodicity());
        vel_sol  .FillBoundary(geom.periodicity());

        for (int m=0; m<5; ++m) {

            Vector<Real> L2(6,0.);
            for (int i=0; i<6; i++)
                L2[i]=0.;

            for ( MFIter mfi(laplacian,false); mfi.isValid(); ++mfi ) {

                const Box& bx = mfi.validbox();
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                const Array4<Real>& lap = laplacian.array(mfi);

                for (auto n=0; n<6; ++n) {
                for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {

                    L2[n] += lap(i,j,k,n)*lap(i,j,k,n);

                }
                }
                }
                }

            } // end MFIter

            ParallelDescriptor::ReduceRealSum(L2.dataPtr(),6);
            amrex::Long totpts =  domain.numPts();
            L2[0] = sqrt(L2[0]/totpts);
            L2[1] = sqrt(L2[1]/totpts);
            L2[2] = sqrt(L2[2]/totpts);
            L2[3] = sqrt(L2[3]/totpts);
            L2[4] = sqrt(L2[4]/totpts);
            L2[5] = sqrt(L2[5]/totpts);
            Print() << "L2 norm of Laplacian (solenoidal) to power " << m << " is " << L2[0]
                    << " "  << L2[1] << " "  << L2[2] << " " << std::endl;
            Print() << "L2 norm of Laplacian (dilational) to power " << m << " is " << L2[3]
                    << " "  << L2[4] << " "  << L2[5] << " " << std::endl;

            for ( MFIter mfi(laplacian,false); mfi.isValid(); ++mfi ) {

                const Box& bx = mfi.validbox();
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                const Array4<Real>& lap = laplacian.array(mfi);

                for (auto n=0; n<6; ++n) {
                for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {

                    lap(i,j,k,n) = lap(i,j,k,n)/L2[n];

                }
                }
                }
                }

            } // end MFIter

            Vector<Real> bins_sol(nbins+1,0.);
            Vector<Real> bins_dil(nbins+1,0.);

            int halfbin = nbins/2;
            Real hbinwidth = range/nbins;
            Real binwidth = 2.*range/nbins;
            amrex::Long count_sol=0;
            amrex::Long totbin_sol=0;
            amrex::Long count_dil=0;
            amrex::Long totbin_dil=0;
            for (int ind=0 ; ind < nbins+1; ind++) bins_sol[ind]=0;
            for (int ind=0 ; ind < nbins+1; ind++) bins_dil[ind]=0;

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
                        bins_sol[index] += 1;
                        totbin_sol++;
                    }

                    count_sol++;

                }
                }
                }
                }

                for (auto n=3; n<6; ++n) {
                for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {

                    int index = floor((lap(i,j,k,n) + hbinwidth)/binwidth);
                    index += halfbin;

                    if( index >=0 && index <= nbins) {
                        bins_dil[index] += 1;
                        totbin_dil++;
                    }

                    count_dil++;

                }
                }
                }
                }

            } // end MFIter

            ParallelDescriptor::ReduceRealSum(bins_sol.dataPtr(),nbins+1);
            ParallelDescriptor::ReduceLongSum(count_sol);
            ParallelDescriptor::ReduceLongSum(totbin_sol);
            ParallelDescriptor::ReduceRealSum(bins_dil.dataPtr(),nbins+1);
            ParallelDescriptor::ReduceLongSum(count_dil);
            ParallelDescriptor::ReduceLongSum(totbin_dil);
            Print() << "Points outside of range (solenoidal) "<< count_sol - totbin_sol << " " <<
                       (double)(count_sol-totbin_sol)/count_sol << std::endl;
            Print() << "Points outside of range (dilational) "<< count_dil - totbin_dil << " " <<
                       (double)(count_dil-totbin_dil)/count_dil << std::endl;

            // print out contents of bins to the screen
            for (int i=0; i<nbins+1; ++i) {
                Print() << "(solenoidal) For  m= "<< m << " " <<  (i-halfbin)*binwidth << " "
                        << bins_sol[i]/(count_sol*binwidth) << std::endl;
            }
            for (int i=0; i<nbins+1; ++i) {
                Print() << "(dilational) For  m= "<< m << " " <<  (i-halfbin)*binwidth << " "
                        << bins_dil[i]/(count_dil*binwidth) << std::endl;
            }
            if (ParallelDescriptor::IOProcessor()) {
                std::ofstream outfile;
                outfile.open(Lap_out_sol[m]);
                for (int i=0; i<nbins+1; ++i) {
                    outfile << (i-halfbin)*binwidth << " " << bins_sol[i]/(count_sol*binwidth) << std::endl;
                }
                outfile.close();
                outfile.open(Lap_out_dil[m]);
                for (int i=0; i<nbins+1; ++i) {
                    outfile << (i-halfbin)*binwidth << " " << bins_dil[i]/(count_dil*binwidth) << std::endl;
                }
                outfile.close();
            }

            for ( MFIter mfi(vel_grown,false); mfi.isValid(); ++mfi ) {

                const Box& bx = mfi.validbox();
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                const Array4<Real>& vel = vel_grown.array(mfi);
                const Array4<Real>& lap = laplacian.array(mfi);

                for (auto n=0; n<6; ++n) {
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
            Copy(vel_grown,laplacian,0,0,6,0);

            // fill ghost cells of vel_grown
            vel_grown.FillBoundary(geom.periodicity());

        } // end loop
        ////////////////////////////////////////////////////////////////////////
        ////////////// velocity Laplacian PDFs //////////// ////////////////////
        ////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////
        ///////////////////////// scalar  PDFs /////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        MultiFab scalar(ba,dmap,4,0);    // vort_mag, div, vort_x, vort_y, vort_z
        scalar.setVal(0.0);
        Copy(scalar,mf,div_ind,0,1,0);

        // Compute vorticity components and store in scalar
        for ( MFIter mfi(vel_sol,false); mfi.isValid(); ++mfi ) {

            const Box& bx = mfi.validbox();
            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);

            Array4<Real const> const& sol  = vel_sol   .array(mfi);
            Array4<Real>       const& sca  = scalar    .array(mfi);

            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
            for (auto i = lo.x; i <= hi.x; ++i) {
                // dw/dy - dv/dz
                sca(i,j,k,1) =
                    (sol(i,j+1,k,velz_sol_ind) - sol(i,j-1,k,velz_sol_ind)) / (2.*dx[1]) -
                    (sol(i,j,k+1,vely_sol_ind) - sol(i,j,k-1,vely_sol_ind)) / (2.*dx[2]);

                // dv/dx - du/dy
                sca(i,j,k,2) =
                    (sol(i+1,j,k,vely_sol_ind) - sol(i-1,j,k,vely_sol_ind)) / (2.*dx[0]) -
                    (sol(i,j+1,k,velx_sol_ind) - sol(i,j-1,k,velx_sol_ind)) / (2.*dx[1]);

                // du/dz - dw/dx
                sca(i,j,k,3) =
                    (sol(i,j,k+1,velx_sol_ind) - sol(i,j,k-1,velx_sol_ind)) / (2.*dx[2]) -
                    (sol(i+1,j,k,velz_sol_ind) - sol(i-1,j,k,velz_sol_ind)) / (2.*dx[0]);

            }
            }
            }
        }

        // compute spatial mean
        Real mean_div     = scalar.sum(0) / (ncells);
        Real mean_vortx   = scalar.sum(1) / (ncells);
        Real mean_vorty   = scalar.sum(2) / (ncells);
        Real mean_vortz   = scalar.sum(3) / (ncells);

        // get fluctuations
        scalar.plus(-1.0*mean_div,     0, 1);
        scalar.plus(-1.0*mean_vortx,   1, 1);
        scalar.plus(-1.0*mean_vorty,   2, 1);
        scalar.plus(-1.0*mean_vortz,   3, 1);

        // get rms
        Real rms_div     = scalar.norm2(0) / sqrt(ncells);
        Real rms_vortx   = scalar.norm2(1) / sqrt(ncells);
        Real rms_vorty   = scalar.norm2(2) / sqrt(ncells);
        Real rms_vortz   = scalar.norm2(3) / sqrt(ncells);

        // scale by rms
        scalar.mult(1.0/rms_div,     0, 1);
        scalar.mult(1.0/rms_vortx,   1, 1);
        scalar.mult(1.0/rms_vorty,   2, 1);
        scalar.mult(1.0/rms_vortz,   3, 1);

        // now compute pdfs
        for (int m = 0; m < 4; ++m) {

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

        // total vorticity PDF
        {
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

              for (auto n = 1;    n < 4;     ++n) {
              for (auto k = lo.z; k <= hi.z; ++k) {
              for (auto j = lo.y; j <= hi.y; ++j) {
              for (auto i = lo.x; i <= hi.x; ++i) {

                  int index = floor((sca(i,j,k,n) + hbinwidth)/binwidth);
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
              Print() << "For scalar m = "<< 4 << " " <<  (i-halfbin)*binwidth << " "
                      << bins[i]/(count*binwidth) << std::endl;
          }
          if (ParallelDescriptor::IOProcessor()) {
              std::ofstream outfile;
              outfile.open(scalar_out[4]);
              for (int i=0; i<nbins+1; ++i) {
                  outfile << (i-halfbin)*binwidth << " " << bins[i]/(count*binwidth) << std::endl;
              }
              outfile.close();
          }

        }

        // solenoidal  and dilataional velocity PDF
        MultiFab vel_decomp(ba,dmap,6,0);

        Copy(vel_decomp,mf,velx_sol_ind,0,1,0); // sol
        Copy(vel_decomp,mf,vely_sol_ind,1,1,0); // sol
        Copy(vel_decomp,mf,velz_sol_ind,2,1,0); // sol
        Copy(vel_decomp,mf,velx_dil_ind,3,1,0); // dil
        Copy(vel_decomp,mf,vely_dil_ind,4,1,0); // dil
        Copy(vel_decomp,mf,velz_dil_ind,5,1,0); // dil

        // compute spatial mean
        Real mean_solx   = vel_decomp.sum(0) / (ncells);
        Real mean_soly   = vel_decomp.sum(1) / (ncells);
        Real mean_solz   = vel_decomp.sum(2) / (ncells);
        Real mean_dilx   = vel_decomp.sum(3) / (ncells);
        Real mean_dily   = vel_decomp.sum(4) / (ncells);
        Real mean_dilz   = vel_decomp.sum(5) / (ncells);

        // get fluctuations
        vel_decomp.plus(-1.0*mean_solx,     0, 1);
        vel_decomp.plus(-1.0*mean_soly,     1, 1);
        vel_decomp.plus(-1.0*mean_solz,     2, 1);
        vel_decomp.plus(-1.0*mean_dilx,     3, 1);
        vel_decomp.plus(-1.0*mean_dily,     4, 1);
        vel_decomp.plus(-1.0*mean_dilz,     5, 1);

        // get rms
        Real rms_solx   = vel_decomp.norm2(0) / sqrt(ncells);
        Real rms_soly   = vel_decomp.norm2(1) / sqrt(ncells);
        Real rms_solz   = vel_decomp.norm2(2) / sqrt(ncells);
        Real rms_dilx   = vel_decomp.norm2(3) / sqrt(ncells);
        Real rms_dily   = vel_decomp.norm2(4) / sqrt(ncells);
        Real rms_dilz   = vel_decomp.norm2(5) / sqrt(ncells);

        // scale by rms
        vel_decomp.mult(1.0/rms_solx,   0, 1);
        vel_decomp.mult(1.0/rms_soly,   1, 1);
        vel_decomp.mult(1.0/rms_solz,   2, 1);
        vel_decomp.mult(1.0/rms_dilx,   3, 1);
        vel_decomp.mult(1.0/rms_dily,   4, 1);
        vel_decomp.mult(1.0/rms_dilz,   5, 1);

        // solenoidal
        {
          Vector<Real> bins(nbins+1,0.);

          int halfbin = nbins/2;
          Real hbinwidth = range/nbins;
          Real binwidth = 2.*range/nbins;
          amrex::Long count=0;
          amrex::Long totbin=0;
          for (int ind=0 ; ind < nbins+1; ind++) bins[ind]=0;

          for ( MFIter mfi(vel_decomp,false); mfi.isValid(); ++mfi ) {

              const Box& bx = mfi.validbox();
              const auto lo = amrex::lbound(bx);
              const auto hi = amrex::ubound(bx);

              const Array4<Real>& vel = vel_decomp.array(mfi);

              for (auto n = 0;    n < 3;     ++n) {
              for (auto k = lo.z; k <= hi.z; ++k) {
              for (auto j = lo.y; j <= hi.y; ++j) {
              for (auto i = lo.x; i <= hi.x; ++i) {

                  int index = floor((vel(i,j,k,n) + hbinwidth)/binwidth);
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
              Print() << "For solenoid. vel. " <<  (i-halfbin)*binwidth << " "
                      << bins[i]/(count*binwidth) << std::endl;
          }
          if (ParallelDescriptor::IOProcessor()) {
              std::ofstream outfile;
              outfile.open(amrex::Concatenate("solenoidal_pdf",step,9));
              for (int i=0; i<nbins+1; ++i) {
                  outfile << (i-halfbin)*binwidth << " " << bins[i]/(count*binwidth) << std::endl;
              }
              outfile.close();
          }

        }

        // dilatational
        {
          Vector<Real> bins(nbins+1,0.);

          int halfbin = nbins/2;
          Real hbinwidth = range/nbins;
          Real binwidth = 2.*range/nbins;
          amrex::Long count=0;
          amrex::Long totbin=0;
          for (int ind=0 ; ind < nbins+1; ind++) bins[ind]=0;

          for ( MFIter mfi(vel_decomp,false); mfi.isValid(); ++mfi ) {

              const Box& bx = mfi.validbox();
              const auto lo = amrex::lbound(bx);
              const auto hi = amrex::ubound(bx);

              const Array4<Real>& vel = vel_decomp.array(mfi);

              for (auto n = 3;    n < 6;     ++n) {
              for (auto k = lo.z; k <= hi.z; ++k) {
              for (auto j = lo.y; j <= hi.y; ++j) {
              for (auto i = lo.x; i <= hi.x; ++i) {

                  int index = floor((vel(i,j,k,n) + hbinwidth)/binwidth);
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
              Print() << "For dilation. vel. " <<  (i-halfbin)*binwidth << " "
                      << bins[i]/(count*binwidth) << std::endl;
          }
          if (ParallelDescriptor::IOProcessor()) {
              std::ofstream outfile;
              outfile.open(amrex::Concatenate("dilatational_pdf",step,9));
              for (int i=0; i<nbins+1; ++i) {
                  outfile << (i-halfbin)*binwidth << " " << bins[i]/(count*binwidth) << std::endl;
              }
              outfile.close();
          }

        }

    }

    amrex::Finalize();

}


