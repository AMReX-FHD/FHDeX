#include "common_functions.H"

void PrintMF (const MultiFab& MF, int comp_lo, int comp_hi)
{
    BL_PROFILE_VAR("PrintMF()",PrintMF);

    const BoxArray& ba = MF.boxArray();
    const DistributionMapping& dm = MF.DistributionMap();
    const int myProc = ParallelDescriptor::MyProc();

    for (int i = 0; i < ba.size(); ++i) {
        if (dm[i] == myProc) {
            // we want all processors to write, not just the IOProcessor
            std::cout << "Grid #" << i << std::endl;
            std::cout << "Processor #" << myProc << std::endl;

            const Box& validBox = ba[i];

            auto lo = validBox.loVect3d();
            auto hi = validBox.hiVect3d();

            std::cout << "valid box ";
            for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                std::cout << "(" << lo[n] << ", " << hi[n] << ")  ";
            }
            std::cout << std::endl;

            const Array4<const Real> MF_arr = MF.array(i);
            int comp_top;
            if(comp_hi < 0)
            {
                comp_top = MF.nComp()-1;
            }else
            {
                comp_top = comp_hi;
            }

            for (auto k = lo[2]; k <= hi[2]; ++k) {
            for (auto j = lo[1]; j <= hi[1]; ++j) {
            for (auto ii = lo[0]; ii <= hi[0]; ++ii) {
            for (auto comp = comp_lo; comp <= comp_top; ++comp) {
#if (AMREX_SPACEDIM == 2)
                std::cout << "i, j, comp" << " "
                          << ii << " " << j << " " << comp
                          << " " << MF_arr(ii, j, k, comp)
                          << std::endl;
#else
                std::cout << "i, j, k, comp" << " "
                          << ii << " " << j << " " << k << " "
                          << comp << " "
                          << MF_arr(ii, j, k, comp)
                          << std::endl;
#endif
            }
            }
            }
            }
        }

        // add this barrier so only one grid gets printed out at a time
        ParallelDescriptor::Barrier();

    } // end loop over boxes
}

void outputMFAscii2(const MultiFab& MF, int comp_lo, int comp_hi, std::string filename)
{
    BL_PROFILE_VAR("PrintMF()",PrintMF);

    const BoxArray& ba = MF.boxArray();
    const DistributionMapping& dm = MF.DistributionMap();
    const int myProc = ParallelDescriptor::MyProc();

    for (int i = 0; i < ba.size(); ++i) {
        if (dm[i] == myProc) {
            std::ofstream ofs(filename, std::ofstream::app);
            const Box& validBox = ba[i];

            auto lo = validBox.loVect3d();
            auto hi = validBox.hiVect3d();

            const Array4<const Real> MF_arr = MF.array(i);
            int comp_top;
            if(comp_hi < 0)
            {
                comp_top = MF.nComp()-1;
            }else
            {
                comp_top = comp_hi;
            }

            for (auto k = lo[2]; k <= hi[2]; ++k) {
            for (auto j = lo[1]; j <= hi[1]; ++j) {
            for (auto ii = lo[0]; ii <= hi[0]; ++ii) {
            for (auto comp = comp_lo; comp <= comp_top; ++comp) {
#if (AMREX_SPACEDIM == 2)
                ofs << ii+1 << ", " << j+1 << ", " << MF_arr(ii, j, k, comp) << std::endl;
#else
                ofs << ii+1 << ", " << j+1 << ", " << k+1 << ", " << MF_arr(ii, j, k, comp)/18.015 << std::endl;
#endif
            }
            }
            }
            }
            ofs.close();
        }



        // add this barrier so only one grid gets printed out at a time
        ParallelDescriptor::Barrier();

    } // end loop over boxes
}

void outputMFAsciiReduce(const MultiFab& MF, int comp, std::string filename, int rat)
{
    BL_PROFILE_VAR("PrintMF()",PrintMF);

    const BoxArray& ba = MF.boxArray();
    const DistributionMapping& dm = MF.DistributionMap();
    const int myProc = ParallelDescriptor::MyProc();
    BoxArray coarsened_BA = ba;
    coarsened_BA.coarsen(rat);

    MultiFab coarse;
    coarse.define(coarsened_BA,dm,3,MF.nGrow());

    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);

    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    IntVect dom_hi_coarse(AMREX_D_DECL(n_cells[0]/rat-1, n_cells[1]/rat-1, n_cells[2]/rat-1));
    Box domain_coarse(dom_lo, dom_hi_coarse);

    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1) {
            is_periodic[i] = 1;
        }
    }

    Geometry geom(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    Geometry geom_coarse(domain_coarse,&real_box,CoordSys::cartesian,is_periodic.data());

    const IntVect ratio(AMREX_D_DECL(rat,rat,rat));
    amrex::average_down(MF,coarse, geom, geom_coarse,0,3,ratio);

    //CellFillCoarse(MF, geom, coarse, geom_coarse);


    for (int i = 0; i < coarsened_BA.size(); ++i) {
        if (dm[i] == myProc) {
            std::ofstream ofs(filename, std::ofstream::app);
            const Box& validBox = coarsened_BA[i];

            auto lo = validBox.loVect3d();
            auto hi = validBox.hiVect3d();

            const Array4<const Real> MF_arr = coarse.array(i);

            for (auto k = lo[2]; k <= hi[2]; ++k) {
            for (auto j = lo[1]; j <= hi[1]; ++j) {
            for (auto ii = lo[0]; ii <= hi[0]; ++ii) {

#if (AMREX_SPACEDIM == 2)
                ofs << ii+1 << ", " << j+1 << ", " << MF_arr(ii, j, k, comp) << std::endl;
#else
                ofs << ii+1 << ", " << j+1 << ", " << k+1 << ", " << MF_arr(ii, j, k, comp)/18.015 << std::endl;
#endif

            }
            }
            }
            ofs.close();
        }



        // add this barrier so only one grid gets printed out at a time
        ParallelDescriptor::Barrier();

    } // end loop over boxes
}

void outputMFAscii(const MultiFab& output, std::string filename)
{
    BL_PROFILE_VAR("outputMFAscii()",outputMFAscii);

    std::string plotfilename = std::to_string(ParallelDescriptor::MyProc()) + "_" + filename;

    std::ofstream ofs(plotfilename, std::ofstream::out);

    for (MFIter mfi(output); mfi.isValid(); ++mfi) {
        ofs<<std::setprecision(16)<< (output[mfi])<<std::endl;
    }

    ofs.close();

}
