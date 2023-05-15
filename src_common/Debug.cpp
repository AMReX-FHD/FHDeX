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
            
            for (auto k = lo[2]-MF.nGrow(2); k <= hi[2]+MF.nGrow(2); ++k) {
            for (auto j = lo[1]-MF.nGrow(1); j <= hi[1]+MF.nGrow(1); ++j) {
            for (auto ii = lo[0]-MF.nGrow(0); ii <= hi[0]+MF.nGrow(0); ++ii) {
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
