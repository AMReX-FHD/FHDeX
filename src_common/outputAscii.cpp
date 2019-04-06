#include "common_functions.H"
#include "common_functions_F.H"


void outputMFAscii(const MultiFab& output, std::string filename)
{

    std::string plotfilename = std::to_string(ParallelDescriptor::MyProc()) + "_" + filename;

    std::ofstream ofs(plotfilename, std::ofstream::out);
    
    for (MFIter mfi(output); mfi.isValid(); ++mfi) 
    {
        //std::Print(ofs)<<(output[mfi])<<std::endl;
          ofs<<(output[mfi])<<std::endl;                                              
    }

}



