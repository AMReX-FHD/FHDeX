#include "common_functions.H"
#include "common_functions_F.H"


void outputMFAscii(const MultiFab& output, std::string filename)
{
    std::ofstream ofs(filename, std::ofstream::out);
    
    for (MFIter mfi(output); mfi.isValid(); ++mfi) 
    {
        amrex::Print(ofs)<<(output[mfi])<<std::endl;                       
    }

}



