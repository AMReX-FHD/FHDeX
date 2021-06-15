#include "AMReX_PlotFileUtil.H"
#include "AMReX_PlotFileDataImpl.H"
#include "DsmcParticleContainer.H"
#include <sys/stat.h>

#include "common_functions.H"
#include "common_namespace.H"

using namespace common;

namespace {
	void GotoNextLine (std::istream& is) {
		constexpr std::streamsize bl_ignore_max { 100000 };
		is.ignore(bl_ignore_max, '\n');
	}
}

void FhdParticleContainer::WriteCheckPoint(int step, int statsCount, const Real time,
                     const amrex::Geometry geom,
                     const amrex::MultiFab& cuInst,   const amrex::MultiFab& cuMean, 
                     const amrex::MultiFab& cuDel,    const amrex::MultiFab& covar){//,
                //     const amrex::MultiFab& mfselect, const amrex::MultiFab& mfphi,
                  //   const amrex::MultiFab& mfvrmax  ) {
	// timer for profiling
	BL_PROFILE_VAR("WriteCheckPoint()",WriteCheckPoint);

	// checkpoint file name, e.g., chk0000010
	const std::string& checkpointname = amrex::Concatenate("chk",step,9);

	amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

	BoxArray ba = cuInst.boxArray();

	// single level problem
	int nlevels = 1;

	amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);
	VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

	// write Header file
	if (ParallelDescriptor::IOProcessor()) {
		std::ofstream HeaderFile;
		HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
		std::string HeaderFileName(checkpointname + "/Header");
		HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                      std::ofstream::trunc |
                      std::ofstream::binary);

		if( !HeaderFile.good()) { amrex::FileOpenFailed(HeaderFileName); }

		HeaderFile.precision(17);
		HeaderFile << "Checkpoint file for FHDeX/dsmc\n"; // write out title line
		HeaderFile << step << "\n"; // write out the time step number
		HeaderFile << time << "\n"; // write out time
      HeaderFile << statsCount << "\n"; // write out statsCount
      ba.writeOn(HeaderFile); // write the BoxArray (fluid)
		HeaderFile << '\n';
	}

	// write the MultiFab data to, e.g., chk00010/Level_0/

	// Stat MFs
	VisMF::Write(cuInst,
					 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuInst"));
	VisMF::Write(cuMean,
	             amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuMean"));
	VisMF::Write(cuDel,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuDel"));
	VisMF::Write(covar,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "covar"));

	// Collision related MFs (defined in DSMC particle container)
	/*VisMF::Write(mfselect,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "mfselect"));
	VisMF::Write(mfphi,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "mfphi"));
	VisMF::Write(mfvrmax,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "mfvrmax"));*/

}

void FhdParticleContainer::ReadCheckPoint(int &step, int& statsCount, amrex::Real& time,
                    amrex::Geometry geom,
                    amrex::MultiFab& cuInst,   amrex::MultiFab& cuMean, 
                    amrex::MultiFab& cuDel,    amrex::MultiFab& covar){//,
      //              amrex::MultiFab& mfselect, amrex::MultiFab& mfphi,
        //            amrex::MultiFab& mfvrmax) {
    // timer for profiling
    BL_PROFILE_VAR("ReadCheckPoint()",ReadCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate("chk",restart,9);

    amrex::Print() << "Restart from checkpoint " << checkpointname << "\n";

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    std::string line, word;

    // Header
    {
        std::string File(checkpointname + "/Header");
        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);

        // read in title line
        std::getline(is, line);

        // read in time step number
        is >> step;
        ++step;
        GotoNextLine(is);

        // read in time
        is >> statsCount;
        GotoNextLine(is);

        // read in statsCount
        is >> time;
        GotoNextLine(is);

        // read in BoxArray (fluid) from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };
        
        // cu, cuMeans, cuVars
        // replace nvars
        cuInst.define(ba,dm,nvars,0);
        cuMean.define(ba,dm,nvars,0);
        cuDel.define(ba,dm,nvars,0);
		  covar.define(ba,dm,nvars,0);
		  
        // prim, primMeans, primVars
        //mfselect.define(ba,dm,nprimvars,0);
        //mfphi.define(ba,dm,nprimvars,0);
        //mfvrmax.define(ba,dm,nprimvars,0);
    }

	// Stat MFs
	VisMF::Read(cuInst,
					 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuInst"));
	VisMF::Read(cuMean,
	             amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuMean"));
	VisMF::Read(cuDel,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuDel"));
	VisMF::Read(covar,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "covar"));

	// Collision related MFs
	/*VisMF::Read(mfselect,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "mfselect"));
	VisMF::Read(mfphi,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "mfphi"));
	VisMF::Read(mfvrmax,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "mfvrmax"));*/
}
