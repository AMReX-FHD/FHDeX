#include "LocalFunctions.H"
#include "common_functions.H"
#include "species.H"
#include "DsmcParticleContainer.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;
using namespace std;

void writePlotFile(const MultiFab& mfcuInst,
    const MultiFab& mfcuMeans,
    const MultiFab& mfcuVars,
    FhdParticleContainer& particles,
    const Geometry& geom,
    Real time,
    int step) {
    BL_PROFILE_VAR("writePlotFile()",writePlotFile);

    int ncon    = 5;

    amrex::BoxArray ba = mfcuInst.boxArray();
    amrex::DistributionMapping dmap = mfcuInst.DistributionMap();

    std::string pltcu    = amrex::Concatenate("pltcu",step,12);
    std::string pltpart   = amrex::Concatenate("plt",step,12);

    amrex::MultiFab mfcuplt(ba, dmap, 2*ncon, 0);
    amrex::Vector<std::string> cuNames(2*ncon);


    int cnt = 0;
    // Instant values
    cuNames[cnt++] = "densityInstant";
    cuNames[cnt++] = "energyFluxDensityXInstant";
    cuNames[cnt++] = "energyFluxDensityYInstant";
    cuNames[cnt++] = "energyFluxDensityZInstant";
    cuNames[cnt++] = "energyDensity";

    MultiFab::Copy(mfcuplt, mfcuInst, 0, ncon*0, ncon, 0);

    // Mean values
    cuNames[cnt++] = "densityMean";
    cuNames[cnt++] = "energyFluxDensityXMean";
    cuNames[cnt++] = "energyFluxDensityYMean";
    cuNames[cnt++] = "energyFluxDensityZMean";
    cuNames[cnt++] = "energyDensityMean";

    MultiFab::Copy(mfcuplt, mfcuMeans, 0, ncon*1, ncon, 0);
    WriteSingleLevelPlotfile(pltcu, mfcuplt, cuNames, geom, time, step);



    // particle in cplt file
    Vector<std::string> real_comp_names = FHD_realData::names();
    Vector<std::string>  int_comp_names = FHD_intData::names();

    Vector<int> write_real_comp(real_comp_names.size(),1);
    Vector<int> write_int_comp(int_comp_names.size(),1);

    particles.WritePlotFile(pltpart, "particles",
                            write_real_comp, write_int_comp, real_comp_names, int_comp_names);

}