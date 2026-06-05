#include "FPU.H"
#include "common_functions.H"
#include "rng_functions.H"
#include <AMReX_PlotFileDataImpl.H>
#include "AMReX_ParmParse.H"
#include <AMReX_VisMF.H>
#include <fstream>
#include <iomanip>
#include <sstream>

AMREX_GPU_MANAGED int FPU::enable_fluctuations;
AMREX_GPU_MANAGED int FPU::diag_int;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3*3> FPU::A;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3*3> FPU::D;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3>   FPU::B;

//void ComputePhiFromState(MultiFab& phi) {
//
//    // convert state to state-state_eq
//    phi.plus(-FPU::r0,0,1,0);
//    phi.plus(-FPU::p0,1,1,0);
//    phi.plus(-FPU::e0,2,1,0);
//
//    // phi = R(state-state_eq)
//    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
//
//        const Box& bx = mfi.tilebox();
//
//        const Array4<Real>& phi_fab = phi.array(mfi);
//
//        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//        {
//            Real R0 = R_00*phi_fab(i,j,k,0) + R_01*phi_fab(i,j,k,1) + R_02*phi_fab(i,j,k,2);
//            Real R1 = R_10*phi_fab(i,j,k,0) + R_11*phi_fab(i,j,k,1) + R_12*phi_fab(i,j,k,2);
//            Real R2 = R_20*phi_fab(i,j,k,0) + R_21*phi_fab(i,j,k,1) + R_22*phi_fab(i,j,k,2);
//
//            phi_fab(i,j,k,0) = R0;
//            phi_fab(i,j,k,1) = R1;
//            phi_fab(i,j,k,2) = R2;
//
//        });
//    }
//
//}
//
//void ComputeCalphaalpha(MultiFab& C_alphaalpha,
//                        const MultiFab& phi,
//                        const MultiFab& phi0) {
//
//    for (MFIter mfi(C_alphaalpha); mfi.isValid(); ++mfi) {
//
//        const Box& bx = mfi.tilebox();
//
//        const Array4<      Real>& C_alphaalpha_fab = C_alphaalpha.array(mfi);
//        const Array4<const Real>& phi_fab          = phi.array(mfi);
//        const Array4<const Real>& phi0_fab         = phi0.array(mfi);
//
//        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//        {
//            C_alphaalpha_fab(i,j,k,0) = phi_fab(i,j,k,0) * phi0_fab(0,j,k,0);
//            C_alphaalpha_fab(i,j,k,1) = phi_fab(i,j,k,1) * phi0_fab(0,j,k,1);
//            C_alphaalpha_fab(i,j,k,2) = phi_fab(i,j,k,2) * phi0_fab(0,j,k,2);
//        });
//    }
//
//}

void InitializeNamespace() {

        ParmParse pp;

        pp.get("enable_fluctuations",FPU::enable_fluctuations);
        pp.get("diag_int",FPU::diag_int);
        Vector<Real> A_tmp(9);
        Vector<Real> D_tmp(9);
        Vector<Real> B_tmp(3);
        pp.getarr("A", A_tmp, 0, 9);
        pp.getarr("D", D_tmp, 0, 9);
        pp.getarr("B", B_tmp, 0, 3);
        for (int i = 0; i < 9; ++i) {
            FPU::A[i] = A_tmp[i];
            FPU::D[i] = D_tmp[i];
        }
        for (int i = 0; i < 3; ++i) {
            FPU::B[i] = B_tmp[i];
        }
}


void InitConsVarStag(MultiFab& cu,
                     MultiFab& cumom,
                     Geometry& geom)
{
    cu.setVal(0.0, 0, 3, cu.nGrowVect());
    cumom.setVal(0.0, 0, 1, cumom.nGrowVect());

    MultiFabFillRandomNormal(cu, 0, 1, 0.0, 1.0, geom, true, true);
    MultiFabFillRandomNormal(cu, 2, 1, 0.0, 1.0, geom, true, true);
    MultiFabFillRandomNormal(cumom, 0, 1, 0.0, 1.0, geom, true, true);

    cumom.FillBoundary(geom.periodicity());

    for (MFIter mfi(cu, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& u = cu.array(mfi);
        Array4<Real const> const& mx = cumom.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            u(i,j,k,1) = 0.5 * (mx(i,j,k,0) + mx(i+1,j,k,0));
        });
    }

    cu.FillBoundary(geom.periodicity());
}

void WriteCheckPoint(int step,
                     const Real time,
                     int statsCount,
                     const Geometry& /*geom*/,
                     const MultiFab& cu,
                     const MultiFab& cuMeans,
                     const MultiFab& cuVars,
                     const MultiFab& cumom,
                     const MultiFab& cumomMeans,
                     const MultiFab& cumomVars,
                     const MultiFab& coVars)
{
    BL_PROFILE_VAR("WriteCheckPoint()", WriteCheckPoint);

    const std::string checkpointname = amrex::Concatenate("chk", step, 9);
    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    constexpr int nlevels = 1;
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    if (ParallelDescriptor::IOProcessor()) {
        const std::string header_name = checkpointname + "/Header";
        std::ofstream header;
        header.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        header.open(header_name.c_str(), std::ofstream::out |
                    std::ofstream::trunc | std::ofstream::binary);

        if (!header.good()) {
            amrex::FileOpenFailed(header_name);
        }

        header.precision(17);
        header << "Checkpoint file for FHDeX/FPU_FHD\n";
        header << step << "\n";
        header << time << "\n";
        header << statsCount << "\n";
        cu.boxArray().writeOn(header);
        header << '\n';
    }

    const int my_rank = ParallelDescriptor::MyProc();
    const int n_ranks = ParallelDescriptor::NProcs();

    for (int rank = 0; rank < n_ranks; ++rank) {
        if (my_rank == rank) {
            const std::string rng_base = checkpointname + "/rng";
            const std::string rng_name = amrex::Concatenate(rng_base, my_rank, 7);

            std::ofstream rng_file;
            rng_file.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            rng_file.open(rng_name.c_str(), std::ofstream::out |
                          std::ofstream::trunc | std::ofstream::binary);

            if (!rng_file.good()) {
                amrex::FileOpenFailed(rng_name);
            }

            amrex::SaveRandomState(rng_file);
        }
        ParallelDescriptor::Barrier();
    }

    VisMF::Write(cu,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cu"));
    VisMF::Write(cuMeans,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuMeans"));
    VisMF::Write(cuVars,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuVars"));

    VisMF::Write(cumom,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumom"));
    VisMF::Write(cumomMeans,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomMeans"));
    VisMF::Write(cumomVars,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomVars"));

    VisMF::Write(coVars,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "coVars"));
}

void ReadCheckPoint(int& step,
                    Real& time,
                    int& statsCount,
                    Geometry& geom,
                    const amrex::Box& domain,
                    MultiFab& cu,
                    MultiFab& cuMeans,
                    MultiFab& cuVars,
                    MultiFab& cumom,
                    MultiFab& cumomMeans,
                    MultiFab& cumomVars,
                    MultiFab& coVars,
                    BoxArray& ba,
                    DistributionMapping& dmap)
{
    BL_PROFILE_VAR("ReadCheckPoint()", ReadCheckPoint);

    const std::string checkpointname = amrex::Concatenate("chk", restart, 9);
    amrex::Print() << "Restart from checkpoint " << checkpointname << "\n";

    std::string line;
    BoxArray ba_old;

    {
        const std::string header_name = checkpointname + "/Header";
        Vector<char> file_chars;
        ParallelDescriptor::ReadAndBcastFile(header_name, file_chars);
        std::string header_string(file_chars.dataPtr());
        std::istringstream is(header_string, std::istringstream::in);

        std::getline(is, line);

        is >> step;
        ++step;
        GotoNextLine(is);

        is >> time;
        GotoNextLine(is);

        is >> statsCount;
        GotoNextLine(is);

        ba_old.readFrom(is);
        GotoNextLine(is);
    }

    ba.define(domain);
    ba.maxSize(IntVect(max_grid_size));
    dmap.define(ba, ParallelDescriptor::NProcs());

    cu.define(ba, dmap, 3, ngc);
    cuMeans.define(ba, dmap, 3, ngc);
    cuVars.define(ba, dmap, 3, ngc);

    cumom.define(convert(ba, nodal_flag_x), dmap, 1, ngc);
    cumomMeans.define(convert(ba, nodal_flag_x), dmap, 1, 0);
    cumomVars.define(convert(ba, nodal_flag_x), dmap, 1, 0);

    coVars.define(ba, dmap, 3, 0);

    if (common::seed == -1) {
#ifdef AMREX_USE_CUDA
        Abort("Restart with negative seed not supported on GPU");
#endif
        const int my_rank = ParallelDescriptor::MyProc();
        const int n_ranks = ParallelDescriptor::NProcs();

        for (int rank = 0; rank < n_ranks; ++rank) {
            if (my_rank == rank) {
                const std::string rng_base = checkpointname + "/rng";
                const std::string rng_name = amrex::Concatenate(rng_base, my_rank, 7);
                std::ifstream rng_file(rng_name.c_str(), std::ios::in | std::ios::binary);

                if (!rng_file.good()) {
                    amrex::FileOpenFailed(rng_name);
                }

                amrex::RestoreRandomState(rng_file, 1, 0);
            }
            ParallelDescriptor::Barrier();
        }
    }

    VisMF::Read(cu,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cu"));
    VisMF::Read(cumom,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumom"));

    if (reset_stats == 1) {
        statsCount = 1;
        cuMeans.setVal(0.0);
        cuVars.setVal(0.0);
        cumomMeans.setVal(0.0);
        cumomVars.setVal(0.0);
        coVars.setVal(0.0);
    } else {
        VisMF::Read(cuMeans,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuMeans"));
        VisMF::Read(cuVars,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuVars"));
        VisMF::Read(cumomMeans,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomMeans"));
        VisMF::Read(cumomVars,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomVars"));
        VisMF::Read(coVars,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "coVars"));
    }

    cu.FillBoundary(geom.periodicity());
    cuMeans.FillBoundary(geom.periodicity());
    cuVars.FillBoundary(geom.periodicity());
    cumom.FillBoundary(geom.periodicity());
}

void WritePlotFile(int step,
                       const Real time,
                       const Geometry& geom,
                       const MultiFab& cu,
                       const MultiFab& cuMeans,
                       const MultiFab& cuVars,
                       const MultiFab& cumom,
                       const MultiFab& cumomMeans,
                       const MultiFab& cumomVars,
                       const MultiFab& coVars)
{
    BL_PROFILE_VAR("WritePlotFile()", WritePlotFile);

    constexpr int nplot = 15;
    int cnt = 0;

    MultiFab plotfile(cu.boxArray(), cu.DistributionMap(), nplot, 0);
    Vector<std::string> varNames(nplot);

    MultiFab::Copy(plotfile, cu, 0, cnt, 3, 0);
    varNames[cnt++] = "stretch";
    varNames[cnt++] = "mom";
    varNames[cnt++] = "energy";

    ShiftFaceToCC(cumom, 0, plotfile, cnt, 1);
    varNames[cnt++] = "mom_face_shift";

    MultiFab::Copy(plotfile, cuMeans, 0, cnt, 3, 0);
    varNames[cnt++] = "stretchMean";
    varNames[cnt++] = "momMean";
    varNames[cnt++] = "energyMean";

    ShiftFaceToCC(cumomMeans, 0, plotfile, cnt, 1);
    varNames[cnt++] = "momFaceMean";

    MultiFab::Copy(plotfile, cuVars, 0, cnt, 3, 0);
    varNames[cnt++] = "stretchVar";
    varNames[cnt++] = "momVar";
    varNames[cnt++] = "energyVar";

    ShiftFaceToCC(cumomVars, 0, plotfile, cnt, 1);
    varNames[cnt++] = "momFaceVar";

    MultiFab::Copy(plotfile, coVars, 0, cnt, 3, 0);
    varNames[cnt++] = "cov_stretch_mom";
    varNames[cnt++] = "cov_stretch_energy";
    varNames[cnt++] = "cov_mom_energy";

    const std::string plotfilename = amrex::Concatenate("plt", step, 9);
    WriteSingleLevelPlotfile(plotfilename, plotfile, varNames, geom, time, step);
}

void evaluateStats(const MultiFab& cons,
                         MultiFab& consMean,
                         MultiFab& consVar,
                         const MultiFab& cumom,
                         MultiFab& cumomMean,
                         MultiFab& cumomVar,
                         MultiFab& coVar,
                         const int steps,
                         const Geometry& geom)
{
    BL_PROFILE_VAR("evaluateStatsStag1D()", evaluateStatsStag1D);

    const Real stepsminusone = Real(steps - 1);
    const Real stepsinv = Real(1.0) / Real(steps);

    for (MFIter mfi(cons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        Array4<Real const> const& cu = cons.const_array(mfi);
        Array4<Real> const& mean = consMean.array(mfi);
        Array4<Real> const& var = consVar.array(mfi);
        Array4<Real> const& cov = coVar.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real x0 = cu(i,j,k,0);
            const Real x1 = cu(i,j,k,1);
            const Real x2 = cu(i,j,k,2);

            const Real m0_old = mean(i,j,k,0);
            const Real m1_old = mean(i,j,k,1);
            const Real m2_old = mean(i,j,k,2);

            const Real m0_new = m0_old + (x0 - m0_old) * stepsinv;
            const Real m1_new = m1_old + (x1 - m1_old) * stepsinv;
            const Real m2_new = m2_old + (x2 - m2_old) * stepsinv;

            var(i,j,k,0) = (stepsminusone * var(i,j,k,0)
                          + (x0 - m0_old) * (x0 - m0_new)) * stepsinv;
            var(i,j,k,1) = (stepsminusone * var(i,j,k,1)
                          + (x1 - m1_old) * (x1 - m1_new)) * stepsinv;
            var(i,j,k,2) = (stepsminusone * var(i,j,k,2)
                          + (x2 - m2_old) * (x2 - m2_new)) * stepsinv;

            cov(i,j,k,0) = (stepsminusone * cov(i,j,k,0)
                          + (x0 - m0_old) * (x1 - m1_new)) * stepsinv;
            cov(i,j,k,1) = (stepsminusone * cov(i,j,k,1)
                          + (x0 - m0_old) * (x2 - m2_new)) * stepsinv;
            cov(i,j,k,2) = (stepsminusone * cov(i,j,k,2)
                          + (x1 - m1_old) * (x2 - m2_new)) * stepsinv;

            mean(i,j,k,0) = m0_new;
            mean(i,j,k,1) = m1_new;
            mean(i,j,k,2) = m2_new;
        });
    }

    for (MFIter mfi(cumom, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& xbx = mfi.nodaltilebox(0);

        Array4<Real const> const& mx = cumom.const_array(mfi);
        Array4<Real> const& mx_mean = cumomMean.array(mfi);
        Array4<Real> const& mx_var = cumomVar.array(mfi);

        amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real x = mx(i,j,k,0);
            const Real m_old = mx_mean(i,j,k,0);
            const Real m_new = m_old + (x - m_old) * stepsinv;

            mx_var(i,j,k,0) = (stepsminusone * mx_var(i,j,k,0)
                             + (x - m_old) * (x - m_new)) * stepsinv;
            mx_mean(i,j,k,0) = m_new;
        });
    }

    consMean.FillBoundary(geom.periodicity());
    consVar.FillBoundary(geom.periodicity());
    cumomMean.FillBoundary(geom.periodicity());
    cumomVar.FillBoundary(geom.periodicity());
}

void GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max {100000};
    is.ignore(bl_ignore_max, '\n');
}

void WritePlotFilesSF_1D(const amrex::MultiFab& mag, const amrex::MultiFab& realimag,
                         const int step, const Real time, const amrex::Vector< std::string >& names, std::string plotfile_base) {

    // Magnitude of the Structure Factor
    std::string name = plotfile_base;
    name += "_mag";
    const std::string plotfilename1 = amrex::Concatenate(name,step,9);

    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // only needed to make a plotfile

    Geometry geom;

    Box domain_pencil = mag.boxArray().minimalBox();

    Vector<Real> projected_lo(AMREX_SPACEDIM);
    Vector<Real> projected_hi(AMREX_SPACEDIM);

    projected_lo[0] = -domain_pencil.length(0)/2 - 0.5;
    projected_hi[0] = domain_pencil.length(0)/2 - 1 + 0.5;

    projected_lo[1] = projected_lo[2] = -0.5;
    projected_hi[1] = projected_hi[2] =  0.5;

    RealBox real_box_pencil({AMREX_D_DECL(projected_lo[0],projected_lo[1],projected_lo[2])},
                            {AMREX_D_DECL(projected_hi[0],projected_hi[1],projected_hi[2])});

    geom.define(domain_pencil,&real_box_pencil,CoordSys::cartesian,is_periodic.data());

    Vector<std::string> varNames;
    varNames.resize(names.size());
    for (int n=0; n<names.size(); n++) {
        varNames[n] = names[n];
    }

    WriteSingleLevelPlotfile(plotfilename1,mag,varNames,geom,time,step);

    // Components of the Structure Factor
    name = plotfile_base;
    name += "_real_imag";
    const std::string plotfilename2 = amrex::Concatenate(name,step,9);

    varNames.resize(2*names.size());
    int cnt = 0; // keep a counter for plotfile variables
    for (int n=0; n<names.size(); n++) {
        varNames[cnt] = names[cnt];
        varNames[cnt] += "_real";
        cnt++;
    }

    int index = 0;
    for (int n=0; n<names.size(); n++) {
        varNames[cnt] = names[index];
        varNames[cnt] += "_imag";
        index++;
        cnt++;
    }

    WriteSingleLevelPlotfile(plotfilename2,realimag,varNames,geom,time,step);
}

//void refreshCellCenteredMomentum(MultiFab& cu,
//                                  const MultiFab& cumom,
//                                  const Geometry& geom)
//{
//    for (MFIter mfi(cu, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
//        const Box& bx = mfi.tilebox();
//        Array4<Real> const& u = cu.array(mfi);
//        Array4<Real const> const& mx = cumom.const_array(mfi);
//
//        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//        {
//            u(i,j,k,1) = 0.5 * (mx(i,j,k,0) + mx(i+1,j,k,0));
//        });
//    }
//    cu.FillBoundary(geom.periodicity());
//}
//
//void initializeState(MultiFab& cu,
//                     MultiFab& cumom,
//                     const Geometry& geom)
//{
//    cu.setVal(0.0, 0, 3, cu.nGrowVect());
//    cumom.setVal(0.0, 0, 1, cumom.nGrowVect());
//
//    MultiFabFillRandomNormal(cu, 0, 1, 0.0, 1.0, geom, true, true);
//    MultiFabFillRandomNormal(cu, 2, 1, 0.0, 1.0, geom, true, true);
//    MultiFabFillRandomNormal(cumom, 0, 1, 0.0, 1.0, geom, true, true);
//
//    cumom.FillBoundary(geom.periodicity());
//    refreshCellCenteredMomentum(cu, cumom, geom);
//    cu.FillBoundary(geom.periodicity());
//}

