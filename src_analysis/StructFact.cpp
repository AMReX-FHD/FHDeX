#include "common_functions.H"
#include "StructFact.H"

#include <AMReX_MultiFabUtil.H>
#include "AMReX_PlotFileUtil.H"
#include "AMReX_BoxArray.H"

#include <AMReX_FFT.H>

// blank constructor
StructFact::StructFact()
{}

// this constructor takes in var_names, which contains the names of all variables under consideration
// we will compute the covariances of all possible pairs of variables
// var_scaling must be sized to match the total number of pairs of variables
StructFact::StructFact(const BoxArray& ba_in,
                       const DistributionMapping& dmap_in,
                       const Vector< std::string >& var_names,
                       const Vector< Real >& var_scaling_in,
                       const int& verbosity_in) {

    this->define(ba_in, dmap_in, var_names, var_scaling_in, verbosity_in);

}

// this constructor takes in var_names, which contains the names of all variables under consideration
// we will compute the covariances of the pairs of variables defined in s_pairA/B_in
// var_scaling must be sized to match the total number of pairs of variables
StructFact::StructFact(const BoxArray& ba_in,
                       const DistributionMapping& dmap_in,
                       const Vector< std::string >& var_names,
                       const Vector< Real >& var_scaling_in,
                       const Vector< int >& s_pairA_in,
                       const Vector< int >& s_pairB_in,
                       const int& verbosity_in) {

    this->define(ba_in, dmap_in, var_names, var_scaling_in, s_pairA_in, s_pairB_in, verbosity_in);

}

// this define takes in var_names, which contains the names of all variables under consideration
// we will compute the covariances of all possible pairs of variables
// var_scaling must be sized to match the total number of pairs of variables
void StructFact::define(const BoxArray& ba_in,
                       const DistributionMapping& dmap_in,
                       const Vector< std::string >& var_names,
                       const Vector< Real >& var_scaling_in,
                       const int& verbosity_in) {

    NVAR = var_names.size();

    Vector<int> l_s_pairA(NVAR * (NVAR + 1) / 2);
    Vector<int> l_s_pairB(NVAR * (NVAR + 1) / 2);

    int counter = 0;
    for (int i = 0; i < NVAR; ++i) {
        for (int j = i; j < NVAR; ++j) {
            l_s_pairA[counter] = j;
            l_s_pairB[counter] = i;
            ++counter;
        }
    }

    define(ba_in, dmap_in, var_names, var_scaling_in, l_s_pairA, l_s_pairB, verbosity_in);

}

// this define takes in var_names, which contains the names of all variables under consideration
// we will compute the covariances of the pairs of variables defined in s_pairA/B_in
// var_scaling must be sized to match the total number of pairs of variables
void StructFact::define(const BoxArray& ba_in,
                       const DistributionMapping& dmap_in,
                       const Vector< std::string >& var_names,
                       const Vector< Real >& var_scaling_in,
                       const Vector< int >& s_pairA_in,
                       const Vector< int >& s_pairB_in,
                       const int& verbosity_in) {

    BL_PROFILE_VAR("StructFact::define()", StructFactDefine);

    verbosity = verbosity_in;

    if (s_pairA_in.size() != s_pairB_in.size())
        amrex::Error("StructFact::define() - Must have an equal number of components");

    NVAR = var_names.size();

    //////////////////////////////////////////////////////
    // Reorder pair selecting vectors & error checking

    NCOV = s_pairA_in.size();

    if (NCOV != var_scaling_in.size())
        amrex::Error("StructFact::define() - NCOV != var_scaling_in.size()");

    scaling.resize(NCOV);
    for (int n = 0; n < NCOV; n++) {
        scaling[n] = 1.0 / var_scaling_in[n];
    }

    s_pairA.resize(NCOV);
    s_pairB.resize(NCOV);

    // Set vectors identifying covariance pairs
    for (int n = 0; n < NCOV; n++) {
        s_pairA[n] = s_pairA_in[n];
        s_pairB[n] = s_pairB_in[n];
    }

    // Create vector of unique variable indices to select which to take the FFT
    NVARU = 2 * NCOV;    // temporary before selecting unique variables
    amrex::Vector<int> varu_temp(NVARU);

    // "Stack" on vectors A and B
    int indx = 0;
    for (int n = 0; n < NCOV; n++) {
        varu_temp[indx] = s_pairA_in[n];
        indx++;
    }
    for (int n = 0; n < NCOV; n++) {
        varu_temp[indx] = s_pairB_in[n];
        indx++;
    }

    // Reorder: ascending order using "bubble sort"
    int loop = 1;
    int tot_iters = 0;
    int x_temp;
    while (loop == 1) {
        loop = 0;
        for (int n = 1; n < NVARU; n++) {
            if (varu_temp[n] < varu_temp[n - 1]) {
                x_temp = varu_temp[n];
                varu_temp[n] = varu_temp[n - 1];
                varu_temp[n - 1] = x_temp;
                loop = 1;
            }
        }
        tot_iters++;
        if (tot_iters > 2 * NVARU) {
            loop = 0;
            amrex::Error("StructFact::StructFact() - Bubble sort failed to converge");
        }
    }

    // Identify number of repeats
    int N_dup = 0;
    for (int n = 1; n < NVARU; n++) {
        if (varu_temp[n] == varu_temp[n - 1]) {
            N_dup = N_dup + 1;
        }
    }

    // Resize based on number of unique variables
    int N_u = NVARU - N_dup;
    var_u.resize(N_u);

    // Only select unique pairs
    indx = 0;
    var_u[indx] = varu_temp[0];
    for (int n = 1; n < NVARU; n++) {
        if (varu_temp[n] != varu_temp[n - 1]) {
            indx = indx + 1;
            var_u[indx] = varu_temp[n];
        }
    }

    // Update number of variables according to unique vars
    NVARU = N_u;

    for (int n = 0; n < NCOV; n++) {
        Print() << "SF pairs (" << n << ") = " << s_pairA[n] << ", " << s_pairB[n] << std::endl;
    }
    Print() << "SF numPairs = " << NCOV << std::endl;

    for (int n = 0; n < NCOV; n++) {
        if (s_pairA[n] < 0 || s_pairA[n] >= NVAR || s_pairB[n] < 0 || s_pairB[n] >= NVAR)
            amrex::Error("StructFact::StructFact() - Invalid pair select values: must be between 0 and (num of varibles - 1)");
    }
    //////////////////////////////////////////////////////

    // Note that we are defining with NO ghost cells

    cov_real.define(ba_in, dmap_in, NCOV, 0);
    cov_imag.define(ba_in, dmap_in, NCOV, 0);
    cov_mag.define(ba_in, dmap_in, NCOV, 0);
    cov_real.setVal(0.0);
    cov_imag.setVal(0.0);
    cov_mag.setVal(0.0);

    cov_names.resize(NCOV);
    std::string x;
    int cnt = 0;
    for (int n = 0; n < NCOV; n++) {
        x = "struct_fact";
        x += '_';
        x += var_names[s_pairB[n]];
        x += '_';
        x += var_names[s_pairA[n]];
        cov_names[cnt] = x;
        cnt++;
    }

    Box domain = ba_in.minimalBox();

    Vector<Real> kspace_lo(AMREX_SPACEDIM);
    Vector<Real> kspace_hi(AMREX_SPACEDIM);

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (domain.length(d) % 2 == 0) {
            // even number of cells
            kspace_lo[d] = -domain.length(d) / 2. - 0.5;
            kspace_hi[d] = domain.length(d) / 2. - 0.5;
        } else {
            // odd number of cells
            kspace_lo[d] = -domain.length(d) / 2.;
            kspace_hi[d] = domain.length(d) / 2.;
        }
    }

    RealBox kspace({AMREX_D_DECL(kspace_lo[0], kspace_lo[1], kspace_lo[2])},
                   {AMREX_D_DECL(kspace_hi[0], kspace_hi[1], kspace_hi[2])});

    // required only to define geom object
    Vector<int> is_periodic(AMREX_SPACEDIM, 1);

    geom_sf.define(domain, &kspace, CoordSys::cartesian, is_periodic.data());
}

void StructFact::FortStructure(const MultiFab& variables,
                               const int& reset)
{
    BL_PROFILE_VAR("StructFact::FortStructure()", FortStructure);

    const BoxArray& ba = variables.boxArray();
    const DistributionMapping& dm = variables.DistributionMap();

    MultiFab variables_dft_real(ba, dm, NVAR, 0);
    MultiFab variables_dft_imag(ba, dm, NVAR, 0);

    ComputeFFT(variables, variables_dft_real, variables_dft_imag);

    // temporary storage built on BoxArray and DistributionMapping of "variables"
    // One case where "variables" and "cov_real/imag/mag" may have different DistributionMappings
    // is for flattened MFs with one grid newly built flattened MFs may be on a different
    // processor than the flattened MF used to build cov_real/imag/mag
    // or in general, problems that are not perfectly load balanced
    MultiFab cov_temp;
    cov_temp.define(ba, dm, 1, 0);

    // temporary storage built on BoxArray and DistributionMapping of "cov_real/imag/mag"
    MultiFab cov_temp2;
    cov_temp2.define(cov_real.boxArray(), cov_real.DistributionMap(), 1, 0);

    int index = 0;
    for (int n = 0; n < NCOV; n++) {
        int i = s_pairA[n];
        int j = s_pairB[n];

        // Compute temporary real and imaginary components of covariance

        // Real component of covariance
        cov_temp.setVal(0.0);
        MultiFab::AddProduct(cov_temp, variables_dft_real, i, variables_dft_real, j, 0, 1, 0);
        MultiFab::AddProduct(cov_temp, variables_dft_imag, i, variables_dft_imag, j, 0, 1, 0);

        // copy into a MF with same ba and dm as cov_real/imag/mag
        cov_temp2.ParallelCopy(cov_temp, 0, 0, 1);

        if (reset == 1) {
            MultiFab::Copy(cov_real, cov_temp2, 0, index, 1, 0);
        } else {
            MultiFab::Add(cov_real, cov_temp2, 0, index, 1, 0);
        }

        // Imaginary component of covariance
        cov_temp.setVal(0.0);
        MultiFab::AddProduct(cov_temp, variables_dft_imag, i, variables_dft_real, j, 0, 1, 0);
        cov_temp.mult(-1.0, 0);
        MultiFab::AddProduct(cov_temp, variables_dft_real, i, variables_dft_imag, j, 0, 1, 0);

        // copy into a MF with same ba and dm as cov_real/imag/mag
        cov_temp2.ParallelCopy(cov_temp, 0, 0, 1);

        if (reset == 1) {
            MultiFab::Copy(cov_imag, cov_temp2, 0, index, 1, 0);
        } else {
            MultiFab::Add(cov_imag, cov_temp2, 0, index, 1, 0);
        }

        index++;
    }

    bool write_data = false;
    if (write_data) {
        std::string plotname;
        plotname = "a_VAR";
        VisMF::Write(variables, plotname);
        plotname = "a_COV_REAL";
        VisMF::Write(cov_real, plotname);
        plotname = "a_COV_IMAG";
        VisMF::Write(cov_imag, plotname);
    }

    if (reset == 1) {
        nsamples = 1;
    } else {
        nsamples++;
    }
}

void StructFact::Reset() {

    BL_PROFILE_VAR("StructFact::Reset()", StructFactReset);

    cov_real.setVal(0.);
    cov_imag.setVal(0.);
    nsamples = 0;
}

void StructFact::ComputeFFT(const MultiFab& variables,
                           MultiFab& variables_dft_real,
                           MultiFab& variables_dft_imag,
                           bool unpack)
{
    BL_PROFILE_VAR("StructFact::ComputeFFT()", ComputeFFT);

    Box domain = variables.boxArray().minimalBox();
    bool chopped_in_x = false;
    bool chopped_in_y = false;
    bool chopped_in_z = false;

    // figure out which direction the spectral box will be chopped
    if (domain.length(0) > 1) {
        chopped_in_x = true;
    } else if (domain.length(1) > 1) {
        chopped_in_y = true;
#if (AMREX_SPACEDIM == 3)
    } else if (domain.length(2) > 1) {
        chopped_in_z = true;
#endif
    } else {
        Abort("Calling ComputeFFT for a MultiFab with only 1 cell");
    }

    // compute number of points in the domain and the square root
    long npts = (AMREX_SPACEDIM == 2) ? (domain.length(0) * domain.length(1))
                                    : (domain.length(0) * domain.length(1) * domain.length(2));
    Real sqrtnpts = std::sqrt(npts);

    // extract BoxArray and DistributionMapping from variables
    BoxArray ba = variables.boxArray();
    DistributionMapping dm = variables.DistributionMap();

    // create storage for one component of variables
    MultiFab phi(ba, dm, 1, 0);

    // Initialize the boxarray "ba_onegrid" from the single box "domain"
    // Initilize a DistributionMapping for one grid
    BoxArray ba_onegrid(domain);
    DistributionMapping dm_onegrid(ba_onegrid);

    // create amrex::FFT object
    amrex::FFT::R2C my_fft(domain);

    // create storage for the FFT (distributed and single-grid)
    auto const& [ba_fft, dm_fft] = my_fft.getSpectralDataLayout();
    FabArray<BaseFab<GpuComplex<amrex::Real>>> phi_fft(ba_fft, dm_fft, 1, 0);

    Box domain_fft = ba_fft.minimalBox();
    BoxArray ba_fft_onegrid(domain_fft);
    FabArray<BaseFab<GpuComplex<amrex::Real>>> phi_fft_onegrid(ba_fft_onegrid, dm_onegrid, 1, 0);

    MultiFab variables_dft_real_onegrid(ba_onegrid, dm_onegrid, 1, 0);
    MultiFab variables_dft_imag_onegrid(ba_onegrid, dm_onegrid, 1, 0);

    // we will take one FFT at a time and copy the answer into the
    // corresponding component of variables_dft_real/imag
    for (int comp = 0; comp < NVAR; comp++) {
        bool comp_fft = false;
        for (int i = 0; i < NVARU; i++) {
            if (comp == var_u[i]) {
                comp_fft = true;
                break;
            }
        }

        if (comp_fft == false) continue;

        // copy component "comp" into a MultiFab with one component
        MultiFab::Copy(phi, variables, comp, 0, 1, 0);

        // ForwardTransform
        my_fft.forward(phi, phi_fft);

        // copy my_fft into a single-grid MultiFab
        phi_fft_onegrid.ParallelCopy(phi_fft, 0, 0, 1);

        // copy data to a full-sized MultiFab
        // this involves copying the complex conjugate from the half-sized field
        // into the appropriate place in the full MultiFab
        for (MFIter mfi(variables_dft_real_onegrid); mfi.isValid(); ++mfi) {
            Box bx = mfi.fabbox();

            Array4<GpuComplex<Real>> spectral = phi_fft_onegrid.array(mfi);
            Array4<Real> const& realpart = variables_dft_real_onegrid.array(mfi);
            Array4<Real> const& imagpart = variables_dft_imag_onegrid.array(mfi);

            /*
              Unpacking rules:

              For domains from (0,0,0) to (Nx-1,Ny-1,Nz-1) and chopped_in_x (i.e., Nx > 1)

              For any cells with i index > Nx/2, these values are complex conjugates of the corresponding
              entry where (Nx-i,Ny-j,Nz-k) UNLESS that index is zero, in which case you use 0.

              e.g. for an 8^3 domain, any cell with i index

              Cell (6,2,3) is complex conjugate of (2,6,5)

              Cell (4,1,0) is complex conjugate of (4,7,0)  (note that the FFT is computed for 0 <= i <= Nx/2)

              The analogy extends for the chopped_in_y and z directions
            */

            if (chopped_in_x) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i <= bx.length(0) / 2) {
                        // copy value
                        realpart(i, j, k) = spectral(i, j, k).real();
                        imagpart(i, j, k) = spectral(i, j, k).imag();
                    } else {
                        // copy complex conjugate
                        int iloc = bx.length(0) - i;
                        int jloc = (j == 0) ? 0 : bx.length(1) - j;
#if (AMREX_SPACEDIM == 2)
                        int kloc = 0;
#elif (AMREX_SPACEDIM == 3)
                        int kloc = (k == 0) ? 0 : bx.length(2) - k;
#endif
                        if (unpack) {
                            realpart(i, j, k) = spectral(iloc, jloc, kloc).real();
                            imagpart(i, j, k) = -spectral(iloc, jloc, kloc).imag();
                        }
                        else {
                            realpart(i, j, k) = 0.0;
                            imagpart(i, j, k) = 0.0;
                        }
                    }

                    realpart(i, j, k) /= sqrtnpts;
                    imagpart(i, j, k) /= sqrtnpts;
                });
            }

            if (chopped_in_y) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j <= bx.length(1) / 2) {
                        // copy value
                        realpart(i, j, k) = spectral(i, j, k).real();
                        imagpart(i, j, k) = spectral(i, j, k).imag();
                    } else {
                        // copy complex conjugate
                        int iloc = (i == 0) ? 0 : bx.length(0) - i;
                        int jloc = bx.length(1) - j;
#if (AMREX_SPACEDIM == 2)
                        int kloc = 0;
#elif (AMREX_SPACEDIM == 3)
                        int kloc = (k == 0) ? 0 : bx.length(2) - k;
#endif
                        if (unpack) {
                            realpart(i, j, k) = spectral(iloc, jloc, kloc).real();
                            imagpart(i, j, k) = -spectral(iloc, jloc, kloc).imag();
                        }
                        else {
                            realpart(i, j, k) = 0.0;
                            imagpart(i, j, k) = 0.0;
                        }
                    }

                    realpart(i, j, k) /= sqrtnpts;
                    imagpart(i, j, k) /= sqrtnpts;
                });
            }

            if (chopped_in_z) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k <= bx.length(2) / 2) {
                        realpart(i, j, k) = spectral(i, j, k).real();
                        imagpart(i, j, k) = spectral(i, j, k).imag();
                    } else {
                        int iloc = (i == 0) ? 0 : bx.length(0) - i;
                        int jloc = (j == 0) ? 0 : bx.length(1) - j;
                        int kloc = bx.length(2) - k;

                        if (unpack) {
                            realpart(i, j, k) = spectral(iloc, jloc, kloc).real();
                            imagpart(i, j, k) = -spectral(iloc, jloc, kloc).imag();
                        }
                        else {
                            realpart(i, j, k) = 0.0;
                            imagpart(i, j, k) = 0.0;
                        }
                    }

                    realpart(i, j, k) /= sqrtnpts;
                    imagpart(i, j, k) /= sqrtnpts;
                });
            }
        } // end MFIter

        variables_dft_real.ParallelCopy(variables_dft_real_onegrid, 0, comp, 1);
        variables_dft_imag.ParallelCopy(variables_dft_imag_onegrid, 0, comp, 1);
    }
}

void StructFact::WritePlotFile(const int step, const Real time,
                               std::string plotfile_base,
                               const int& zero_avg) {

    BL_PROFILE_VAR("StructFact::WritePlotFile()", StructFactWritePlotFile);

    MultiFab plotfile;
    Vector<std::string> varNames;
    int nPlot = 1;

    // Build temp real & imag components
    const BoxArray& ba = cov_mag.boxArray();
    const DistributionMapping& dm = cov_mag.DistributionMap();

    MultiFab cov_real_temp(ba, dm, NCOV, 0);
    MultiFab cov_imag_temp(ba, dm, NCOV, 0);
    MultiFab::Copy(cov_real_temp, cov_real, 0, 0, NCOV, 0);
    MultiFab::Copy(cov_imag_temp, cov_imag, 0, 0, NCOV, 0);

    // Finalize covariances - scale & compute magnitude
    Finalize(cov_real_temp, cov_imag_temp, zero_avg);

    //////////////////////////////////////////////////////////////////////////////////
    // Write out structure factor magnitude to plot file
    //////////////////////////////////////////////////////////////////////////////////

    std::string name = plotfile_base;
    name += "_mag";

    const std::string plotfilename1 = amrex::Concatenate(name, step, 9);
    nPlot = NCOV;
    plotfile.define(cov_mag.boxArray(), cov_mag.DistributionMap(), nPlot, 0);
    varNames.resize(nPlot);

    for (int n = 0; n < NCOV; n++) {
        varNames[n] = cov_names[n];
    }

    MultiFab::Copy(plotfile, cov_mag, 0, 0, NCOV, 0); // copy structure factor into plotfile

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename1, plotfile, varNames, geom_sf, time, step);

    //////////////////////////////////////////////////////////////////////////////////
    // Write out real and imaginary components of structure factor to plot file
    //////////////////////////////////////////////////////////////////////////////////

    name = plotfile_base;
    name += "_real_imag";

    const std::string plotfilename2 = amrex::Concatenate(name, step, 9);
    nPlot = 2 * NCOV;
    plotfile.define(cov_mag.boxArray(), cov_mag.DistributionMap(), nPlot, 0);
    varNames.resize(nPlot);

    int cnt = 0; // keep a counter for plotfile variables
    for (int n = 0; n < NCOV; n++) {
        varNames[cnt] = cov_names[cnt];
        varNames[cnt] += "_real";
        cnt++;
    }

    int index = 0;
    for (int n = 0; n < NCOV; n++) {
        varNames[cnt] = cov_names[index];
        varNames[cnt] += "_imag";
        index++;
        cnt++;
    }

    MultiFab::Copy(plotfile, cov_real_temp, 0, 0, NCOV, 0);
    MultiFab::Copy(plotfile, cov_imag_temp, 0, NCOV, NCOV, 0);

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename2, plotfile, varNames, geom_sf, time, step);
}

void StructFact::Finalize(MultiFab& cov_real_in, MultiFab& cov_imag_in,
                          const int& zero_avg) {

    BL_PROFILE_VAR("StructFact::Finalize()", StructFactFinalize);

    Real nsamples_inv = 1.0/(Real)nsamples;

    ShiftFFT(cov_real_in, zero_avg);
    ShiftFFT(cov_imag_in, zero_avg);

    cov_real_in.mult(nsamples_inv);
    for (int d = 0; d < NCOV; d++) {
        cov_real_in.mult(scaling[d], d, 1);
    }

    cov_imag_in.mult(nsamples_inv);
    for (int d = 0; d < NCOV; d++) {
        cov_imag_in.mult(scaling[d], d, 1);
    }

    cov_mag.setVal(0.0);
    MultiFab::AddProduct(cov_mag, cov_real_in, 0, cov_real_in, 0, 0, NCOV, 0);
    MultiFab::AddProduct(cov_mag, cov_imag_in, 0, cov_imag_in, 0, 0, NCOV, 0);

    SqrtMF(cov_mag);
}

// Finalize covariances - scale & compute magnitude
void StructFact::CallFinalize(const int& zero_avg)
{
    BL_PROFILE_VAR("CallFinalize()", CallFinalize);

    // Build temp real & imag components
    const BoxArray& ba = cov_mag.boxArray();
    const DistributionMapping& dm = cov_mag.DistributionMap();

    MultiFab cov_real_temp(ba, dm, NCOV, 0);
    MultiFab cov_imag_temp(ba, dm, NCOV, 0);
    MultiFab::Copy(cov_real_temp, cov_real, 0, 0, NCOV, 0);
    MultiFab::Copy(cov_imag_temp, cov_imag, 0, 0, NCOV, 0);

    // Finalize covariances - scale & compute magnitude
    Finalize(cov_real_temp, cov_imag_temp, zero_avg);
}

void StructFact::ShiftFFT(MultiFab& dft_out, const int& zero_avg) {

    BL_PROFILE_VAR("StructFact::ShiftFFT()", ShiftFFT);

    /*
        Shifting rules:

        For each direction d, shift the data in the +d direction by n_cell[d]/2 and then modulo by n_cell[d].

        e.g. for an 8^3 domain

        Cell (7,2,3) is shifted to ( (7+4)%8, (2+4)%8, (3+4)%8 ) =  (3,6,7)

    */
    BoxArray ba_onegrid;
    {
        Box domain = dft_out.boxArray().minimalBox();

        // Initialize the boxarray "ba" from the single box "bx"
        ba_onegrid.define(domain);
    }

    DistributionMapping dmap_onegrid(ba_onegrid);

    MultiFab dft_onegrid;
    MultiFab dft_onegrid_temp;
    dft_onegrid.define(ba_onegrid, dmap_onegrid, 1, 0);
    dft_onegrid_temp.define(ba_onegrid, dmap_onegrid, 1, 0);

    for (int d = 0; d < NCOV; d++) {
        dft_onegrid_temp.ParallelCopy(dft_out, d, 0, 1);

        if (zero_avg == 1) {
            for (MFIter mfi(dft_onegrid); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& dft_temp = dft_onegrid_temp.array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i == 0 && j == 0 && k == 0) {
                        dft_temp(i, j, k) = 0.;
                    }
                });
            }
        }

        // Shift DFT by N/2+1 (pi)
        for (MFIter mfi(dft_onegrid); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            const Array4<Real>& dft = dft_onegrid.array(mfi);
            const Array4<Real>& dft_temp = dft_onegrid_temp.array(mfi);

            int nx = bx.length(0);
            int nxh = nx / 2;
            int ny = bx.length(1);
            int nyh = ny / 2;
#if (AMREX_SPACEDIM == 3)
            int nz = bx.length(2);
            int nzh = nz / 2;
#endif

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int ip = (i + nxh) % nx;
                int jp = (j + nyh) % ny;
                int kp = 0;
#if (AMREX_SPACEDIM == 3)
                kp = (k + nzh) % nz;
#endif
                dft(ip, jp, kp) = dft_temp(i, j, k);
            });
        }

        dft_out.ParallelCopy(dft_onegrid, 0, d, 1);
    }
}

// integrate cov_mag over k shells
void StructFact::IntegratekShells(const int& step, const std::string& name) {

    BL_PROFILE_VAR("StructFact::IntegratekShells", IntegratekShells);

    GpuArray<int, AMREX_SPACEDIM> center;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        center[d] = n_cells[d] / 2;
    }

    int npts = n_cells[0] / 2;
    //int npts_sq = npts*npts;

    Gpu::DeviceVector<Real> phisum_device(npts);
    Gpu::DeviceVector<int> phicnt_device(npts);

    Gpu::HostVector<Real> phisum_host(npts);

    Real* phisum_ptr = phisum_device.dataPtr();  // pointer to data
    int* phicnt_ptr = phicnt_device.dataPtr();  // pointer to data

    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        phisum_ptr[d] = 0.;
        phicnt_ptr[d] = 0;
    });


//    Gpu::DeviceVector<Real> phisum_device_large(npts_sq);
//    Gpu::DeviceVector<int>  phicnt_device_large(npts_sq);

//    for (int d=0; d<npts_sq; ++d) {
//        phisum_device_large[d] = 0.;
//        phicnt_device_large[d] = 0;
//    }

  //  Real* phisum_large_gpu = phisum_device_large.dataPtr();  // pointer to data
 //   int*  phicnt_large_gpu = phicnt_device_large.dataPtr();  // pointer to data

    // only consider cells that are within 15k of the center point

    for (MFIter mfi(cov_mag, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        const Array4<Real>& cov = cov_mag.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ilen = amrex::Math::abs(i - center[0]);
            int jlen = amrex::Math::abs(j - center[1]);
            int klen = (AMREX_SPACEDIM == 3) ? amrex::Math::abs(k - center[2]) : 0;

            Real dist = (ilen*ilen + jlen*jlen + klen*klen);
 //           int idist = (ilen*ilen + jlen*jlen + klen*klen);
            dist = std::sqrt(dist);

            if (dist <= center[0] - 0.5) {
                dist = dist + 0.5;
                int cell = int(dist);
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    amrex::HostDevice::Atomic::Add(&(phisum_ptr[cell]), cov(i, j, k, d));
//                    phisum_large_gpu[idist]  += cov(i,j,k,d);
                }
                amrex::HostDevice::Atomic::Add(&(phicnt_ptr[cell]), 1);
 //               ++phicnt_large_gpu[idist];
            }
        });
    }

    for (int d = 1; d < npts; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_device[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_device[d]);
    }

#if 0
    for (int d=1; d<npts_sq; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_device_large[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_device_large[d]);
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb_disc;
        std::string turbNamedisc = "turb_disc";
        turbNamedisc += std::to_string(step);
        turbNamedisc += ".txt";

        turb_disc.open(turbNamedisc);
        for (int d=1; d<npts_sq; ++d) {
            if(phicnt_device_large[d]>0) {
                Real dreal = d;
                turb_disc << sqrt(dreal) << " " << 4.*M_PI*d*phisum_device_large[d]/phicnt_device_large[d] << std::endl;
            }
        }
    }

    Real dk = 1.;
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb_alt;
        std::string turbNamealt = "turb_alt";
        turbNamealt += std::to_string(step);
        turbNamealt += ".txt";

        turb_alt.open(turbNamealt);
        for (int d=1; d<npts; ++d) {
            turb_alt << d << " " << phisum_device[d] << std::endl;
        }
    }
#endif

    Real dk = 1.;

#if (AMREX_SPACEDIM == 2)
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0) {
            // phisum_ptr[d] *= 2.*M_PI*d*dk*dk/phicnt_ptr[d];
            phisum_ptr[d] *= 2.*M_PI*(d*dk+.5*dk*dk)/phicnt_ptr[d];
        }
    });
#else
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0) {
            // phisum_ptr[d] *= 4.*M_PI*(d*d)*dk*dk*dk/phicnt_ptr[d];
            // phisum_ptr[d] *= 4.*M_PI*(d*d*dk+d*dk*dk+dk*dk*dk/3.)/phicnt_ptr[d];
            phisum_ptr[d] *= 4. * M_PI * (d * d * dk + dk * dk * dk / 12.) / phicnt_ptr[d];
        }
    });
#endif

    Gpu::copy(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb;
        std::string turbBaseName = "turb_";
        turbBaseName += name;
        //turbName += std::to_string(step);
        std::string turbName = Concatenate(turbBaseName, step, 7);
        turbName += ".txt";

        turb.open(turbName);
        for (int d = 1; d < npts; ++d) {
            turb << d << " " << phisum_host[d] << std::endl;
        }
    }
}

// integrate cov_mag over k shells for scalar qtys
void StructFact::IntegratekShellsScalar(const int& step,
                                        const amrex::Vector<std::string>& names) {

    BL_PROFILE_VAR("StructFact::IntegratekShellsMisc", IntegratekShellsMisc);

    int turbvars = NVAR;

    if (names.size() != turbvars) amrex::Error("StructFact::IntegratekShellsMisc() requires names to be of the same length as NVAR");

    GpuArray<int, AMREX_SPACEDIM> center;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        center[d] = n_cells[d] / 2;
    }

    int npts = n_cells[0] / 2;

    Gpu::DeviceVector<Real> phisum_device(npts);
    Gpu::DeviceVector<int> phicnt_device(npts);

    Gpu::HostVector<Real> phisum_host(npts);

    Real* phisum_ptr = phisum_device.dataPtr();  // pointer to data
    int* phicnt_ptr = phicnt_device.dataPtr();  // pointer to data

    for (int var_ind = 0; var_ind < turbvars; var_ind++) {

        amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
        {
            phisum_ptr[d] = 0.;
            phicnt_ptr[d] = 0;
        });

        // only consider cells that are within 15k of the center point

        for (MFIter mfi(cov_mag, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            const Array4<Real>& cov = cov_mag.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int ilen = amrex::Math::abs(i - center[0]);
                int jlen = amrex::Math::abs(j - center[1]);
                int klen = (AMREX_SPACEDIM == 3) ? amrex::Math::abs(k - center[2]) : 0;

                Real dist = (ilen*ilen + jlen*jlen + klen*klen);
                dist = std::sqrt(dist);

                if (dist <= center[0] - 0.5) {
                    dist = dist + 0.5;
                    int cell = int(dist);
                    amrex::HostDevice::Atomic::Add(&(phisum_ptr[cell]), cov(i, j, k, var_ind));
                    amrex::HostDevice::Atomic::Add(&(phicnt_ptr[cell]), 1);
                }
            });
        }

        for (int d = 1; d < npts; ++d) {
            ParallelDescriptor::ReduceRealSum(phisum_device[d]);
            ParallelDescriptor::ReduceIntSum(phicnt_device[d]);
        }

        Real dk = 1.;

#if (AMREX_SPACEDIM == 2)
#if 0
        for (int d=1; d<npts; ++d) {
            //  phisum_vect[d] *= 2.*M_PI*d*dk*dk/phicnt_vect[d];
            phisum_vect[d] *= 2.*M_PI*(d*dk+.5*dk*dk)/phicnt_vect[d];
        }
#endif
        amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
        {
            if (d != 0) {
                // phisum_ptr[d] *= 2.*M_PI*d*dk*dk/phicnt_ptr[d];
                phisum_ptr[d] *= 2.*M_PI*(d*dk+.5*dk*dk)/phicnt_ptr[d];
            }
        });
#else
#if 0
        for (int d=1; d<npts; ++d) {
            //  phisum_vect[d] *= 4.*M_PI*(d*d)*dk*dk*dk/phicnt_vect[d];
            //  phisum_vect[d] *= 4.*M_PI*(d*d*dk+d*dk*dk+dk*dk*dk/3.)/phicnt_vect[d];
            phisum_vect[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_vect[d];
        }
#endif
        amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
        {
            if (d != 0) {
            // phisum_ptr[d] *= 4.*M_PI*(d*d)*dk*dk*dk/phicnt_ptr[d];
            // phisum_ptr[d] *= 4.*M_PI*(d*d*dk+d*dk*dk+dk*dk*dk/3.)/phicnt_ptr[d];
                phisum_ptr[d] *= 4. * M_PI * (d * d * dk + dk * dk * dk / 12.) / phicnt_ptr[d];
            }
        });
#endif
        Gpu::copy(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());

        if (ParallelDescriptor::IOProcessor()) {
            std::ofstream turb;
            std::string turbBaseName = "turb_" + names[var_ind];
            std::string turbName = Concatenate(turbBaseName, step, 7);
            turbName += ".txt";

            turb.open(turbName);
            for (int d = 1; d < npts; ++d) {
                turb << d << " " << phisum_host[d] << std::endl;
            }
        }
    }
}

void StructFact::AddToExternal(MultiFab& x_mag, MultiFab& x_realimag, const int& zero_avg) {

    BL_PROFILE_VAR("StructFact::AddToExternal", AddToExternal);

    MultiFab plotfile;
    int nPlot = 1;

    // Build temp real & imag components
    const BoxArray& ba = cov_mag.boxArray();
    const DistributionMapping& dm = cov_mag.DistributionMap();

    MultiFab cov_real_temp(ba, dm, NCOV, 0);
    MultiFab cov_imag_temp(ba, dm, NCOV, 0);
    MultiFab::Copy(cov_real_temp, cov_real, 0, 0, NCOV, 0);
    MultiFab::Copy(cov_imag_temp, cov_imag, 0, 0, NCOV, 0);

    // Finalize covariances - scale & compute magnitude
    Finalize(cov_real_temp, cov_imag_temp, zero_avg);

    nPlot = NCOV;
    plotfile.define(cov_mag.boxArray(), cov_mag.DistributionMap(), nPlot, 0);
    MultiFab::Copy(plotfile, cov_mag, 0, 0, NCOV, 0); // copy structure factor into plotfile
    MultiFab::Add(x_mag, plotfile, 0, 0, NCOV, 0);

    nPlot = 2 * NCOV;
    plotfile.define(cov_mag.boxArray(), cov_mag.DistributionMap(), nPlot, 0);
    MultiFab::Copy(plotfile, cov_real_temp, 0, 0, NCOV, 0);
    MultiFab::Copy(plotfile, cov_imag_temp, 0, NCOV, NCOV, 0);
    MultiFab::Add(x_realimag, plotfile, 0, 0, 2 * NCOV, 0);

}

void StructFact::WriteCheckPoint(const int& step,
                                 std::string checkfile_base)
{
    // checkpoint file name, e.g., chk_SF0000010 (digits is how many digits...)
    const std::string& checkpointname = amrex::Concatenate(checkfile_base, step, 9);

    amrex::Print() << "Writing structure factor checkpoint " << checkpointname << "\n";

    BoxArray ba = cov_real.boxArray();

    // single level problem
    int nlevels = 1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
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

        if (!HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        // write out title line
        HeaderFile << "Structure factor checkpoint file\n";

        // write out misc structure factor member objects
        HeaderFile << NVAR << "\n";
        HeaderFile << NVARU << "\n";
        HeaderFile << NCOV << "\n";
        HeaderFile << nsamples << "\n";
        for (int i = 0; i < NCOV; ++i) {
            HeaderFile << scaling[i] << "\n";
        }
        for (int i = 0; i < NCOV; ++i) {
            HeaderFile << cov_names[i] << "\n";
        }
        for (int i = 0; i < NCOV; ++i) {
            HeaderFile << s_pairA[i] << "\n";
        }
        for (int i = 0; i < NCOV; ++i) {
            HeaderFile << s_pairB[i] << "\n";
        }
        for (int i = 0; i < NVARU; ++i) {
            HeaderFile << var_u[i] << "\n";
        }

        /*
        // write the BoxArray
        ba.writeOn(HeaderFile);
        HeaderFile << '\n';
        */
    }

    // write the MultiFab data to, e.g., chk_SF00010/Level_0/
    VisMF::Write(cov_real,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cov_real"));
    VisMF::Write(cov_imag,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cov_imag"));
    VisMF::Write(cov_mag,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cov_mag"));
}

namespace {
    void GotoNextLine(std::istream& is)
    {
        constexpr std::streamsize bl_ignore_max{100000};
        is.ignore(bl_ignore_max, '\n');
    }
}

void StructFact::ReadCheckPoint(std::string checkfile_base,
                                BoxArray& ba_in,
                                DistributionMapping& dmap_in)
{
    const std::string& checkpointname = amrex::Concatenate(checkfile_base, restart, 9);

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

        // write out misc structure factor member objects
        is >> NVAR;
        GotoNextLine(is);
        is >> NVARU;
        GotoNextLine(is);
        is >> NCOV;
        GotoNextLine(is);
        is >> nsamples;
        GotoNextLine(is);

        scaling.resize(NCOV);
        cov_names.resize(NCOV);
        s_pairA.resize(NCOV);
        s_pairB.resize(NCOV);
        var_u.resize(NVARU);

        for (int i = 0; i < NCOV; ++i) {
            is >> scaling[i];
            GotoNextLine(is);
        }
        for (int i = 0; i < NCOV; ++i) {
            is >> cov_names[i];
            GotoNextLine(is);
        }
        for (int i = 0; i < NCOV; ++i) {
            is >> s_pairA[i];
            GotoNextLine(is);
        }
        for (int i = 0; i < NCOV; ++i) {
            is >> s_pairB[i];
            GotoNextLine(is);
        }
        for (int i = 0; i < NVARU; ++i) {
            is >> var_u[i];
            GotoNextLine(is);
        }

        // read in level 'lev' BoxArray from Header
        /*
        ba_in.readFrom(is);
        GotoNextLine(is);
        */

        // build MultiFab data
        cov_real.define(ba_in, dmap_in, NCOV, 0);
        cov_imag.define(ba_in, dmap_in, NCOV, 0);
        cov_mag.define(ba_in, dmap_in, NCOV, 0);
    }

    // read in the MultiFab data
    VisMF::Read(cov_real,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cov_real"));
    VisMF::Read(cov_imag,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cov_imag"));
    VisMF::Read(cov_mag,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cov_mag"));
}
