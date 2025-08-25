#include "hydro_functions.H"
#include "rng_functions.H"
#include "gmres_functions.H"

// used to add noise to an initial momentum field
void addMomFluctuations(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                        const MultiFab& rhotot, const MultiFab& Temp,
                        const amrex::Real& variance,
                        const Geometry& geom) {

    BL_PROFILE_VAR("addMomFluctuations()",addMomFluctuations);

    std::array< MultiFab, AMREX_SPACEDIM > m_old;
    std::array< MultiFab, AMREX_SPACEDIM > rhotot_fc;
    std::array< MultiFab, AMREX_SPACEDIM > Temp_fc;

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        m_old[d].define(     umac[d].boxArray(), umac[d].DistributionMap(), 1, 0);
        rhotot_fc[d].define( umac[d].boxArray(), umac[d].DistributionMap(), 1, 0);
        Temp_fc[d].define(   umac[d].boxArray(), umac[d].DistributionMap(), 1, 0);
    }

    // NOTE: these only operate on valid cells
    AverageCCToFace(rhotot, rhotot_fc, 0, 1, RHO_BC_COMP, geom);
    AverageCCToFace(Temp,   Temp_fc,   0, 1, TEMP_BC_COMP, geom);

    // Convert umac to momenta, rho*umac
    ConvertMToUmac(rhotot_fc,umac,m_old,0);

    addMomFluctuations_stag(m_old, rhotot_fc, Temp_fc, variance, geom);

    // Convert momenta to umac, (1/rho)*momentum
    ConvertMToUmac(rhotot_fc,umac,m_old,1);
}

// used to add noise to an initial momentum field
// called by addMomFluctuations
void addMomFluctuations_stag(std::array< MultiFab, AMREX_SPACEDIM >& m_old,
                             const std::array< MultiFab, AMREX_SPACEDIM >& rhotot_fc,
                             const std::array< MultiFab, AMREX_SPACEDIM >& Temp_fc,
                             const amrex::Real& variance,
                             const Geometry& geom) {

    BL_PROFILE_VAR("addMomFluctuations_stag()",addMomFluctuations_stag);

    const Real* dx = geom.CellSize();
    Real dVol = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2];

    // Initialize variances
    Real variance_mom = amrex::Math::abs(variance)*k_B/dVol;

    std::array<MultiFab, AMREX_SPACEDIM> variance_mfab;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        variance_mfab[d].define(m_old[d].boxArray(), m_old[d].DistributionMap(),1,0);
    }

    std::array< MultiFab, AMREX_SPACEDIM > mac_temp;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        mac_temp[d].define(m_old[d].boxArray(), m_old[d].DistributionMap(),1,0);
    }

    // Fill momentum multifab with random numbers, scaled by equilibrium variances
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // Set variance multifab to sqrt(rho*temp)
        MultiFab::Copy(     variance_mfab[d],rhotot_fc[d],0,0,1,0);
        MultiFab::Multiply( variance_mfab[d],Temp_fc[d],  0,0,1,0);
        SqrtMF(variance_mfab[d]);

        // Fill momentum with random numbers, scaled by sqrt(var*k_B/dV)
        MultiFabFillRandom(mac_temp[d],0,variance_mom,geom);

        // Scale random momenta further by factor of sqrt(rho*temp)
        MultiFab::Multiply(mac_temp[d],variance_mfab[d],0,0,1,0);

        MultiFab::Saxpy(m_old[d], 1.0, mac_temp[d],0,0,1,0);

        // For safety, although called by MultiFabFillRandom()
        m_old[d].OverrideSync(geom.periodicity());
        m_old[d].FillBoundary(geom.periodicity());
    }

    if (variance < 0.0) {
        // Ensure zero total momentum
        Vector<Real> av_mom;
        // take staggered sum & divide by number of cells
        SumStag(m_old,av_mom,true);
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            // subtract off average
            m_old[d].plus(-av_mom[d],1);
            m_old[d].OverrideSync(geom.periodicity());
            m_old[d].FillBoundary(geom.periodicity());
        }
    }

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        m_old[i].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(m_old[i], geom,i);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(m_old[i], geom, i, is_inhomogeneous);
    }
}
