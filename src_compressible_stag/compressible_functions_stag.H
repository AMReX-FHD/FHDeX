#ifndef _compressible_functions_stag_H_
#define _compressible_functions_stag_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "common_functions.H"
#include "compressible_functions.H"

using namespace amrex;
using namespace compressible;

void WritePlotFileStag(int step,
                       const Real time,
                       const Geometry geom,
                       const MultiFab& cu,
                       const MultiFab& cuMeans,
                       const MultiFab& cuVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                       const MultiFab& prim,
                       const MultiFab& primMeans,
                       const MultiFab& primVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velVars,
                       const MultiFab& coVars,
                       const MultiFab& eta, 
                       const MultiFab& kappa);

void WriteSpatialCross3D(const Vector<Real>& spatialCross, int step, const Geometry& geom, const int ncross);

void WriteSpatialCross1D(const MultiFab& spatialCross    , int step, const Geometry& geom, const int ncross); 

void WritePlotFilesSF_2D(const amrex::MultiFab& mag, const amrex::MultiFab& realimag, const amrex::Geometry& geom,
                         const int step, const Real time, const amrex::Vector< std::string >& names, std::string plotfile_base);

void conservedToPrimitiveStag(MultiFab& prim_in, std::array<MultiFab, AMREX_SPACEDIM>& velStag_in,
                              MultiFab& cons_in, const std::array<MultiFab, AMREX_SPACEDIM>& momStag_in);

void InitConsVarStag(MultiFab& cons,
                     std::array< MultiFab, AMREX_SPACEDIM >& momStag,
                     const amrex::Geometry geom);

void SetupBCStag();

void SetupCWallStag();

void setBCStag(MultiFab& prim_in, MultiFab& cons_in,
                 std::array< MultiFab, AMREX_SPACEDIM >& cumom_in,
                 std::array< MultiFab, AMREX_SPACEDIM >& vel_in,
                 const amrex::Geometry geom);

void BCWallReservoirFluxStag(std::array< MultiFab, AMREX_SPACEDIM >& faceflux,
                             std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in,
                             const amrex::Geometry geom);

void BCMem(MultiFab& prim_in, MultiFab& cons_in,
           std::array< MultiFab, AMREX_SPACEDIM >& cumom_in,
           std::array< MultiFab, AMREX_SPACEDIM >& vel_in,
           const amrex::Geometry geom);

void BCMassTempPress(MultiFab& prim_in,MultiFab& cons_in,const amrex::Geometry geom,int dim);

void BCMomNormal(MultiFab& mom_in, MultiFab& vel_in, MultiFab& cons_in,
                 const amrex::Geometry geom, int dim);

void BCMomTrans(MultiFab& mom_in, MultiFab& vel_in,
                 const amrex::Geometry geom, int dim);

void BCRhoRhoE(MultiFab& cons_in, MultiFab& prim_in, 
               std::array< MultiFab, AMREX_SPACEDIM >& cumom_in, 
               const amrex::Geometry geom);

void StochFluxStag(std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in, std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in, 
                   std::array< MultiFab, 2 >& edgeflux_x_in, std::array< MultiFab, 2 >& edgeflux_y_in, 
                   std::array< MultiFab, 2 >& edgeflux_z_in, const amrex::Geometry geom);

void StochFluxMem(std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in, std::array< MultiFab, 2 >& edgeflux_x_in,
                   std::array< MultiFab, 2 >& edgeflux_y_in, std::array< MultiFab, 2 >& edgeflux_z_in);

void RK3stepStag(MultiFab& cu, 
                 std::array< MultiFab, AMREX_SPACEDIM >& cumom,
                 MultiFab& prim, std::array< MultiFab, AMREX_SPACEDIM >& facevel,
                 MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 MultiFab& chi, MultiFab& D,
                 std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real dt, const int step);

void calculateFluxStag(const MultiFab& cons_in, const std::array< MultiFab, AMREX_SPACEDIM >& momStag_in, 
                       const MultiFab& prim_in, const std::array< MultiFab, AMREX_SPACEDIM >& velStag_in,
                       const MultiFab& eta_in, const MultiFab& zeta_in, const MultiFab& kappa_in,
                       const MultiFab& chi_in, const MultiFab& D_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
                       std::array< MultiFab, 2 >& edgeflux_x_in,
                       std::array< MultiFab, 2 >& edgeflux_y_in,
                       std::array< MultiFab, 2 >& edgeflux_z_in,
                       std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in,
                       std::array< MultiFab, AMREX_SPACEDIM>& stochface_in,
                       std::array< MultiFab, 2 >& stochedge_x_in,
                       std::array< MultiFab, 2 >& stochedge_y_in,
                       std::array< MultiFab, 2 >& stochedge_z_in,
                       std::array< MultiFab, AMREX_SPACEDIM>& stochcen_in,
                       const amrex::Geometry geom,
		                   const amrex::Vector< amrex::Real >& stoch_weights,
		                   const amrex::Real dt);

void doMembraneStag(MultiFab& cons, 
                    std::array< MultiFab, AMREX_SPACEDIM >& cumom,
                    MultiFab& prim, 
                    std::array< MultiFab, AMREX_SPACEDIM >& vel,
                    std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                    const amrex::Geometry geom, const amrex::Real dt);

void doLangevin(MultiFab& cons_in, MultiFab& prim_in,
                std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                const amrex::Geometry geom,
                const amrex::Real dt);

void applyEffusion(std::array<MultiFab, AMREX_SPACEDIM>& faceflux, MultiFab& cons);

void WriteCheckPoint3D(int step,
                       const amrex::Real time,
                       int statsCount,
                       const amrex::Geometry geom,
                       const amrex::MultiFab& cu,
                       const amrex::MultiFab& cuMeans,
                       const amrex::MultiFab& cuVars,
                       const amrex::MultiFab& prim,
                       const amrex::MultiFab& primMeans,
                       const amrex::MultiFab& primVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velVars,
                       const amrex::MultiFab& coVars,
                       const Vector<Real>& spatialCross, int ncross);

void ReadCheckPoint3D(int& step,
                      amrex::Real& time,
                      int& statsCount,
                      amrex::Geometry geom,
                      const amrex::Box& domain,
                      amrex::MultiFab& cu,
                      amrex::MultiFab& cuMeans,
                      amrex::MultiFab& cuVars,
                      amrex::MultiFab& prim,
                      amrex::MultiFab& primMeans,
                      amrex::MultiFab& primVars,
                      std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                      std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                      std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                      std::array<MultiFab, AMREX_SPACEDIM>& vel,
                      std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                      std::array<MultiFab, AMREX_SPACEDIM>& velVars,
                      amrex::MultiFab& coVars,
                      Vector<Real>& spatialCross, int ncross,
                      BoxArray& ba, DistributionMapping& dmap);

void WriteCheckPoint2D(int step,
                       const amrex::Real time,
                       int statsCount,
                       const amrex::Geometry geom,
                       const amrex::MultiFab& cu,
                       const amrex::MultiFab& cuMeans,
                       const amrex::MultiFab& cuVars,
                       const amrex::MultiFab& prim,
                       const amrex::MultiFab& primMeans,
                       const amrex::MultiFab& primVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velVars,
                       const amrex::MultiFab& coVars,
                       const amrex::MultiFab& spatialCross, int ncross);

void ReadCheckPoint2D(int& step,
                       amrex::Real& time,
                       int& statsCount,
                       amrex::Geometry geom,
                       const amrex::Box& domain,
                       amrex::MultiFab& cu,
                       amrex::MultiFab& cuMeans,
                       amrex::MultiFab& cuVars,
                       amrex::MultiFab& prim,
                       amrex::MultiFab& primMeans,
                       amrex::MultiFab& primVars,
                       std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                       std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                       std::array<MultiFab, AMREX_SPACEDIM>& vel,
                       std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                       std::array<MultiFab, AMREX_SPACEDIM>& velVars,
                       amrex::MultiFab& coVars,
                       amrex::MultiFab& spatialCross, int ncross,
                       BoxArray& ba, DistributionMapping& dmap);

void WriteCheckPoint1D(int step,
                       const amrex::Real time,
                       int statsCount,
                       const amrex::Geometry geom,
                       const amrex::MultiFab& cu,
                       const amrex::MultiFab& cuMeans,
                       const amrex::MultiFab& cuVars,
                       const amrex::MultiFab& prim,
                       const amrex::MultiFab& primMeans,
                       const amrex::MultiFab& primVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velVars,
                       const amrex::MultiFab& coVars,
                       const amrex::MultiFab& spatialCross, int ncross);

void ReadCheckPoint1D(int& step,
                       amrex::Real& time,
                       int& statsCount,
                       amrex::Geometry geom,
                       const amrex::Box& domain,
                       amrex::MultiFab& cu,
                       amrex::MultiFab& cuMeans,
                       amrex::MultiFab& cuVars,
                       amrex::MultiFab& prim,
                       amrex::MultiFab& primMeans,
                       amrex::MultiFab& primVars,
                       std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                       std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                       std::array<MultiFab, AMREX_SPACEDIM>& vel,
                       std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                       std::array<MultiFab, AMREX_SPACEDIM>& velVars,
                       amrex::MultiFab& coVars,
                       amrex::MultiFab& spatialCross, int ncross,
                       BoxArray& ba, DistributionMapping& dmap);

void ReadFile(const std::string& filename, Vector<char>& charBuf,
              bool bExitOnError=true);

void Read_Copy_MF_Checkpoint(amrex::MultiFab& mf, std::string mf_name, const std::string& checkpointname, 
                             BoxArray& ba_old, DistributionMapping& dmap_old,
                             int NVARS, int NGC, int nodal_flag=-1);

void PrintFluxes(std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in, std::array< MultiFab, 2 >& edgeflux_x_in,
                 std::array< MultiFab, 2 >& edgeflux_y_in, std::array< MultiFab, 2 >& edgeflux_z_in,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in, std::string prefix); 

void evaluateStatsStag3D(MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                         MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                         const std::array<MultiFab, AMREX_SPACEDIM>& vel, 
                         std::array<MultiFab, AMREX_SPACEDIM>& velMean, 
                         std::array<MultiFab, AMREX_SPACEDIM>& velVar, 
                         const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
                         MultiFab& coVar, 
                         Vector<Real>& dataSliceMeans_xcross,
                         Vector<Real>& spatialCross3D, const int ncross,
                         const amrex::Box& domain,
                         const int steps);

void evaluateStatsStag2D(MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                         MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                         const std::array<MultiFab, AMREX_SPACEDIM>& vel, 
                         std::array<MultiFab, AMREX_SPACEDIM>& velMean, 
                         std::array<MultiFab, AMREX_SPACEDIM>& velVar, 
                         const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
                         MultiFab& coVar, 
                         MultiFab& spatialCross2D, const int ncross,
                         const int steps);

void evaluateStatsStag1D(MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                         MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                         const std::array<MultiFab, AMREX_SPACEDIM>& vel, 
                         std::array<MultiFab, AMREX_SPACEDIM>& velMean, 
                         std::array<MultiFab, AMREX_SPACEDIM>& velVar, 
                         const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
                         MultiFab& coVar, 
                         MultiFab& spatialCross1D, const int ncross,
                         const int steps);

void EvaluateStatsMeans(MultiFab& cons, MultiFab& consMean,
                        MultiFab& prim_in, MultiFab& primMean,
                        const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                        std::array<MultiFab, AMREX_SPACEDIM>& velMean,
                        const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                        std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                        const int steps);

void EvaluateVarsCoVars(const MultiFab& cons, const MultiFab& consMean, MultiFab& consVar,
                        const MultiFab& prim_in, const MultiFab& primMean, MultiFab& primVar,
                        const std::array<MultiFab, AMREX_SPACEDIM>& velMean,
                        std::array<MultiFab, AMREX_SPACEDIM>& velVar,
                        const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                        const std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                        std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
                        MultiFab& coVar,
                        const int steps);

void GetSliceAverageCross(Vector<Real>& dataAvMeans_x,
                         Vector<Real>& dataAvMeans_xcross,
                         const MultiFab& consMean,
                         const MultiFab& primMean,
                         const MultiFab& prim_in,
                         const MultiFab& cons,
                         const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                         const std::array<MultiFab, AMREX_SPACEDIM>& velMean,
                         const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                         const std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                         const int nstats);

void GetPencilCross(amrex::Gpu::DeviceVector<Real>& data_xcross,
                    const MultiFab& consMean,
                    const MultiFab& primMean,
                    const MultiFab& prim_in,
                    const MultiFab& cons,
                    const int nstats,
                    const int x_star);

void EvaluateSpatialCorrelations3D(Vector<Real>& spatialCross,
                                   Vector<Real>& data_xcross,
                                   amrex::Gpu::HostVector<Real>& cu_avg, 
                                   amrex::Gpu::HostVector<Real>& cumeans_avg, 
                                   amrex::Gpu::HostVector<Real>& prim_avg, 
                                   amrex::Gpu::HostVector<Real>& primmeans_avg, 
                                   const int steps,
                                   const int nstats,
                                   const int ncross);

void EvaluateSpatialCorrelations1D(MultiFab& spatialCross1D,
                                   amrex::Gpu::DeviceVector<Real>& data_xcross,
                                   const MultiFab& consMean,
                                   const MultiFab& primMean,
                                   const MultiFab& prim_in,
                                   const MultiFab& cons,
                                   const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                                   const std::array<MultiFab, AMREX_SPACEDIM>& velMean,
                                   const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                                   const std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                                   const int steps,
                                   const int nstats,
                                   const int ncross,
                                   const int star_index);

#endif
