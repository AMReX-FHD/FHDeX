#ifndef _compressible_functions_H_
#define _compressible_functions_H_

#define LOHI 2

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "compressible_functions_F.H"
#include "compressible_namespace.H"

#include "common_functions.H"

using namespace compressible;
using namespace amrex;

// copy contents of compressible_params_module to C++ compressible namespace
void InitializeCompressibleNamespace();

void calculateFlux(const MultiFab& cons, const MultiFab& prim,
                   const MultiFab& eta, const MultiFab& zeta, const MultiFab& kappa,
		   const MultiFab& chi, const MultiFab& D,
                   std::array<MultiFab, AMREX_SPACEDIM>& flux,
                   std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
                   std::array<MultiFab, AMREX_SPACEDIM>& cornx,
                   std::array<MultiFab, AMREX_SPACEDIM>& corny,
                   std::array<MultiFab, AMREX_SPACEDIM>& cornz,
                   MultiFab& visccorn,
                   MultiFab& rancorn,
                   const Geometry geom,
		   const Vector< Real >& stoch_weights,
	           const Real* dx, const Real dt);

void calculateTransportCoeffs(const MultiFab& prim_in,
			      MultiFab& eta_in, MultiFab& zeta_in, MultiFab& kappa_in,
			      MultiFab& chi_in, MultiFab& Dij_in);

void RK3step(MultiFab& cu, MultiFab& cup, MultiFab& cup2, MultiFab& cup3,
             MultiFab& prim, MultiFab& source,
             MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
             MultiFab& chi, MultiFab& D,
             std::array<MultiFab, AMREX_SPACEDIM>& flux,
             std::array<MultiFab, AMREX_SPACEDIM>& stochFlux, 
             std::array<MultiFab, AMREX_SPACEDIM>& cornx,
             std::array<MultiFab, AMREX_SPACEDIM>& corny,
             std::array<MultiFab, AMREX_SPACEDIM>& cornz,
             MultiFab& visccorn, MultiFab& rancorn,
             const Geometry geom, const Real* dx, const Real dt);

void conservedToPrimitive(MultiFab& prim_in, const MultiFab& cons_in);

void primitiveToConserved(const MultiFab& prim, MultiFab& cons);

void setBC(MultiFab& prim, MultiFab& cons);

void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                   const MultiFab& prim, MultiFab& primMean, MultiFab& primVar,
                   MultiFab& spatialCross,
                   Real* delHolder1, Real* delHolder2, Real* delHolder3,
                   Real* delHolder4, Real* delHolder5, Real* delHolder6,
                   const int steps, const Real* dx);

void doMembrane(MultiFab& cons, MultiFab& prim, std::array<MultiFab, AMREX_SPACEDIM>& flux,
                const Geometry geom, const Real* dx, const Real dt);

void WritePlotFile(int step,
                   const Real time,
                   const Geometry geom,
	           const MultiFab& cu,
	           const MultiFab& cuMeans,
	           const MultiFab& cuVars,
	           const MultiFab& prim,
	           const MultiFab& primMeans,
	           const MultiFab& primVars,
                   const MultiFab& spatialCross,
                   const MultiFab& eta, 
		   const MultiFab& kappa);

void yzAverage(const MultiFab& consMean, 
               const MultiFab& consVar, 
               const MultiFab& primMean, 
               const MultiFab& primVar, 
               const MultiFab& spatialCross, 
               MultiFab& consMeanAv, 
               MultiFab& consVarAv, 
               MultiFab& primMeanAv, 
               MultiFab& primVarAv, 
               MultiFab& spatialCrossAv);

void WriteCheckPoint(int step,
                     const Real time,
                     int statsCount,
                     const Geometry geom,
                     const MultiFab& cu,
                     const MultiFab& cuMeans,
                     const MultiFab& cuVars,
                     const MultiFab& prim,
                     const MultiFab& primMeans,
                     const MultiFab& primVars,
                     const MultiFab& spatialCross, 
                     const MultiFab& eta, 
                     const MultiFab& kappa);

void ReadCheckPoint(int& step,
                     Real& time,
                     int& statsCount,
                     Geometry geom,
                     MultiFab& cu,
                     MultiFab& cuMeans,
                     MultiFab& cuVars,
                     MultiFab& prim,
                     MultiFab& primMeans,
                     MultiFab& primVars,
                     MultiFab& spatialCross,
                     MultiFab& eta,
                     MultiFab& kappa);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetMolfrac ( int i, int j, int k,
                  const GpuArray<Real,MAX_SPECIES>& Yk,
                  GpuArray<Real,MAX_SPECIES>& Xk,
                  int const nspecies,
                  const GpuArray<Real,MAX_SPECIES>& molmass)
{
    Real molmix = 0.;

    for (int n=0; n<nspecies; ++n) {
        molmix += Yk[n]/molmass[n];
    }
    molmix = 1./molmix;
    for (int n=0; n<nspecies; ++n) {
        Xk[n] = Yk[n]*molmix/molmass[n];
    }

};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetTemperature ( int i, int j, int k,
                      Real const energy,
                      const GpuArray<Real,MAX_SPECIES>& Yk,
                      Real& temp,
                      int const nspecies,
                      const GpuArray<Real,MAX_SPECIES>& hcv)
{
    Real cvmix = 0.;
    Real e0 = 0.;

    for (int n=0; n<nspecies; ++n) {
        cvmix += Yk[n]*hcv[n];
    }

    temp = (energy-e0)/cvmix;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetPressureGas ( int i, int j, int k,
                      Real& pressure,
                      const GpuArray<Real,MAX_SPECIES>& Yk,
                      Real const density,
                      Real const temp,
                      int const nspecies,
                      Real const Runiv,
                      const GpuArray<Real,MAX_SPECIES>& molmass)
{
    Real molmix = 0.;

    for (int n=0; n<nspecies; ++n) {
        molmix += Yk[n]/molmass[n];
    }
    molmix = 1./molmix;

    pressure = density*(Runiv/molmix)*temp;
};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Decomp ( int const neq,
              GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& a,
              GpuArray<int,MAX_SPECIES> ip,
              int const nspecies )
{
    Real amult, apm, apn;
    int nm1, ip1;
    int i, j, k, l, n;
    Real eps  = 1.e-30;
    Real eps2 = 1.e-60;
      
    nm1=neq-1;
    for (n=0; n<neq; ++n) {
        ip[n]=n+1;
    }
    for (i=0; i<nm1; ++i) {
        ip1=i+2;
        k=i+1;
        apm = std::abs(a[ip[i]*nspecies-1+i]);
        for (j=ip1-1; j<neq; ++j) {
            apn=std::abs(a[(ip[j]-1)*nspecies+i]);
            if (apm >= apn) continue;
            apm=apn;
            k=j+1;
        }
        j=ip[k-1];
        ip[k-1]=ip[i];
        ip[i]=j;
        if (apm < eps) {
            std::printf("%e",a[(ip[neq-1]-1)*nspecies+neq-1]);
            Abort("singular matrix in decomp");
        }
        for (l=ip1-1; l<neq; ++l) {
            n=ip[l];
            amult=a[(n-1)*nspecies+i]/a[(j-1)*nspecies+i];
            a[(n-1)*nspecies+i]=amult;
            if (std::abs(amult) < eps2) continue;
            for (k=ip1-1; k<neq; ++k) {
                a[(n-1)*nspecies+k]=a[(n-1)*nspecies+k]-amult*a[(j-1)*nspecies+k];
            }
        }
    }
    if(std::abs(a[(ip[neq-1]-1)*nspecies+neq-1]) < eps) {
            std::printf("%e",a[(ip[neq-1]-1)*nspecies+neq-1]);
            Abort("singular matrix in decomp");
    }
            
}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Solve ( int const neq,
              GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& a,
              GpuArray<Real,MAX_SPECIES>& b,
              GpuArray<int,MAX_SPECIES> ip,
              int const nspecies )
{
    Real scr;
    int j,k,l,n,nm1,lm1,jp1;

    nm1 = neq-1;
    for (l=1; l<neq; ++l) {
        n=ip[l];
        lm1=l;
        for (k=0; k<lm1; ++k) {
            b[n-1]=b[n-1]-a[(n-1)*nspecies+k]*b[ip[k]-1];
        }
    }
    b[ip[neq-1]-1]=b[ip[neq-1]-1]/a[ip[neq-1]-1,neq-1];
    for (l=0; l<nm1; ++l) {
        j=neq-(l+1);
        jp1=j+1;
        n=ip[j-1];
        for (k=jp1-1; k<neq; ++k) {
            b[n-1]=b[n-1]-a[(n-1)*nspecies+k]*b[ip[k]-1];
        }
        b[n-1]=b[n-1]/a[(n-1)*nspecies+j-1];
    }
    for (n=0; n<neq; ++n) {
        while (ip[n] != n+1) {
            j=ip[n];
            scr=b[j-1];
            ip[n]=ip[j-1];
            b[j-1]=b[ip[j-1]-1];
            b[ip[j-1]-1]=scr;
            ip[j-1]=j;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ViscLin ( const GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& omega11,
               const GpuArray<Real,MAX_SPECIES>& Ykp,
               const Real T,
               const Real rho,
               const GpuArray<Real,MAX_SPECIES>& mk,
               Real& etaMix,
               const int nspecies,
               Real const k_B )
{

    GpuArray<Real,MAX_SPECIES*MAX_SPECIES> QoR;
    GpuArray<Real,MAX_SPECIES> nk, bSonine, diag;
    GpuArray<int,MAX_SPECIES> ip;
    
    Real sum1;
    
    // find number density
    for (int ii=0; ii<nspecies; ++ii) {
        nk[ii] = rho*Ykp[ii]/mk[ii];
    }

    // print*, "Hack (visc_lin): T = ", T
    // stop

    for (int ii=0; ii<nspecies; ++ii) {
        diag[ii] = 0.;
        for (int kk=0; kk<nspecies; ++kk) {
            diag[ii] = diag[ii] + nk[kk]*mk[kk]/pow(mk[ii]+mk[kk],2.) *
               (5.*mk[ii]*omega11[ii*MAX_SPECIES+kk] +
                3.*mk[kk]*omega11[ii*MAX_SPECIES+kk]);
        }
    }

    // print*, "Hack (visc_lin): ", mk, nk, Ykp
    // stop

    for (int ii=0; ii<nspecies; ++ii) {
        for (int jj=0; jj<nspecies; ++jj) {
            sum1 = -nk[jj]*mk[jj]/pow(mk[ii]+mk[jj],2.) *
                2.*mk[jj]*omega11[ii*MAX_SPECIES+jj];

            if (ii==jj) {
                sum1 = sum1 + diag[ii];
            }

            QoR[ii*MAX_SPECIES+jj] = -(16./15.)*(mk[ii]/mk[jj])*sum1;
            // HCB 7.4-62,63
        }
        bSonine[ii] = -1.;
    }

    // print*, "Hack (visc_lin): original QoR = ", QoR

    // QoR = 1.0d24*QoR
    // bSonine = 1.0d24*bSonine

    // print*, "Hack (visc_lin): scaled QoR = ", QoR
        
    Decomp(nspecies,QoR,ip,nspecies);
    Solve(nspecies,QoR,bSonine,ip,nspecies);

    // print*, "Hack (visc_lin): factored QoR = ", QoR
    // stop

    sum1 = 0.;
    for (int ii=0; ii<nspecies; ++ii) {
        sum1 = sum1 + nk[ii]*bSonine[ii];
        // HCB 7.4-56
    }

    etaMix = 0.5*k_B*T*sum1;
    
    // print*, "Hack (visc_lin): etaMix = ", etaMix, T
    // stop
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void IdealMixtureTransport ( int i, int j, int k,
                             Real const density,
                             Real const temperature,
                             Real const pressure,
                             const GpuArray<Real,MAX_SPECIES>& Yk,
                             const GpuArray<Real,MAX_SPECIES>& Xk,
                             Real& eta,
                             Real& kappa,
                             Real& zeta,
                             const Array4<Real>& diff_ij,
                             const Array4<Real>& chitil,
                             int const nspecies,
                             const GpuArray<Real,MAX_SPECIES>& molmass,
                             const GpuArray<Real,MAX_SPECIES>& diameter,
                             Real const k_B,
                             Real const Runiv)
{

    GpuArray<Real,MAX_SPECIES*MAX_SPECIES> Dbin;
    GpuArray<Real,MAX_SPECIES*MAX_SPECIES> omega11;
    GpuArray<Real,MAX_SPECIES*MAX_SPECIES> sigma11;
    GpuArray<Real,MAX_SPECIES*MAX_SPECIES> a_ij1;
    GpuArray<Real,MAX_SPECIES*MAX_SPECIES> a_ij2;
    GpuArray<Real,MAX_SPECIES*MAX_SPECIES> alphabar;
    
    GpuArray<Real,MAX_SPECIES> xxtr;
    GpuArray<Real,MAX_SPECIES> yytr;
    GpuArray<Real,MAX_SPECIES> molecular_mass;

    Real sigma11bar, diamat, mu, Fijstar, Fij, MWmix, sqrtT;
    Real pi = 3.1415926535897932;

    // compute molecular_mass by dividing molmass by Avogadro's
    for (int n=0; n<nspecies; ++n) {
        molecular_mass[n] = molmass[n]*(k_B/Runiv);;
    }
            
    // mole fracions correction - EGLIB
    for (int ii=0; ii<nspecies; ++ii) {
        Real sumXk = 0.;
        for (int jj=0; jj<nspecies; ++jj) {
            sumXk += Xk[jj];
        }
        // GM: Why this factor of 1E-15???
        xxtr[ii] = Xk[ii] + 1.e-15*(sumXk/nspecies - Xk[ii]);
    }
    
    // molecular weight of mixture - EGLIB
    MWmix = 0.;
    for (int ii=0; ii<nspecies; ++ii) {
        MWmix += xxtr[ii]*molmass[ii];
    }

    // mass fractions correction - EGLIB
    for (int ii=0; ii<nspecies; ++ii) {
        yytr[ii] = molmass[ii]/MWmix*xxtr[ii];
    }

    Fijstar = -27./5.;

    // find binary diffusion coefficients  
    // HCB 8.2-9   
    sqrtT = sqrt(temperature);

    for (int i=0; i<nspecies; ++i) {
        for (int j=0; j<nspecies; ++j) {

            // These matrices are computed in init_chemistry in FluctHydro code

            diamat = 0.5*(diameter[i] + diameter[j]);
            mu = molecular_mass[i]*molecular_mass[j]/(molecular_mass[i] + molecular_mass[j]);
            Fij = (6.0*molecular_mass[i]*molecular_mass[i] + 13.0/5.0*molecular_mass[j]*molecular_mass[j] +
                   16.0/5.0*molecular_mass[i]*molecular_mass[j])/(pow(molecular_mass[i]+molecular_mass[j],2.));
                
            sigma11bar = sqrt( k_B/(2.0*pi*mu) )*pi*pow(diamat,2.);
          
            alphabar[i*nspecies+j] = 8.0/(3.0*k_B)*mu*mu*(-.5*sigma11bar);

            Dbin[i*nspecies+j] = 3.0/16.0*
                sqrt(2.0*pi*pow(k_B,3.)
                     *(molecular_mass[i]+molecular_mass[j])/molecular_mass[i]/molecular_mass[j])/(pi*pow(diamat,2.))
                *temperature*sqrtT/pressure;
            omega11[i*nspecies+j] = sqrt(pi*k_B/(2.0*mu))*pow(diamat,2.)*sqrtT;
            sigma11[i*nspecies+j] = sigma11bar*sqrtT;
            a_ij1[i*nspecies+j] = (5.0/(k_B)*molecular_mass[i]*molecular_mass[j]/
                                      (molecular_mass[i]+molecular_mass[j])*Fij*sigma11bar) /sqrtT;
            a_ij2[i*nspecies+j] = (5.0/(k_B)*molecular_mass[i]*molecular_mass[j]*molecular_mass[i]*molecular_mass[j]/
                                      (pow(molecular_mass[i]+molecular_mass[j],3.))*Fijstar*sigma11bar) /sqrtT;

        }
    }

    ViscLin(omega11,yytr,temperature,density,molecular_mass,eta,nspecies,k_B);
    
              
};
                      
#endif
