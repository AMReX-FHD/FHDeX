#ifndef _compressible_functions_H_
#define _compressible_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "common_functions.H"
#include "compressible_namespace.H"

using namespace amrex;
using namespace compressible;

// copy contents of compressible_params_module to C++ compressible namespace
void InitializeCompressibleNamespace();

// if gas heat capacities in the namelist are negative, calculate them using using dofs
void GetHcGas();

void InitConsVar(MultiFab& cons,
                 const amrex::Geometry& geom);

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
                   const Geometry& geom,
                   const Vector< Real >& stoch_weights,
                   const Real dt);

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
             MultiFab& visccorn, MultiFab& rancorn, MultiFab& ranchem,
             const Geometry& geom, const Real dt);

void conservedToPrimitive(MultiFab& prim_in, const MultiFab& cons_in);

void primitiveToConserved(const MultiFab& prim, MultiFab& cons);

/////////////////////////////////////
// in boundary.cpp

void SetupBC();

void SetupCWall();

void setBC(MultiFab& prim, MultiFab& cons);

void BCWallSpeciesFlux(std::array< MultiFab, AMREX_SPACEDIM >& flux,
                       const amrex::Geometry& geom);

void StochFlux(std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
               const amrex::Geometry& geom);

void MembraneFlux(std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
                  const amrex::Geometry& geom);
/////////////////////////////////////

void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                   const MultiFab& prim, MultiFab& primMean, MultiFab& primVar,
                   MultiFab& spatialCross, MultiFab& miscStats, Real* miscVals,
                   const int steps, const Real* dx);

void doMembrane(MultiFab& cons, MultiFab& prim, std::array<MultiFab, AMREX_SPACEDIM>& flux,
                const Geometry& geom, const Real* dx, const Real dt);

void WritePlotFile(int step,
                   const Real time,
                   const Geometry& geom,
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
                     const Geometry& geom,
                     const MultiFab& cu,
                     const MultiFab& cuMeans,
                     const MultiFab& cuVars,
                     const MultiFab& prim,
                     const MultiFab& primMeans,
                     const MultiFab& primVars,
                     const MultiFab& spatialCross,
                     const MultiFab& miscStats,
                     const MultiFab& eta,
                     const MultiFab& kappa);

void ReadCheckPoint(int& step,
                     Real& time,
                     int& statsCount,
                     Geometry& geom,
                     MultiFab& cu,
                     MultiFab& cuMeans,
                     MultiFab& cuVars,
                     MultiFab& prim,
                     MultiFab& primMeans,
                     MultiFab& primVars,
                     MultiFab& spatialCross,
                     MultiFab& miscStats,
                     MultiFab& eta,
                     MultiFab& kappa);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetMolfrac (const GpuArray<Real,MAX_SPECIES>& Yk,
                 GpuArray<Real,MAX_SPECIES>& Xk)
{
    Real molmix = 0.;

    for (int n=0; n<nspecies; ++n) {
        molmix += Yk[n]/molmass[n];
    }
    molmix = 1./molmix;
    for (int n=0; n<nspecies; ++n) {
        Xk[n] = Yk[n]*molmix/molmass[n];
    }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetMassfrac (const GpuArray<Real,MAX_SPECIES>& Xk,
                 GpuArray<Real,MAX_SPECIES>& Yk)
{
    Real molmix = 0.;

    for (int n=0; n<nspecies; ++n) {
        molmix += Xk[n]*molmass[n];
    }
    for (int n=0; n<nspecies; ++n) {
        Yk[n] = Xk[n]*molmass[n]/molmix;
    }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetTemperature ( Real const energy,
                      const GpuArray<Real,MAX_SPECIES>& Yk,
                      Real& temp)
{
    Real cvmix = 0.;
    Real e0mix = 0.;

    for (int n=0; n<nspecies; ++n) {
        cvmix += Yk[n]*hcv[n];
        e0mix += Yk[n]*e0[n];
    }

    temp = (energy-e0mix)/cvmix;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetPressureGas ( Real& pressure,
                      const GpuArray<Real,MAX_SPECIES>& Yk,
                      Real const density,
                      Real const temp)
{
    Real molmix = 0.;

    for (int n=0; n<nspecies; ++n) {
        molmix += Yk[n]/molmass[n];
    }
    molmix = 1./molmix;

    pressure = density*(Runiv/molmix)*temp;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetSoundSpeed (Real& c,
                    const GpuArray<Real,MAX_SPECIES>& Yk,
                    Real const temp)
{
    Real molmix = 0.;
    Real cvmix = 0.;
    Real cpmix = 0.;

    for (int n=0; n<nspecies; ++n) {
        molmix += Yk[n]/molmass[n];
        cvmix += Yk[n]*hcv[n];
        cpmix += Yk[n]*hcp[n];
    }
    molmix = 1./molmix;
    Real gamma = cpmix/cvmix;

    c = sqrt(gamma*(Runiv/molmix)*temp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Decomp ( int const neq,
              GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& a,
              GpuArray<int,MAX_SPECIES>& ip,
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
        apm = std::abs(a[(ip[i]-1)*nspecies+i]);
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
           // std::printf("%e",a[(ip[neq-1]-1)*nspecies+neq-1]);
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
            //std::printf("%e",a[(ip[neq-1]-1)*nspecies+neq-1]);
            Abort("singular matrix in decomp");
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Decomp2 ( int const neq,
               GpuArray<Real,4*MAX_SPECIES*MAX_SPECIES>& a,
               GpuArray<int,2*MAX_SPECIES>& ip,
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
        apm = std::abs(a[(ip[i]-1)*2*nspecies+i]);
        for (j=ip1-1; j<neq; ++j) {
            apn=std::abs(a[(ip[j]-1)*2*nspecies+i]);
            if (apm >= apn) continue;
            apm=apn;
            k=j+1;
        }
        j=ip[k-1];
        ip[k-1]=ip[i];
        ip[i]=j;
        if (apm < eps) {
            //std::printf("%e",a[(ip[neq-1]-1)*2*nspecies+neq-1]);
            Abort("singular matrix in decomp");
        }
        for (l=ip1-1; l<neq; ++l) {
            n=ip[l];
            amult=a[(n-1)*2*nspecies+i]/a[(j-1)*2*nspecies+i];
            a[(n-1)*2*nspecies+i]=amult;
            if (std::abs(amult) < eps2) continue;
            for (k=ip1-1; k<neq; ++k) {
                a[(n-1)*2*nspecies+k]=a[(n-1)*2*nspecies+k]-amult*a[(j-1)*2*nspecies+k];
            }
        }
    }
    if(std::abs(a[(ip[neq-1]-1)*2*nspecies+neq-1]) < eps) {
            //std::printf("%e",a[(ip[neq-1]-1)*2*nspecies+neq-1]);
            Abort("singular matrix in decomp");
    }

}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Solve ( int const neq,
              GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& a,
              GpuArray<Real,MAX_SPECIES>& b,
              GpuArray<int,MAX_SPECIES>& ip,
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
    b[ip[neq-1]-1]=b[ip[neq-1]-1]/a[(ip[neq-1]-1)*nspecies+neq-1];
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
void Solve2 ( int const neq,
              GpuArray<Real,4*MAX_SPECIES*MAX_SPECIES>& a,
              GpuArray<Real,2*MAX_SPECIES>& b,
              GpuArray<int,2*MAX_SPECIES>& ip,
              int const nspecies )
{
    Real scr;
    int j,k,l,n,nm1,lm1,jp1;

    nm1 = neq-1;
    for (l=1; l<neq; ++l) {
        n=ip[l];
        lm1=l;
        for (k=0; k<lm1; ++k) {
            b[n-1]=b[n-1]-a[(n-1)*2*nspecies+k]*b[ip[k]-1];
        }
    }
    b[ip[neq-1]-1]=b[ip[neq-1]-1]/a[(ip[neq-1]-1)*2*nspecies+neq-1];
    for (l=0; l<nm1; ++l) {
        j=neq-(l+1);
        jp1=j+1;
        n=ip[j-1];
        for (k=jp1-1; k<neq; ++k) {
            b[n-1]=b[n-1]-a[(n-1)*2*nspecies+k]*b[ip[k]-1];
        }
        b[n-1]=b[n-1]/a[(n-1)*2*nspecies+j-1];
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
void Decompnp ( int const neq,
                int const ndim,
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& a)
{
    Real eps = 1.e-30;
    Real eps2 = 1.e-60;

    Real amult;

    int i,k,l,nm1,ip1;

    nm1=neq-1;
    for (i=0; i<nm1; ++i) {
        ip1 = i+2;
        for (l=ip1-1; l<neq; ++l) {
            amult=a[l*ndim+i]/a[i*ndim+i];
            a[l*ndim+i]=amult;
            if (std::abs(amult) < eps2) continue;
            for (k=ip1-1; k<neq; ++k) {
                a[l*ndim+k]=a[l*ndim+k]-amult*a[i*ndim+k];
            }
        }
    }
    if (std::abs(a[(neq-1)*ndim+neq-1])<eps) {
        //std::printf("%e",a[(neq-1)*ndim+neq-1]);
        Abort("singular matrix in decomp");
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Solvenp ( int const neq,
               int const ndim,
               GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& a,
               GpuArray<Real,MAX_SPECIES>& b)
{
    int j,k,l,nm1,lm1,jp1;

    nm1 = neq-1;
    for (l=1; l<neq; ++l) {
        lm1=l;
        for (k=0; k<lm1; ++k) {
            b[l]=b[l]-a[l*ndim+k]*b[k];
        }
    }
    b[neq-1]=b[neq-1]/a[(neq-1)*ndim+neq-1];
    for (l=0; l<nm1; ++l) {
        j=neq-(l+1);
        jp1=j+1;
        for (k=jp1-1; k<neq; ++k) {
            b[j-1]=b[j-1]-a[(j-1)*ndim+k]*b[k];
        }
        b[j-1]=b[j-1]/a[(j-1)*ndim+j-1];
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void CholeskyDecomp( GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& a,
                     int const nspecies,
                     GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& sqda)
{
    GpuArray<Real,MAX_SPECIES> p;

    Real sum1 = 0.;
    int ising;

    for (int i=0; i<nspecies; ++i) {
        ising = 0;

        for (int j=0; j<nspecies; ++j) {
            sum1 = a[i*nspecies+j];

            for (int k=i-1; k>=0; --k) {
                sum1 = sum1 - a[i*nspecies+k]*a[j*nspecies+k];
            }

            if (i == j) {
                if (sum1 <= 0.) {
                    p[i] = 0.;
                    ising = 1;
                } else {
                    p[i] = sqrt(sum1);
                }
            }
            if (ising == 0) {
                a[j*nspecies+i] = sum1/p[i];
            } else {
                a[j*nspecies+i] = 0.;
            }
        }
    }

    for (int i=0; i<nspecies; ++i) {
        for (int j=0; j<nspecies; ++j) {
            sqda[i*nspecies+j] = 0.;
        }
    }

    for (int i=0; i<nspecies; ++i) {
        for (int j=i-1; j>=0; --j) {
            sqda[i*nspecies+j] = a[i*nspecies+j];
        }
        sqda[i*nspecies+i] = p[i];
    }

}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void D_GIO ( const GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& Dbin,
             const GpuArray<Real,MAX_SPECIES>& Ykp,
             const GpuArray<Real,MAX_SPECIES>& Xkp,
             const Array4<Real>& D_tilde,
             const int nspecies,
             const int iloc,
             const int jloc,
             const int kloc)
{
    int i,j,k,ll;
    Real term1, term2, Deltamat, scr;

    GpuArray<Real,MAX_SPECIES*MAX_SPECIES> Diff_ij, Zmat, Pmat, Jmat, PJ, matrix1, matrix2;
    GpuArray<Real,MAX_SPECIES> Di, Minv, Mmat;

    int jmax = 3;

    // Find Di matrix
    for (i=0; i<nspecies; ++i) {
        term2 = 0.;
        for (j=0; j<nspecies; ++j) {
            if (j != i) {
                term2 = term2 + Xkp[j]/Dbin[i*nspecies+j];
            }
        }
        Di[i] = (1.-Ykp[i])/term2;
    }

    // Compute Mmat and Minv
    for (i=0; i<nspecies; ++i) {
        Mmat[i] = Xkp[i]/Di[i];
        Minv[i] = Di[i]/Xkp[i];
    }

    // Compute P matrix
    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            Pmat[i*nspecies+j] = - Ykp[j];
            if (i == j) {
                Pmat[i*nspecies+j] =  Pmat[i*nspecies+j] + 1.;
            }
       }
    }

    // Compute Deltamat
    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            if (i == j) {
                term1 = 0.;
                for (k=0; k<nspecies; ++k) {
                    if (k != i) {
                        term1 = term1 + Xkp[i]*Xkp[k]/Dbin[i*nspecies+k];
                    }
                }
                Deltamat = term1;
            } else {
                Deltamat = -Xkp[i]*Xkp[j]/Dbin[i*nspecies+j];
            }
            Zmat[i*nspecies+j] = -Deltamat;
        }
    }

    // Compute Zmat
    for (i=0; i<nspecies; ++i) {
        Zmat[i*nspecies+i] = Zmat[i*nspecies+i] + Mmat[i];
    }

    // Compute Jmat
    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            Jmat[i*nspecies+j] = Minv[i]*Zmat[i*nspecies+j];
        }
    }

    // Compute PJ
    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            PJ[i*nspecies+j] = 0.;
            for (k=0; k<nspecies; ++k) {
                PJ[i*nspecies+j] = PJ[i*nspecies+j] + Pmat[i*nspecies+k]*Jmat[k*nspecies+j];
            }
        }
    }

    // Compute P M^-1 Pt; store it in matrix2
    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            scr = 0.;
            for (k=0; k<nspecies; ++k) {
                scr = scr + Pmat[i*nspecies+k]*Minv[k]*Pmat[j*nspecies+k];
                // notice the change in indices for Pmat to represent Pmat^t
            }
            matrix2[i*nspecies+j] = scr;
            Diff_ij[i*nspecies+j] = scr;
        }
    }

    if (jmax > 0) {

        for (ll=0; ll<jmax; ++ll) {

            // matrix1=0
            for (i=0; i<nspecies; ++i) {
                for (j=0; j<nspecies; ++j) {
                    scr = 0.;
                    for (k=0; k<nspecies; ++k) {
                        scr = scr + PJ[i*nspecies+k]*Diff_ij[k*nspecies+j];
                    }
                    matrix1[i*nspecies+j] = scr+matrix2[i*nspecies+j];
                }
            }

            for (i=0; i<nspecies; ++i) {
                for (j=0; j<nspecies; ++j) {
                    Diff_ij[i*nspecies+j]=matrix1[i*nspecies+j];
                }
            }
        }

    }

    // Compute D_tilde
    for (i=0; i<nspecies; ++i) {
       for (j=0; j<nspecies; ++j) {
           // indicies on LHS are switched due to C++/F90 memory layout
           // this is the only location this needs to happen I believe
           D_tilde(iloc,jloc,kloc,j*nspecies+i) = Diff_ij[i*nspecies+j]*Ykp[i];
       }
    }

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void D_GIO1 ( const Array2D<Real, 0, MAX_SPECIES-1, 0, MAX_SPECIES-1>& Dbin,
             const GpuArray<Real,MAX_SPECIES>& Ykp,
             const GpuArray<Real,MAX_SPECIES>& Xkp,
             const Array4<Real>& D_tilde,
             const int nspecies,
             const int iloc,
             const int jloc,
             const int kloc)
{
    int i,j,k,ll;
    Real term1, term2, Deltamat, scr;

    GpuArray<Real,MAX_SPECIES*MAX_SPECIES> Diff_ij, Zmat, Pmat, Jmat, PJ, matrix1, matrix2;
    GpuArray<Real,MAX_SPECIES> Di, Minv, Mmat;

    int jmax = 3;

    // Find Di matrix
    for (i=0; i<nspecies; ++i) {
        term2 = 0.;
        for (j=0; j<nspecies; ++j) {
            if (j != i) {
                term2 = term2 + Xkp[j]/Dbin(i,j);
            }
        }
        Di[i] = (1.-Ykp[i])/term2;
    }

    // Compute Mmat and Minv
    for (i=0; i<nspecies; ++i) {
        Mmat[i] = Xkp[i]/Di[i];
        Minv[i] = Di[i]/Xkp[i];
    }

    // Compute P matrix
    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            Pmat[i*nspecies+j] = - Ykp[j];
            if (i == j) {
                Pmat[i*nspecies+j] =  Pmat[i*nspecies+j] + 1.;
            }
       }
    }

    // Compute Deltamat
    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            if (i == j) {
                term1 = 0.;
                for (k=0; k<nspecies; ++k) {
                    if (k != i) {
                        term1 = term1 + Xkp[i]*Xkp[k]/Dbin(i,k);
                    }
                }
                Deltamat = term1;
            } else {
                Deltamat = -Xkp[i]*Xkp[j]/Dbin(i,j);
            }
            Zmat[i*nspecies+j] = -Deltamat;
        }
    }

    // Compute Zmat
    for (i=0; i<nspecies; ++i) {
        Zmat[i*nspecies+i] = Zmat[i*nspecies+i] + Mmat[i];
    }

    // Compute Jmat
    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            Jmat[i*nspecies+j] = Minv[i]*Zmat[i*nspecies+j];
        }
    }

    // Compute PJ
    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            PJ[i*nspecies+j] = 0.;
            for (k=0; k<nspecies; ++k) {
                PJ[i*nspecies+j] = PJ[i*nspecies+j] + Pmat[i*nspecies+k]*Jmat[k*nspecies+j];
            }
        }
    }

    // Compute P M^-1 Pt; store it in matrix2
    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            scr = 0.;
            for (k=0; k<nspecies; ++k) {
                scr = scr + Pmat[i*nspecies+k]*Minv[k]*Pmat[j*nspecies+k];
                // notice the change in indices for Pmat to represent Pmat^t
            }
            matrix2[i*nspecies+j] = scr;
            Diff_ij[i*nspecies+j] = scr;
        }
    }

    if (jmax > 0) {

        for (ll=0; ll<jmax; ++ll) {

            // matrix1=0
            for (i=0; i<nspecies; ++i) {
                for (j=0; j<nspecies; ++j) {
                    scr = 0.;
                    for (k=0; k<nspecies; ++k) {
                        scr = scr + PJ[i*nspecies+k]*Diff_ij[k*nspecies+j];
                    }
                    matrix1[i*nspecies+j] = scr+matrix2[i*nspecies+j];
                }
            }

            for (i=0; i<nspecies; ++i) {
                for (j=0; j<nspecies; ++j) {
                    Diff_ij[i*nspecies+j]=matrix1[i*nspecies+j];
                }
            }
        }

    }

    // Compute D_tilde
    for (i=0; i<nspecies; ++i) {
       for (j=0; j<nspecies; ++j) {
           // indicies on LHS are switched due to C++/F90 memory layout
           // this is the only location this needs to happen I believe
           D_tilde(iloc,jloc,kloc,j*nspecies+i) = Diff_ij[i*nspecies+j]*Ykp[i];
       }
    }

}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void LambdaLin ( const GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& Dbin,
                 const GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& omega11,
                 const GpuArray<Real,MAX_SPECIES>& Ykp,
                 const Real T,
                 const Real rho,
                 const GpuArray<Real,MAX_SPECIES>& mk,
                 Real& lammix,
                 const int nspecies,
                 Real const k_B )
{

    GpuArray<Real,MAX_SPECIES> nk, beta, diag1, diag2, diag3, diag4, D_T;
    GpuArray<Real,4*MAX_SPECIES*MAX_SPECIES> QQ;
    GpuArray<Real,2*MAX_SPECIES> aSonine;
    GpuArray<int ,2*MAX_SPECIES> ip;

    int i,j,k;
    Real sum00, sum01, sum11, sum1, lamdaprime, ntotal, ratm, sqratm, scr;


    // find number density
    for (i=0; i<nspecies; ++i) {
        nk[i] = rho*Ykp[i]/mk[i];
        diag1[i] = 0.;
        diag2[i] = 0.;
        diag3[i] = 0.;
        diag4[i] = 0.;
    }

    for (i=0; i<nspecies; ++i) {
        for (k=0; k<nspecies; ++k) {

            diag1[i] =diag1[i]+ nk[k]*mk[k]/(mk[i]+mk[k]) *
                nk[i]*mk[i]*omega11[i*nspecies+k];
            if (k != i) {
                diag2[i] = diag2[i]+nk[k]*mk[k]/(mk[i]+mk[k])*omega11[i*nspecies+k];
            }
            diag3[i] = diag3[i] + nk[i]*nk[k]*pow(mk[k],2.)/pow(mk[i]+mk[k],2.)
                *  .5*omega11[i*nspecies+k];

            diag4[i] = diag4[i]+ nk[i]*nk[k]*mk[k]/pow(mk[i]+mk[k],3.) *
                ((5./4.*(6.*pow(mk[i],2.)+5.*pow(mk[k],2.))*omega11[i*nspecies+k]
                  - 3.*pow(mk[k],2.)*omega11[i*nspecies+k] )
                 +  4.*mk[i]*mk[k]*omega11[i*nspecies+k] );
       }
    }

    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {

            sum00 = -nk[j]*mk[j]*diag2[i] - nk[j]*mk[j]/(mk[i]+mk[j])*nk[i]*mk[i]
                *omega11[i*nspecies+j];
            sum01 = - nk[i]*nk[j]*pow(mk[j],2.)/pow(mk[i]+mk[j],2.) *
                .5*omega11[i*nspecies+j];

            sum11 =  nk[i]*nk[j]*mk[j]/pow(mk[i]+mk[j],3.) *
                ( -(5./4.*(6.*pow(mk[j],2.)+5.*pow(mk[j],2.))*omega11[i*nspecies+j]
                    - 3.*pow(mk[j],2)*omega11[i*nspecies+j] )
                  +  4.*mk[j]*mk[j]*omega11[i*nspecies+j] );

            if (i == j) {
                sum00 = sum00+diag1[i];
                sum01 = sum01+diag3[i];
                sum11 = sum11+diag4[i];
            }

            ratm = mk[i]/mk[j];
            sqratm = sqrt(ratm);
            scr = ratm*sqratm;
            QQ[i*2*nspecies+j] =  8.*sqratm/mk[i] * sum00;         // HCB 7.4-50
            QQ[i*2*nspecies+j+nspecies] = -8.*scr * sum01;    // HCB 7.4-51
            QQ[(i+nspecies)*2*nspecies+j] = -8.*sqratm*sum01;
            QQ[(i+nspecies)*2*nspecies+j+nspecies] = 8.*ratm*sqratm* sum11;  // HCB 7.4-53
       }
    }

    // Build vector rhs ; see HCB 7.4-54
    for (i=0; i<nspecies; ++i) {
        aSonine[i] = 0.;
        beta[i]= sqrt(2.*k_B*T/mk[i]);
        aSonine[i+nspecies] = -(15./4.)*nk[i]*beta[i];
    }

    // print*, "Hack (lambda_lin) predecomp: ", QQ

    // NOTE: the minus sign below; see HCB p 488
    Decomp2(2*nspecies,QQ,ip,nspecies);
    Solve2(2*nspecies,QQ,aSonine,ip,nspecies);

    // print*, "Hack (lambda_lin) postdecomp: ", QQ
    // stop

    for (i=0; i<nspecies; ++i) {
       // HCB 7.4-9
        D_T[i] = 0.5*nk[i]*beta[i]*mk[i]*aSonine[i];
    }

    sum1=0;
    // HCB 7.4-33
    for (i=0; i<nspecies; ++i) {
        sum1 = sum1 + nk[i]*beta[i]* aSonine[nspecies+i];
    }
    lamdaprime = -5./4. * k_B * sum1;


    sum1 = 0;
    // HCB 7.4-65
    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
          sum1 = sum1 + nk[i]*nk[j]/Dbin[i*nspecies+j] *
              pow((D_T[i]/(nk[i]*mk[i]) - D_T[j]/(nk[j]*mk[j])),2.);
       }
    }

    ntotal = 0;
    for (i=0; i<nspecies; ++i) {
        ntotal += nk[i];
    }

    lammix = lamdaprime - 0.5*k_B/ntotal * sum1;  // HCB 7.4-65

    // print*, "Hack (lambda_lin): ", D_T, nk, beta, mk, aSonine
    // print*, "Hack (lambda_lin): ", sum1, nk, Dbin, D_T, mk
    // print*, "Hack (lambda_lin): ", lammix, lamdaprime, ntotal, sum1
    // print*, "Hack (lambda_lin) inputs: ", Dbin,omega11,Ykp,T,rho,mk,lammix
    // stop



}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ThermalDiff ( int iloc, int jloc, int kloc,
                   const GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& /*sigma11*/,
                   const GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& a_ij1,
                   const GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& a_ij2,
                   const GpuArray<Real,MAX_SPECIES*MAX_SPECIES>& alphabar,
                   const GpuArray<Real,MAX_SPECIES>& Xkp,
                   const Real sqrtT,
                   const GpuArray<Real,MAX_SPECIES>& mk,
                   const Array4<Real>& kT,
                   const int nspecies)
{

    GpuArray<Real,MAX_SPECIES*MAX_SPECIES> Aij, alphaij;
    GpuArray<Real,MAX_SPECIES> AA;

    int i,j,k;
    Real sumTemp;

    // Based on Valk 1963 (Waldmann)

    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            Aij[i*nspecies+j] = Xkp[j]*a_ij2[i*nspecies+j];
        }
        AA[i] = 15./4.;
    }

    for (i=0; i<nspecies; ++i) {
        sumTemp = 0.;
        for (k=0; k<nspecies; ++k) {
            sumTemp = sumTemp + Xkp[k]*a_ij1[i*nspecies+k];
        }
        Aij[i*nspecies+i] = Aij[i*nspecies+i] + sumTemp;
    }

    Decompnp(nspecies,nspecies,Aij);
    Solvenp(nspecies,nspecies,Aij,AA);

    for (i=0; i<nspecies; ++i) {
        for (j=0; j<nspecies; ++j) {
            alphaij[i*nspecies+j] = alphabar[i*nspecies+j]*(AA[i]/mk[i] - AA[j]/mk[j])/sqrtT;
        }
    }

    for (k=0; k<nspecies; ++k) {
        sumTemp = 0.;
        for (j=0; j<nspecies; ++j) {
            sumTemp = sumTemp + Xkp[j]*alphaij[k*nspecies+j];
        }
        kT(iloc,jloc,kloc,k) = sumTemp;
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
               (5.*mk[ii]*omega11[ii*nspecies+kk] +
                3.*mk[kk]*omega11[ii*nspecies+kk]);
        }
    }

    // print*, "Hack (visc_lin): ", mk, nk, Ykp
    // stop

    for (int ii=0; ii<nspecies; ++ii) {
        for (int jj=0; jj<nspecies; ++jj) {
            sum1 = -nk[jj]*mk[jj]/pow(mk[ii]+mk[jj],2.) *
                2.*mk[jj]*omega11[ii*nspecies+jj];

            if (ii==jj) {
                sum1 = sum1 + diag[ii];
            }

            QoR[ii*nspecies+jj] = -(16./15.)*(mk[ii]/mk[jj])*sum1;
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
// Transport Coefficients from Valk/Waldmann
void IdealMixtureTransportVW ( int iloc, int jloc, int kloc,
                               Real const density,
                               Real const temperature,
                               Real const pressure,
                               const GpuArray<Real,MAX_SPECIES>& /*Yk*/,
                               const GpuArray<Real,MAX_SPECIES>& Xk,
                               Real& eta,
                               Real& kappa,
                               Real& zeta,
                               const Array4<Real>& diff_ij,
                               const Array4<Real>& chitil)
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
    // eta = 0.0d0

    // // GCM: Why hard-coded in? from original LLNS
    zeta = 0.;

    LambdaLin(Dbin,omega11,yytr,temperature,density,molecular_mass,kappa,nspecies,k_B);
    // kappa = 0.d0

    D_GIO(Dbin,yytr,xxtr,diff_ij,nspecies,iloc,jloc,kloc);
    // diff_ij = 0.0d0

    ThermalDiff(iloc,jloc,kloc,sigma11,a_ij1,a_ij2,alphabar,xxtr,sqrtT,molecular_mass,
                chitil,nspecies);
    // chitil = 0.0d0
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
// Transport Coefficients from Giovangigli
void IdealMixtureTransportGIO (int iloc, int jloc, int kloc,
                               Real const /*density*/,
                               Real const temperature,
                               Real const pressure,
                               const GpuArray<Real,MAX_SPECIES>& Yk,
                               Real& eta,
                               Real& kappa,
                               Real& zeta,
                               const Array4<Real>& diff_ij,
                               const Array4<Real>& chitil)
{
    Array2D<Real, 0, MAX_SPECIES-1, 0, MAX_SPECIES-1> mu;
    Array2D<Real, 0, MAX_SPECIES-1, 0, MAX_SPECIES-1> diam;
    Array2D<Real, 0, MAX_SPECIES-1, 0, MAX_SPECIES-1> Dbin;
    Array2D<Real, 0, MAX_SPECIES-1, 0, MAX_SPECIES-1> amat;

    GpuArray<Real,MAX_SPECIES> molecular_mass;
    GpuArray<Real,MAX_SPECIES> Xk;
    GpuArray<Real,MAX_SPECIES> Xkp;
    GpuArray<Real,MAX_SPECIES> Ykp;
    GpuArray<Real,MAX_SPECIES> eta1;
    GpuArray<Real,MAX_SPECIES> rhs;
    GpuArray<Real,MAX_SPECIES> Bi;

    Real mbar, Xksum, kT, AKL, BKL, CKL, fact1, sum1;
    Real pi = 3.1415926535897932;

    // compute molecular_mass by dividing molmass by Avogadro's
    for (int n=0; n<nspecies; ++n) {
        molecular_mass[n] = molmass[n]*(k_B/Runiv);;
    }

    // mixture molecular weight
    mbar = 0.0;
    for (int n=0; n<nspecies; ++n) {
        mbar += Yk[n]/molecular_mass[n];
    }
    mbar = 1.0/mbar;

    // mole fracions
    Xksum = 0.0;
    for (int n=0; n<nspecies; ++n) {
        Xk[n] = Yk[n]*mbar/molecular_mass[n];
        Xksum += Xk[n];
    }

    // Perturb Xk for zero or negative values -- store in Xkp
    for (int n=0; n<nspecies; ++n) {
        Xkp[n] = Xk[n] + 1.0e-16*(Xksum/double(nspecies) - Xk[n]);
    }

    // Perturbed mass fractions Ykp
    for (int n=0; n<nspecies; ++n) {
        Ykp[n] = Xkp[n]*molecular_mass[n]/mbar;
    }

    for (int i=0; i<nspecies; ++i) {
        for (int j=0; j<nspecies; ++j) {
            diam(i,j) = 0.5*(diameter[i] + diameter[j]);
            mu(i,j)   = molecular_mass[i]*molecular_mass[j]/(molecular_mass[i] + molecular_mass[j]);
        }
    }

    // Binary Diffusion Coefficients
    kT = k_B*temperature;
    for (int i=0; i<nspecies; ++i) {
        for (int j=0; j<nspecies; ++j) {
            Dbin(i,j) = (3.0/16.0)*sqrt(2.0*pi*kT*kT*kT/mu(i,j))/(pressure*pi*diam(i,j)*diam(i,j));
        }
    }

    AKL = 1.0;
    BKL = 1.0;
    CKL = 1.0;

    // Viscosity
    for (int i=0; i<nspecies; ++i) {
        eta1[i] = 5.0/(16.0*diam(i,i)*diam(i,i))*sqrt(molecular_mass[i]*k_B*temperature/pi);
    }

    for (int i=0; i<nspecies; ++i) {
        rhs[i] = Xkp[i];
        Bi[i] = rhs[i];
        for (int j=0; j<nspecies; ++j) {
            if (i==j) {
                amat(i,i) = Xkp[i]*Xkp[i]/eta1[i];
                for (int k=0; k<nspecies; ++k) {
                    if (k!=i) {
                        amat(i,i) += 2.0*kT*Xkp[i]*Xkp[k]/(pressure*Dbin(i,k)*(molecular_mass[i]+molecular_mass[k]))*
                                     (1.0+0.6*molecular_mass[k]/molecular_mass[i]*AKL);
                    }
                }
            }
            else {
                amat(i,j) = 2.0*kT*Xkp[i]*Xkp[j]/(pressure*Dbin(i,j)*(molecular_mass[i]+molecular_mass[j]))*(0.6*AKL-1.0);
            }
        }
    }

    for (int j=0; j<nspecies-1; ++j) {
        for (int i=j+1; i<nspecies; ++i) {
            fact1 = amat(i,j)/amat(j,j);
            for (int k=0; k<nspecies; ++k) {
                amat(i,k) -= fact1*amat(j,k);
            }
            Bi[i] -= fact1*Bi[j];
        }
    }
    Bi[nspecies-1] /= amat(nspecies-1,nspecies-1);

    for (int i=nspecies-2; i>-1; --i) {
        sum1 = 0.0;
        for (int j=nspecies-1; j>i; --j) {
            sum1 += amat(i,j)*Bi[j];
        }
        Bi[i] = (Bi[i] - sum1)/amat(i,i);
    }

    eta = 0.0;
    for (int i=0; i<nspecies; ++i) {
        eta += Xkp[i]*Bi[i];
    }

    // zero bulk viscosity for monoatomic HS
    zeta = 0.;

    // Thermal Conductivity
    for (int i=0; i<nspecies; ++i) {
        rhs[i] = 2.5*Xkp[i];
        Bi[i] = rhs[i];
        for (int j=0; j<nspecies; ++j) {
            if (i==j) {
                amat(i,i) = Xkp[i]*Xkp[i]/Dbin(i,j)*2.0*AKL;
                for (int k=0; k<nspecies; ++k) {
                    if (k!=i) {
                        amat(i,i) += Xkp[i]*Xkp[k]/Dbin(i,k)*mu(i,k)/(molecular_mass[i]+molecular_mass[k])*
                                     (7.5*molecular_mass[i]/molecular_mass[k]+6.25*molecular_mass[k]/molecular_mass[i]
                                      -3.0*molecular_mass[k]/molecular_mass[i]*BKL+4.0*AKL);
                    }
                }
            }
            else {
                amat(i,j) = -1.0*Xkp[i]*Xkp[j]/Dbin(i,j)*mu(i,j)/(molecular_mass[i]+molecular_mass[j])*(55.0/4.0-3.0*BKL-4.0*AKL);
            }
        }
    }

    for (int j=0; j<nspecies-1; ++j) {
        for (int i=j+1; i<nspecies; ++i) {
            fact1 = amat(i,j)/amat(j,j);
            for (int k=0; k<nspecies; ++k) {
                amat(i,k) -= fact1*amat(j,k);
            }
            Bi[i] -= fact1*Bi[j];
        }
    }
    Bi[nspecies-1] /= amat(nspecies-1,nspecies-1);

    for (int i=nspecies-2; i>-1; --i) {
        sum1 = 0.0;
        for (int j=nspecies-1; j>i; --j) {
            sum1 += amat(i,j)*Bi[j];
        }
        Bi[i] = (Bi[i] - sum1)/amat(i,i);
    }

    kappa = 0.0;
    for (int i=0; i<nspecies; ++i) {
        kappa += rhs[i]*Bi[i];
    }
    kappa *= pressure/temperature;

    // Thermal Diffusion Ratios
    for (int i=0; i<nspecies; ++i) {
        chitil(iloc,jloc,kloc,i) = 0.0;
        for (int j=0; j<nspecies; ++j) {
            if (i==j) {
                amat(i,i) = 0.0;
                for (int k=0; k<nspecies; ++k) {
                    if (k!=i) {
                        amat(i,i) -= 0.5*Xkp[k]/Dbin(i,k)*molecular_mass[k]/(molecular_mass[i]+molecular_mass[k])*(6.0*CKL-5.0);
                    }
                }
            }
            else {
                amat(i,j) = 0.5*Xkp[j]/Dbin(i,j)*molecular_mass[i]/(molecular_mass[i]+molecular_mass[j])*(6.0*CKL-5.0);
            }
        }
    }

    for (int i=0; i<nspecies; ++i) {
        for (int j=0; j<nspecies; ++j) {
            chitil(iloc,jloc,kloc,i) += amat(i,j)*Bi[j];
        }
    }

    D_GIO1(Dbin,Ykp,Xkp,diff_ij,nspecies,iloc,jloc,kloc);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
// Transport Coefficients from Hirschfelder-Curtiss-Bird for Binary Mixtures
void IdealMixtureTransportHCBBin ( int iloc, int jloc, int kloc,
                                   Real const density,
                                   Real const temperature,
                                   Real const /*pressure*/,
                                   const GpuArray<Real,MAX_SPECIES>& /*Yk*/,
                                   const GpuArray<Real,MAX_SPECIES>& Xk,
                                   Real& eta,
                                   Real& kappa,
                                   Real& zeta,
                                   const Array4<Real>& diff_ij,
                                   const Array4<Real>& chitil)
{

    // This routine calculates transport coefficients for a
    // binary mixtures of gases of monoatomic hard spheres
    // The components 0 and 1 of Yk and Xk are the only relevant
    // ones.

    GpuArray<Real,2> xxtr;
    GpuArray<Real,2> yytr;
    GpuArray<Real,2> molecular_mass;

    Real pi = 3.1415926535897932;

    // compute molecular_mass by dividing molmass by Avogadro's
    for (int n=0; n<2; ++n) {
        molecular_mass[n] = molmass[n]*(k_B/Runiv);;
    }

    // mole fracions correction - EGLIB
    for (int ii=0; ii<2; ++ii) {
        Real sumXk = 0.;
        for (int jj=0; jj<2; ++jj) {
            sumXk += Xk[jj];
        }
        // GM: Why this factor of 1E-15???
        xxtr[ii] = Xk[ii] + 1.e-15*(sumXk/2.0 - Xk[ii]);
    }

    // molecular weight of mixture - EGLIB
    Real MWmix = 0.;
    for (int ii=0; ii<2; ++ii) {
        MWmix += xxtr[ii]*molmass[ii];
    }

    // mass fractions correction - EGLIB
    for (int ii=0; ii<2; ++ii) {
        yytr[ii] = molmass[ii]/MWmix*xxtr[ii];
    }

    // Mass, diameter, concentration
    Real m0 = molecular_mass[0];
    Real m1 = molecular_mass[1];

    Real d0 = diameter[0];
    Real d1 = diameter[1];

    Real mx = 2.0*m0*m1/(m0 + m1);
    Real dx = 0.5*(d0 + d1);

    Real conc = xxtr[1]; // c'
    Real onc  = xxtr[0]; // 1-c'

    // Viscosity
    Real eta0 = (5.0/(16.0*d0*d0)) * sqrt(m0*k_B*temperature/pi);
    Real eta1 = (5.0/(16.0*d1*d1)) * sqrt(m1*k_B*temperature/pi);
    Real etax = (5.0/(16.0*dx*dx)) * sqrt(mx*k_B*temperature/pi);

    Real Xeta = (onc*onc/eta0) + (2.0*conc*onc/etax) + (conc*conc/eta1);
    Real Yeta = (3.0/5.0) * ( (onc*onc*m0)/(eta0*m1) +
                              (2.0*conc*onc*etax*(m0+m1)*(m0+m1))/(4.0*eta0*eta1*m0*m1) +
                              (conc*conc*m1)/(eta1*m0) );
    Real Zeta = (3.0/5.0) * ( (onc*onc*m0/m1) + 2.0*conc*onc*( ((m0+m1)*(m0+m1)/(4.0*m0*m1))*((etax/eta0) + (etax/eta1)) - 1.0 ) + conc*conc*m1/m0 );

    eta = 1.016*(1.0 + Zeta)/(Xeta + Yeta);

    // Thermal conductivity
    Real kappa0 = 15.0*k_B*eta0/(4.0*m0);
    Real kappa1 = 15.0*k_B*eta1/(4.0*m1);
    Real kappax = 15.0*k_B*etax/(4.0*mx);

    Real U0 = 4.0/15.0 - 17.0*m0/(60.0*m1) + (m0-m1)*(m0-m1)/(2.0*m0*m1);
    Real U1 = 4.0/15.0 - 17.0*m1/(60.0*m0) + (m0-m1)*(m0-m1)/(2.0*m0*m1);
    Real UY = 4.0*(m0+m1)*(m0+m1)*kappax*kappax/(15.0*4.0*m0*m1*kappa0*kappa1) - 17.0/60.0 + 13.0*(m0-m1)*(m0-m1)/(32.0*m0*m1);
    Real UZ = (4.0/15.0)* ( ((m0+m1)*(m0+m1)/(4.0*m0*m1))*((kappax/kappa0) + (kappax/kappa1)) - 1.0) - 17.0/60.0;

    Real Xkappa = onc*onc/kappa0 + 2.0*conc*onc/kappax + conc*conc/kappa1;
    Real Ykappa = onc*onc*U0/kappa0 + 2.0*conc*onc*UY/kappax + conc*conc*U1/kappa1;
    Real Zkappa = onc*onc*U0 + 2.0*conc*onc*UZ + conc*conc*U1;

    kappa = 1.025*(1.0 + Zkappa)/(Xkappa + Ykappa);

    zeta = 0.0;

    // Diffusion Coefficient
    Real diff = 6.0*1.133*etax/(5.0*density);

    // Thermal diffusion ratio (k_T)
    Real S0 = (m0+m1)*kappax/(2.0*m1*kappa0) - 15*(m1-m0)/(8.0*m0) - 1.0;
    Real S1 = (m0+m1)*kappax/(2.0*m0*kappa1) - 15*(m0-m1)/(8.0*m1) - 1.0;

    Real thermal_diff = 1.299*conc*onc*(S1*conc - S0*onc)/(6.0*kappax*(Xkappa + Ykappa));

    // Barodiffusion
    // Real baro_diff    = (m0 - m1)*yytr[1]*yytr[0]*(yytr[0]/m0 + yytr[1]/m1);

    // Fill the Diffusion Matrix
    Real mR = m0/m1;
    Real diffO_tilde = mR*diff*yytr[1]*yytr[1]*yytr[1]/(yytr[0]*xxtr[1]*xxtr[1]);
    Real diff1_tilde = mR*diff*yytr[1]*yytr[1]/(xxtr[1]*xxtr[1]);

    diff_ij(iloc,jloc,kloc,0) =  1.0*yytr[0]*diffO_tilde;
    diff_ij(iloc,jloc,kloc,1) = -1.0*yytr[0]*diffO_tilde;
    diff_ij(iloc,jloc,kloc,2) = -1.0*yytr[0]*diff1_tilde;
    diff_ij(iloc,jloc,kloc,3) =  1.0*yytr[0]*diff1_tilde;

    // Fil the thermal diffusion ratio vector
    thermal_diff = thermal_diff*mR*yytr[1]*yytr[1]/(xxtr[1]*xxtr[1]);
    chitil(iloc,jloc,kloc,0) = -1.0*xxtr[1]*xxtr[1]*thermal_diff/(mR*xxtr[0]*yytr[1]*yytr[1]);
    chitil(iloc,jloc,kloc,1) =  1.0*xxtr[1]*thermal_diff/(mR*yytr[1]*yytr[1]);

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetDensity ( const Real pressure,
                  Real& density,
                  const Real temp,
                  const GpuArray<Real,MAX_SPECIES>& massfrac)
{
    Real molmix = 0.;
    for (int i=0; i<nspecies; ++i) {
        molmix += massfrac[i]/molmass[i];
    }

    molmix = 1./molmix;

    density = pressure/(Runiv/molmix)/temp;

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetEnergy ( Real& energy,
                 const GpuArray<Real,MAX_SPECIES>& massvec,
                 const Real temp)
{
    Real cvmix = 0.;
    Real e0mix = 0.;

    for (int i=0; i<nspecies; ++i) {
        cvmix += massvec[i]*hcv[i];
        e0mix += massvec[i]*e0[i];
    }

    energy = e0mix + temp*cvmix;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetEnthalpies (const Real& T,
                    GpuArray<Real,MAX_SPECIES>& hk)
{
    for (int i=0; i<nspecies; ++i) {
        hk[i] = e0[i]+hcp[i]*T;
    }
}

#endif
