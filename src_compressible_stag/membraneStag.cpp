#include "compressible_functions_stag.H"
#include "compressible_functions.H"
#include "common_functions.H"
#include "rng_functions.H"
#include <math.h>

void doMembraneStag(MultiFab& cons,
                    std::array< MultiFab, AMREX_SPACEDIM >& cumom,
                    MultiFab& prim,
                    std::array< MultiFab, AMREX_SPACEDIM >& vel,
                    std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                    const amrex::Geometry& geom, const amrex::Real dt)
{
    BL_PROFILE_VAR("doMembraneStag()",doMembraneStag);

#ifdef AMREX_USE_GPU
    // doLangevin/applyEffusion loop on the host over Array4 views of device memory and draw
    // host random numbers; they have not been ported to GPU.
    amrex::Abort("doMembraneStag: membrane routines are not ported to GPU");
#endif

    faceflux[0].setVal(0.0); // set mass flux to zero

    doLangevin(cons,prim,faceflux,geom,dt);

    faceflux[0].OverrideSync(geom.periodicity());

    applyEffusion(faceflux,cons);

    cons.FillBoundary(geom.periodicity()); // need correct periodicity for density to calculate velocities below
    conservedToPrimitiveStag(prim, vel, cons, cumom);
    cons.FillBoundary(geom.periodicity()); // cell-centered momentum is filled in the line above
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        vel[d].FillBoundary(geom.periodicity());
    }
    prim.FillBoundary(geom.periodicity());
}

void doLangevin(MultiFab& cons_in, MultiFab& prim_in,
                std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                const amrex::Geometry& geom,
                const amrex::Real dt)
{

    BL_PROFILE_VAR("doLangevin()",doLangevin);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // The effusion coefficients involve k_B^2.5 (~1e-40) and products like dt*area*mass (~1e-45),
    // which underflow single precision. This loop runs on the host, so evaluate it entirely in double
    // and cast to Real only when storing the flux.
    double vol = double(dx[0])*double(dx[1])*double(dx[2]);
    double area = double(dx[1])*double(dx[2]);
    double dtd = double(dt);

    GpuArray<double,MAX_SPECIES> mass;
    GpuArray<double,MAX_SPECIES> fac1;
    GpuArray<double,MAX_SPECIES> fac3;
    GpuArray<double,MAX_SPECIES> fac5;
    for (int l=0;l<nspecies;++l) {
        mass[l] = double(molmass[l])/6.02e23;
        fac5[l] = double(transmission[l])*std::pow(double(k_B),2.5)*6.0/std::sqrt(2*mass[l]*3.142);
        fac3[l] = double(transmission[l])*std::pow(double(k_B),1.5)*2.0/std::sqrt(2*mass[l]*3.142);
        fac1[l] = double(transmission[l])*std::sqrt(double(k_B))*1.0/std::sqrt(2*mass[l]*3.142);
    }

    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<Real>& prim = prim_in.array(mfi);
        const Array4<Real>& cons = cons_in.array(mfi);
        const Array4<Real>& xflux = faceflux[0].array(mfi);

        if ((bx.smallEnd(0) == membrane_cell) or (bx.bigEnd(0) == membrane_cell - 1)) {

            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {

                double TL = prim(membrane_cell-1,j,k,4);
                double TR = prim(membrane_cell,j,k,4);
                double sqrtTL = std::sqrt(TL);
                double sqrtTR = std::sqrt(TR);

                GpuArray<double,MAX_SPECIES> rhoL;
                GpuArray<double,MAX_SPECIES> rhoR;
                GpuArray<double,MAX_SPECIES> delUmean;
                GpuArray<double,MAX_SPECIES> delNmean;
                GpuArray<double,MAX_SPECIES> delUvar;
                GpuArray<double,MAX_SPECIES> delNvar;
                GpuArray<double,MAX_SPECIES> cross;
                GpuArray<double,MAX_SPECIES> corr;

                for (int l=0;l<nspecies;++l) {

                    rhoL[l] = cons(membrane_cell-1,j,k,5+l);
                    rhoR[l] = cons(membrane_cell,j,k,5+l);

                    delUmean[l] = fac3[l]*(sqrtTL*TL*rhoL[l] - sqrtTR*TR*rhoR[l]);
                    delNmean[l] = fac1[l]*(sqrtTL*rhoL[l] - sqrtTR*rhoR[l]);

                    delUvar[l] = fac5[l]*(sqrtTL*TL*TL*rhoL[l] + sqrtTR*TR*TR*rhoR[l]);
                    delNvar[l] = fac1[l]*(sqrtTL*rhoL[l] + sqrtTR*rhoR[l]);

                    cross[l] = fac3[l]*(sqrtTL*TL*rhoL[l] + sqrtTR*TR*rhoR[l]);
                    corr[l] = cross[l]/(std::sqrt(delUvar[l])*std::sqrt(delNvar[l]));

                    double rn1, rn2, rn3;
                    if (stoch_stress_form == 1) {
                        //rn1 = get_fhd_normal_func();
                        //rn2 = get_fhd_normal_func();
                        rn1 = amrex::RandomNormal(0.,1.);
                        rn2 = amrex::RandomNormal(0.,1.);
                        // |corr| <= 1 analytically; clamp so rounding cannot make the sqrt argument negative
                        rn3 = rn1*corr[l] + std::sqrt(amrex::max(0.0, 1-(corr[l]*corr[l])))*rn2;
                    }
                    else {
                        rn1 = 0.;
                        rn2 = 0.;
                        rn3 = 0.;
                    }

                    xflux(membrane_cell,j,k,5+l) = Real((dtd*area*delNmean[l] + std::sqrt(dtd*area*mass[l]*delNvar[l])*rn1)/vol);
                    xflux(membrane_cell,j,k,0) +=  xflux(membrane_cell,j,k,5+l);
                    xflux(membrane_cell,j,k,4) +=  Real((dtd*area*delUmean[l] + std::sqrt(dtd*area*mass[l]*delUvar[l])*rn3)/(vol*mass[l]));

                }

            }
            }

        }

    }

}

void applyEffusion(std::array<MultiFab, AMREX_SPACEDIM>& faceflux, MultiFab& cons_in)
{

    BL_PROFILE_VAR("applyEffusion()",applyEffusion);

    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<Real>& cons = cons_in.array(mfi);
        const Array4<Real>& xflux = faceflux[0].array(mfi);

        if (bx.smallEnd(0) == membrane_cell) { // membrane at the left end (cell to the right of the membrane)

            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {

                cons(membrane_cell,j,k,0) += xflux(membrane_cell,j,k,0);
                cons(membrane_cell,j,k,4) += xflux(membrane_cell,j,k,4);

                for (int l=0;l<nspecies;++l) {
                    cons(membrane_cell,j,k,5+l) += xflux(membrane_cell,j,k,5+l);
                }
            }
            }
        }

        else if (bx.bigEnd(0) == membrane_cell - 1) { // membrane at the right end (cell to the left of the membrane)

            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {

                cons(membrane_cell-1,j,k,0) -= xflux(membrane_cell,j,k,0);
                cons(membrane_cell-1,j,k,4) -= xflux(membrane_cell,j,k,4);

                for (int l=0;l<nspecies;++l) {
                    cons(membrane_cell-1,j,k,5+l) -= xflux(membrane_cell,j,k,5+l);
                }
            }
            }

        }

    }

}
