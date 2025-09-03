#include "compressible_functions.H"
#include "common_functions.H"

void doMembrane(MultiFab& cons, MultiFab& prim, std::array<MultiFab, AMREX_SPACEDIM>& flux,
                const amrex::Geometry& geom, const amrex::Real* dx, const amrex::Real dt)
{
    BL_PROFILE_VAR("doMembrane()",doMembrane);

    AMREX_D_TERM(flux[0].setVal(0.0);,
                 flux[1].setVal(0.0);,
                 flux[2].setVal(0.0););

    Real vol = (AMREX_SPACEDIM == 2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2];
    Real area = (AMREX_SPACEDIM == 2) ? dx[1]*cell_depth : dx[1]*dx[2];

    Real mm = molmass[0] / 6.02e23;
    // assumes membrane is on x-faces
    // only works if membrane is on grid boundaries

#if 0
    // SSA
    Real rr = k_B / mm;
    Real hole = transmission*area;

    for ( MFIter mfi(flux[0]); mfi.isValid(); ++mfi) {

        // note bx is nodal in x
        const Box& bx = mfi.validbox();

        const Array4<Real>& xflux = flux[0].array(mfi);
        const Array4<Real>& pr = prim.array(mfi);
        const Array4<Real>& cu = cons.array(mfi);

        if (bx.smallEnd(0) == membrane_cell && bx.bigEnd(0) == membrane_cell) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == membrane_cell) {

                    Real tm = 0.;
                    Real massflux = 0.;
                    Real energyflux = 0.;

                    Real nl = cu(i-1,j,k,0)/mm;
                    Real nr = cu(i,j,k,0)/mm;

                    Real taul = (1./(hole*nl))*sqrt(2.*3.142/(rr*pr(i-1,j,k,4)));
                    Real taur = (1./(hole*nr))*sqrt(2.*3.142/(rr*pr(i,j,k,4)));

                    Real ratel = 1./taul;
                    Real rater = 1./taur;

                    Real rn1 = amrex::Random();

                    Real theta = -std::log(rn1)/(ratel + rater);
                    tm = tm + theta;

                    while (tm < dt) {

                        Real rn1 = amrex::Random();
                        Real rn2 = amrex::Random();
                        Real rn3 = amrex::Random();
                        Real rn4 = amrex::Random();

                        Real energy;

                        if(ratel/(ratel+rater) > rn3) {
                            // left to right
                            energy = -k_B*pr(i-1,j,k,4)*std::log(rn1*rn2);

                            massflux = massflux + mm;
                            energyflux = energyflux + energy;
                        } else {
                            // right to left
                            energy = -k_B*pr(i,j,k,4)*std::log(rn1*rn2);

                            massflux = massflux - mm;
                            energyflux = energyflux - energy;
                        }

                        Real theta = -std::log(rn4)/(ratel + rater);

                        tm = tm + theta;
                    }

                    xflux(i,j,k,0) = massflux/vol;
                    xflux(i,j,k,4) = energyflux/vol;
                }
            });
        }
    }
#endif

    // langevin
    Real fac5 = transmission*(std::pow(k_B,2.5))*6.0/std::sqrt(2*mm*3.142);;
    Real fac3 = transmission*(std::pow(k_B,1.5))*2.0/std::sqrt(2*mm*3.142);
    Real fac1 = transmission*(std::pow(k_B,0.5))/std::sqrt(2*mm*3.142);

    for ( MFIter mfi(flux[0]); mfi.isValid(); ++mfi) {

        // note bx is nodal in x
        const Box& bx = mfi.validbox();

        const Array4<Real>& xflux = flux[0].array(mfi);
        const Array4<Real>& pr = prim.array(mfi);
        const Array4<Real>& cu = cons.array(mfi);

        if (bx.smallEnd(0) == membrane_cell && bx.bigEnd(0) == membrane_cell) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == membrane_cell) {

                    Real rhol = cu(i-1,j,k,0);
                    Real rhor = cu(i,j,k,0);

                    Real tl = pr(i-1,j,k,4);
                    Real tr = pr(i,j,k,4);

                    Real sqrttl = std::sqrt(tl);
                    Real sqrttr = std::sqrt(tr);

                    Real um = fac3*(sqrttl*tl*rhol-sqrttr*tr*rhor);
                    Real nm = fac1*(sqrttl*rhol-sqrttr*rhor);

                    Real uv = fac5*(sqrttl*tl*tl*rhol+sqrttr*tr*tr*rhor);
                    Real nv = fac1*(sqrttl*rhol+sqrttr*rhor);

                    Real cross = fac3*(sqrttl*tl*rhol+sqrttr*tr*rhor);

                    Real corr = cross/(std::sqrt(uv)*std::sqrt(nv));

                    Real rn1 = amrex::RandomNormal(0.,1.);
                    Real rn2 = amrex::RandomNormal(0.,1.);
                    Real rn3 = rn1*corr + std::sqrt(1-pow(corr,2.))*rn2;

                    xflux(i,j,k,0) = (dt*area*nm + std::sqrt(dt*area*mm*nv)*rn1)/vol;
                    xflux(i,j,k,4) = (dt*area*um + std::sqrt(dt*area*mm*uv)*rn3)/(vol*mm);

                }
            });
        }
    }

    flux[0].OverrideSync(geom.periodicity());

    for ( MFIter mfi(cons); mfi.isValid(); ++mfi) {

        // note bx is cell-centered
        const Box& bx = mfi.validbox();

        const Array4<Real>& xflux = flux[0].array(mfi);
        const Array4<Real>& pr = prim.array(mfi);
        const Array4<Real>& cu = cons.array(mfi);

        if (bx.smallEnd(0) == membrane_cell) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == membrane_cell) {

                    cu(i,j,k,0) = cu(i,j,k,0) + xflux(i,j,k,0);
                    cu(i,j,k,4) = cu(i,j,k,4) + xflux(i,j,k,4);

                    if((cu(i,j,k,0) < 0.) || (cu(i,j,k,4) < 0.)) {
                        Print() << "Negative effusion removed";
                        cu(i,j,k,0) = cu(i,j,k,0) - xflux(i,j,k,0);
                        cu(i,j,k,4) = cu(i,j,k,4) - xflux(i,j,k,4);
                    }
                }
            });
        }

        if (bx.bigEnd(0) == membrane_cell-1) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == membrane_cell-1) {

                    cu(i,j,k,0) = cu(i,j,k,0) - xflux(membrane_cell,j,k,0);
                    cu(i,j,k,4) = cu(i,j,k,4) - xflux(membrane_cell,j,k,4);

                    if((cu(i,j,k,0) < 0.) || (cu(i,j,k,4) < 0.)) {
                        Print() << "Negative effusion removed";
                        cu(i,j,k,0) = cu(i,j,k,0) + xflux(membrane_cell,j,k,0);
                        cu(i,j,k,4) = cu(i,j,k,4) + xflux(membrane_cell,j,k,4);
                    }
                }
            });
        }
    }

    conservedToPrimitive(prim, cons);
    cons.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());
}
