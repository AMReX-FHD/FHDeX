#include "multispec_test_functions.H"

#include "multispec_namespace.H"

using namespace amrex;
using namespace multispec;

void InitRhoUmac(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                 MultiFab& rho_in,
                 const Geometry& geom)
{

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    GpuArray<Real,AMREX_SPACEDIM> center;
    GpuArray<Real,AMREX_SPACEDIM> L; // domain length

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        center[d] = 0.5*(prob_hi[d]+prob_lo[d]);
        L[d] = prob_hi[d] - prob_lo[d];
    }

    BoxArray ba = rho_in.boxArray();
    DistributionMapping dmap = rho_in.DistributionMap();

    MultiFab conc(ba,dmap,nspecies,0);

    // set velocity to zero; overwrite below if needed
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac[d].setVal(0.);
    }

    for (MFIter mfi(rho_in); mfi.isValid(); ++mfi ) {

        Box bx = mfi.tilebox();

        const Array4<Real>& c = conc.array(mfi);

        if (prob_type == 1) {

            /*
              bubble with radius = 1/4 of domain in x
              c=c_init_1(:) inside, c=c_init_2(:) outside
              can be discontinous or smooth depending on smoothing_width
            */
            Real rad = L[0] / 4.;
            rad = radius_cyl;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x,y,z;
                AMREX_D_TERM(x = prob_lo[0] + (i+0.5)*dx[0] - center[0];,
                             y = prob_lo[1] + (j+0.5)*dx[1] - center[1];,
                             z = prob_lo[2] + (k+0.5)*dx[2] - center[2];);

                Real r = (AMREX_SPACEDIM == 1) ? std::sqrt(x*x+y*y) : std::sqrt(x*x+y*y+z*z);

                if (smoothing_width == 0.) {

                    // discontinuous interface
                    if (r < rad) {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = c_init_1[n];
                        }
                    } else {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = c_init_2[n];
                        }
                    }

                } else {
                    // smooth interface
                    for (int n=0; n<nspecies; ++n) {
                        c(i,j,k,n) = c_init_1[n] + (c_init_2[n]-c_init_1[n]) *
                            0.5*(1. + std::tanh((r-rad)/(smoothing_width*dx[0])));
                    }
                }

            });

        } else if (prob_type == 5) {

            /*
              bubble with radius = 1/4 of domain in x
              c=c_init_1(:) inside, c=c_init_2(:) outside
              can be discontinous or smooth depending on smoothing_width
            */

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x,y,z;
                AMREX_D_TERM(x = prob_lo[0] + (i+0.5)*dx[0] - center[0];,
                             y = prob_lo[1] + (j+0.5)*dx[1] - center[1];,
                             z = prob_lo[2] + (k+0.5)*dx[2] - center[2];);


                if (smoothing_width == 0.) {
                    for (int n=0; n<nspecies; ++n) {
                        c(i,j,k,n) = 0.01;
                    }

                    // discontinuous interface
                    if (x< 0. && y < 0.) {
                        c(i,j,k,0) = 0.97;
                    } else if ( x > 0. && y < 0.){
                        c(i,j,k,1) = 0.97;
                    } else if ( x < 0. && y > 0.){
                        c(i,j,k,2) = 0.97;
                    } else {
                        c(i,j,k,3) = 0.97;
                    }


//
//                } else {
//                    // smooth interface
//                    for (int n=0; n<nspecies; ++n) {
//                        c(i,j,k,n) = c_init_1[n] + (c_init_2[n]-c_init_1[n]) *
//                            0.5*(1. + std::tanh((r-rad)/(smoothing_width*dx[0])));
//                    }
                }

            });
        } else if (prob_type == 6) {

            /*
              cylinder with radius = 1/4 of domain in x
              c=c_init_1(:) inside, c=c_init_2(:) outside
              can be discontinous or smooth depending on smoothing_width
            */
            //Real rad = L[0] / 8.;
            Real rad = L[1]/8.;
            amrex::Print() << "smoothing width " << smoothing_width << " radius " << rad << std::endl;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x,y,z;
                AMREX_D_TERM(x = prob_lo[0] + (i+0.5)*dx[0] - center[0];,
                             y = prob_lo[1] + (j+0.5)*dx[1] - center[1];,
                             z = prob_lo[2] + (k+0.5)*dx[2] - center[2];);

                Real r = (AMREX_SPACEDIM == 2) ? std::sqrt(x*x+y*y) : std::sqrt(x*x+y*y+z*z);

                if (smoothing_width == 0.) {

                    // discontinuous interface
                    if (r < rad) {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = c_init_1[n];
                        }
                    } else {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = c_init_2[n];
                        }
                    }

                } else {
                    // smooth interface
                    //    c(i,j,k,1) = c_init_1[0] + (c_init_2[0]-c_init_1[0]) *
                    //        0.5*(1. - std::tanh((r-rad)/(smoothing_width*dx[0])));
                    c(i,j,k,0) = c_init_1[0] + (c_init_2[0]-c_init_1[0]) *
                        (1.- std::tanh(y/(smoothing_width*dx[0])))
                      * 0.25*(1. + std::tanh((r-rad)/(smoothing_width*dx[0])));
                    c(i,j,k,2) = c_init_1[0] + (c_init_2[0]-c_init_1[0]) *
                        (1.+ std::tanh(y/(smoothing_width*dx[0])))
                      * 0.25*(1. + std::tanh((r-rad)/(smoothing_width*dx[0])));
                    //    c(i,j,k,2) = 1.-c(i,j,k,0)-c(i,j,k,1);
                    c(i,j,k,1) = 1.-c(i,j,k,0)-c(i,j,k,2);

                }

            });
        } else if (prob_type == 7) {

            /*
              cylinder with radius = 1/4 of domain in x
              c=c_init_1(:) inside, c=c_init_2(:) outside
              can be discontinous or smooth depending on smoothing_width
            */
            //Real rad = L[0] / 8.;
            Real rad = radius_cyl;
            amrex::Print() << "smoothing width " << smoothing_width << " radius " << rad << std::endl;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x,y,z;
                AMREX_D_TERM(x = prob_lo[0] + (i+0.5)*dx[0] - center[0];,
                             y = prob_lo[1] + (j+0.5)*dx[1] - center[1];,
                             z = prob_lo[2] + (k+0.5)*dx[2] - center[2];);

                Real r = (AMREX_SPACEDIM == 2) ? std::sqrt(y*y) : std::sqrt(y*y+z*z);

                if (smoothing_width == 0.) {

                    // discontinuous interface
                    if (r < rad) {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = c_init_1[n];
                        }
                    } else {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = c_init_2[n];
                        }
                    }

                } else {
                    // smooth interface
                    for (int n=0; n<nspecies; ++n) {
                        c(i,j,k,n) = c_init_1[n] + (c_init_2[n]-c_init_1[n]) *
                            0.5*(1. + std::tanh((r-rad)/(smoothing_width*dx[0])));
                    }
                }

            });
        } else if (prob_type == 8) {

            /*
              cylinder with radius = 1/4 of domain in x
              c=c_init_1(:) inside, c=c_init_2(:) outside
              can be discontinous or smooth depending on smoothing_width
            */
            //Real rad = L[0] / 8.;
            Real rad = radius_cyl;
            int nsub = 10;
            Real factor = nsub;
            Real dxsub = dx[0]/factor;
            Real dysub = dx[1]/factor;
            Real x,y,z;
            amrex::Print() << "smoothing width " << smoothing_width << " radius " << rad << std::endl;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               for (int n=0; n<nspecies; ++n) {
                   c(i,j,k,n) = 0.;
               }
               Real x,y,z;

               for(int i1=0; i1<nsub; ++i1) {
               for(int j1=0; j1<nsub; ++j1) {

                   AMREX_D_TERM(x = prob_lo[0] + i*dx[0] + (i1+0.5)*dxsub - center[0];,
                                y = prob_lo[1] + j*dx[1] + (j1+0.5)*dysub - center[1];,
                                z = prob_lo[2] + (k+0.5)*dx[2] - center[2];);

                   Real r = (AMREX_SPACEDIM == 2) ? std::sqrt(x*x+y*y) : std::sqrt(x*x+y*y);

                   if (smoothing_width == 0.) {

                       // discontinuous interface
                       if (r < rad) {
                           for (int n=0; n<nspecies; ++n) {
                               c(i,j,k,n) += c_init_1[n];
                           }
                       } else {
                           for (int n=0; n<nspecies; ++n) {
                               c(i,j,k,n) += c_init_2[n];
                           }
                       }

                   } else {
                       // smooth interface
                       for (int n=0; n<nspecies; ++n) {
                           c(i,j,k,n) += c_init_1[n] + (c_init_2[n]-c_init_1[n]) *
                               0.5*(1. + std::tanh((r-rad)/(smoothing_width*dx[0])));
                       }
                   }
                }
                }
               for (int n=0; n<nspecies; ++n) {
                   c(i,j,k,n) = c(i,j,k,n)/(factor*factor);
               }
            });

#if (AMREX_SPACEDIM == 3)


            const Array4<Real> & wmac = (umac[2]).array(mfi);
            Box bx_wmac = mfi.tilebox(nodal_flag_z);

           //  Real veljet = 4082.e0;
            Real veljet = 0.e0;

            amrex::ParallelFor(bx_wmac, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x,y,z;

                for(int i1=0; i1<nsub; ++i1) {
                for(int j1=0; j1<nsub; ++j1) {

                AMREX_D_TERM(x = prob_lo[0] + i*dx[0] + (i1+0.5)*dxsub - center[0];,
                             y = prob_lo[1] + j*dx[1] + (j1+0.5)*dysub - center[1];,
                             z = prob_lo[2] + (k+0.5)*dx[2] - center[2];);

                Real r = (AMREX_SPACEDIM == 2) ? std::sqrt(x*x+y*y) : std::sqrt(x*x+y*y);

                if (smoothing_width == 0.) {

                    // discontinuous interface
                    if (r < rad) {
                            wmac(i,j,k) += veljet;
                    } else {
                            wmac(i,j,k) -= veljet;
                    }

                } else {
                    // smooth interface
                        wmac(i,j,k) += veljet -2.e0*veljet *
                            0.5*(1. + std::tanh((r-rad)/(smoothing_width*dx[0])));
                }
             }
             }
                   wmac(i,j,k) = wmac(i,j,k)/(factor*factor);

            });

#endif

        } else if (prob_type == 9) {

            /*
               torus
            */
            //Real rad = L[0] / 8.;
            Real router = 1.5*5.73e-6;
            Real rad = radius_cyl;
            int nsub = 10;
            Real factor = nsub;
            Real dxsub = dx[0]/factor;
            Real dysub = dx[1]/factor;
            Real dzsub = dx[2]/factor;
            Real x,y,z;
            amrex::Print() << "smoothing width " << smoothing_width << " radius " << rad << std::endl;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               for (int n=0; n<nspecies; ++n) {
                   c(i,j,k,n) = 0.;
               }
               Real x,y,z;

               for(int i1=0; i1<nsub; ++i1) {
               for(int j1=0; j1<nsub; ++j1) {
               for(int k1=0; k1<nsub; ++k1) {

                   AMREX_D_TERM(x = prob_lo[0] + i*dx[0] + (i1+0.5)*dxsub - center[0];,
                                y = prob_lo[1] + j*dx[1] + (j1+0.5)*dysub - center[1];,
                                z = prob_lo[2] + k*dx[2] + (k1+0.5)*dzsub - center[2];);

                   Real r_ring = (AMREX_SPACEDIM == 2) ? std::sqrt(x*x+y*y) : std::sqrt(x*x+y*y);
                   Real r = std::sqrt((r_ring-router)*(r_ring-router)+ z*z);

                   if (smoothing_width == 0.) {

                       // discontinuous interface
                       if (r < rad) {
                           for (int n=0; n<nspecies; ++n) {
                               c(i,j,k,n) += c_init_1[n];
                           }
                       } else {
                           for (int n=0; n<nspecies; ++n) {
                               c(i,j,k,n) += c_init_2[n];
                           }
                       }

                   } else {
                       // smooth interface
                       for (int n=0; n<nspecies; ++n) {
                           c(i,j,k,n) += c_init_1[n] + (c_init_2[n]-c_init_1[n]) *
                               0.5*(1. + std::tanh((r-rad)/(smoothing_width*dx[0])));
                       }
                   }
                }
                }
                }
               for (int n=0; n<nspecies; ++n) {
                   c(i,j,k,n) = c(i,j,k,n)/(factor*factor*factor);
               }
            });

        } else if (prob_type == 16) {

            /*
               thin film
            */
            //Real rad = L[0] / 8.;
            int nsub = 10;
            Real factor = nsub;
            Real dxsub = dx[0]/factor;
            Real dysub = dx[1]/factor;
            Real dzsub = dx[2]/factor;
            Real x,y,z;
            amrex::Print() << "smoothing width " << smoothing_width << " film_thickness " << film_thickness << std::endl;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               for (int n=0; n<nspecies; ++n) {
                   c(i,j,k,n) = 0.;
               }
               Real x,y,z;

               for(int j1=0; j1<nsub; ++j1) {

                               y = prob_lo[1] + j*dx[1] + (j1+0.5)*dysub ;


                   if (smoothing_width == 0.) {

                       // discontinuous interface
                       if (y < film_thickness) {
                           for (int n=0; n<nspecies; ++n) {
                               c(i,j,k,n) += c_init_1[n];
                           }
                       } else {
                           for (int n=0; n<nspecies; ++n) {
                               c(i,j,k,n) += c_init_2[n];
                           }
                       }

                   } else {
                       // smooth interface
                       for (int n=0; n<nspecies; ++n) {
                           c(i,j,k,n) += c_init_1[n] + (c_init_2[n]-c_init_1[n]) *
                               0.5*(1. + std::tanh((y-film_thickness)/(smoothing_width*dx[0])));
                       }
                   }
                }
               for (int n=0; n<nspecies; ++n) {
                   c(i,j,k,n) = c(i,j,k,n)/factor;
               }
            });


        } else if (prob_type == 20) {

            /*
              cylinder with radius = 1/4 of domain in x
              c=c_init_1(:) inside, c=c_init_2(:) outside
              can be discontinous or smooth depending on smoothing_width
            */
            //Real rad = L[0] / 8.;
            int nsub = 100;
            Real factor = nsub;
            Real dxsub = dx[0]/factor;
            Real dysub = dx[1]/factor;
            Real x,y,z;
            amrex::Print() << "smoothing width " << smoothing_width << " film_thickness " << film_thickness << std::endl;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               for (int n=0; n<nspecies; ++n) {
                   c(i,j,k,n) = 0.;
               }
               Real x,y,z;

               for(int i1=0; i1<nsub; ++i1) {
               for(int j1=0; j1<nsub; ++j1) {

                   AMREX_D_TERM(x = prob_lo[0] + i*dx[0] + (i1+0.5)*dxsub ;,
                                y = prob_lo[1] + j*dx[1] + (j1+0.5)*dysub ;,
                                z = prob_lo[2] + (k+0.5)*dx[2] - center[2];);

                   Real height = film_thickness + 0.4*film_thickness*amrex::Math::cospi(2.*x/L[0]);
                   //height = film_thickness ;

   //               amrex::Print{} << "x,thickness,height " << x << " " << film_thickness << " " << height << std::endl;

                   if (smoothing_width == 0.) {

                       // discontinuous interface
                       if (y < height) {
                           for (int n=0; n<nspecies; ++n) {
                               c(i,j,k,n) += c_init_1[n];
                           }
                       } else {
                           for (int n=0; n<nspecies; ++n) {
                               c(i,j,k,n) += c_init_2[n];
                           }
                       }

                   } else {
                       // smooth interface
                       for (int n=0; n<nspecies; ++n) {
                           c(i,j,k,n) += c_init_1[n] + (c_init_2[n]-c_init_1[n]) *
                               0.5*(1. + std::tanh((y-height)/4.64475e-8));
                               //0.5*(1. + std::tanh((y-height)/9.2895e-8));
                               //0.5*(1. + std::tanh((y-height)/(smoothing_width*dx[0])));
                       }
                   }
                }
                }
               for (int n=0; n<nspecies; ++n) {
                   c(i,j,k,n) = c(i,j,k,n)/(factor*factor);
               }
            });



        } else if (prob_type == 3) {

            Real rad1 = 0.0119;
            Real rad2 = 0.0022853;
            rad1 = 0.01;
            Real shift1 = 1.05*rad1;
            Real shift2 = 2.1*rad1 + rad2;
            Real velbub = 100.;


            Real bub1[3];
            Real bub2[3];
            Real back[3];

//            bub1[0]=0.98;
//            bub1[1]=0.01;
//            bub1[2]=0.01;
//            bub2[0]=0.01;
//            bub2[1]=0.98;
//            bub2[2]=0.01;
//            back[0]=0.01;
//            back[1]=0.01;
//            back[2]=0.98;

            bub1[0]=1.0;
            bub1[1]=0.0;
            bub1[2]=0.0;
            bub2[0]=0.0;
            bub2[1]=1.0;
            bub2[2]=0.0;
            back[0]=0.0;
            back[1]=0.0;
            back[2]=1.0;

/*
            bub1[0]=.8;
            bub1[1]=.1;
            bub1[2]=.1;
            bub2[0]=.1;
            bub2[1]=.8;
            bub2[2]=.1;
            back[0]=.1;
            back[1]=.1;
            back[2]=.8;
*/

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x,y,z;
                Real x2,y2,z2;
                AMREX_D_TERM(x = prob_lo[0] + (i+0.5)*dx[0] - center[0];,
                             y = prob_lo[1] + (j+0.5)*dx[1] - center[1];,
                             z = prob_lo[2] + (k+0.5)*dx[2] -shift1;);
                AMREX_D_TERM(x2 = prob_lo[0] + (i+0.5)*dx[0] - center[0];,
                             y2 = prob_lo[1] + (j+0.5)*dx[1] - center[1];,
                             z2 = prob_lo[2] + (k+0.5)*dx[2] -shift2;);

                Real r1 = (AMREX_SPACEDIM == 2) ? std::sqrt(x*x+y*y) : std::sqrt(x*x+y*y+z*z);
                Real r2 = (AMREX_SPACEDIM == 2) ? std::sqrt(x2*x2+y2*y2) : std::sqrt(x2*x2+y2*y2+z2*z2);

                if (smoothing_width == 0.) {

                    // discontinuous interface
                    if (r1 < rad1) {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = bub1[n];
                        }
                    } else if (r2 < rad2) {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = bub2[n];
                        }
                    } else {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = back[n];
                        }
                    }

                } else {
                    // smooth interface
                    // not coded for this
                    c(i,j,k,0) = bub1[0] + (1.-bub2[0]-2.*bub1[0]) *
                        0.5*(1. + std::tanh((r1-rad1)/(smoothing_width*dx[0])));
                    c(i,j,k,1) = bub2[1] + (1.-bub1[1]-2.*bub2[1]) *
                        0.5*(1. + std::tanh((r2-rad2)/(smoothing_width*dx[0])));
                    c(i,j,k,2) = 1.-c(i,j,k,0)-c(i,j,k,1);
                   // for (int n=0; n<nspecies; ++n) {
                   //     c(i,j,k,n) = c_init_1[n] + (c_init_2[n]-c_init_1[n]) *
                   //        0.5*(1. + std::tanh((r1-rad1)/smoothing_width*dx[0]));
                   // }
                }

            });

            const Array4<Real> & wmac = (umac[2]).array(mfi);
            Box bx_wmac = mfi.tilebox(nodal_flag_z);

            amrex::ParallelFor(bx_wmac, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x2,y2,z2;
                AMREX_D_TERM(x2 = prob_lo[0] + (i+0.5)*dx[0] - center[0];,
                             y2 = prob_lo[1] + (j+0.5)*dx[1] - center[1];,
                             z2 = prob_lo[2] + (k)*dx[2] - shift2;);

                Real r2 = (AMREX_SPACEDIM == 2) ? std::sqrt(x2*x2+y2*y2) : std::sqrt(x2*x2+y2*y2+z2*z2);
                if (r2 < rad2) {
                    for (int n=0; n<nspecies; ++n) {
                        wmac(i,j,k) = -velbub;
                    }
                }
            });

        } else if (prob_type == 2) {

            /*
              constant concentration gradient along y
              c=c_init_1(:) on bottom, c=c_init_2(:) on top
            */
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x,y,z;
                AMREX_D_TERM(x = prob_lo[0] + (i+0.5)*dx[0];,
                             y = prob_lo[1] + (j+0.5)*dx[1];,
                             z = prob_lo[2] + (k+0.5)*dx[2];);

                for (int n=0; n<nspecies; ++n) {
                    c(i,j,k,n) = c_init_1[n] + (c_init_2[n]-c_init_1[n]) *
                        (y-prob_lo[1]) / L[1];
                }


            });
        } else if (prob_type == 4) {

            /*
              bubble with radius = 1/4 of domain in x
              c=c_init_1(:) inside, c=c_init_2(:) outside
              can be discontinous or smooth depending on smoothing_width
            */
            Real l1 = L[1] / 4.;
            Real l2 = 3.*l1;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x,y,z;
                AMREX_D_TERM(x = prob_lo[0] + (i+0.5)*dx[0];,
                             y = prob_lo[1] + (j+0.5)*dx[1];,
                             z = prob_lo[2] + (k+0.5)*dx[2];);

                if (smoothing_width == 0) {

                    // discontinuous interface
                    if (y < l1) {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = c_init_1[n];
                        }
                    } else if (y < l2) {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = c_init_2[n];
                        }
                    } else {
                        for (int n=0; n<nspecies; ++n) {
                            c(i,j,k,n) = c_init_1[n];
                        }
                    }
                } else {

                    // smooth interface
                    for (int n=0; n<nspecies; ++n) {
                        c(i,j,k,n) = c_init_1[n] +
                            (c_init_2[n] - c_init_1[n])*
                            (1./(1.+std::exp(-smoothing_width*(y-l1))) - 1./(1.+std::exp(-smoothing_width*(y-l2))));
                    }
                }
            });
        } else if (prob_type == 12) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                /*
                  Gaussian bubble centered in domain
                  c=c_init_1(:) inside; c=c_init_2(:) outside
                  lo- and hi-y walls move with prescribed velocity,
                  see inhomogeneous_bc_val.f90
                */
                Real x,y,z;
                AMREX_D_TERM(x = prob_lo[0] + (i+0.5)*dx[0] - center[0];,
                             y = prob_lo[1] + (j+0.5)*dx[1] - center[1];,
                             z = prob_lo[2] + (k+0.5)*dx[2] - center[2];);

                Real r = (AMREX_SPACEDIM == 2) ? std::sqrt(x*x+y*y) : std::sqrt(x*x+y*y+z*z);

                for (int n=0; n<nspecies; ++n) {
                    c(i,j,k,n) = c_init_1[n]*std::exp(-75.*r*r);
                }
            });

        } else if (prob_type == 15 || prob_type == -15) {

            /*
              case +/-15: mostly for testing electrodiffusion
              Discontinuous band in central 1/2 (case 15)
              c=c_init_1(:) inside; c=c_init_2(:) outside
            */

            // first quarter of domain
            Real y1 = (3.*prob_lo[1] + prob_hi[1]) / 4.;
            Real x1 = (3.*prob_lo[0] + prob_hi[0]) / 4.;

            // last quarter of domain
            Real y2 = (prob_lo[1] + 3*prob_hi[1]) / 4.;
            Real x2 = (prob_lo[0] + 3*prob_hi[0]) / 4.;

            amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {

                Real x,y,z;
                AMREX_D_TERM(x = prob_lo[0] + (i+0.5)*dx[0] - x1;,
                             y = prob_lo[1] + (j+0.5)*dx[1] - y1;,
                             z = prob_lo[2] + (k+0.5)*dx[2];);

                // tanh smoothing in y
                Real coeff=0.5*(std::tanh(y/(smoothing_width*dx[1]))+1.)*0.5*(std::tanh((-y+y2-y1)/(smoothing_width*dx[1]))+1.);

                // Donev: prob_type = -15: a special case for doing ternary diffusion NaCl + KCl
                // Here the last two species have a tanh profile in both x and y (species are Na+,Cl-,K+,water)
                // Aadd a tanh smoothing for central 50% of domain in x for second-to-last species
                if ( (prob_type==-15) && (n==nspecies-2) ) {
                    coeff=0.5*(std::tanh(x/(smoothing_width*dx[0]))+1.)*0.5*(std::tanh((-x+x2-x1)/(smoothing_width*dx[0]))+1.)*coeff;
                }

                // smooth between c_init_1(:) and c_init_2(:)
                Real c_loc = c_init_2[n] + (c_init_1[n]-c_init_2[n])*coeff;
                c(i,j,k,n) = c_loc;

                // for 4-species test, need to add Cl to central square to balance the K
                if ( (prob_type==-15) && (nspecies == 4) && (n == nspecies-2) ) {
                    c(i,j,k,1) = c(i,j,k,1) - charge_per_mass[2]/charge_per_mass[1]*c_loc;
                }
            });
        } else {
            Abort("Init.cpp: Invalid prob_type");
        }
    }

    for (MFIter mfi(rho_in); mfi.isValid(); ++mfi ) {

        Box bx = mfi.tilebox();

        const Array4<Real>& c = conc.array(mfi);
        const Array4<Real>& rho = rho_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // set final c_i such that sumtot(c_i) = 1 to within roundoff
            Real sumtot = 0.;
            for (int n=0; n<nspecies-1; ++n) {
                sumtot = sumtot + c(i,j,k,n);
            }

            c(i,j,k,nspecies-1) = 1. - sumtot;

            Real rho_total;

            // calculate rho_total from eos
            if (algorithm_type == 6) {
                rho_total = rho0;
            } else {
                sumtot = 0.;
                for (int n=0; n<nspecies; ++n) {
                    // sumtot represents rhoinv
                    sumtot = sumtot + c(i,j,k,n)/rhobar[n];
                }
                rho_total = 1./sumtot;
            }

            // calculate rho_i
            for (int n=0; n<nspecies; ++n) {
                rho(i,j,k,n) = rho_total*c(i,j,k,n);
            }

        });
     }
}
