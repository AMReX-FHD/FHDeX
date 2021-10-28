#include "compressible_functions.H"
#include "compressible_functions_stag.H"

void InitConsVarStag(MultiFab& cons,
                     std::array< MultiFab, AMREX_SPACEDIM >& momStag,
                     const amrex::Geometry geom) {

    const Real* dx_host = geom.CellSize();
    const RealBox& realDomain = geom.ProbDomain();

    GpuArray<Real,AMREX_SPACEDIM> dx;
    GpuArray<Real,AMREX_SPACEDIM> reallo;
    GpuArray<Real,AMREX_SPACEDIM> realhi;
    GpuArray<Real,AMREX_SPACEDIM> center;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        dx[d] = dx_host[d];
        reallo[d] = realDomain.lo(d);
        realhi[d] = realDomain.hi(d);
        center[d] = ( realhi[d] - reallo[d] ) / 2.;
    }
    
    Real t_lo_y = t_lo[1];
    Real t_hi_y = t_hi[1];

    // local variables
    Real mach = 0.3;
    Real velscale = 30565.2*mach;

    Real hy = ( prob_hi[1] - prob_lo[1] ) / 3.;
    Real pi = acos(-1.);
    Real Lf = realhi[0] - reallo[0];

    // compute some values and overwrite based on prob_type
        
    // compute internal energy
    Real intEnergy;
    GpuArray<Real,MAX_SPECIES  > massvec;
    for(int i=0;i<nspecies;i++) {
        massvec[i] = rhobar[i];
    }
    GetEnergy(intEnergy, massvec, T_init[0]);

    cons.setVal(0.0,0,nvars,ngc);
    cons.setVal(rho0,0,1,ngc);           // density
    cons.setVal(0,1,3,ngc);              // x/y/z momentum
    cons.setVal(rho0*intEnergy,4,1,ngc); // total energy
    for(int i=0;i<nspecies;i++) {
        cons.setVal(rho0*rhobar[i],5+i,1,ngc); // mass densities
    }

    for (int d=0; d<AMREX_SPACEDIM; d++) { // staggered momentum & velocities
        momStag[d].setVal(0.);
    }
    
    for ( MFIter mfi(cons); mfi.isValid(); ++mfi ) {
        const Array4<      Real> cu = cons.array(mfi);

        const Box& bx = mfi.tilebox();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            GpuArray<Real,AMREX_SPACEDIM> itVec;
            GpuArray<Real,AMREX_SPACEDIM> pos;
            GpuArray<Real,AMREX_SPACEDIM> relpos;

            GpuArray<Real,MAX_SPECIES> massvec;
            GpuArray<Real,MAX_SPECIES> Yk;

            AMREX_D_TERM(itVec[0] = (i+0.5)*dx[0]; ,
                         itVec[1] = (j+0.5)*dx[1]; ,
                         itVec[2] = (k+0.5)*dx[2]);
            
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                pos[d] = reallo[d] + itVec[d];
                relpos[d] = pos[d] - center[d];
            }

            // Total density must be pre-set
            
            if (prob_type == 2) { // Rayleigh-Taylor
                
                if (relpos[2] >= 0.) {
                    massvec[0] = 0.4;
                    massvec[1] = 0.4;
                    massvec[2] = 0.1;
                    massvec[3] = 0.1;
                } else {
                    massvec[0] = 0.1;
                    massvec[1] = 0.1;
                    massvec[2] = 0.4;
                    massvec[3] = 0.4;
                }

                Real pamb;
                GetPressureGas(pamb, massvec, cu(i,j,k,0), T_init[0]);
                
                Real molmix = 0.;

                for (int l=0; l<nspecies; ++l) {
                    molmix = molmix + massvec[l]/molmass[l];
                }
                molmix = 1.0/molmix;
                Real rgasmix = Runiv/molmix;
                Real alpha = grav[2]/(rgasmix*T_init[0]);

                // rho = exponential in z-dir to init @ hydrostatic eqm.
                // must satisfy system: dP/dz = -rho*g & P = rhogasmix*rho*T
                // Assumes temp=const
                cu(i,j,k,0) = pamb*exp(alpha*pos[2])/(rgasmix*T_init[0]);
                
                for (int l=0; l<nspecies; ++l) {
                    cu(i,j,k,5+l) = cu(i,j,k,0)*massvec[l];
                }

                Real intEnergy;
                GetEnergy(intEnergy, massvec, T_init[0]);

                cu(i,j,k,4) = cu(i,j,k,0)*intEnergy + 0.5*(cu(i,j,k,1)*cu(i,j,k,1) +
                                                           cu(i,j,k,2)*cu(i,j,k,2) +
                                                           cu(i,j,k,3)*cu(i,j,k,3)) / cu(i,j,k,0);
            } else if (prob_type == 3) { // diffusion barrier

                for (int l=0; l<nspecies; ++l) {
                    Real Ygrad = (bc_Yk_y_hi[l] - bc_Yk_y_lo[l])/(realhi[1] - reallo[1]);
                    massvec[l] = Ygrad*pos[1] + bc_Yk_y_lo[l];
                    cu(i,j,k,5+l) = cu(i,j,k,0)*massvec[l];
                }

                Real intEnergy;
                GetEnergy(intEnergy, massvec, T_init[0]);
                cu(i,j,k,4) = cu(i,j,k,0)*intEnergy + 0.5*(cu(i,j,k,1)*cu(i,j,k,1) +
                                                           cu(i,j,k,2)*cu(i,j,k,2) +
                                                           cu(i,j,k,3)*cu(i,j,k,3)) / cu(i,j,k,0);
            } else if (prob_type == 4) { // Taylor Green Vortex

                Real x=itVec[0];
                Real y=itVec[1];
                Real z=itVec[2];

                cu(i,j,k,0) = 1.784e-3;
                cu(i,j,k,1) =  velscale*cu(i,j,k,0)*sin(2.*pi*x/Lf)*cos(2.*pi*y/Lf)*cos(2.*pi*z/Lf);
                cu(i,j,k,2) = -velscale*cu(i,j,k,0)*cos(2.*pi*x/Lf)*sin(2.*pi*y/Lf)*cos(2.*pi*z/Lf);
                cu(i,j,k,3) = 0.;
                Real pres = 1.01325e6+cu(i,j,k,0)*velscale*velscale*cos(2.*pi*x/Lf)*cos(4.*pi*y/Lf)*(cos(4.*pi*z/Lf)+2.);
                cu(i,j,k,4) = pres/(5./3.-1.) + 0.5*(cu(i,j,k,1)*cu(i,j,k,1) +
                                                     cu(i,j,k,2)*cu(i,j,k,2) +
                                                     cu(i,j,k,3)*cu(i,j,k,3)) / cu(i,j,k,0);
                cu(i,j,k,5) = cu(i,j,k,0);
                cu(i,j,k,6) = 0.;
                
            } else if (prob_type == 5) { // Taylor Green Vortex

                Real intEnergy;
                
                cu(i,j,k,0) = rho0;
                cu(i,j,k,1) = 0;
                cu(i,j,k,2) = 0;
                cu(i,j,k,3) = 0;
                if((prob_lo[1] + itVec[1]) < hy) {
                    massvec[0] = bc_Yk_x_lo[0];
                    massvec[1] = bc_Yk_x_lo[1];
                    GetEnergy(intEnergy, massvec, t_lo_y);
                    cu(i,j,k,4) = cu(i,j,k,0)*intEnergy;
                    cu(i,j,k,5) = cu(i,j,k,0)*bc_Yk_x_lo[0];
                    cu(i,j,k,6) = cu(i,j,k,0)*bc_Yk_x_lo[1];
                } else if ((prob_lo[1] + itVec[1]) < 2*hy) {
                    massvec[0] = bc_Yk_x_hi[0];
                    massvec[1] = bc_Yk_x_hi[1];
                    GetEnergy(intEnergy, massvec, t_hi_y);
                    cu(i,j,k,4) = cu(i,j,k,0)*intEnergy;
                    cu(i,j,k,5) = cu(i,j,k,0)*bc_Yk_x_hi[0];
                    cu(i,j,k,6) = cu(i,j,k,0)*bc_Yk_x_hi[1];
                } else {
                    massvec[0] = bc_Yk_x_lo[0];
                    massvec[1] = bc_Yk_x_lo[1];
                    GetEnergy(intEnergy, massvec, t_lo_y);
                    cu(i,j,k,4) = cu(i,j,k,0)*intEnergy;
                    cu(i,j,k,5) = cu(i,j,k,0)*bc_Yk_x_lo[0];
                    cu(i,j,k,6) = cu(i,j,k,0)*bc_Yk_x_lo[1];
                }
            } else if (prob_type == 100) { // sinusoidal density variation

                    Real y = itVec[1];
                    Real Ly = realhi[1] - reallo[1];
                    for (int l=0;l<nspecies;l++) {
                        Yk[l] = cu(i,j,k,5+l)/cu(i,j,k,0);
                    }
                    cu(i,j,k,0) = rho0 + 0.1*rho0*sin(2.*pi*y/Ly);
                    for (int l=0;l<nspecies;l++) {
                        cu(i,j,k,5+l) = cu(i,j,k,0)*Yk[l];
                    }
            }
            else if (prob_type == 101) { // sinusoidal temperature variation (constant pressure)

                   Real y = itVec[1];
                   Real Ly = realhi[1] - reallo[0];

                   for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];

                   Real pressure;
                   GetPressureGas(pressure,massvec,rho0,T_init[0]);

                   Real temperature;
                   temperature = T_init[0] + 0.1*T_init[0]*sin(2.*pi*y/Ly);

                   Real density;
                   GetDensity(pressure,density,temperature,massvec);
                   cu(i,j,k,0) = density;
                   for (int ns=0;ns<nspecies;++ns) cu(i,j,k,5+ns) = density*massvec[ns];

                   Real intEnergy;
                   GetEnergy(intEnergy, massvec, temperature);
                   cu(i,j,k,4) = density*intEnergy;
            }
            else if (prob_type == 102) { // two temperature across membrane
                    
                    Real intEnergy;
                    cu(i,j,k,0) = rho0;
                    massvec[0] = 0.25; massvec[1] = 0.25; massvec[2] = 0.25; massvec[3] = 0.25;
                    if (i < membrane_cell) {
                        GetEnergy(intEnergy, massvec, t_lo[0]);
                        cu(i,j,k,4) = cu(i,j,k,0)*intEnergy;
                    }
                    else {
                        GetEnergy(intEnergy, massvec, t_hi[0]);
                        cu(i,j,k,4) = cu(i,j,k,0)*intEnergy;
                    }
            } // prob type
            else if (prob_type == 103) { // double the pressure in other half
                
                    for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
                    Real pressure;
                    GetPressureGas(pressure,massvec,rho0,T_init[0]);

                    if (relpos[0] > 0.0) {
                        Real pressure_new = 2.0*pressure;
                        Real temperature;
                        temperature = T_init[0];

                        Real density;
                        GetDensity(pressure_new,density,temperature,massvec);
                        cu(i,j,k,0) = density;
                        for (int ns=0;ns<nspecies;++ns) cu(i,j,k,5+ns) = density*massvec[ns];

                        Real intEnergy;
                        GetEnergy(intEnergy, massvec, temperature);
                        cu(i,j,k,4) = density*intEnergy;
                    }
                    
            } // prob type
            else if (prob_type == 104) { // temperature discontinuity 
                
                    if (i==15) 
                    {
                        for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
                        Real pamb;
                        GetPressureGas(pamb,massvec,rho0,T_init[0]);

                        Real density;
                        GetDensity(pamb,density,400.0,massvec);
                        cu(i,j,k,0) = density;
                        for (int ns=0;ns<nspecies;++ns) cu(i,j,k,5+ns) = density*massvec[ns];

                        Real intEnergy;
                        GetEnergy(intEnergy, massvec, 400.0);
                        cu(i,j,k,4) = density*intEnergy;
                    }
                    
            }
            else if (prob_type == 105) { // pressure discontinuity 
                
                    if (i==15) 
                    {
                        for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
                        Real pamb;
                        GetPressureGas(pamb,massvec,rho0,T_init[0]);

                        Real density;
                        GetDensity(pamb*1.5,density,T_init[0],massvec);
                        cu(i,j,k,0) = density;
                        for (int ns=0;ns<nspecies;++ns) cu(i,j,k,5+ns) = density*massvec[ns];

                        Real intEnergy;
                        GetEnergy(intEnergy, massvec, T_init[0]);
                        cu(i,j,k,4) = density*intEnergy;
                    }
            }
            else if (prob_type == 106) { // concentration discontinuity 
                
                    if (i==15) 
                    {
                        for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
                        Real pamb;
                        GetPressureGas(pamb,massvec,rho0,T_init[0]);

                        massvec[0] = 0.4;
                        massvec[1] = 0.6;
                        Real density;
                        GetDensity(pamb,density,T_init[0],massvec);
                        cu(i,j,k,0) = density;
                        for (int ns=0;ns<nspecies;++ns) cu(i,j,k,5+ns) = density*massvec[ns];

                        Real intEnergy;
                        GetEnergy(intEnergy, massvec, T_init[0]);
                        cu(i,j,k,4) = density*intEnergy;
                    }
                    
            } // prob type

            else if (prob_type == 107) { // two-fluid Rayleigh Taylor. 0: lighter species; 1: heavier
                
                if (relpos[2] >= 0.) {
                    massvec[0] = 0.1;
                    massvec[1] = 0.9;
                } else {
                    massvec[0] = 0.9;
                    massvec[1] = 0.1;
                }

                Real pamb;
                GpuArray<Real,MAX_SPECIES  > massvec_bot;
                massvec_bot[0] = 0.9; massvec_bot[1] = 0.1;
                GetPressureGas(pamb, massvec_bot, rho0, T_init[0]);
                
                Real molmix = 0.;

                for (int l=0; l<nspecies; ++l) {
                    molmix = molmix + massvec[l]/molmass[l];
                }
                molmix = 1.0/molmix;
                Real rgasmix = Runiv/molmix;
                Real alpha = grav[2]/(rgasmix*T_init[0]);

                // rho = exponential in z-dir to init @ hydrostatic eqm.
                // must satisfy system: dP/dz = -rho*g & P = rhogasmix*rho*T
                // Assumes temp=const
                Real press = pamb*exp(alpha*pos[2]);
                cu(i,j,k,0) = pamb*exp(alpha*pos[2])/(rgasmix*T_init[0]);
                
                for (int l=0; l<nspecies; ++l) {
                    cu(i,j,k,5+l) = cu(i,j,k,0)*massvec[l];
                }

                Real intEnergy;
                GetEnergy(intEnergy, massvec, T_init[0]);

                cu(i,j,k,4) = cu(i,j,k,0)*intEnergy + 0.5*(cu(i,j,k,1)*cu(i,j,k,1) +
                                                           cu(i,j,k,2)*cu(i,j,k,2) +
                                                           cu(i,j,k,3)*cu(i,j,k,3)) / cu(i,j,k,0);
            }


        });
    } // end MFIter
}

void PrintFluxes(std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in, std::array< MultiFab, 2 >& edgeflux_x_in,
                 std::array< MultiFab, 2 >& edgeflux_y_in, std::array< MultiFab, 2 >& edgeflux_z_in,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in, std::string prefix)
{
    std::string xflux = prefix + "xface";
    std::string yflux = prefix + "yface";
    std::string zflux = prefix + "zface";

    std::string edgex_v = prefix + "edgexv";
    std::string edgex_w = prefix + "edgexw";
    std::string edgey_u = prefix + "edgeyu";
    std::string edgey_w = prefix + "edgeyw";
    std::string edgez_u = prefix + "edgezu";
    std::string edgez_v = prefix + "edgezv";

    std::string cenx = prefix + "cenx";
    std::string ceny = prefix + "ceny";
    std::string cenz = prefix + "cenz";

    outputMFAscii(faceflux_in[0],xflux);
    outputMFAscii(faceflux_in[1],yflux);
    outputMFAscii(faceflux_in[2],zflux);

    outputMFAscii(edgeflux_x_in[0],edgex_v);
    outputMFAscii(edgeflux_x_in[1],edgex_w);
    outputMFAscii(edgeflux_y_in[0],edgey_u);
    outputMFAscii(edgeflux_y_in[1],edgey_w);
    outputMFAscii(edgeflux_z_in[0],edgez_u);
    outputMFAscii(edgeflux_z_in[1],edgez_v);

    outputMFAscii(cenflux_in[0],cenx);
    outputMFAscii(cenflux_in[1],ceny);
    outputMFAscii(cenflux_in[2],cenz);
}
