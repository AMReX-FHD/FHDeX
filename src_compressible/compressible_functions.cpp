#include "compressible_functions.H"
#include "AMReX_ParmParse.H"

AMREX_GPU_MANAGED int compressible::transport_type;
AMREX_GPU_MANAGED int compressible::membrane_cell;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> compressible::transmission;
AMREX_GPU_MANAGED int compressible::do_1D;
AMREX_GPU_MANAGED int compressible::do_2D;
AMREX_GPU_MANAGED int compressible::all_correl;

void InitializeCompressibleNamespace()
{
    // extract inputs parameters
    ParmParse pp;

    int temp_max = std::max(3,MAX_SPECIES*MAX_SPECIES);    
    amrex::Vector<amrex::Real> temp    (temp_max,0.);
    amrex::Vector<int>         temp_int(temp_max,0 );

    // get transport type (1: Giovangigli; 2: Valk/Waldmann; 2: HCB)
    transport_type = 1; // Giovangigli (default)
    pp.query("transport_type",transport_type);
    switch (transport_type) {
        case 1:
            amrex::Print() << "Giovangigli transport model selected" << "\n";
            break;
        case 2:
            amrex::Print() << "Valk/Waldmann transport model selected" << "\n";
            break;
        case 3:
            amrex::Print() << "HCB binary transport model selected" << "\n";
            break;
    }

    // get membrane cell 
    membrane_cell = -1; // location of membrane (default)
    pp.query("membrane_cell",membrane_cell);
    if (membrane_cell >= 0) amrex::Print() << "Membrane cell is: " << membrane_cell << "\n";

    // get membrane transmission
    for (int i=0; i<MAX_SPECIES; ++i) {
        transmission[i] = 0.0;
    }
    if (pp.queryarr("transmission",temp,0,nspecies)) {
        amrex::Print() << "Membrane cell transmissions are: ";
        for (int i=0; i<nspecies; ++i) {
            transmission[i] = temp[i];
            amrex::Print() << transmission[i] << " ";
        }
        amrex::Print() << "\n";
    }

    // 1D simulation toggle
    do_1D = 0;
    pp.query("do_1D",do_1D);

    // 2D simulation toggle
    do_2D = 0;
    pp.query("do_2D",do_2D);
    // options for spatial correlations at multiple x*
    all_correl = 0;
    pp.query("all_correl",all_correl);


    return;
}


void GetHcGas() {
    for (int i=0; i<nspecies; ++i) {
        if (hcv[i] < 0.) {
            hcv[i] = 0.5*dof[i]*Runiv/molmass[i];
        }
        if (hcp[i] < 0.) {
            hcp[i] = 0.5*(2.+dof[i])*Runiv/molmass[i];
        }
    }
}


void InitConsVar(MultiFab& cons, const MultiFab& prim,
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
    
    // from namelist
    Real t_lo_y = t_lo[1];
    Real t_hi_y = t_hi[1];

    // local variables
    Real mach = 0.3;
    Real velscale = 30565.2*mach;

    Real hy = ( prob_hi[1] - prob_lo[1] ) / 3.;
    Real pi = acos(-1.);
    Real Lf = realhi[0] - reallo[0];

    for ( MFIter mfi(cons); mfi.isValid(); ++mfi ) {
        const Array4<const Real> pu = prim.array(mfi);
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
                GetPressureGas(pamb, massvec, cu(i,j,k,0), pu(i,j,k,4));
                
                Real molmix = 0.;

                for (int l=0; l<nspecies; ++l) {
                    molmix = molmix + massvec[l]/molmass[l];
                }
                molmix = 1.0/molmix;
                Real rgasmix = Runiv/molmix;
                Real alpha = grav[2]/(rgasmix*pu(i,j,k,4));

                // rho = exponential in z-dir to init @ hydrostatic eqm.
                // must satisfy system: dP/dz = -rho*g & P = rhogasmix*rho*T
                // Assumes temp=const
                cu(i,j,k,0) = pamb*exp(alpha*pos[2])/(rgasmix*pu(i,j,k,4));
                
                for (int l=0; l<nspecies; ++l) {
                    cu(i,j,k,5+l) = cu(i,j,k,0)*massvec[l];
                }

                Real intEnergy;
                GetEnergy(intEnergy, massvec, pu(i,j,k,4));

                cu(i,j,k,4) = cu(i,j,k,0)*intEnergy + 0.5*cu(i,j,k,0)*(pu(i,j,k,1)*pu(i,j,k,1) +
                                                                       pu(i,j,k,2)*pu(i,j,k,2) +
                                                                       pu(i,j,k,3)*pu(i,j,k,3));
            } else if (prob_type == 3) { // diffusion barrier

                for (int l=0; l<nspecies; ++l) {
                    Real Ygrad = (bc_Yk_y_hi[l] - bc_Yk_y_lo[l])/(realhi[1] - reallo[1]);
                    massvec[l] = Ygrad*pos[1] + bc_Yk_y_lo[l];
                    cu(i,j,k,5+l) = cu(i,j,k,0)*massvec[l];
                }

                Real intEnergy;
                GetEnergy(intEnergy, massvec, pu(i,j,k,4));
                cu(i,j,k,4) = cu(i,j,k,0)*intEnergy + 0.5*cu(i,j,k,0)*(pu(i,j,k,1)*pu(i,j,k,1) +
                                                                       pu(i,j,k,2)*pu(i,j,k,2) +
                                                                       pu(i,j,k,3)*pu(i,j,k,3));
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
                   Real Ly = realhi[1] - reallo[0];
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

        });
    } // end MFIter

}
