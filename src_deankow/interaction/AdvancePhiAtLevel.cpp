#include <AmrCoreAdv.H>
#include <Kernels.H>
#include <myfunc.H>

using namespace amrex;

// Advance a single level for a single time step, updates flux registers
void
AmrCoreAdv::AdvancePhiAtLevel (int lev, Real /*time*/, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{
    constexpr double  PI = 3.14159265358979323846264338327950288;
    if (lev != 0) {
        amrex::Abort("FFT convolution is only initialized for level 0.");
    }

    Array<MultiFab, AMREX_SPACEDIM> fluxes;
    Array<MultiFab, AMREX_SPACEDIM> stochFluxes;
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        BoxArray ba = grids[lev];
        ba.surroundingNodes(i);
        fluxes[i].define(ba, dmap[lev], phi_new[lev].nComp(), 0);
        stochFluxes[i].define(ba, dmap[lev], phi_new[lev].nComp(), 0);
        stochFluxes[i].setVal(0.);
    }

//  compute convolution of U with phi to define interaction

        r2c_forward->forward(U, Uhat);
        r2c_forward->forward(phi_old[lev], Phihat);

       // Spectral multiply: Chat = Uhat * Phihat
        using Complex = GpuComplex<Real>;
        for (MFIter mfi(Uhat); mfi.isValid(); ++mfi) {
              auto const& ub = Uhat[mfi].box();
              auto const& u = Uhat[mfi].array();
              auto const& p = Phihat[mfi].array();
              auto const& c = Chat[mfi].array();
              ParallelFor(ub, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                  c(i,j,k) = u(i,j,k) * p(i,j,k);
//                  amrex::Print() << " i j fft " << i << " " << j << " " << p(i,j,k) << " " << u(i,j,k) << " " << c(i,j,k) << std::endl;
              });
          }


/*
    phi_old[lev].FillBoundary(Geom(lev).periodicity());


          MultiFab hack;
          BoxArray ba = grids[lev];
          hack.define(ba,dmap[0],1,0);
          hack.setVal(0.);
          Real dx = geom[0].CellSize(0);
          Real dy = geom[0].CellSize(1);
          
        for (MFIter mfi(hack); mfi.isValid(); ++mfi) {
              auto const& ub = hack[mfi].box();
              auto const& rho_arr = phi_old[lev][mfi].array();
              auto const& hack_arr = hack[mfi].array();
              ParallelFor(ub, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                  Real xi = (i+.5)*dx;
                  Real yj = (j+.5)*dy;
                  for(int ip=0 ; ip<32; ip++){
                  for(int jp=0 ; jp<32; jp++){
                     Real xip = (ip+.5)*dx;
                     Real yjp = (jp+.5)*dy;
                     amrex::Real eps, R, alpha;
                     eps = 0.0333;
                     R = 0.1;
                     alpha = 3.;

//    shiftinging definit of U by 1/2 for fft
                     amrex::Real x,y;
                     x = xi -xip;
                     y = yj -yjp;
                     Real r = std::sqrt(x*x+y*y);

                     Real utemp = eps * std::exp( - std::pow(r/R,3));

                     if(i == 15 && j == 15){
                          amrex::Print() << " ip jp " << ip << " " << jp << " " << utemp << std::endl;
                     }


                     hack_arr(i,j,k) += utemp * rho_arr(i,j,k)*dx*dy;
                     
                  }
                  }
                    // hack_arr(i,j,k) = 512.*std::sin(2.*PI*i*dx);
                  if(i > 10 && i < 20 && j > 10 && j < 20){
                  for(int ip=-9 ; ip<10; ip++){
                  for(int jp=-9 ; jp<10; jp++){

//                     amrex::Real ifact =(ip == 0) ? 1. : .5;
//                     amrex::Real jfact =(jp == 0) ? 1. : .5;
                     
                     Real r = std::sqrt(ip*ip*dx*dx+jp*jp*dy*dy);
                     amrex::Real eps, R, alpha;
                     eps = 0.0333;
                     R = 0.1;
                     alpha = 3.;

                     Real utemp = eps * std::exp( - std::pow(r/R,3));
                     if(i == 15 && j == 15){
                          amrex::Print() << " ip jp " << ip << " " << jp << " " << utemp << std::endl;
                     }
                     hack_arr(i,j,k) += utemp*rho_arr(i+ip,j+jp,k)*dx*dy;
                     
//                     hack_arr(i,j,k) += ifact*jfact*rho_arr(i+ip,j+jp,k)*dx*dy;
                  }
                  }
                  }
              });
          }
*/


          // Inverse FFT
          r2c_backward->backward(Chat, C);

          // Scale by 1/N (FFTW convention)
          Real mesh_scale =  AMREX_D_TERM(   geom[0].CellSize(0),
                              * geom[0].CellSize(1),
                              * geom[0].CellSize(2));
          const Real scaling = mesh_scale / geom[lev].Domain().d_numPts();
          C.mult(scaling, 0, 1);

/*

        for (MFIter mfi(hack); mfi.isValid(); ++mfi) {
              auto const& ub = hack[mfi].box();
              auto const& hack_arr = hack[mfi].array();
              auto const& C_arr = C[mfi].array();
              ParallelFor(ub, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                   amrex::Print() << "i j " << i << " " << j << " " << C_arr(i,j,k) << " " << hack_arr(i,j,k) << std::endl;
              });
         }

*/

    C.FillBoundary(Geom(lev).periodicity());

    phi_old[lev].FillBoundary(Geom(lev).periodicity());

    // We do this here so we can print the FABs for debugging
    phi_new[lev].setVal(0.0);

    advance_phi(phi_old[lev], phi_new[lev], fluxes, stochFluxes, C, dt_lev, num_part,  diff_coeff, dorand,  geom[lev], bcs);

    // Increment or decrement the flux registers by area and time-weighted fluxes
    // Note that the fluxes have already been scaled by dt and area
    // In this example we are solving phi_t = -div(+F)
    // The fluxes contain, e.g., F_{i+1/2,j} = (phi*u)_{i+1/2,j}
    // Keep this in mind when considering the different sign convention for updating
    // the flux registers from the coarse or fine grid perspective
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    if (do_reflux) {
        if (flux_reg[lev+1]) {
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                // update the lev+1/lev flux register (index lev+1)
                flux_reg[lev+1]->CrseInit(fluxes[i],i,0,0,fluxes[i].nComp(),1.0);
            }
        }
        if (flux_reg[lev]) {
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                // update the lev/lev-1 flux register (index lev)
                flux_reg[lev]->FineAdd(fluxes[i],i,0,0,fluxes[i].nComp(),-1.0);
            }
        }
    }
}
