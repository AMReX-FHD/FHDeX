#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Print.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_phi (int i, int j, int k,
               amrex::Array4<amrex::Real> const& phi,
               amrex::Array4<amrex::Real> const& det,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
               amrex::Real surf_area, int pure_part, int Ncomp, int ext_pot)
{
    using amrex::Real;;
    Real cellvol = dx[0]*dx[1];
    constexpr double  TWOPI = 2.*3.14159265358979323846264338327950288;

#if(AMREX_SPACEDIM > 2)
    cellvol *=dx[2];
#endif

    if(pure_part == 0){

    //phi(i,j,k,0) = cellvol * det(i,j,k,1) /surf_area;
//    phi(i,j,k,0) = 1. /surf_area;
//   Real x = prob_lo[0] + (i+Real(0.5)) * dx[0] -0.5*TWOPI ;
   Real y = prob_lo[1] + (j+Real(0.5)) * dx[1] -0.5*TWOPI ;
   Real x = prob_lo[0] + (i+Real(0.5)) * dx[0]  ;

   if(ext_pot == 1){

//    if(4.*(x-0.1*TWOPI)*(x-0.1*TWOPI)+.25*y*y < .01*TWOPI*TWOPI){
    if(2.*(x-0.15*TWOPI)*(x-0.15*TWOPI)+.5*(y-.1*TWOPI)*(y-.1*TWOPI) < .01*TWOPI*TWOPI){
        phi(i,j,k,0) = 1000.;
    } else {
        phi(i,j,k) = 0.;
    }

    } else {

       phi(i,j,k,0) = 1./surf_area;

    }

    } else {

   Real x = prob_lo[0] + (i+Real(0.5)) * dx[0] -0.5*TWOPI ;
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1] -0.5*TWOPI ;
//
#if 0
    if(x*x+y*y < .01*TWOPI*TWOPI){
       phi(i,j,k,0) = 806.1;
    } else {
       phi(i,j,k,0) = 0.;
    }
#else
      phi(i,j,k,0) = 10.1;
#endif
     

    }

    //amrex::Print() << "vol,dat,valu" << cellvol << " " << det(i,j,k,1)  << " " << surf_area << " " << phi(i,j,k,0) << std::endl;


#if 0
    Real rad = Real(0.1);
    Real xcent = Real(0.5);
    Real ycent = Real(0.5);
    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0] - xcent;
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1] - ycent;

#if (AMREX_SPACEDIM > 2)
    Real zcent = Real(0.5);
    Real z = prob_lo[2] + (k+Real(0.5)) * dx[2] - zcent;
    Real r2 = (x*x + y*y + z*z) / (rad*rad);
#else
    Real z = Real(0.);
    Real r2 = (x*x + y*y) / (rad*rad);
#endif

    // Gaussian
    phi(i,j,k,0) = Real(1.) + std::exp(-r2);
#endif

#if 0
    // Jump in x-direction
    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0];
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1];
    if ( x > 0.5) {
        phi(i,j,k,0) = 1.5*npts_scale;
    } else {
        phi(i,j,k,0) = 0.5*npts_scale;
    }
#endif

#if 0
    // Delta function
    if (i==4 && j ==4)
    {
       phi(i,j,k,0) = 1./(dx[0]*dx[1]);
    }
    else
    {
       phi(i,j,k,0) = 0;;
    }
#endif

#if 0
    // 2D Test for multilevel
    Real rad = Real(0.2);
    Real xcent = Real(0.5);
    Real ycent = Real(0.5);
    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0] - xcent;
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1] - ycent;
    Real r2 = (x*x + y*y);
//    r2 = 4.*(x-y)*(x-y)+ .25*(x+y)*(x+y);
#if (AMREX_SPACEDIM == 3)
    Real zcent = Real(0.75);
    Real z = prob_lo[2] + (k+Real(0.5)) * dx[2] - zcent;
    r2 = (x*x + y*y+z*z);
#endif

//    if ( r2 < rad*rad) {
    if (i == 50 && j == 50) {
        phi(i,j,k,0) = 10.1/cellvol;
    } else {
        phi(i,j,k,0) = 0./cellvol;
    }
#endif

#if 0
    // 2D Test for multilevel
    Real rad = Real(0.4);
    Real radi = Real(0.1);
    Real xcent = Real(0.5);
    Real ycent = Real(0.5);
    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0] - xcent;
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1] - ycent;
    Real r2 = (x*x + y*y);
    r2 = 2.*(x-y)*(x-y)+ .5*(x+y)*(x+y);
    Real r3 = 0.5*(x-y)*(x-y)+ 2.*(x+y)*(x+y);
#if (AMREX_SPACEDIM == 3)
    Real zcent = Real(0.5);
    Real z = prob_lo[2] + (k+Real(0.5)) * dx[2] - zcent;
    r3 = (x*x + y*y+z*z);
    r2 = (x-.1)*(x-.1)+(y-.2)*(y-.2)+z*z;
#endif

    if ( r3 > rad*rad) {
        //phi(i,j,k,0) = 50./cellvol;
        phi(i,j,k,0) = 30./cellvol;
    } else if (r2 < radi*radi) {
        //phi(i,j,k,0) = 20./cellvol;
        phi(i,j,k,0) = 15./cellvol;
    } else {
        phi(i,j,k,0) = 0./cellvol;
    }
    //if(std::abs(x) <= .1 && y > 0.){
    if(x >= .3 && y > .3){
       phi(i,j,k,0) = 0.;
    } else {
       phi(i,j,k,0) = 30./cellvol;
    }

//
#endif
#if 0
//    int ioff = 25;
//    int joff = 25;
    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0] ;
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1] ;

    Real value = 200./(.01 * .01);

//    if( i-ioff >= 0 && i - ioff < 10 && j-joff >= 0 && j - joff < 50 ) {
    if( x >= .25 && x < .35 && y >= .25 && y < .75 ) {
    //if( x >= .2 && x < .3 && y >= .25 && y < .75 ) {
       //phi(i,j,k,0) = 200./cellvol;
       phi(i,j,k,0) = value;
    } else {
       phi(i,j,k,0) = 0./cellvol;
    }

//

#endif
#if 0

        int ext_pot = 1;
        amrex::Real alpha = .25;
        amrex::Real beta = .75;
        amrex::Real gamma = 5.e-4;
        gamma = 1.1e-3;
        Real x = prob_lo[0] + (i+Real(0.5)) * dx[0] ;
        Real y = prob_lo[1] + (j+Real(0.5)) * dx[1] ;

    if( x >= .2 && x < .4 && y >= .25 && y < .75 ) {
        amrex::Real pot = (x-alpha)*(x-alpha)*(x-beta)*(x-beta);
        pot += 0.5* (y-.5)*(y-.5)*(y-.5)*(y-.5);
        pot /= gamma;
        Real value = 20000.* 200.* std::exp(-pot);
        phi(i,j,k,0) = value;
    } else {
       phi(i,j,k,0) = 0./cellvol;
    }



#endif

    if (Ncomp == 2){
        phi(i,j,k,1) = phi(i,j,k,0);
    }
}

#if 0
void init_metric(int i, int j, int k,
               amrex::Array4<amrex::Real> const& gmet,
               amrex::Array4<amrex::Real> const& sqrmet,
               amrex::Array4<amrex::Real> const& detg,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo)
{
    constexpr double  TWOPI = 2.*3.14159265358979323846264338327950288;
    constexpr double  PI = 3.14159265358979323846264338327950288;

    amrex::Real xloc = prob_lo[0] + (i+0.5) * dx[0] ;
    amrex::Real yloc = prob_lo[1] + (j+0.5) * dx[1] ;

    amrex::Real amp = 4.0;

    amrex::Real caps, capsy, capr, caprx;

    caps  = 0.;
    capsy = 0.;

    capr  = std::pow(std::cos(0.5*(xloc-PI)),4);
    caprx = -2.*std::sin(0.5*(xloc-PI))*std::pow(std::cos(0.5*(xloc-PI)),3);

    for (int iterm  = 0; iterm < 3; iterm++)
    {
       amrex::Real shift = PI/3. + 2.*iterm*PI/3.;
       caps += std::pow(std::cos(0.5*(yloc-shift)),16);
       capsy += -8*std::sin(0.5*(yloc-shift))*std::pow(std::cos(0.5*(yloc-shift)),15);
    }

//    amrex::Real sinx,siny,cosx,cosy,det,sig11,sig12,sig21,sig22;

    amrex::Real p,q;
    p = amp*caps*caprx;
    q = amp*capsy*capr;


//    sinx = std::sin(xloc);
//    siny = std::sin(yloc);
//    cosx = std::cos(xloc);
//    cosy = std::cos(yloc);

//    det = 1. + amp*amp*(cosx*cosx*siny*siny+sinx*sinx*cosy*cosy);

    amrex::Real det = 1. + p*p + q*q;

    detg(i,j,k,0) = det;
    detg(i,j,k,1) = std::sqrt(det);
    detg(i,j,k,2) = std::sqrt(detg(i,j,k,1));

    gmet(i,j,k,0) = (1.+ q*q)/det;
    gmet(i,j,k,1) = -p*q/det;
    gmet(i,j,k,2) = -p*q/det;
    gmet(i,j,k,3) = (1.+ p*p)/det;

    amrex::Real factor = 1./(detg(i,j,k,1)*(1.+detg(i,j,k,1)));

    sqrmet(i,j,k,0) = (1.+q*q+detg(i,j,k,1))*factor;
    sqrmet(i,j,k,1) =  - p*q*factor;
    sqrmet(i,j,k,2) =  - p*q*factor;
    sqrmet(i,j,k,3) = (1.+p*p+detg(i,j,k,1))*factor;

}
#endif


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_x (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxx,
                     amrex::Array4<amrex::Real> const& stochfluxx,
                     amrex::Array4<amrex::Real> const& stochfluxy,
                     amrex::Array4<amrex::Real const> const& phi,
                     amrex::Array4<amrex::Real const> const& gmet,
                     amrex::Array4<amrex::Real const> const& sqrgmet,
                     amrex::Array4<amrex::Real const> const& det,
                     amrex::Real num_part,
                     amrex::Real dxinv,
                     amrex::Real dyinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi,int Ncomp, int num_flux, int ext_pot)
{
    constexpr double  TWOPI = 2.*3.14159265358979323846264338327950288;
    constexpr double  PI = 3.14159265358979323846264338327950288;

    int coeff_index = Ncomp-1;
    amrex::Real phibdry = 0.;
    amrex::Real stofact = 1.;

//    int ext_pot = 0;
#if 0
    amrex::Real alpha = .25;
    amrex::Real beta = .75;
    amrex::Real gamma = 5.e-4;
    gamma = 2.e-1;
#else
    amrex::Real alpha = .3;
    amrex::Real beta = 0.;
    amrex::Real delta = 2.;
    amrex::Real gamma = 1.e-1;
#endif

    amrex::Real phip,phim, phipn, phimn;

/*
    if (lo == dom_lo &&
        (bc_lo == amrex::BCType::foextrap ||
         bc_lo == amrex::BCType::ext_dir))
    {
        if(i == lo)
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv / 0.5;
        }
        else
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv;
        }
    }
    else if (hi == dom_hi &&
             (bc_hi == amrex::BCType::foextrap ||
              bc_hi == amrex::BCType::ext_dir))
    {
        if(i == hi+1)
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv / 0.5;
        }
        else
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv;
        }
    }
*/
    //amrex::Print() << "BC info " << i << " " << dom_lo << " " << bc_lo << " " << dom_hi << " " << bc_hi << std::endl;
#if 0
    if (i == dom_lo && (bc_lo == amrex::BCType::foextrap || bc_lo == amrex::BCType::ext_dir)) {

        stofact = std::sqrt(2.);

        phim = 2.*phibdry - phi(i,j,k,0);
        phip = phi(i,j,k,0);

        //phimn = 2.*phibdry - phi(i,j,k,coeff_index);
        //phipn = phi(i,j,k,coeff_index);
        phimn = phibdry;
        phipn = phibdry;
        //amrex::Print() << " A " << std::endl;

    } else if (i == dom_hi+1 && (bc_hi == amrex::BCType::foextrap || bc_hi == amrex::BCType::ext_dir)) {

        stofact = std::sqrt(2.);
        phip = 2.*phibdry - phi(i-1,j,k,0);
        phim = phi(i-1,j,k,0);
        //phipn = 2.*phibdry - phi(i-1,j,k,coeff_index);
        //phimn = phi(i-1,j,k,coeff_index);
        phimn = phibdry;
        phipn = phibdry;
        //amrex::Print() << " B " << std::endl;
    } else {

        phip = phi(i,j,k,0);
        phim = phi(i-1,j,k,0);
        phipn = std::abs(phi(i,j,k,coeff_index));
        phimn = std::abs(phi(i-1,j,k,coeff_index));
        stofact = 1.;
        //amrex::Print() << " C " << std::endl;

    }
#endif

//	amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i-1,j,k))/(phi(i,j,k)+phi(i-1,j,k)+1.e-16);
	//amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
	//amrex::Real phiavg = 0.5*(phip+phim);

	// HACK TO COMMENT OUT STOCHASTIC TERMS
	// phiavg = 0.;

        //fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv
    	//                    + stochfluxx(i,j,k) * std::sqrt(phiavg);

//          if(ext_pot == 1) {
//
//             amrex::Real xloc;
//             xloc = i / dxinv;
//           amrex::Real phiavg = 0.5*(phim+phip);

//               amrex::Real Vsubx;

//             Vsubx = 2.*(xloc - beta) * (xloc - alpha)* (2.*xloc - alpha - beta) / gamma;

 //            fluxx(i,j,k,0) += phiavg * Vsubx;

//          }
          //amrex::Print() << " point " << i << " " << j << " " <<  dxinv << std::endl;
          //amrex::Print() << "phip " << phip << std::endl;
          //amrex::Print() << "phib " << phibdry << std::endl;
          //amrex::Print() << "phi " << phi(i,j,k,0) << std::endl;
          //amrex::Print() << "phim " << phim << std::endl;

          amrex::Real diffphixm =( phi(i,j,k,0) - phi(i-1,j,k,0))*dxinv;
          amrex::Real diffphiym =( phi(i,j,k,0) - phi(i,j-1,k,0))*dyinv;

          fluxx(i,j,k,0) = det(i,j,k,1)*(gmet(i,j,k,0)*diffphixm + gmet(i,j,k,1)*diffphiym);

          if(num_flux ==4){
            // Four one-sided gradients on the 2x2 stencil around the face (i,j)
            // -- at (i,j)
            amrex::Real gx_mm = ( phi(i,j,k,0)   - phi(i-1,j,k,0) )*dxinv;
            amrex::Real gy_mm = ( phi(i,j,k,0)   - phi(i,j-1,k,0) )*dyinv;
            // ++ at (i-1,j-1)
            //amrex::Real gx_pp = ( phi(i,j-1,k,0) - phi(i-1,j-1,k,0) )*dxinv;
            //amrex::Real gy_pp = ( phi(i-1,j,k,0) - phi(i-1,j-1,k,0) )*dyinv;
            // +- at (i-1,j)
            //amrex::Real gx_pm = ( phi(i,j,k,0)   - phi(i-1,j,k,0) )*dxinv;
            //amrex::Real gy_pm = ( phi(i-1,j,k,0) - phi(i-1,j-1,k,0) )*dyinv;
            // -+ at (i,j-1)
            //amrex::Real gx_mp = ( phi(i,j-1,k,0) - phi(i-1,j-1,k,0) )*dxinv;
            //amrex::Real gy_mp = ( phi(i,j,k,0)   - phi(i,j-1,k,0) )*dyinv;

            //fluxx(i,j,k,0) = det(i,j,k,1)*(gmet(i,j,k,0)*gx_mm + gmet(i,j,k,1)*gy_mm);
            //fluxx(i,j,k,1) = det(i-1,j-1,k,1)*(gmet(i-1,j-1,k,0)*gx_pp + gmet(i-1,j-1,k,1)*gy_pp);
            //fluxx(i,j,k,2) = det(i-1,j,k,1)*(gmet(i-1,j,k,0)*gx_pm + gmet(i-1,j,k,1)*gy_pm);
            //fluxx(i,j,k,3) = det(i,j-1,k,1)*(gmet(i,j-1,k,0)*gx_mp + gmet(i,j-1,k,1)*gy_mp);

            amrex::Real gx_pp = gx_mm;
            amrex::Real gy_pp = ( phi(i-1,j+1,k,0) - phi(i-1,j,k,0) )*dyinv;

            amrex::Real gx_pm = ( phi(i,j,k,0)   - phi(i-1,j,k,0) )*dxinv;
            amrex::Real gy_pm = ( phi(i,j+1,k,0) - phi(i,j,k,0) )*dyinv;

            amrex::Real gx_mp = ( phi(i,j,k,0) - phi(i-1,j,k,0) )*dxinv;
            amrex::Real gy_mp = ( phi(i-1,j,k,0)   - phi(i-1,j-1,k,0) )*dyinv;

            fluxx(i,j,k,0) = det(i,j,k,1)*(gmet(i,j,k,0)*gx_mm + gmet(i,j,k,1)*gy_mm);
            fluxx(i,j,k,1) = det(i-1,j,k,1)*(gmet(i-1,j,k,0)*gx_pp + gmet(i-1,j,k,1)*gy_pp);
            fluxx(i,j,k,2) = det(i,j,k,1)*(gmet(i,j,k,0)*gx_pm + gmet(i,j,k,1)*gy_pm);
            fluxx(i,j,k,3) = det(i-1,j,k,1)*(gmet(i-1,j,k,0)*gx_mp + gmet(i-1,j,k,1)*gy_mp);


          }

          if(ext_pot == 1) {


             amrex::Real xloc,yloc;
             xloc = (i+0.5) / dxinv;
             yloc = (j+0.5) / dyinv;


             amrex::Real Vsubx;
             amrex::Real Vsuby;
#if 1

             amrex::Real caps, capsy, capr, caprx;

          
          if(num_flux == 1){

             caps = 0;
             capsy = 0;
             capr  = (delta-alpha*std::pow(std::sin(xloc),2))*(1.+std::pow(std::cos(xloc/2),6));
             caprx = -2.*alpha*(1.+std::pow(std::cos(xloc/2),6))*std::cos(xloc)*std::sin(xloc)
                     -3.*std::pow(std::cos(xloc/2),5)*std::sin(xloc/2)*(delta-alpha*std::pow(std::sin(xloc),2));


             for (int iterm  = 0; iterm < 3; iterm++)
             {
                amrex::Real shift = PI/3. + 2.*iterm*PI/3.;
                caps += std::pow(std::cos(0.5*(yloc-shift)),16);
                capsy += -8*std::sin(0.5*(yloc-shift))*std::pow(std::cos(0.5*(yloc-shift)),15);
             }


             Vsubx = caprx/gamma;
             Vsuby = capsy/gamma;


#else
             amrex::Real sinxpy = std::sin(0.5*(xloc+yloc));
             amrex::Real cosxpy = std::cos(0.5*(xloc+yloc));
             amrex::Real sinxmy = std::sin(0.5*(xloc-yloc));
             amrex::Real cosxmy = std::cos(0.5*(xloc-yloc));

             Vsubx = (sinxpy*cosxpy*sinxmy*sinxmy + sinxpy*sinxpy*sinxmy*cosxmy)/gamma;
             Vsuby = (sinxpy*cosxpy*sinxmy*sinxmy - sinxpy*sinxpy*sinxmy*cosxmy)/gamma;
elseif
             amrex::Real sinx,cosx,siny,cosy;

             sinx = std::sin(xloc);
             cosx = std::cos(xloc);
             siny = std::sin(yloc);
             cosy = std::cos(yloc);
             Vsubx = 2.*sinx*cosx*siny*siny/gamma;
             Vsuby = 2.*siny*cosy*sinx*sinx/gamma;
#endif

             fluxx(i,j,k,0) +=  phi(i,j,k,0)*det(i,j,k,1)*(gmet(i,j,k,0)*Vsubx + gmet(i,j,k,1)*Vsuby);

          } else {

             amrex::GpuArray<amrex::Real, 4>  Vsub_x;
             amrex::GpuArray<amrex::Real, 4>  Vsub_y;

             for(int ii = 0; ii < 2; ii++){
             for(int jj = 0; jj < 2; jj++){

                xloc = (i+ii-0.5) / dxinv;
                yloc = (j+jj-0.5) / dyinv;

/*
                caps = 0;
                capsy = 0;
                capr  = (delta-alpha*std::pow(std::sin(xloc),2))*(1.+std::pow(std::cos(xloc/2),6));
                caprx = -2.*alpha*(1.+std::pow(std::cos(xloc/2),6))*std::cos(xloc)*std::sin(xloc)
                     -3.*std::pow(std::cos(xloc/2),5)*std::sin(xloc/2)*(delta-alpha*std::pow(std::sin(xloc),2));


                for (int iterm  = 0; iterm < 3; iterm++)
                {
                   amrex::Real shift = PI/3. + 2.*iterm*PI/3.;
                   caps += std::pow(std::cos(0.5*(yloc-shift)),16);
                   capsy += -8*std::sin(0.5*(yloc-shift))*std::pow(std::cos(0.5*(yloc-shift)),15);
                }
**/

                capsy = -3*std::cos(yloc/2. - 0.41887902047863906*std::cos(xloc))*std::pow(std::sin(yloc/2. - 0.41887902047863906*std::cos(xloc)),5);
                caprx = -0.27890625*std::sin(3*xloc) - 0.05625*std::sin(4*xloc) - 0.01171875*std::sin(5*xloc) + 
                        std::sin(xloc)*(-0.9046875 - 2.175*std::cos(xloc) - 2.5132741228718345*std::cos(0.5*yloc - 0.41887902047863906*std::cos(xloc))*
                        std::pow(std::sin(0.5*yloc - 0.41887902047863906*std::cos(xloc)),5));

                Vsub_x[ii+2*jj] = caprx/gamma;
                Vsub_y[ii+2*jj] = capsy/gamma;

             }
             }

//             fluxx(i,j,k,0) +=  phi(i,j,k,0)*det(i,j,k,1)*(gmet(i,j,k,0)*Vsub_x[3] + gmet(i,j,k,1)*Vsub_y[3]);
//             fluxx(i,j,k,1) +=  phi(i-1,j-1,k,0)*det(i-1,j-1,k,1)*(gmet(i-1,j-1,k,0)*Vsub_x[0] + gmet(i-1,j-1,k,1)*Vsub_y[0]);
//             fluxx(i,j,k,2) +=  phi(i,j-1,k,0)*det(i,j-1,k,1)*(gmet(i,j-1,k,0)*Vsub_x[1] + gmet(i,j-1,k,1)*Vsub_y[1]);
//             fluxx(i,j,k,3) +=  phi(i-1,j,k,0)*det(i-1,j,k,1)*(gmet(i-1,j,k,0)*Vsub_x[2] + gmet(i-1,j,k,1)*Vsub_y[2]);

             fluxx(i,j,k,0) +=  phi(i,j,k,0)*det(i,j,k,1)*(gmet(i,j,k,0)*Vsub_x[3] + gmet(i,j,k,1)*Vsub_y[3]);
             fluxx(i,j,k,1) +=  phi(i-1,j,k,0)*det(i-1,j,k,1)*(gmet(i-1,j,k,0)*Vsub_x[2] + gmet(i-1,j,k,1)*Vsub_y[2]);
             fluxx(i,j,k,2) +=  phi(i,j,k,0)*det(i,j,k,1)*(gmet(i,j,k,0)*Vsub_x[3] + gmet(i,j,k,1)*Vsub_y[3]);
             fluxx(i,j,k,3) +=  phi(i-1,j,k,0)*det(i-1,j,k,1)*(gmet(i-1,j,k,0)*Vsub_x[2] + gmet(i-1,j,k,1)*Vsub_y[2]);
         
          }

          }

        if(Ncomp == 2){
           fluxx(i,j,k,1) = fluxx(i,j,k,0);
        }

//          amrex::Print() << "computing flux " << i << " " << j << " " << fluxx(i,j,k,0) << std::endl;
          if(num_flux == 1){

          amrex::Real phi_loc = phi(i,j,k,0);
          amrex::Real sqrphi = std::sqrt(std::max(phi_loc,0.));
          fluxx(i,j,k,0) +=  std::sqrt(2./num_part)*det(i,j,k,2)*sqrphi*(sqrgmet(i,j,k,0)*stochfluxx(i,j,k,0)+sqrgmet(i,j,k,1)*stochfluxy(i,j,k,0));

          } else if(num_flux == 4){

 //            phi_loc = phi(i-1,j-1,k,0);
 //            sqrphi = std::sqrt(std::max(phi_loc,0.));
 //            fluxx(i,j,k,1) +=  std::sqrt(2./num_part)*det(i-1,j-1,k,2)*sqrphi*(sqrgmet(i-1,j-1,k,0)*stochfluxx(i-1,j-1,k,1)+sqrgmet(i-1,j-1,k,1)*stochfluxy(i-1,j-1,k,1));

 //            phi_loc = phi(i-1,j,k,0);
 //            sqrphi = std::sqrt(std::max(phi_loc,0.));
 //            fluxx(i,j,k,2) +=  std::sqrt(2./num_part)*det(i-1,j,k,2)*sqrphi*(sqrgmet(i-1,j,k,0)*stochfluxx(i-1,j,k,2)+sqrgmet(i-1,j,k,1)*stochfluxy(i-1,j,k,2));

 //            phi_loc = phi(i,j-1,k,0);
 //            sqrphi = std::sqrt(std::max(phi_loc,0.));
 //            fluxx(i,j,k,3) +=  std::sqrt(2./num_part)*det(i,j-1,k,2)*sqrphi*(sqrgmet(i,j-1,k,0)*stochfluxx(i,j-1,k,3)+sqrgmet(i,j-1,k,1)*stochfluxy(i,j-1,k,3));

             amrex::Real phi_loc = phi(i,j,k,0);
             amrex::Real sqrphi = std::sqrt(std::max(phi_loc,0.));
             fluxx(i,j,k,0) += 2.* std::sqrt(2./num_part)*det(i,j,k,2)*sqrphi*(sqrgmet(i,j,k,0)*stochfluxx(i,j,k,0)+sqrgmet(i,j,k,1)*stochfluxy(i,j,k,0));

             phi_loc = phi(i-1,j,k,0);
             sqrphi = std::sqrt(std::max(phi_loc,0.));
             fluxx(i,j,k,1) +=  2.*std::sqrt(2./num_part)*det(i-1,j,k,2)*sqrphi*(sqrgmet(i-1,j,k,0)*stochfluxx(i-1,j,k,1)+sqrgmet(i-1,j,k,1)*stochfluxy(i-1,j,k,1));

             phi_loc = phi(i-1,j,k,0);
             sqrphi = std::sqrt(std::max(phi_loc,0.));
             fluxx(i,j,k,2) += 2.* std::sqrt(2./num_part)*det(i-1,j,k,2)*sqrphi*(sqrgmet(i-1,j,k,0)*stochfluxx(i-1,j,k,2)+sqrgmet(i-1,j,k,1)*stochfluxy(i-1,j,k,2));

             phi_loc = phi(i,j,k,0);
             sqrphi = std::sqrt(std::max(phi_loc,0.));
             fluxx(i,j,k,3) += 2.* std::sqrt(2./num_part)*det(i,j,k,2)*sqrphi*(sqrgmet(i,j,k,0)*stochfluxx(i,j,k,3)+sqrgmet(i,j,k,1)*stochfluxy(i,j,k,3));
          }
//          amrex::Print() << " x flux " << i << " " << j << " " << fluxx(i,j,k,0) << std::endl;
}



AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_y (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxy,
                     amrex::Array4<amrex::Real> const& stochfluxx,
                     amrex::Array4<amrex::Real> const& stochfluxy,
                     amrex::Array4<amrex::Real const> const& phi,
                     amrex::Array4<amrex::Real const> const& gmet,
                     amrex::Array4<amrex::Real const> const& sqrgmet,
                     amrex::Array4<amrex::Real const> const& det,
                     amrex::Real num_part,
                     amrex::Real dxinv,
                     amrex::Real dyinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi,int Ncomp, int num_flux, int ext_pot)
{
    constexpr double  TWOPI = 2.*3.14159265358979323846264338327950288;
    constexpr double  PI = 3.14159265358979323846264338327950288;

    int coeff_index = Ncomp-1;
    amrex::Real phibdry = 0.;
    amrex::Real stofact = 1.;

//    int ext_pot = 0;
#if 0
    amrex::Real alpha = .25;
    amrex::Real beta = .75;
    amrex::Real gamma = 5.e-4;
    gamma = 2.e-1;
#else
    amrex::Real alpha = .3;
    amrex::Real beta = 0.;
    amrex::Real delta = 2.;
    amrex::Real gamma = 1.e-1;
#endif

    amrex::Real phip,phim, phipn, phimn;

/*
    if (lo == dom_lo &&
        (bc_lo == amrex::BCType::foextrap ||
         bc_lo == amrex::BCType::ext_dir))
    {
        if(i == lo)
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv / 0.5;
        }
        else
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv;
        }
    }
    else if (hi == dom_hi &&
             (bc_hi == amrex::BCType::foextrap ||
              bc_hi == amrex::BCType::ext_dir))
    {
        if(i == hi+1)
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv / 0.5;
        }
        else
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv;
        }
    }
*/
    //amrex::Print() << "BC info " << i << " " << dom_lo << " " << bc_lo << " " << dom_hi << " " << bc_hi << std::endl;
#if 0
    if (i == dom_lo && (bc_lo == amrex::BCType::foextrap || bc_lo == amrex::BCType::ext_dir)) {

        stofact = std::sqrt(2.);

        phim = 2.*phibdry - phi(i,j,k,0);
        phip = phi(i,j,k,0);

        //phimn = 2.*phibdry - phi(i,j,k,coeff_index);
        //phipn = phi(i,j,k,coeff_index);
        phimn = phibdry;
        phipn = phibdry;
        //amrex::Print() << " A " << std::endl;

    } else if (i == dom_hi+1 && (bc_hi == amrex::BCType::foextrap || bc_hi == amrex::BCType::ext_dir)) {

        stofact = std::sqrt(2.);
        phip = 2.*phibdry - phi(i-1,j,k,0);
        phim = phi(i-1,j,k,0);
        //phipn = 2.*phibdry - phi(i-1,j,k,coeff_index);
        //phimn = phi(i-1,j,k,coeff_index);
        phimn = phibdry;
        phipn = phibdry;
        //amrex::Print() << " B " << std::endl;
    } else {

        phip = phi(i,j,k,0);
        phim = phi(i-1,j,k,0);
        phipn = std::abs(phi(i,j,k,coeff_index));
        phimn = std::abs(phi(i-1,j,k,coeff_index));
        stofact = 1.;
        //amrex::Print() << " C " << std::endl;

    }
#endif

//	amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i-1,j,k))/(phi(i,j,k)+phi(i-1,j,k)+1.e-16);
	//amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
	//amrex::Real phiavg = 0.5*(phip+phim);

	// HACK TO COMMENT OUT STOCHASTIC TERMS
	// phiavg = 0.;

        //fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv
    	//                    + stochfluxx(i,j,k) * std::sqrt(phiavg);

//          if(ext_pot == 1) {
//
//             amrex::Real yloc;
//             yloc = j / dyinv;
//	     amrex::Real phiavg = 0.5*(phip+phim);

 //            amrex::Real Vsuby;

  //           Vsuby =0.5* 4.*(yloc - .5)*(yloc - .5)*(yloc - .5) / gamma;

   //          fluxy(i,j,k,0) += phiavg * Vsuby;

    //      }
          //amrex::Print() << " point " << i << " " << j << " " <<  dxinv << std::endl;
          //amrex::Print() << "phip " << phip << std::endl;
          //amrex::Print() << "phib " << phibdry << std::endl;
          //amrex::Print() << "phi " << phi(i,j,k,0) << std::endl;
          //amrex::Print() << "phim " << phim << std::endl;

          amrex::Real diffphixm =( phi(i,j,k,0) - phi(i-1,j,k,0))*dxinv;
          amrex::Real diffphiym =( phi(i,j,k,0) - phi(i,j-1,k,0))*dyinv;

          fluxy(i,j,k,0) = det(i,j,k,1)*(gmet(i,j,k,2)*diffphixm + gmet(i,j,k,3)*diffphiym);

          if(num_flux ==4){
             // Four one-sided gradients on the 2x2 stencil around the face (i,j)
             // -- at (i,j)
/*             amrex::Real gx_mm = ( phi(i,j,k,0)   - phi(i-1,j,k,0) )*dxinv;
             amrex::Real gy_mm = ( phi(i,j,k,0)   - phi(i,j-1,k,0) )*dyinv;
             // ++ at (i-1,j-1)
             amrex::Real gx_pp = ( phi(i,j-1,k,0) - phi(i-1,j-1,k,0) )*dxinv;
             amrex::Real gy_pp = ( phi(i-1,j,k,0) - phi(i-1,j-1,k,0) )*dyinv;
             // +- at (i-1,j)
             amrex::Real gx_pm = ( phi(i,j,k,0)   - phi(i-1,j,k,0) )*dxinv;
             amrex::Real gy_pm = ( phi(i-1,j,k,0) - phi(i-1,j-1,k,0) )*dyinv;
             // -+ at (i,j-1)
             amrex::Real gx_mp = ( phi(i,j-1,k,0) - phi(i-1,j-1,k,0) )*dxinv;
             amrex::Real gy_mp = ( phi(i,j,k,0)   - phi(i,j-1,k,0) )*dyinv;

             fluxy(i,j,k,0) = det(i,j,k,1)*(gmet(i,j,k,2)*gx_mm + gmet(i,j,k,3)*gy_mm);
             fluxy(i,j,k,1) = det(i-1,j-1,k,1)*(gmet(i-1,j-1,k,2)*gx_pp + gmet(i-1,j-1,k,3)*gy_pp);
             fluxy(i,j,k,2) = det(i-1,j,k,1)*(gmet(i-1,j,k,2)*gx_pm + gmet(i-1,j,k,3)*gy_pm);
             fluxy(i,j,k,3) = det(i,j-1,k,1)*(gmet(i,j-1,k,2)*gx_mp + gmet(i,j-1,k,3)*gy_mp);
*/
             amrex::Real gx_mm = ( phi(i,j,k,0)   - phi(i-1,j,k,0) )*dxinv;
             amrex::Real gy_mm = ( phi(i,j,k,0)   - phi(i,j-1,k,0) )*dyinv;
             // ++ at (i-1,j-1)
             amrex::Real gx_pp = ( phi(i+1,j-1,k,0) - phi(i,j-1,k,0) )*dxinv;
             amrex::Real gy_pp = gy_mm;
             // +- at (i-1,j)
             amrex::Real gx_pm = ( phi(i,j-1,k,0)   - phi(i-1,j-1,k,0) )*dxinv;
             amrex::Real gy_pm = gy_mm;
             // -+ at (i,j-1)
             amrex::Real gx_mp = ( phi(i+1,j,k,0) - phi(i,j,k,0) )*dxinv;
             amrex::Real gy_mp = gy_mm;

             fluxy(i,j,k,0) = det(i,j,k,1)*(gmet(i,j,k,2)*gx_mm + gmet(i,j,k,3)*gy_mm);
             fluxy(i,j,k,1) = det(i,j-1,k,1)*(gmet(i,j-1,k,2)*gx_pp + gmet(i,j-1,k,3)*gy_pp);
             fluxy(i,j,k,2) = det(i,j-1,k,1)*(gmet(i,j-1,k,2)*gx_pm + gmet(i,j-1,k,3)*gy_pm);
             fluxy(i,j,k,3) = det(i,j,k,1)*(gmet(i,j,k,2)*gx_mp + gmet(i,j,k,3)*gy_mp);

          }

          if(ext_pot == 1) {

             amrex::Real xloc,yloc;


             amrex::Real Vsubx;
             amrex::Real Vsuby;

#if 1

             amrex::Real caps, capsy, capr, caprx;

          if(num_flux == 1){

             xloc = (i+0.5) / dxinv;
             yloc = (j+0.5) / dyinv;

             caps = 0;
             capsy = 0;
             capr  = (delta-alpha*std::pow(std::sin(xloc),2))*(1.+std::pow(std::cos(xloc/2),6));
             caprx = -2.*alpha*(1.+std::pow(std::cos(xloc/2),6))*std::cos(xloc)*std::sin(xloc)
                     -3.*std::pow(std::cos(xloc/2),5)*std::sin(xloc/2)*(delta-alpha*std::pow(std::sin(xloc),2));


             for (int iterm  = 0; iterm < 3; iterm++)
             {
                amrex::Real shift = PI/3. + 2.*iterm*PI/3.;
                caps += std::pow(std::cos(0.5*(yloc-shift)),16);
                capsy += -8*std::sin(0.5*(yloc-shift))*std::pow(std::cos(0.5*(yloc-shift)),15);
             }


             Vsubx = caprx/gamma;
             Vsuby = capsy/gamma;

#else

             amrex::Real sinxpy = std::sin(0.5*(xloc+yloc));
             amrex::Real cosxpy = std::cos(0.5*(xloc+yloc));
             amrex::Real sinxmy = std::sin(0.5*(xloc-yloc));
             amrex::Real cosxmy = std::cos(0.5*(xloc-yloc));

             Vsubx = (sinxpy*cosxpy*sinxmy*sinxmy + sinxpy*sinxpy*sinxmy*cosxmy)/gamma;
             Vsuby = (sinxpy*cosxpy*sinxmy*sinxmy - sinxpy*sinxpy*sinxmy*cosxmy)/gamma;
elseif
             amrex::Real sinx,cosx,siny,cosy;

             sinx = std::sin(xloc);
             cosx = std::cos(xloc);
             siny = std::sin(yloc);
             cosy = std::cos(yloc);
             Vsubx = 2.*sinx*cosx*siny*siny/gamma;
             Vsuby = 2.*siny*cosy*sinx*sinx/gamma;
#endif


             fluxy(i,j,k,0) +=  phi(i,j,k,0)*det(i,j,k,1)*(gmet(i,j,k,2)*Vsubx + gmet(i,j,k,3)*Vsuby);

          } else if (num_flux == 4){

             amrex::GpuArray<amrex::Real, 4>  Vsub_x;
             amrex::GpuArray<amrex::Real, 4>  Vsub_y;

             for(int ii = 0; ii < 2; ii++){
             for(int jj = 0; jj < 2; jj++){

                xloc = (i+ii-0.5) / dxinv;
                yloc = (j+jj-0.5) / dyinv;
/*
                caps = 0;
                capsy = 0;
                capr  = (delta-alpha*std::pow(std::sin(xloc),2))*(1.+std::pow(std::cos(xloc/2),6));
                caprx = -2.*alpha*(1.+std::pow(std::cos(xloc/2),6))*std::cos(xloc)*std::sin(xloc)
                     -3.*std::pow(std::cos(xloc/2),5)*std::sin(xloc/2)*(delta-alpha*std::pow(std::sin(xloc),2));


                for (int iterm  = 0; iterm < 3; iterm++)
                {
                   amrex::Real shift = PI/3. + 2.*iterm*PI/3.;
                   caps += std::pow(std::cos(0.5*(yloc-shift)),16);
                   capsy += -8*std::sin(0.5*(yloc-shift))*std::pow(std::cos(0.5*(yloc-shift)),15);
                }

*/

                capsy = -3*std::cos(yloc/2. - 0.41887902047863906*std::cos(xloc))*std::pow(std::sin(yloc/2. - 0.41887902047863906*std::cos(xloc)),5);
                caprx = -0.27890625*std::sin(3*xloc) - 0.05625*std::sin(4*xloc) - 0.01171875*std::sin(5*xloc) + 
                        std::sin(xloc)*(-0.9046875 - 2.175*std::cos(xloc) - 2.5132741228718345*std::cos(0.5*yloc - 0.41887902047863906*std::cos(xloc))*
                        std::pow(std::sin(0.5*yloc - 0.41887902047863906*std::cos(xloc)),5));

                Vsub_x[ii+2*jj] = caprx/gamma;
                Vsub_y[ii+2*jj] = capsy/gamma;

             }
             }

             fluxy(i,j,k,0) +=  phi(i,j,k,0)*det(i,j,k,1)*(gmet(i,j,k,2)*Vsub_x[3] + gmet(i,j,k,3)*Vsub_y[3]);
             fluxy(i,j,k,1) +=  phi(i,j-1,k,0)*det(i,j-1,k,1)*(gmet(i,j-1,k,2)*Vsub_x[1] + gmet(i,j-1,k,3)*Vsub_y[1]);
             fluxy(i,j,k,2) +=  phi(i,j-1,k,0)*det(i,j-1,k,1)*(gmet(i,j-1,k,2)*Vsub_x[1] + gmet(i,j-1,k,3)*Vsub_y[1]);
             fluxy(i,j,k,3) +=  phi(i,j,k,0)*det(i,j,k,1)*(gmet(i,j,k,2)*Vsub_x[3] + gmet(i,j,k,3)*Vsub_y[3]);

          }

          }

        if(Ncomp == 2){
           fluxy(i,j,k,1) = fluxy(i,j,k,0);
        }

          if(num_flux == 1){

          amrex::Real phi_loc = phi(i,j,k,0);
          amrex::Real sqrphi = std::sqrt(std::max(phi_loc,0.));
          fluxy(i,j,k,0) +=  std::sqrt(2./num_part)*det(i,j,k,2)*sqrphi*(sqrgmet(i,j,k,2)*stochfluxx(i,j,k,0)+sqrgmet(i,j,k,3)*stochfluxy(i,j,k,0));

          } else if(num_flux == 4){

//.             phi_loc = phi(i-1,j-1,k,0);
//.             sqrphi = std::sqrt(std::max(phi_loc,0.));
//.             fluxy(i,j,k,1) +=  std::sqrt(2./num_part)*det(i-1,j-1,k,2)*sqrphi*(sqrgmet(i-1,j-1,k,2)*stochfluxx(i-1,j-1,k,1)+sqrgmet(i-1,j-1,k,3)*stochfluxy(i-1,j-1,k,1));

//.             phi_loc = phi(i-1,j,k,0);
//.             sqrphi = std::sqrt(std::max(phi_loc,0.));
//.             fluxy(i,j,k,2) +=  std::sqrt(2./num_part)*det(i-1,j,k,2)*sqrphi*(sqrgmet(i-1,j,k,2)*stochfluxx(i-1,j,k,2)+sqrgmet(i-1,j,k,3)*stochfluxy(i-1,j,k,2));

//.             phi_loc = phi(i,j-1,k,0);
//.             sqrphi = std::sqrt(std::max(phi_loc,0.));
//.             fluxy(i,j,k,3) +=  std::sqrt(2./num_part)*det(i,j-1,k,2)*sqrphi*(sqrgmet(i,j-1,k,2)*stochfluxx(i,j-1,k,3)+sqrgmet(i,j-1,k,3)*stochfluxy(i,j-1,k,3));

             amrex::Real phi_loc = phi(i,j,k,0);
             amrex::Real sqrphi = std::sqrt(std::max(phi_loc,0.));
             fluxy(i,j,k,0) += 2.* std::sqrt(2./num_part)*det(i,j,k,2)*sqrphi*(sqrgmet(i,j,k,2)*stochfluxx(i,j,k,0)+sqrgmet(i,j,k,3)*stochfluxy(i,j,k,0));

             phi_loc = phi(i,j-1,k,0);
             sqrphi = std::sqrt(std::max(phi_loc,0.));
             fluxy(i,j,k,1) += 2.* std::sqrt(2./num_part)*det(i,j-1,k,2)*sqrphi*(sqrgmet(i,j-1,k,2)*stochfluxx(i,j-1,k,1)+sqrgmet(i,j-1,k,3)*stochfluxy(i,j-1,k,1));

             phi_loc = phi(i,j,k,0);
             sqrphi = std::sqrt(std::max(phi_loc,0.));
             fluxy(i,j,k,2) += 2.* std::sqrt(2./num_part)*det(i,j,k,2)*sqrphi*(sqrgmet(i,j,k,2)*stochfluxx(i,j,k,2)+sqrgmet(i,j,k,3)*stochfluxy(i,j,k,2));

             phi_loc = phi(i,j-1,k,0);
             sqrphi = std::sqrt(std::max(phi_loc,0.));
             fluxy(i,j,k,3) += 2.* std::sqrt(2./num_part)*det(i,j-1,k,2)*sqrphi*(sqrgmet(i,j-1,k,2)*stochfluxx(i,j-1,k,3)+sqrgmet(i,j-1,k,3)*stochfluxy(i,j-1,k,3));

          }
//          amrex::Print() << " x flux " << i << " " << j << " " << fluxx(i,j,k,0) << std::endl;
}

#if 0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_y (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxy,
                     amrex::Array4<amrex::Real> const& stochfluxy,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dyinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi, int Ncomp)
{
    int coeff_index = Ncomp-1;
    amrex::Real phibdry = 0.;
    amrex::Real stofact = 1.;

    int ext_pot = 1;
    amrex::Real alpha = .25;
    amrex::Real beta = .75;
    amrex::Real gamma = 5.e-4;
    gamma = 1.1e-3;

    amrex::Real phip,phim, phipn, phimn;

    if (j == dom_lo && (bc_lo == amrex::BCType::foextrap || bc_lo == amrex::BCType::ext_dir))
    {

        stofact = std::sqrt(2.);
        phim = 2.*phibdry - phi(i,j,k,0);
        phip = phi(i,j,k,0);
        //phimn = 2.*phibdry - std::abs(phi(i,j,k,coeff_index));
        //phipn = std::abs( phi(i,j,k,coeff_index));
        phipn = phibdry;
        phimn = phibdry;

    } else if (j == dom_hi+1 && (bc_hi == amrex::BCType::foextrap || bc_hi == amrex::BCType::ext_dir)) {

        stofact = std::sqrt(2.);
        phip = 2.*phibdry - phi(i,j-1,k,0);
        phim = phi(i,j-1,k,0);
        //phipn = 2.*phibdry - std::abs(phi(i,j-1,k,coeff_index));
        //phimn = std::abs(phi(i,j-1,k,coeff_index));
        phipn = phibdry;
        phimn = phibdry;

     } else {

        phip = phi(i,j,k,0);
        phim = phi(i,j-1,k,0);
        phipn = std::abs(phi(i,j,k,coeff_index));
        phimn = std::abs(phi(i,j-1,k,coeff_index));
        stofact = 1.;

    }
//	amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i,j-1,k))/(phi(i,j,k)+phi(i,j-1,k)+1.e-16);
	//amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
	//amrex::Real phiavg = 0.5*(phip+phim);

	// HACK TO COMMENT OUT STOCHASTIC TERMS
	// phiavg = 0.;

        //fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv
	//+ stochfluxy(i,j,k)* std::sqrt(phiavg);

        fluxy(i,j,k,0) = 0.5*(phip-phim) * dyinv;

          if(ext_pot == 1) {

             amrex::Real yloc;
             yloc = j / dyinv;
	     amrex::Real phiavg = 0.5*(phip+phim);

             amrex::Real Vsuby;

             Vsuby =0.5* 4.*(yloc - .5)*(yloc - .5)*(yloc - .5) / gamma;

             fluxy(i,j,k,0) += phiavg * Vsuby;

          }

        if(Ncomp == 2){
           fluxy(i,j,k,1) = fluxy(i,j,k,0);
        }

	amrex::Real sqphiavg = 0.5*(std::sqrt(phipn)+std::sqrt(phimn));

	fluxy(i,j,k,0) += stofact * stochfluxy(i,j,k)* sqphiavg;
//        amrex::Print() << " y flux " << i << " " << j << " " << fluxy(i,j,k,0) << std::endl;
}


#endif
#if (AMREX_SPACEDIM > 2)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_z (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxz,
                     amrex::Array4<amrex::Real> const& stochfluxz,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dzinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi, int Ncomp)
{
    int coeff_index = Ncomp-1;
    if (lo == dom_lo &&
        (bc_lo == amrex::BCType::foextrap ||
         bc_lo == amrex::BCType::ext_dir))
    {
        if(k == lo)
        {
            fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv / 0.5;
        }
        else
        {
            fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv;
        }
    }
    else if (hi == dom_hi &&
             (bc_hi == amrex::BCType::foextrap ||
              bc_hi == amrex::BCType::ext_dir))
    {
        if(k == hi+1)
        {
            fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv / 0.5;
        }
        else
        {
            fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv;
        }
    }
    else
    {
//	amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i,j,k-1))/(phi(i,j,k)+phi(i,j,k-1)+1.e-16);
        amrex::Real phip = std::max(phi(i,j,k,coeff_index),0.);
        amrex::Real phim = std::max(phi(i,j,k-1,coeff_index),0.);
	//amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
        //fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv
	// +stochfluxz(i,j,k) * std::sqrt(phiavg);
	amrex::Real sqphiavg = 0.5*(std::sqrt(phip)+std::sqrt(phim));
        fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv
	 +stochfluxz(i,j,k) * sqphiavg;
        if(Ncomp == 2){
           fluxz(i,j,k,1) = 0.5*(phi(i,j,k,1)-phi(i,j,k-1,1)) * dzinv;
        }
    }
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void update_phi (int i, int j, int k,
                 amrex::Array4<amrex::Real const> const& phiold,
                 amrex::Array4<amrex::Real      > const& phinew,
                 amrex::Array4<amrex::Real const> const& det,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fluxx,
                              amrex::Array4<amrex::Real const> const& fluxy,
                              amrex::Array4<amrex::Real const> const& fluxz),
                 amrex::Real dt,
                 AMREX_D_DECL(amrex::Real dxinv,
                              amrex::Real dyinv,
                              amrex::Real dzinv),
                 int Ncomp,
                 int num_flux)
{
 // amrex::Print() << "Ncomp " << Ncomp << " " << dt << " " << dxinv << " " << dyinv << std::endl;
 // amrex::Print() << " i,j " << i << " " << j << " " << fluxx(i+1,j,k,0) << " " << fluxx(i,j,k,0) << " " << fluxy(i,j+1,k,0)<< " " << fluxy(i,j,k,0) << std::endl;
 // amrex::Print() << "det " << det(i,j,k,1) << std::endl;
    if(num_flux == 1){

       for (int n = 0; n < Ncomp; n++){
            phinew(i,j,k,n) = phiold(i,j,k,n)
                + ( dt * dxinv * (fluxx(i+1,j  ,k,n  ) - fluxx(i,j,k,n))
                +  dt * dyinv * (fluxy(i  ,j+1,k,n  ) - fluxy(i,j,k,n))
#if (AMREX_SPACEDIM > 2)
                +  dt * dzinv * (fluxz(i  ,j  ,k+1,n) - fluxz(i,j,k,n)))/det(i,j,k,1);
#else
                )/det(i,j,k,1);
#endif
       }
    } else if (num_flux == 4 ){

       phinew(i,j,k,0) = phiold(i,j,k,0)
           + 0.25* dt* (  dxinv * (fluxx(i+1,j  ,k,0  ) - fluxx(i,j,k,0))
                       +  dxinv * (fluxx(i+1,j  ,k,1  ) - fluxx(i,j,k,1))
                       +  dxinv * (fluxx(i+1,j  ,k,2  ) - fluxx(i,j,k,2))
                       +  dxinv * (fluxx(i+1,j  ,k,3  ) - fluxx(i,j,k,3))
                       +  dyinv * (fluxy(i  ,j+1,k,0  ) - fluxy(i,j,k,0))
                       +  dyinv * (fluxy(i  ,j+1,k,1  ) - fluxy(i,j,k,1))
                       +  dyinv * (fluxy(i  ,j+1,k,2  ) - fluxy(i,j,k,2))
                       +  dyinv * (fluxy(i  ,j+1,k,3  ) - fluxy(i,j,k,3)))/det(i,j,k,1);


    }
}

#endif
