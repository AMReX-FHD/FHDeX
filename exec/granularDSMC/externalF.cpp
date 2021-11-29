#include "DsmcParticleContainer.H"
#include <common_namespace.H>
#include <math.h>

void FhdParticleContainer::externalForce(Real dt)
{
  int lev = 0;

  const Real* dx = Geom(lev).CellSize();
  const Real* plo = Geom(lev).ProbLo();
  const Real* phi = Geom(lev).ProbHi();

	// drag = -mu*(v_p - U_gas) + mass*stochastic_term
	// mu = mu0*R(\phi) where mu_0 is solvent viscosity
	// mu_0 = 3*pi*mu_air*dp
	// stochastic_term = sqrt(2*mu*Text)
	Real Fb, Ff, Fs; // mean drag, drag due to friciton, fluid forces due neighbors
	
	Real rho;

	Real R_diss, mu_air;
	Real sphi = 0.;
	R_diss = 1. + 3.*sqrt(sphi*0.5);
	
	Real eta; // for stochastic term
	eta = sqrt(2.*mu_air*tbath); // from above, needs a sqrt(dp)

	Real Ugx, Ugy, Ugz; // need as inputs?
	Ugx = 0.;
	Ugy = 0.; // for length 2L, v = - (x-L)^2 + L^2 = 2xL - x^2
	Ugz = 0.;

  for (FhdParIter pti(* this, lev); pti.isValid(); ++pti)
  {
    const int grid_id = pti.index();
    const int tile_id = pti.LocalTileIndex();
    const Box& tile_box  = pti.tilebox();

    auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
    auto& particles = particle_tile.GetArrayOfStructs();
    const long np = particles.numParticles();

    for (int i = 0; i < np; ++ i)
    {
      ParticleType & p = particles[i];
      Real mass = p.rdata(FHD_realData::mass);
			Real u    = p.rdata(FHD_realData::velx);
			Real v    = p.rdata(FHD_realData::vely);
			Real w    = p.rdata(FHD_realData::velz);
			Real dp   = (p.rdata(FHD_realData::radius)*2.);
			Real mu_p = mu_air*dp;
			// Weiner Process
			Real fDragx = amrex::RandomNormal(0.,1.)*eta*sqrt(2.*eta*tbath/mass)*sqrt(dt)*sqrt(dp);
			Real fDragy = amrex::RandomNormal(0.,1.)*eta*sqrt(2.*eta*tbath/mass)*sqrt(dt)*sqrt(dp);
			Real fDragz = amrex::RandomNormal(0.,1.)*eta*sqrt(2.*eta*tbath/mass)*sqrt(dt)*sqrt(dp);
			Real mDragx = -mu_p*(u - Ugx);
			Real mDragy = -mu_p*(v - Ugy);
			Real mDragz = -mu_p*(w - Ugz);
			// TEMPORARY
			mDragx = 0.;
			mDragy = 0.;
			mDragz = 0.;
			fDragx = 0.;
			fDragy = 0.;
			fDragz = 0.;
			u = u + (grav[0]+mDragx)*dt+fDragx;
			v = v + (grav[1]+mDragy)*dt+fDragy;
			w = w + (grav[2]+mDragz)*dt+fDragz;
			p.rdata(FHD_realData::velx) = u;
			p.rdata(FHD_realData::vely) = v;
			p.rdata(FHD_realData::velz) = w;
    }
  }
}
