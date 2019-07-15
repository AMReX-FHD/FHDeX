#include "common_functions.H"
#include "common_functions_F.H"

#include "common_namespace.H"

using namespace common;

int greatest_common_divisor(int,int);
void factor(int,int*,int);

//Computes divergence at cell centres from velcocities at cell faces
void ComputeVerticalAverage(const MultiFab& mf, MultiFab& mf_avg, 
			    const Geometry& geom, const int dir, 
			    const int incomp, const int outcomp, const int ncomp)
{

  MultiFab mf_flattened, mf_pencil;

  BoxArray ba_in = mf.boxArray();
  BoxArray ba_flattened;
  BoxArray ba_pencil;
  BoxArray ba_flat;

  int nbox = ba_in.size();
  Box domain(geom.Domain());
  const DistributionMapping& dmap = mf.DistributionMap();

  Vector<int> max_grid_size_flattened(AMREX_SPACEDIM);
  Vector<int> max_grid_size_pencil(AMREX_SPACEDIM);

  int nx[3], nbx[3];
  int mx[2], mbx[2];

  int nxprod, indlo, indhi, a, b, c;

  // We assume that all grids have the same size hence 
  // we have the same nx,ny,nz on all ranks
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    nx[d]  = ba_in[0].size()[d];
    nbx[d] = domain.length(d) / nx[d];
  }
  if (nbx[0]*nbx[1]*nbx[2] != ba_in.size()) 
    amrex::Error("ALL GRIDS DO NOT HAVE SAME SIZE");

  nxprod = nx[0]*nx[1]*nx[2]/nx[dir];
  if (nxprod%nbx[dir] != 0) {
    amrex::Error("CURRENT PENCIL REFACTORING DOESN'T WORK");
  } else {
    nxprod /= nbx[dir];
  }

  indlo = (dir-1+AMREX_SPACEDIM)%AMREX_SPACEDIM;
  indhi = (dir+1+AMREX_SPACEDIM)%AMREX_SPACEDIM;

  IntVect dom_lo(domain.loVect3d());
  IntVect dom_hi(domain.hiVect3d());
  dom_hi[dir] = nbx[dir];
  Box domain_flattened(dom_lo, dom_hi);
  dom_hi[dir] = 1;
  Box domain_flat(dom_lo, dom_hi);

  a = greatest_common_divisor( nxprod,domain.length(indlo) );
  b = greatest_common_divisor( nxprod,domain.length(indhi) );
  c = (a*b)/nxprod;
  factor(c, mx, 2);
  mx[0] = a/mx[0];
  mx[1] = b/mx[1];

  // Print() << "Hack: " << a << "," << b << "," << c << std::endl;
  // Print() << "Hack: " << nxprod << ", grids: " << mx[0] << "x" << mx[1] << std::endl;

  if (mx[0]*mx[1] != nxprod) 
    amrex::Error("FACTORING NOT POSSIBLE DUE TO UNCOMMON PRIME FACTOR");

  mbx[0] = domain.length(indlo)/mx[0];
  mbx[1] = domain.length(indhi)/mx[1];

  Print() << "2D redist: " << mbx[0] << "x" << mbx[1] << ", grids: " << mx[0] << "x" << mx[1] << std::endl;

  max_grid_size_pencil[indlo] = mx[0];
  max_grid_size_pencil[indhi] = mx[1];
  max_grid_size_pencil[dir]   = 1;      // nx[dir]

  max_grid_size_flattened      = max_grid_size;
  max_grid_size_flattened[dir] = 1;      // nx[dir]

  ba_flattened.define(domain_flattened);
  ba_flattened.maxSize(IntVect(max_grid_size_flattened));
  mf_flattened.define(ba_flattened,dmap,ncomp,0);

  ba_pencil.define(domain_flattened);
  ba_pencil.maxSize(IntVect(max_grid_size_pencil));
  mf_pencil.define(ba_pencil,dmap,ncomp,0);

  ba_flat.define(domain_flat);
  ba_flat.maxSize(IntVect(max_grid_size_pencil));
  mf_avg.define(ba_flat,dmap,ncomp,0);

  // exit(0);

  for ( MFIter mfi(mf); mfi.isValid(); ++mfi ) {
    const Box& bx = mfi.validbox();
    compute_vert_average(BL_TO_FORTRAN_FAB(mf[mfi]),BL_TO_FORTRAN_FAB(mf_flattened[mfi]), 
			 &dir,
			 &incomp, &outcomp, &ncomp);
  }

  // Copy/redistrubute to pencils
  // ...

  // for ( MFIter mfi(mf); mfi.isValid(); ++mfi ) {
  //   const Box& bx = mfi.validbox();
  //   compute_vert_average(BL_TO_FORTRAN_FAB(mf_pencil[mfi]),BL_TO_FORTRAN_FAB(mf_avg[mfi]), 
  // 			 &dir,
  // 			 &incomp, &outcomp, &ncomp);
  // }

}

int greatest_common_divisor(int a, int b) {
  return b == 0 ? a : greatest_common_divisor(b, a % b);
}

void factor(int num, int* factors, int nf) {
  
  int n = num;
  int cnt = 0;
  
  for (int i=0; i<nf; i++)  
    factors[i]=1;

  // factor out powers of 2 
  while (n%2 == 0)  
    {  
      // cout << 2 << " ";  
      n /= 2; 
      factors[cnt%nf]*=2;
      cnt++;
    }  
  for (int i=3; i<=sqrt(n); i+=2)  
    {  
      // find other odd factors
      while (n % i == 0)  
	{  
	  // cout << i << " ";  
	  n /= i;  
	  factors[cnt%nf]*=i;
	  cnt++;
	}  
    }  
  
  // check if n is a prime number greater than 2  
  if (n > 2)  {
    // cout << n << " ";  
    amrex::Error("CANNOT FACTOR PRIME NUMBER");
  }
  
  n = 1;
  for (int i=0; i<nf; i++)  
    n *= factors[i];

  if (n != num)  {
    amrex::Error("ERROR IN FACTORING");
  }
}  
