#include "common_functions.H"
#include "common_functions_F.H"

#include "common_namespace.H"

#include "AMReX_PlotFileUtil.H"

using namespace common;

int greatest_common_factor(int,int);
void factor(int,int*,int);

void ComputeVerticalAverage(const MultiFab& mf, MultiFab& mf_avg,
			    const Geometry& geom, const int dir,
			    const int incomp, const int ncomp,
			    const int findredist)
{
    BL_PROFILE_VAR("ComputVerticalAverage()",ComputeVerticalAverage);
  
    bool write_data = false;
    std::string plotname;

    int outcomp = 0;
    int inputcomp = 0;

    MultiFab mf_pencil;

    BoxArray ba_in = mf.boxArray();
    BoxArray ba_pencil;
    BoxArray ba_flat;

    int nbox = ba_in.size();
    Box domain(geom.Domain());
    const DistributionMapping& dmap = mf.DistributionMap();

    Vector<int> max_grid_size_pencil(AMREX_SPACEDIM);
    Vector<int> max_grid_size_flat(AMREX_SPACEDIM);

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

    indlo = (dir-1+AMREX_SPACEDIM)%AMREX_SPACEDIM;
    indhi = (dir+1+AMREX_SPACEDIM)%AMREX_SPACEDIM;

    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());
    dom_hi[dir] = 0;
    Box domain_flat(dom_lo, dom_hi);

    if (findredist == 1) {

        nxprod = nx[0]*nx[1]*nx[2]/nx[dir];
        if (nxprod%nbx[dir] != 0) {
            amrex::Error("CURRENT PENCIL REFACTORING DOESN'T WORK");
        } else {
            nxprod /= nbx[dir];
        }

        // Find a,b,&c such that (a*b)/c = (nx*ny)/pz, with c as a common factor to a & b
        a = greatest_common_factor( nxprod,domain.length(indlo) );
        b = greatest_common_factor( nxprod,domain.length(indhi) );
        c = (a*b)/nxprod; // c is a factor of both a & b
        factor(c, mx, 2); // factor c into two numbers
        mx[0] = a/mx[0];
        mx[1] = b/mx[1];

        if (mx[0]*mx[1] != nxprod)
            amrex::Error("FACTORING NOT POSSIBLE DUE TO UNCOMMON PRIME FACTOR");

    } else {

        mx[0] = max_grid_projection[0];
        mx[1] = max_grid_projection[1];

    }

    mbx[0] = domain.length(indlo)/mx[0];
    mbx[1] = domain.length(indhi)/mx[1];

    Print() << "2D redist: " << mbx[0] << "x" << mbx[1] << ", grids: " << mx[0]
            << "x" << mx[1] << std::endl;

    // Print() << "HACK: dmap = " << dmap << std::endl;

    max_grid_size_pencil[indlo] = mx[0];
    max_grid_size_pencil[indhi] = mx[1];
    max_grid_size_pencil[dir]   = domain.length(dir);
    ba_pencil.define(domain);
    ba_pencil.maxSize(IntVect(max_grid_size_pencil));
    DistributionMapping dmap_pencil(ba_pencil);
    mf_pencil.define(ba_pencil,dmap_pencil,ncomp,0);

    max_grid_size_flat      = max_grid_size_pencil;
    max_grid_size_flat[dir] = 1;
    ba_flat.define(domain_flat);
    ba_flat.maxSize(IntVect(max_grid_size_flat));
    mf_avg.define(ba_flat,dmap_pencil,ncomp,0);

    // copy/redistrubute to pencils

    mf_pencil.ParallelCopy(mf, incomp, 0, ncomp);

    if (write_data) {
        plotname = "mf_pencil";
        VisMF::Write(mf_pencil,plotname);
    }

    inputcomp = 0;
    for ( MFIter mfi(mf_pencil); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();
        compute_vert_average(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_FAB(mf_pencil[mfi]),BL_TO_FORTRAN_FAB(mf_avg[mfi]),
                             &dir, &inputcomp, &outcomp, &ncomp);
    }

    if (write_data) {
        plotname = "mf_avg";
        VisMF::Write(mf_avg,plotname);
    }

}

// Functions for computing automatic refactorization of 3D -> 2D grids

int greatest_common_factor(int a, int b) {
    return b == 0 ? a : greatest_common_factor(b, a % b);
}

void factor(int num, int* factors, int nf) {
    int n = num;
    int cnt = 0;

    for (int i=0; i<nf; i++)
        factors[i]=1;

    while (n%2 == 0) { // factor out powers of 2
        n /= 2;
        factors[cnt%nf]*=2;
        cnt++;
    }
    for (int i=3; i<=sqrt(n); i+=2) { // find other odd factors
        while (n % i == 0) {
            n /= i;
            factors[cnt%nf]*=i;
            cnt++;
        }
    }
    if (n > 2)  { // check if n is a prime number greater than 2
        amrex::Error("CANNOT FACTOR PRIME NUMBER");
    }

    // check:
    n = 1;
    for (int i=0; i<nf; i++)
        n *= factors[i];

    if (n != num)  {
        amrex::Error("ERROR IN FACTORING");
    }
}
