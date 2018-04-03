#ifndef FHD_H_
#define FHD_H_

class FHD
{
public:

    /*
      public member functions
    */

    // constructor and destructor
    FHD ();
    virtual ~FHD ();

    void Init ();

    void Evolve ();

private:

    /*
      private data members
    */

    // these are defined in Maestro.cpp
    static amrex::IntVect nodal_flag;
    static amrex::IntVect nodal_flag_x;
    static amrex::IntVect nodal_flag_y;
    static amrex::IntVect nodal_flag_z;

    // keep track of old time, new time, and time step at each level
    amrex::Real t_old;
    amrex::Real t_new;
    amrex::Real dt;

    amrex::MultiFab sold;

    // stores domain boundary conditions
    // these muse be vectors (rather than arrays) so we can ParmParse them
    amrex::Vector<int> phys_bc;

    // boundary condition objects needed for FillPatch routines
    // this is essentially an array (over components)
    // of 2*DIM integer arrays storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<amrex::BCRec> bcs_s;  // for sold/snew
    amrex::Vector<amrex::BCRec> bcs_u;  // for uold/unew
    amrex::Vector<amrex::BCRec> bcs_f;  // "first-order extrap" used for S_cc, gpi, dSdt

};

#endif
