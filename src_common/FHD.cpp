
#include <FHD.H>

using namespace amrex;

// helper IntVects used to define face/nodal MultiFabs
#if (AMREX_SPACEDIM == 2)
IntVect FHD::nodal_flag(1,1);
IntVect FHD::nodal_flag_x(1,0);
IntVect FHD::nodal_flag_y(0,1);
#elif (AMREX_SPACEDIM == 3)
IntVect FHD::nodal_flag(1,1,1);
IntVect FHD::nodal_flag_x(1,0,0);
IntVect FHD::nodal_flag_y(0,1,0);
IntVect FHD::nodal_flag_z(0,0,1);
#endif

FHD::FHD ()
{}

FHD::~FHD ()
{}

void
FHD::Init ()
{}

void
FHD::Evolve ()
{}
