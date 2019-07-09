void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
	           const amrex::MultiFab& cu,
	           const amrex::MultiFab& cuMeans,
	           const amrex::MultiFab& cuVars);
