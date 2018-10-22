
void ssaStep(MultiFab& cu, MultiFab& prim, std::array<MultiFab, AMREX_SPACEDIM>& stochFlux, const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt)
{


    AMREX_D_TERM(MultiFABFillRandom(stochFlux[0], 0, 1, geom);
         MultiFABFillRandom(stochFlux[0], 1, 1, geom);
         MultiFABFillRandom(stochFlux[0], 2, 1, geom);
         MultiFABFillRandom(stochFlux[0], 3, 1, geom);,
         MultiFABFillRandom(stochFlux[1], 0, 1, geom);
         MultiFABFillRandom(stochFlux[1], 1, 1, geom);
         MultiFABFillRandom(stochFlux[1], 2, 1, geom);
         MultiFABFillRandom(stochFlux[1], 3, 1, geom);,
         MultiFABFillRandom(stochFlux[2], 0, 1, geom);
         MultiFABFillRandom(stochFlux[2], 1, 1, geom);
         MultiFABFillRandom(stochFlux[2], 2, 1, geom);
         MultiFABFillRandom(stochFlux[2], 3, 1, geom););


        


}
