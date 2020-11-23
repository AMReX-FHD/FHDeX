AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE

#if (AMREX_SPACEDIM == 2)

Real InhomogeneousBCVal ( int comp,
                  Real x,
                  Real y,
                  Real time)
{
    return 0.;
};

#elif (AMREX_SPACEDIM == 3)


Real InhomogeneousBCVal ( int comp,
                  Real x,
                  Real y,
                  Real z,
                  Real time)
{
    return 0.;
};


#endif
