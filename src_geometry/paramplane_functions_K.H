#ifndef _paramplane_functions_K_H_
#define _paramplane_functions_K_H_

#include <AMReX.H>
#include <common_namespace.H>

void pre_check(FhdParticleContainer::ParticleType& part, const amrex_real delt, paramPlane paramplanes, 
                                                const int ns, const amrex_real* phi, const amrex_real* plo)
{


}

void find_inter_gpu(FhdParticleContainer::ParticleType& part, const amrex_real delt, paramPlane paramplanes, const int ns, int* intsurf,
                        amrex_real* inttime, int* intside, const amrex_real* phi, const amrex_real* plo)
{

    int flag = 0;
    *inttime = delt;
    int intSurf = -1;
    
    pre_check(part, delt, paramplanes, ns, phi, plo);

}


#endif
