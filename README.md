# FHDeX

Code libraries and implementations of incompressible and low Mach number
multispecies mixtures with thermal fluctuations.

The algorithms use a fluctuating hydrodynamics approach, where each dissipative flux
is supplemented with a stochastic contribution.

The code uses a staggered grid velocity field and cell-centered scalars.

Also includes modules for:

* implicit viscosity via a projection-preconditioned Stokes solver

* a general multispecies diffusion model based on the Maxwell-Stefan formulation

* stochastic chemistry

* electrodiffusion

The development of the models and algorithms contained within are described
in a series of papers:

1. A. Donev, A. Nonaka, Y. Sun, T. Fai, A. Garcia and J. Bell,
*Low Mach Number Fluctuating Hydrodynamics of Diffusively Mixing Fluids*,
Comm. App. Math. and Comp. Sci., 9, 1, 2014.

2. A. Nonaka, Y. Sun, J. B. Bell, and A. Donev, 
*Low Mach Number Fluctuating Hydrodynamics of Binary Liquid Mixtures*,
Comm. App. Math. and Comp. Sci., 10, 2, 2015.

3. A. Donev, A. Nonaka, A. K. Bhattacharjee, A. L. Garcia, and J. B. Bell,
*Low Mach Number Fluctuating Hydrodynamics of Multispecies Liquid Mixtures*,
Phys. Fluids, 27, 037103, 2015.

4. J. Peraud, A. Nonaka, A. Chaudhri, J. B. Bell, A. Donev, and A. L. Garcia,
*Low Mach Number Fluctuating Hydrodynamics for Electrolytes*,
Phys. Rev. Fluids, 1, 074103, 2016.

5. C. Kim, A. Nonaka, J. B. Bell, A. L. Garcia, and A. Donev,
*Fluctuating Hydrodynamics of Reactive Liquid Mixtures*,
in preparation.

https://github.com/AMReX-Codes/FHDeX.git
