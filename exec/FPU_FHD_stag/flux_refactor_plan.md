# Staggered 1D SPDE Flux Refactor Plan

## Goal

Repurpose `flux.cpp` from the legacy staggered compressible Navier-Stokes flux routine into a compact 1D anomalous transport SPDE flux routine:

```text
partial_t u + partial_x(A u - D partial_x u + B xi) = 0
```

The target layout is staggered:

- stretch is cell-centered in `cu` component 0
- momentum is x-face-centered in `cumom` component 0
- energy is cell-centered in `cu` component 2
- `cu` component 1 stores cell-centered averaged momentum for diagnostics and matrix coupling

## Implementation Summary

- Replace the old `calculateFluxStag` implementation in `flux.cpp` with a new 1D-only `calculateFlux`.
- Add the matching `calculateFlux` declaration to `FPU.H`.
- Compute only x-direction fluxes.
- Store stretch and energy fluxes on x-faces in `faceflux`.
- Store momentum flux at cell centers in `cenflux`.
- Use `A` and `D` as row-major 3x3 matrices from `FPU_namespace.H`.
- Use `B` as a diagonal 3-vector.
- Keep the existing RK3 structure, but update only the true conserved staggered variables:
  stretch and energy through `faceflux`, and momentum through `cenflux`.

## Flux Placement

- `faceflux(i,j,k,0)` is the stretch flux, using row 0 of `A`, `D`, and `B`.
- `cenflux(i,j,k,0)` is the momentum flux, using row 1 of `A`, `D`, and `B`.
- `faceflux(i,j,k,1)` is the energy flux, using row 2 of `A`, `D`, and `B`.

The physical flux sign should match the equation directly:

```text
flux = A u - D partial_x u + B xi
```

This matches the existing RK update form:

```text
u_new = u_old - dt * div(flux)
```

## Supporting Fixes

- Fix stochastic MultiFab setup syntax in `timeStep.cpp`.
- Scale stochastic fields by `1 / sqrt(dt * dV)`.
- Refresh `cu(...,1)` from staggered momentum after each RK stage.
- Add `driver_main.cpp`, `main_driver.cpp`, `timeStep.cpp`, and `flux.cpp` to the build in `Make.package`.
- Remove `main.cpp` from `Make.package` so it is retained in the repository but no longer required or compiled.
- Normalize input files to use array-style matrix parameters:
  - `A = ...` with 9 row-major entries
  - `D = ...` with 9 row-major entries
  - `B = ...` with 3 diagonal entries

## Test Plan

- Build with `make -j`.
- Run a deterministic constant-field test with `enable_fluctuations = 0`; flux divergence should be zero.
- Run a deterministic linear-gradient test to verify signs for `A u - D partial_x u`.
- Run a stochastic test with only one nonzero `B` entry and verify forcing affects only that component.
- Confirm periodic ghost fills occur before flux evaluation and no out-of-bounds face/cell indexing occurs.

## Assumptions

- The target executable is the staggered driver path, not the older standalone `main.cpp` SPDE loop.
- `main.cpp` remains in place as a defunct historical file and is not part of the executable build.
- `B` is diagonal and represented as `GpuArray<Real,3> B`.
- `A` and `D` are row-major 3x3 arrays where row 0 is stretch, row 1 is momentum, and row 2 is energy.
- The existing RK3 structure remains in place.
