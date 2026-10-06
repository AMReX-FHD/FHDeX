# src_deankow/flocking

Dean–Kawasaki model for active Brownian particles (`src_deankow/flocking.pdf`),
in two versions that share parameters and initial conditions so they can be
compared directly. This is the base for a later flocking (alignment)
interaction.

**Particles:** N particles move at constant speed v in the direction θ_i, which
diffuses:

  dX_i = g(θ_i) dt,  dθ_i = √(2D) dW_i,  g(θ) = v (cos θ, sin θ)

**SPDE** (Müller, von Renesse, Zimmer 2025), with the density written φ in the
code (ρ in the write-up):

  dφ = D ∂θθ φ dt − ∇ₓ·(φ g(θ)) dt + N^(−1/2) ∂θ(√(2Dφ) dB)

φ is normalized to integrate to 1 over [0, Lx) × [0, Ly) × [0, 2π). All
directions are periodic.

## Two builds from one directory

| Build | What | Executable |
| --- | --- | --- |
| `make DIM=3` | SPDE on (x, y, θ), with θ the z direction, always [0, 2π) | `main3d.gnu.MPI.ex` |
| `make DIM=2` | particles in (x, y), with θ a particle attribute | `main2d.gnu.MPI.ex` |

`Make.package` selects the sources by `DIM`. The code is built and run in
`exec/dean_kow/flocking`:

```
cd exec/dean_kow/flocking
make DIM=3 -j8 && mpirun -np 4 ./main3d.gnu.MPI.ex inputs_spde
make DIM=2 -j8 && mpirun -np 4 ./main2d.gnu.MPI.ex inputs_particles
python3 ../../../python/flocking_compare.py plt_spde_000180 plt_part_000035 --num-part 4e6
```

The inputs files share `inputs_common` through `FILE = inputs_common`, which is
resolved relative to the directory the code is run from.

## Files

| File | Build | Contents |
| --- | --- | --- |
| `FlockingParams.H` | both | `FlockingParams`, `read_flocking_params`, the initial density `flocking_ic` |
| `main.cpp` | both | driver: the SPDE for DIM = 3, the particles for DIM = 2 |
| `AmrCoreFlock.H/.cpp` | 3 | single-level `AmrCore` class: setup, time loop, time step, plotfiles, checkpoint/restart |
| `AdvancePhiAtLevel.cpp` | 3 | one time step (Euler–Maruyama or Heun) |
| `myfunc.H/.cpp` | 3 | deterministic fluxes, noise flux, flux update |
| `mykernel.H` | 3 | face-flux and update kernels |
| `FlockPC.H/.cpp` | 2 | particle container (NeighborParticleContainer with a θ attribute): sampling, Euler–Maruyama step, θ histogram |
| `FlockParticles.H/.cpp` | 2 | particle driver: grid, time loop, plotfiles |

## Discretization (SPDE)

The fluxes are finite-volume fluxes on cell faces:

- **x, y:** advection φ·v(cos θ_k, sin θ_k), with θ_k the cell-centre θ of the
  plane. The flux is set by `adv_order`:
  - 0: centred;
  - 1: first-order upwind;
  - 2: MUSCL with a `limiter`, 0 (minmod) or 1 (MC).
- **θ:** −D ∂θφ (centred), plus a noise flux √(2Dφ_face/(N·dV·dt))·Z with
  Z ~ N(0, 1) on every θ face. Here φ_face is the average of the two adjacent
  cells, clipped at 0.
- **Turning-rate slot:** an upwinded turning-rate term ω·φ_face is in the
  kernel, with ω = 0. This is where the alignment interaction will enter.

**Time integration** (`time_integrator`):
- 0: Euler–Maruyama.
- 1: Heun for the deterministic fluxes, with the noise built once from φⁿ and
  used in both stages.

**Time step:** dt = `adv.cfl` / (v(1/dx + 1/dy) + 2D/dθ²). This is the
positivity limit of the upwind scheme at cfl = 1.

**Plotfile:** the SPDE plotfile has two variables:
- `phi`: φ(x, y, θ).
- `rho`: the spatial density ρ(x, y) = Σ_k φ(x, y, θ_k)·dθ, the integral over
  θ. It is the same in every θ plane, so it can be viewed alongside `phi` on
  the 3D grid and compared directly with the particle plotfile's `rho`.

**Which advection to use:**

| Scheme | Advection accuracy (blob, L1 error order) | Fluctuations from a uniform start (var / Poisson) | Positivity |
| --- | --- | --- | --- |
| centred + Heun | 2nd (1.99) | φ: 1.00, ρ: 1.00 | not guaranteed |
| upwind + Euler–Maruyama | 1st (0.86) | φ: 0.49, ρ: 0.04 | yes, cfl ≤ 1 |
| MUSCL (MC) + Heun | 2nd (1.90) | φ: 0.42, ρ: 0.07 | not guaranteed with Heun |

The continuum SPDE has Poisson white noise φ₀/N·δ as its stationary
covariance, because the advection is skew-adjoint. Only the centred flux keeps
that property after discretization. Upwind and limited fluxes add numerical
diffusion without matching noise, which damps the fluctuations, the spatial
density ρ above all. So:
- **Stochastic runs:** use centred + Heun, the default in `inputs_spde`.
- **Deterministic mean-field runs with sharp features:** use MUSCL or upwind
  (`inputs_compare_det` uses MUSCL + Heun).
- **Not allowed:** centred advection with forward Euler, which the code rejects.

## Particles

- **Initialization:** each cell receives about N times its share of
  ∫ flocking_ic (midpoint rule in x, y; 64 points in θ). Within a cell, x and y
  are uniform and θ is drawn from f(x, y, ·) by rejection.
- **Step:** each Euler–Maruyama step moves particles with their start-of-step
  θ, adds √(2D dt)·Z to θ, wraps θ into [0, 2π), and redistributes. The time
  step is `particle_dt`, or `adv.cfl`·min(dx, dy)/v by default.
- **Output plotfile** (2D):
  - `phi_t000 … phi_tNNN`: the density in each of the `n_theta_bins` θ bins,
    normalized like the SPDE's φ. The default is the third entry of
    `amr.n_cell`, so it matches the SPDE grid cell for cell.
  - `rho`, `px`, `py`: the moments, ρ = Σ φ_k dθ and p = Σ φ_k (cos, sin)θ_k dθ.
  - With `write_particles = 1`, the particle positions and θ.

## Inputs

| Input | Default | Build | Meaning |
| --- | --- | --- | --- |
| `num_part` | (required) | both | N |
| `diff_coeff` | 1 | both | D, rotational diffusion |
| `flock.speed` | 1 | both | v |
| `flock.ic_type` | 0 | both | Initial density: 0 uniform; 1 1 + ic_amp·cos(2π(ic_kx x/Lx + ic_ky y/Ly))·cos(ic_m(θ − ic_theta0)); 2 Gaussian blob (ic_sigma at ic_x0, ic_y0) × (1 + ic_amp·cos(θ − ic_theta0)) |
| `flock.ic_amp`, `ic_kx`, `ic_ky`, `ic_m`, `ic_theta0`, `ic_sigma`, `ic_x0`, `ic_y0` | 0.5, 1, 0, 1, 0, 0.1, 0.5, 0.5 | both | Initial-condition parameters. \|ic_amp\| ≤ 1 for types 1 and 2. |
| `flock.align_K` | 0 | both | Alignment strength. Not implemented yet; must be 0. |
| `geometry.prob_lo`, `prob_hi` | 0 0, 1 1 | both | x, y extent. A third entry, for θ, must be 0 and 2π. |
| `amr.n_cell` | (required) | both | nx ny ntheta |
| `amr.max_grid_size` | 32 / 64 | both | Grid size |
| `stop_time`, `max_step` | ∞ | both | Run length |
| `amr.plot_dt` | off | both | Write a plotfile at multiples of this time, with dt clipped to hit it, so both codes write at the same times |
| `amr.plot_int` | off | both | Write a plotfile every n steps |
| `amr.plot_file` | `plt_spde_` / `plt_part_` | both | Plotfile prefix |
| `seed` | 0 | both | Random seed; 0 means from the clock |
| `adv.cfl` | 0.5 | both | Time-step fraction (see above). For particles it sets the default `particle_dt`. |
| `dorand` | 1 | 3 | 0 turns the noise off, giving the deterministic mean-field equation |
| `adv_order` | 1 | 3 | 0 centred, 1 upwind, 2 MUSCL (`inputs_spde` uses 0) |
| `limiter` | 0 | 3 | MUSCL limiter: 0 minmod, 1 MC |
| `time_integrator` | 0 | 3 | 0 Euler–Maruyama, 1 Heun (`inputs_spde` uses 1) |
| `diag_int` | 1 | 3 | Print the mass and min φ every n steps |
| `amr.chk_int`, `amr.chk_file`, `amr.restart` | off | 3 | Checkpoint and restart |
| `n_theta_bins` | ntheta | 2 | Number of θ histogram bins |
| `particle_dt` | cfl·min(dx, dy)/v | 2 | Particle time step |
| `write_particles` | 1 | 2 | Also write particle data to the plotfile |

## Checks done

All on 64² × 32 with N = 4·10⁶ and v = D = 1, unless noted.

- **Conservation:** the SPDE conserves mass to round-off.
- **Particles vs deterministic SPDE** (MUSCL + Heun, mode initial condition,
  t = 0.25): the mode amplitudes differ by 0.82σ of particle sampling noise. ρ
  and the full φ agree to the Poisson noise level (3% and 18%).
- **θ diffusion only** (v = 0, m = 2, t = 0.25):
  - SPDE: A(t)/A(0) = 0.37267, against 0.37261 for the discrete Laplacian.
  - Particles: 0.36717, against the exact e⁻¹ = 0.36788.
- **D = 0 advection against the exact solution:** observed orders 1.99
  (centred), 1.90 (MUSCL) and 0.86 (upwind).
- **Fluctuations:** see the table above.

## Alignment interaction, later

The turning rate ω(x, θ) is already a slot in the θ flux (`compute_flux_theta`).
A Vicsek or Kuramoto type alignment would be:
- **SPDE:** ω = −K·Im(e^(−iθ)·(W ⋆ P)(x)), with polarization
  P = Σ_k φ_k e^(iθ_k) dθ. That is one 2D FFT convolution of the θ-moments,
  in the same pattern as the `U ⋆ φ` convolution in `interaction`. `USE_FFT`
  is already on.
- **Particles:** Σ_j W(|x_i − x_j|)·sin(θ_j − θ_i)/N, a pair sum over
  neighbours. `FlockPC` is already a `NeighborParticleContainer` for this.
