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

## Density-dependent speed (motility-induced phase separation)

Particles slow down where the smoothed density is high. Both codes use the same
speed function s, kernel W and grid:

  v(x) = v₀ · s(ρ̃(x)/ρ̄),  ρ̃ = W ⋆ ρ,  ρ = ∫φ dθ,  ρ̄ = 1/(Lx·Ly)

- **Particles:** dXᵢ = v₀·s(ρ̃(Xᵢ)/ρ̄)·(cos θᵢ, sin θᵢ) dt. The θ equation is
  unchanged.
- **SPDE:** the x and y fluxes become φ·v(ρ̃)·(cos θ, sin θ), with ρ̃ averaged
  from the two cells on each side of the face. The θ diffusion and θ noise are
  unchanged: spatial motion is deterministic given ρ̃, so no new noise term
  appears.

**Speed functions** (`flock.speed_type`), with u = ρ̃/ρ̄ and λ = `flock.speed_lambda`:
- 0: s = 1, no interaction (the default).
- 1: s = exp(−λu).
- 2: s = max(0, 1 − λu). Needs λ < 1.

**Kernels** (`flock.kernel_type`, radius R = `flock.kernel_R`):

| Type | W(r) ∝ | Ŵ(k) | Notes |
| --- | --- | --- | --- |
| 0 Gaussian (default) | exp(−r²/2R²) | exp(−k²R²/2) | Ŵ > 0; cut off at 4R for the pair sum |
| 1 top-hat | 1 for r < R | 2J₁(kR)/(kR) | "neighbours within R"; negative lobe about −0.13 |
| 2 parabolic | (1 − r²/R²)₊ | 8J₂(kR)/(kR)² | same shape as `hk` |
| 3 bump | (1 − r²/R²)₊² | 48J₃(kR)/(kR)³ | C¹, compact; best for the pair sum |

On the grid, W is sampled at the cell offsets and normalized so that
Σ W dx dy = 1, which makes the mean of ρ̃ exactly ρ̄. `DensityConv.H/.cpp`
does the periodic FFT convolution, in the same pattern as the `U ⋆ φ`
convolution in `interaction`:
- **SPDE:** W is replicated over θ, so convolving φ gives ρ̃ in every θ plane.
- **Particles:** it convolves the deposited density.

**How the particles get ρ̃** (`flock.density_method`):
- 0, particle-mesh (default): cloud-in-cell deposit onto the grid, FFT
  convolution, then cloud-in-cell interpolation back to each particle.
  O(N + M log M).
- 1, pair sum: ρ̃ᵢ = (1/N)·Σⱼ W(|xᵢ − xⱼ|) over neighbours, plus the particle
  itself, with the continuous kernel. This is the grid-free particle model,
  for small N. `flock.density_check = 1` prints its difference from
  particle-mesh at the initial positions.

**Linear stability.** At startup the code prints d ln s/d ln u at the mean
density, the persistence length v/D, and the growth rate of the first four box
modes in the diffusive limit:

  σ(k) = −D_eff·k²·(1 + Ŵ(k)·d ln s/d ln u),  D_eff = v(ρ̄)²/(2D)

The uniform state is unstable when 1 + Ŵ·d ln s/d ln u < 0; for exponential s
that means λŴ(k) > 1. The formula assumes k·ℓ_p ≪ 1. At k·ℓ_p ≈ 0.6–1 it is
about 8% off. The exact rate for the discrete θ grid is the leading
eigenvalue of the linearized kinetic equation, and the SPDE matches that to 1%
(see the checks below).

**Example:** `exec/dean_kow/flocking/inputs_mips` (SPDE) and
`inputs_mips_particles`:
- 128² × 32 grid;
- v₀ = 4, D = 10, exponential s with λ = 1.5;
- Gaussian kernel with R = 0.03;
- N = 10⁶, from a uniform start.

**Plotfiles:**
- SPDE: `rhot` (ρ̃) alongside `phi` and `rho`.
- Particles: `rhot` on the grid (by particle-mesh), and a `rhot` particle
  attribute.

## Inputs

| Input | Default | Build | Meaning |
| --- | --- | --- | --- |
| `num_part` | (required) | both | N |
| `diff_coeff` | 1 | both | D, rotational diffusion |
| `flock.speed` | 1 | both | v |
| `flock.ic_type` | 0 | both | Initial density: 0 uniform; 1 1 + ic_amp·cos(2π(ic_kx x/Lx + ic_ky y/Ly))·cos(ic_m(θ − ic_theta0)); 2 Gaussian blob (ic_sigma at ic_x0, ic_y0) × (1 + ic_amp·cos(θ − ic_theta0)) |
| `flock.ic_amp`, `ic_kx`, `ic_ky`, `ic_m`, `ic_theta0`, `ic_sigma`, `ic_x0`, `ic_y0` | 0.5, 1, 0, 1, 0, 0.1, 0.5, 0.5 | both | Initial-condition parameters. \|ic_amp\| ≤ 1 for types 1 and 2. |
| `flock.align_K` | 0 | both | Alignment strength. Not implemented yet; must be 0. |
| `flock.speed_type` | 0 | both | Speed function: 0 constant, 1 exp(−λu), 2 max(0, 1 − λu) |
| `flock.speed_lambda` | 1 | both | λ; must be < 1 for speed_type 2 |
| `flock.kernel_type` | 0 | both | Sensing kernel: 0 Gaussian, 1 top-hat, 2 parabolic, 3 bump |
| `flock.kernel_R` | 0.05 | both | Kernel radius R (keep at least about 4 dx) |
| `flock.density_method` | 0 | 2 | ρ̃ at the particles: 0 particle-mesh, 1 pair sum |
| `flock.density_check` | 0 | 2 | With the pair sum: print its difference from particle-mesh at startup |
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

### Density-dependent speed

- **No interaction is unchanged:** with `speed_type = 0`, both codes reproduce
  the earlier results bit for bit (SPDE φ, and all particle histogram fields).
- **Convolution:** ρ̃ has mean exactly ρ̄ and is identical in every θ plane.
  A cos(2πx) density mode is damped by 0.951975, against the exact Gaussian
  Ŵ = 0.951850 at R = 3.2 dx.
- **Particle-mesh vs pair sum:** bump kernel, N = 2·10⁴, 64² grid, same
  particles. The rms relative difference is 1.4% at R = 2.6 dx, 0.15% at
  5.1 dx and 0.017% at 10.2 dx.
- **Linear growth rates:** deterministic centred SPDE, mode k = 2π, v₀ = 4,
  D = 10, R = 0.03, against the leading kinetic eigenvalue:

| λ | SPDE | Kinetic eigenvalue | Diffusive-limit formula |
| --- | --- | --- | --- |
| 0.8 (stable) | −1.448 | −1.462 | −1.365 |
| 1.5 (unstable) | +0.684 | +0.686 | +0.745 |

## Alignment interaction, later

The turning rate ω(x, θ) is already a slot in the θ flux (`compute_flux_theta`).
A Vicsek or Kuramoto type alignment would be:
- **SPDE:** ω = −K·Im(e^(−iθ)·(W ⋆ P)(x)), with polarization
  P = Σ_k φ_k e^(iθ_k) dθ. That is one 2D FFT convolution of the θ-moments,
  in the same pattern as the `U ⋆ φ` convolution in `interaction`. `USE_FFT`
  is already on.
- **Particles:** Σ_j W(|x_i − x_j|)·sin(θ_j − θ_i)/N, a pair sum over
  neighbours. `FlockPC` is already a `NeighborParticleContainer` for this.
