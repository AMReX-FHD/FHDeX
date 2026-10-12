# The flocking model: formulation, stability and kernels

This note describes the model solved by the codes in `src_deankow/flocking`:
- active Brownian particles moving at a speed that can depend on the local
  smoothed density;
- the Dean–Kawasaki SPDE for their empirical density;
- the linear stability of the uniform state;
- the choice of sensing kernel.

`README.md` covers the code layout, inputs and numerical checks in detail.
The density is written φ in the code; it is ρ in `src_deankow/flocking.pdf`.

## 1. Particle model

N particles move in the periodic box [0, Lx) × [0, Ly). Each particle has a
position Xᵢ and a heading θᵢ ∈ [0, 2π):

  dXᵢ = v(Xᵢ) e(θᵢ) dt,  e(θ) = (cos θ, sin θ)
  dθᵢ = √(2D) dWᵢ

The Wᵢ are independent standard Brownian motions and D is the rotational
diffusivity (`diff_coeff`). The speed is

  v(x) = v₀ · s(ρ̃(x)/ρ̄)

with self-propulsion speed v₀ (`flock.speed`) and a dimensionless speed
function s with s(u) ≤ 1. ρ̃ is the density seen through a sensing kernel W
of range R:

  ρ̃(x) = (1/N) Σⱼ W(x − Xⱼ),  ∫ W dA = 1

ρ̄ = 1/(Lx·Ly) is the mean of ρ̃. With s ≡ 1 the particles are independent
active Brownian particles. With s decreasing they slow down in crowded
regions, which can lead to motility-induced phase separation (MIPS).

Positions are moved deterministically given the headings; all the noise is in
θ. The coupling is mean-field: each particle contributes W/N to ρ̃.

**Length and time scales:**
- **Persistence time** 1/D: a particle loses memory of its heading on this
  time.
- **Persistence length** ℓ_p = v/D: the distance a particle travels before its
  heading decorrelates.
- **Long-time diffusivity** D_eff = v²/(2D), in 2D. On scales ≫ ℓ_p, a
  particle at constant speed v performs a random walk with this diffusivity.

**The dimensionless groups** that control the behaviour:
- L/ℓ_p, the box size in persistence lengths;
- R/ℓ_p, the sensing range in persistence lengths;
- the speed-function parameter λ;
- N, which sets the size of fluctuations.

## 2. Phase-space density

Let

  φ(x, θ, t) = (1/N) Σᵢ δ(x − Xᵢ(t)) δ(θ − θᵢ(t))

be the empirical phase-space density, normalized so that ∫∫ φ dθ dx = 1. Its
moments are:
- the spatial density ρ(x) = ∫ φ dθ;
- the polarization p(x) = ∫ e(θ) φ dθ.

The smoothed density is ρ̃ = W ⋆ ρ, the spatial convolution. In the code, φ is
a cell-averaged field on a 3D grid (x, y, θ), and θ is the z direction.

### 2.1 Mean-field (kinetic) equation

As N → ∞, φ obeys the nonlinear Fokker–Planck equation

  ∂ₜφ = D ∂θθφ − ∇ₓ·( v(ρ̃) e(θ) φ ),  v(ρ̃) = v₀ s(ρ̃/ρ̄)

### 2.2 Dean–Kawasaki SPDE

For finite N, Itô's formula applied to the empirical measure gives, formally
(Dean 1996; Müller, von Renesse and Zimmer 2025),

  dφ = D ∂θθφ dt − ∇ₓ·( v(ρ̃) e(θ) φ ) dt + N^(−1/2) ∂θ( √(2Dφ) dB )

where B is space–angle–time white noise. For a test function f(x, θ),

  d⟨φ, f⟩ = ⟨φ, v e·∇ₓf + D ∂θθ f⟩ dt + (1/N) Σᵢ √(2D) ∂θ f(Xᵢ, θᵢ) dWᵢ

The martingale term has quadratic variation (2D/N)·⟨φ, (∂θf)²⟩ dt. That is
exactly the conservative noise N^(−1/2) ∂θ(√(2Dφ) ξ). Two points follow:
- **Noise only in θ:** there is no noise in x, because the spatial motion is
  deterministic given the headings.
- **The interaction adds no noise:** it acts only through v(ρ̃), which is
  itself a functional of φ.

The noise conserves mass at each spatial point: it moves probability in θ
only.

### 2.3 Moment equations and the diffusive limit

Integrating the kinetic equation over θ, and over θ against e(θ), gives

  ∂ₜρ = −∇·( v p )
  ∂ₜp = −D p − ∇·( v ∫ e eᵀ φ dθ )

On length scales ≫ ℓ_p and time scales ≫ 1/D, the polarization relaxes
quickly and the second moment is close to isotropic, ∫ e eᵀ φ dθ ≈ (ρ/2) I.
Then p ≈ −∇(vρ)/(2D), and

  ∂ₜρ = ∇·( (v/2D) ∇(vρ) )

For a local speed v(ρ) this becomes

  ∂ₜρ = ∇·( D_eff(ρ) [1 + d ln v / d ln ρ] ∇ρ ),  D_eff(ρ) = v(ρ)²/(2D)

The effective collective diffusivity turns negative when
d ln v / d ln ρ < −1. That is the spinodal condition for MIPS (Tailleur and
Cates 2008; Cates and Tailleur 2015). Particles accumulate where they are slow,
and if they slow down fast enough with density, the accumulation runs away.

## 3. Linear stability of the uniform state

The uniform state φ₀ = ρ̄/(2π) is a stationary solution for any s and W.
Write u = ρ̃/ρ̄ and define, at the mean density u = 1:

  v₁ = v₀ s(1),  g = d ln s / d ln u |_{u=1} = s′(1)/s(1)

Perturb φ = φ₀ + δφ. Then δρ̃ = W ⋆ δρ, and the speed changes by
δv = v₁ g δρ̃/ρ̄.

### 3.1 Exact (kinetic) dispersion relation

Take a Fourier mode e^{ik·x} with k along x, and write δφ in angular
harmonics, δφ = Σₘ aₘ e^{imθ}, so that δρ = 2π a₀. The linearized kinetic
equation is

  ∂ₜ δφ = D ∂θθ δφ − ik v₁ cos θ [ δφ + g Ŵ(k) a₀ ]

Here Ŵ is the 2D Fourier transform of W. In harmonics:

  ∂ₜ a₀ = −(ik v₁/2) (a₁ + a₋₁)
  ∂ₜ a±1 = −D a±1 − (ik v₁/2) [ (1 + g Ŵ(k)) a₀ + a±2 ]
  ∂ₜ aₘ = −D m² aₘ − (ik v₁/2) (aₘ₋₁ + aₘ₊₁),  |m| ≥ 2

The interaction enters only through the m = 0 → ±1 coupling, multiplied by
1 + gŴ(k). The growth rate σ(k) is the eigenvalue of this tridiagonal system
with the largest real part; it can also be written as a continued fraction.
The SPDE resolves θ with n_θ cells, so the code's exact rate is the leading
eigenvalue of the same operator with θ discretized by the n_θ-point
second-difference Laplacian.

### 3.2 Diffusive limit

For kℓ₁ ≪ 1, with ℓ₁ = v₁/D, the harmonics with |m| ≥ 2 are negligible and
a±1 relax quasi-statically:

  a±1 ≈ −(ik v₁ / 2D) (1 + g Ŵ(k)) a₀

That gives

  σ(k) = −D_eff k² ( 1 + g Ŵ(k) ),  D_eff = v₁²/(2D)

This is the Fourier form of §2.3 with the nonlocal speed. The uniform state is
unstable at wavenumber k when

  **1 + g Ŵ(k) < 0**

**Consequences:**
- **Onset.** Ŵ(0) = 1, so long waves are the first to go unstable. The onset
  is g = −1: d ln v/d ln ρ = −1 at the mean density.
- **Finite unstable band.** For a decreasing Ŵ (Gaussian), the band is
  0 < k < k_c with Ŵ(k_c) = −1/g. The sensing range sets the short-wavelength
  cutoff: k_c ~ 1/R.
- **Fastest mode.** σ is largest where k²(−1 − gŴ(k)) peaks. That sets the
  initial domain spacing, which coarsens later (nonlinear regime).
- **Mass conservation.** σ → 0 as k → 0, and the k = 0 mode (total mass) is
  neutral.
- **Kinetic corrections.** At kℓ₁ ≈ 0.5–1 the diffusive formula overestimates
  the growth rate. Finite persistence lags the response. The startup printout
  reports the diffusive rates and kℓ_p for the first four box modes.

### 3.3 Speed functions

| `flock.speed_type` | s(u) | s(1) | g = d ln s/d ln u at u = 1 | Unstable (long waves) when |
| --- | --- | --- | --- | --- |
| 0 | 1 | 1 | 0 | never: independent particles |
| 1 | exp(−λu) | e^(−λ) | −λ | λ > 1 |
| 2 | max(0, 1 − λu), λ < 1 | 1 − λ | −λ/(1 − λ) | λ > 1/2 |

The exponential form stays positive at any density. Stronger λ also lowers
the speed at the mean density, v₁ = v₀e^(−λ), which reduces D_eff = v₁²/2D
and slows the dynamics. The linear form stops particles completely at
u = 1/λ. It needs λ < 1 so the mean-density speed is positive.

### 3.4 Checks of the code against the theory

The deterministic SPDE (`dorand = 0`, centred advection) was run with
v₀ = 4, D = 10, a Gaussian kernel with R = 0.03, and mode k = 2π. Its measured
growth rates compare as follows:

| λ | SPDE | Leading kinetic eigenvalue (32-point θ) | Diffusive limit |
| --- | --- | --- | --- |
| 0.8 (stable) | −1.448 | −1.462 | −1.365 |
| 1.5 (unstable) | +0.684 | +0.686 | +0.745 |

kℓ₁ is 1.1 and 0.56 for these two runs. The SPDE matches the kinetic theory
to 1%. The diffusive formula is 6–9% off at these kℓ₁.

### 3.5 The non-interacting case

With s ≡ 1 the equation is linear, and every perturbation decays. The θ
harmonics damp at rates Dm². Spatial structure relaxes diffusively with D_eff
on scales ≫ ℓ_p, and ballistically below ℓ_p.

The fluctuations are Poisson. The stationary covariance of the linear SPDE is
⟨δφ(x,θ) δφ(x′,θ′)⟩ = (φ₀/N) δ(x − x′) δ(θ − θ′). The θ diffusion and the θ
noise balance each other, and the advection −e·∇ₓ is skew-adjoint, so it
neither creates nor destroys variance.

The discretization must keep this property. Centred advection does, so the
code reproduces Poisson statistics (variance ratio 1.00). Upwind and limited
(MUSCL) fluxes add numerical diffusion without matching noise. They damp the
spatial density variance to about 5% of the correct value. This is why the
stochastic inputs use centred advection, with the stochastic RK3 integrator.

## 4. Choice of sensing kernel

W is radially symmetric, periodic, and normalized to ∫W dA = 1. By §3.2,
stability depends on W only through Ŵ(k). Four kernels are available
(`flock.kernel_type`, radius R = `flock.kernel_R`):

| Type | W(r) ∝ | Normalization | Ŵ(k) | Smoothness | Support |
| --- | --- | --- | --- | --- | --- |
| 0 Gaussian | exp(−r²/2R²) | 1/(2πR²) | exp(−k²R²/2) | C^∞ | infinite: cut off at 4R for the pair sum |
| 1 top-hat | 1, r < R | 1/(πR²) | 2J₁(kR)/(kR) | discontinuous | compact |
| 2 parabolic | 1 − r²/R², r < R | 2/(πR²) | 8J₂(kR)/(kR)² | continuous, kink at R | compact |
| 3 bump | (1 − r²/R²)², r < R | 3/(πR²) | 48J₃(kR)/(kR)³ | C¹ | compact |

**How the kernel enters the dynamics:**
- **Onset is kernel-independent.** It is set by g at k → 0, because Ŵ(0) = 1
  for every kernel.
- **The unstable band and fastest mode depend on the kernel.** They follow
  from 1 + gŴ(k) < 0. A larger R narrows the band toward long waves (k_c ~ 1/R)
  and gives larger initial domains.
- **Negative lobes.** The compact kernels have Ŵ < 0 in lobes; the top-hat
  goes down to about −0.13 near kR ≈ 5.1. For a decreasing s (g < 0), a
  negative Ŵ makes 1 + gŴ > 1 there, so the lobes are stabilizing and cannot
  create a finite-wavelength instability. They would matter for an increasing
  speed function (g > 0). The Gaussian has no lobes and gives the cleanest
  long-wave instability.
- **Grid resolution.** On the grid, W is sampled at the cell offsets and
  normalized discretely, Σ W dx dy = 1, so that the mean of ρ̃ is exactly ρ̄.
  A cos(kx) mode is then damped by the discrete transform, which matches the
  continuous Ŵ to O((dx/R)²) for smooth kernels. Measured: 0.951975 against
  0.951850 for the Gaussian at R = 3.2 dx. Keep R ≥ about 4 dx; the top-hat,
  with its edge, needs more.
- **Particles: particle-mesh vs pair sum.** Particle-mesh deposits the
  particles with cloud-in-cell weights, convolves with the discrete W, and
  interpolates back with cloud-in-cell weights. This adds a smoothing of
  about one cell. Against the exact pair sum with the continuous kernel (bump,
  same particles), the rms difference is 1.4% at R = 2.6 dx, 0.15% at 5.1 dx
  and 0.017% at 10.2 dx: about (dx/R)³.
- **Cost of the pair sum.** It grows like N·πR²/A per particle, so compact
  kernels (bump or parabolic) keep it bounded. The Gaussian needs its 4R
  cutoff, which costs about 16× the neighbours of a compact kernel of the
  same R.

**Recommendations:**
- **SPDE and particle-mesh runs:** the Gaussian, for its positive, monotone Ŵ
  and its smoothness.
- **Pair-sum particle runs:** the bump, for compact support and C¹
  smoothness.
- **Comparison with the literature:** the top-hat, for "neighbours within R"
  models.

## 5. Numerical formulation

**SPDE** (3D build; `AmrCoreFlock`, `mykernel.H`). A finite-volume scheme on
(x, y, θ) cells, with θ_k the cell-centre heading. The fluxes are:
- **x, y faces:** advection φ·v_face·(cos θ_k, sin θ_k), with
  v_face = v₀ s(ρ̄⁻¹·(ρ̃_L + ρ̃_R)/2). The flux is centred (`adv_order = 0`),
  first-order upwind (1), or limited MUSCL (2).
- **θ faces:** −D ∂θφ (centred), plus the noise √(2Dφ_face/(N ΔV Δt))·Z,
  with Z ~ N(0, 1).

**How ρ̃ is computed in the SPDE** (`DensityConv.cpp`). In the continuum,
ρ̃(x) = ∫ W(x − x′) ρ(x′) dx′ with ρ = ∫ φ dθ. On the grid, with
nx × ny × nθ cells, Δθ = 2π/nθ, and φ normalized so that
Σ φ ΔxΔyΔθ = 1, the steps are:

1. **Integrate over θ.** ρ_ij = Σ_k φ_ijk Δθ, the spatial density on the
   x–y plane (`ReduceToPlaneMF`, times Δθ).
2. **Sample and normalize the kernel.** For each cell offset (a, b), with
   a = −nx/2, …, nx/2 − 1 and b likewise (the nearest periodic image):

     w_ab = W_shape(√((aΔx)² + (bΔy)²)),  W_ab = w_ab / (Σ_a′b′ w_a′b′ ΔxΔy)

   W_shape is the unnormalized shape of `flock.kernel_type`, with the
   Gaussian cut off at 4R. The discrete normalization Σ W_ab ΔxΔy = 1
   replaces the continuous one. The code prints the ratio of the two at
   startup; it is 1 when R ≳ 4Δx. This is done once, at setup, and the 2D FFT
   Ŵ_d of W_ab is stored.
3. **Convolve periodically.**

     ρ̃_ij = Σ_a Σ_b W_ab ρ_(i−a mod nx),(j−b mod ny) ΔxΔy

   This is a circular convolution, done as ρ̃ = IFFT(Ŵ_d · FFT(ρ)) with
   real-to-complex 2D FFTs on the x–y plane. The inverse is scaled by
   ΔxΔy/(nx ny). In Fourier space, ρ̃̂(k) = Ŵ_d(k) ρ̂(k), with Ŵ_d(0) = 1.
4. **Copy into every θ plane.** ρ̃_ijk = ρ̃_ij for all k. The periodic ghost
   cells (one layer) are then filled.
5. **Face values for the speed.** On the x face (i − ½, j, k):

     u = ρ̃_(i−½),j / ρ̄,  ρ̃_(i−½),j = ½(ρ̃_(i−1),j + ρ̃_ij),  ρ̄ = 1/(Lx Ly)

   The face velocity is v₀ s(u) cos θ_k, and the y faces are the same with
   sin θ_k. ρ̃ is averaged before s is applied, not s(ρ̃) after.

**Properties:**
- **Mass:** the discrete normalization makes the mean of ρ̃ exactly the mean
  of ρ, which is ρ̄.
- **No θ dependence:** W does not depend on θ, so only the θ-integral of φ
  enters. A 2D FFT pair per evaluation replaces a 3D one.
- **Kernel range:** the support of W (4R for the Gaussian, R for the others)
  must be below L/2, so that each offset has a single periodic image.

**When ρ̃ is recomputed:** from the current state at every stage of the time
integrator. It is computed from φⁿ for Euler–Maruyama (`time_integrator = 0`);
from φⁿ and then the predictor for Heun (1); and from φⁿ, u₁ and u₂ for the
stochastic SSP-RK3 of Delong et al. (2013) (2). RK3 is stable for centred
advection and third order for the deterministic part. With `speed_type = 0`
(s ≡ 1), ρ̃ is not needed and is skipped during the step. It is still
computed from the current φ for the `rhot` plotfile component.

**Time step:** dt = cfl/(v₀(1/dx + 1/dy) + 2D/dθ²). This is the positivity
limit of upwind with forward Euler; v₀ bounds the speed because s ≤ 1. The
θ-diffusion term grows like n_θ², so with fine θ grids an implicit θ
diffusion would be the next improvement.

**Particles** (2D build; `FlockPC`):
- **Step:** Euler–Maruyama, with the position advanced using the speed and
  heading at the start of the step, and θ += √(2D dt) Z.
- **ρ̃:** from particle-mesh (`density_method = 0`) or the pair sum (1).

**Comparison:** the particle histogram over θ bins on the SPDE grid is
compared cell for cell with the SPDE by `python/flocking_compare.py`.

## 6. Parameter regimes

For MIPS:
- **Box size:** L ≫ ℓ₁ = v₁/D, so the box holds many persistence lengths.
- **Speed function:** g < −1, for example λ > 1 for the exponential.
- **Kernel range:** R small enough that the unstable band contains several box
  modes, but at least 4 dx.
- **Particle count:** N large enough that a kernel area holds many particles,
  N·πR²/A ≫ 1.

The example `inputs_mips` uses:
- v₀ = 4, D = 10, λ = 1.5;
- a Gaussian kernel with R = 0.03;
- a 128² × 32 grid and N = 10⁶.

That gives v₁ ≈ 0.89, ℓ₁ ≈ 0.09 and box modes 1–4 unstable. Starting from
uniform, both codes form dense, slow domains about a third of the box apart.
By t ≈ 3 about 26% of the area has ρ̃ > 1.5ρ̄ and about 55% has ρ̃ < 0.5ρ̄.

## 7. Outlook: alignment

Flocking proper adds a turning rate to the heading equation:

  dθᵢ = ω(Xᵢ, θᵢ) dt + √(2D) dWᵢ

The Vicsek/Kuramoto form is ω = −K·Im( e^(−iθ) (W ⋆ P)(x) ), with the
complex polarization P = ∫ e^(iθ) φ dθ. In the SPDE this is a θ-advection
term −∂θ(ωφ). The kernel already has a slot for it in the θ flux, and ω
needs one more FFT convolution of the first θ-moment. Mean-field theory then
predicts:
- an isotropic-to-polar transition near Kρ̄/2 ≈ D, with the exact threshold
  depending on how K is normalized;
- travelling bands near the transition;
- giant number fluctuations in the ordered phase.

These are a natural next test of whether the SPDE reproduces finite-N
effects.

## References

- D. S. Dean, *Langevin equation for the density of a system of interacting
  Langevin processes*, J. Phys. A 29, L613 (1996).
- J. Tailleur and M. E. Cates, *Statistical mechanics of interacting
  run-and-tumble bacteria*, Phys. Rev. Lett. 100, 218103 (2008).
- M. E. Cates and J. Tailleur, *Motility-induced phase separation*, Annu. Rev.
  Condens. Matter Phys. 6, 219 (2015).
- Müller, von Renesse and Zimmer (2025), as cited in `src_deankow/flocking.pdf`,
  for the Dean–Kawasaki equation of this active-particle model.
