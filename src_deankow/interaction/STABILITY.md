# Interaction Potential Stability

The Dean-Kawasaki solver in `src_deankow/interaction` loses stability by two
unrelated mechanisms depending on the sign of the interaction kernel: an
attractive kernel demixes at `k = 0` once the mean-field strength exceeds the
diffusivity, while a repulsive one forms clusters at a finite wavelength near
`1.2 R` once the diffusivity drops below a much smaller threshold. At the
code's own parameters the two critical diffusivities differ by a factor of 34
in 2D and 63 in 3D.

## Linear stability

The deterministic part of the solver advances a density against a
self-consistent potential. `AdvancePhiAtLevel.cpp` builds `C` by FFT
convolution and `mykernel.H` adds its gradient to the diffusive flux:

$$\partial_t \phi = \nabla \cdot \left( D \nabla \phi + \phi \nabla C \right),
\qquad C = U \star \phi$$

Linearizing about a uniform state and transforming gives a growth rate that
depends on the kernel only through its Fourier transform:

$$\sigma(k) = -k^{2} \left[ D + \phi_{0} \, \hat{U}(k) \right]$$

The uniform state is unstable at wavenumber `k` whenever `phi0 * Uhat(k) < -D`.
Every threshold below is a statement about where `Uhat(k)` is negative and how
negative it gets. Here `phi0 = 1`: `AmrCoreAdv.cpp` renormalizes the initial
condition so `phi` integrates to 1 over the unit box, and `num_part` sets only
the noise amplitude and the particle count, not the mean-field strength.

The kernel is `U(r) = eps * exp(-(r/R)^alpha)` with `alpha = 3` and `R = 0.1`.
Its `k = 0` value has a closed form:

$$\hat{U}(0) = \int_{\mathbb{R}^{d}} U \, d^{d}r
= \epsilon \, \frac{2 \pi^{d/2}}{\Gamma(d/2)} \, \frac{R^{d}}{\alpha} \,
\Gamma\!\left( \frac{d}{\alpha} \right)$$

That evaluates to `0.028361 eps` in 2D and `0.004189 eps` in 3D. Whether
`Uhat` dips negative anywhere else is decided by `alpha` alone: for a bounded
kernel `exp(-(r/sigma)^alpha)`, `alpha <= 2` (the Gaussian core model) gives a
transform that is positive everywhere, while `alpha > 2` gives an oscillating
transform with negative lobes. This is the criterion Likos and co-workers use
to separate clustering from reentrant-melting behaviour, and `alpha = 3` puts
this code on the clustering side.

The first zero crossing sits at `kR = 4.10` in 2D and `4.58` in 3D, and the
lobe bottoms out at `kR = 5.03` and `5.38`. Its depth relative to `Uhat(0)` is
2.90 percent in 2D and 1.59 percent in 3D. The attractive case lives at `k = 0`;
the repulsive case lives entirely inside those few percent.

## Attractive kernel

This is the sign committed in `mykernel.H`: `U = -eps exp(-(r/R)^alpha)`, so the
kernel is a potential well everywhere and `Uhat(0)` is negative. The transform
is most negative at `k = 0`, so the uniform state goes unstable there first, as
soon as

$$\Lambda = \frac{\phi_{0} \, \epsilon \, |\hat{S}(0)|}{D} > 1$$

with `Shat(0)` the kernel integral above, `eps` divided out. Equivalently the
critical diffusivity is `phi0 * eps * |Shat(0)|`, which is 18.89 in 2D and 2.79
in 3D at `eps = 666`.

Because the unstable band runs all the way down to `k = 0`, no wavelength is
selected. Whatever structure appears early is transient: the system coarsens
without limit and ends as one blob whose size is set by the box rather than by
the kernel. This is ordinary mean-field demixing, and it is the regime the
repo's 2D runs sit in.

Growth still vanishes at `k = 0` like `k^2`, so the early pattern is set by
where `k^2 |Uhat(k)|` peaks:

| Case | Lambda | Fastest kR | Wavelength |
| --- | --- | --- | --- |
| 2D, eps = 666, D = 0.6 | 31.5 | 2.27 | 0.277 |
| 3D, eps = 666, D = 0.6 | 4.65 | 1.99 | 0.316 |
| 3D, eps = 4509, D = 0.6 | 31.5 | 2.35 | 0.267 |

The middle row is what happens if the 2D potential is carried into 3D
untouched: still unstable, but nearly seven times weaker and about ten times
slower to develop. The third row is the same potential recalibrated below.

## Repulsive kernel

This is the sign in the current working tree: `U = +eps exp(-(r/R)^alpha)`. Now
`Uhat(0)` is positive, so the effective diffusivity at long wavelength is
`D + 18.9` rather than `D`. Demixing is impossible at any diffusivity:
repulsion actively protects the uniform state against it. It would be easy to
stop there and conclude the repulsive kernel is unconditionally stable, and for
`alpha <= 2` that would be correct.

It is not correct at `alpha = 3`. The negative lobe is still there, and with
`U` positive it is the lobe rather than `k = 0` that carries the instability.
The condition is

$$D < D_{\mathrm{crit}} = \phi_{0} \, \epsilon \,
\left| \min_{k} \hat{S}(k) \right|$$

which evaluates to 0.5484 in 2D and 0.0445 in 3D at `eps = 666`. These are one
and a half orders of magnitude below the attractive thresholds, because the
lobe is only 2.90 percent as deep as `Uhat(0)` in 2D and 1.59 percent in 3D.

What this produces is qualitatively different from demixing. The unstable band
is a finite window around `kR` of about 5, bounded away from `k = 0`, so a
wavelength is genuinely selected: `1.26 R` in 2D and `1.17 R` in 3D. Particles
pile into overlapping clumps on a lattice with that spacing, and the spacing is
fixed by the kernel rather than by the density. This is the cluster-crystal
behaviour of the generalized exponential model, and it does not coarsen.

Two practical consequences follow from the lobe being shallow. The threshold is
demanding, so a repulsive run needs a diffusivity well below anything the repo
currently uses. And because `sigma` vanishes on both edges of the window,
growth is slow near threshold: the number of steps per e-folding is minimized
at `D = D_crit / 2`, and scales as `1 / D_crit^2`, so it falls quadratically as
`eps` rises.

## Where the repo's runs sit

All four thresholds at `eps = 666`, `phi0 = 1`, with the diffusivity each
inputs file actually ends on (later assignments in these files shadow earlier
ones):

| Sign | Dim | D_crit | inputs_fv / inputs_part, D = 0.46 | inputs_fv_3d_det, D = 0.6 |
| --- | --- | --- | --- | --- |
| Attractive | 2D | 18.89 | unstable, Lambda = 41 | unstable, Lambda = 31 |
| Attractive | 3D | 2.79 | unstable, Lambda = 6.1 | unstable, Lambda = 4.7 |
| Repulsive | 2D | 0.5484 | unstable, D/D_crit = 0.84 | stable, D/D_crit = 1.09 |
| Repulsive | 3D | 0.0445 | stable, D/D_crit = 10 | stable, D/D_crit = 13 |

Two entries are worth pausing on. The 2D repulsive row is close: at `D = 0.46`
the repo's 2D runs are inside the cluster regime, and at `D = 0.6` they are 9
percent outside it. A change to the diffusivity that looks inconsequential
crosses that boundary.

The 3D repulsive row is the one that will look like a broken run. At `D = 0.6`
it is thirteen times above threshold, so a 3D repulsive simulation started from
the current inputs will stay uniform indefinitely and give no indication why.

The sign in the working tree is the repulsive one, while the sign committed to
git is the attractive one, so which of these rows applies depends on an
uncommitted edit to `mykernel.H`.

## Calibrating from 2D to 3D

`eps` is a peak height, but nothing in the dynamics responds to a peak height.
Both thresholds are set by integrals of the kernel, and those shrink with
dimension at fixed `R`. Carrying `eps` across unchanged therefore weakens the
model, whichever sign is in use.

For the attractive case the invariant is `Uhat(0)`, and preserving it fixes the
ratio exactly:

$$\frac{\epsilon_{3D}}{\epsilon_{2D}}
= \frac{\Gamma(2/\alpha)}{2 R \, \Gamma(3/\alpha)} = 6.7706
\qquad (\alpha = 3, \; R = 0.1)$$

so `eps` goes from 666 to 4509. For the repulsive case the invariant is the
depth of the negative lobe instead, and preserving that needs `eps = 8208`, a
further factor of 1.8. The two choices disagree because the lobe is relatively
shallower in 3D than `Uhat(0)` is, so one number cannot hold both fixed.

| Target | eps in 3D | D_crit in 3D | D matching the 2D ratio |
| --- | --- | --- | --- |
| Preserve Uhat(0) | 4509 | 0.3013 | 0.253 |
| Preserve the lobe depth | 8208 | 0.5484 | 0.460 |

Which to pick is a modelling decision, not a correctness question. If `eps` and
`R` came from a physical pair potential, keep them and accept that 3D behaves
differently; that is the physics. If `eps` was tuned to put the 2D runs
somewhere interesting, rescale. For a repulsive study the second row is the
natural choice, since the lobe is the quantity of interest and it leaves the
diffusivity where the 2D runs had it.

Whatever is chosen has to be applied to the particle solver as well.
`AdvectParticles` uses `strength = interaction_strength / num_part` with
`U = strength exp(-(r/interaction_scale)^3)`, and moves particles down the
gradient of the sum over neighbours, so `num_part` times the pair potential is
exactly the mesh kernel. The defaults `interaction_strength = 666` and
`interaction_scale = 0.1` are deliberately the same pair as the mesh.
Rescaling one without the other makes the SPDE and particle regions model
different systems.

## HK kernel

`amr.ip_type = hk` selects a truncated paraboloid with compact support:

$$U(r) = \frac{\epsilon R}{2} \left( 1 - \frac{r^{2}}{R^{2}} \right) \quad (r < R),
\qquad U(r) = 0 \quad (r \ge R)$$

The growth rate `sigma(k) = -k^2 [D + phi0 Uhat(k)]` is unchanged; only the
transform differs. Writing `Uhat(k) = eps R^(d+1) Shat(kR)`, the transform of
`(1 - r^2)_+ / 2` at unit radius is a Bessel function:

$$\hat{S}(q) = (2\pi)^{d/2} \, \frac{J_{d/2+1}(q)}{q^{d/2+1}},
\qquad \hat{S}(0) = \frac{\pi}{4} \; (2D), \quad \frac{4\pi}{15} \; (3D)$$

Like the GEM kernel with `alpha > 2`, it oscillates in sign, here because of
the compact support rather than the steepness of the core. The minimum sits at
the first zero of `J_{d/2+2}`:

| Quantity | 2D | 3D |
| --- | --- | --- |
| Shat(0) | 0.7854 | 0.8378 |
| first zero, kR | 5.136 | 5.763 |
| minimum, kR | 6.380 | 6.988 |
| Shat_min | -0.04604 | -0.03444 |
| lobe depth / Shat(0) | 5.86 percent | 4.11 percent |
| selected wavelength | 0.985 R | 0.899 R |

The lobe is roughly twice as deep relative to `Uhat(0)` as the GEM lobe, so the
two thresholds sit closer together: a factor of 17 in 2D and 24 in 3D, against
34 and 63 for GEM.

**Attractive (`eps < 0`).** The instability is at `k = 0` and sets in once

$$\Lambda = \frac{\phi_{0} \, |\epsilon| \, R^{d+1} \, \hat{S}(0)}{D} > 1$$

As with GEM, no wavelength is selected and the system coarsens.

**Repulsive (`eps > 0`).** The instability is carried by the lobe:

$$D < D_{\mathrm{crit}} = \phi_{0} \, \epsilon \, R^{d+1} \, |\hat{S}_{\min}|
= 0.04604 \, \phi_{0} \epsilon R^{3} \; (2D), \quad
0.03444 \, \phi_{0} \epsilon R^{4} \; (3D)$$

The unstable band is bounded away from `k = 0`, and clusters form on a lattice
with spacing just under `R`.

**2D to 3D.** Preserving `Uhat(0)` needs `eps_3D / eps_2D = 15 / (16 R)`
(9.375 at `R = 0.1`). Preserving the lobe depth needs
`0.04604 / (0.03444 R)` (13.37 at `R = 0.1`).

**Numerics.**
- The slope of `U` jumps from `-eps` to 0 at `r = R`, so `Uhat` decays only
  like `k^{-(d+3)/2}` and aliases more than GEM's. The selected wavelength is
  about `R`, so keep `R/dx` at 8 or more. The code warns below 8.
- The slowly decaying positive tail also makes `int_stiff`, the interaction
  term in `EstTimeStep`, relatively larger than for GEM at the same `Uhat(0)`.
- The particle force is a constant-stiffness spring, `-(eps / R)(x_j - x_i)`
  divided by `num_part`, cut off exactly at `R`. When `amr.ip_range` is not
  set, it defaults to `ip_R`.

At startup, `PrintUhatMinMax` prints these continuum values next to the
discrete minimum of `Uhat` over `k != 0`. When the two disagree by more than a
few percent, the grid is too coarse.

## Piecewise parabolic kernel

`amr.ip_type = pp` is the piecewise parabolic potential of Gerber et al.
(arXiv:2510.17629), written in the same convention as HK, `U = -eps R w(r/R)`:

$$w(s) = \begin{cases}
\alpha (s^{2} - a^{2})/2 + \beta (a^{2} - 1)/2 & s \le a \\
\beta (s^{2} - 1)/2 & a < s \le 1 \\
0 & s > 1
\end{cases}$$

with `alpha = amr.ip_pp_alpha`, `beta = amr.ip_pp_beta` and `a = amr.ip_pp_a`.
The paper's `W = gamma ell w(x/ell)` maps to `ip_R = ell`, `ip_eps = -gamma`.
`alpha = beta = 1` is exactly HK (bit for bit), whatever `a` is.

It is the sum of two HK paraboloids, a tail of radius `R` and a core of
radius `aR`, so its transform follows from the HK one:

$$\hat{U}(k) = \epsilon R^{d+1} \left[ \beta \, \hat{S}(kR)
+ (\alpha - \beta) \, a^{d+2} \, \hat{S}(k a R) \right]$$

with `Shat` the HK transform above. `PrintContinuumStability` evaluates this
and scans the box wavenumbers for the minimum and the fastest growing mode.

**Why it exists.** A single paraboloid of range `R` forms clusters about
`2.06 R` apart in 2D, beyond its own range, so neighbouring clusters do not
attract. They merge only when their random walks close the gap, at a rate
`D / (N m)` that slows with the particle count. With a strong core and a weak,
longer tail, the core sets the initial pattern while the tail pulls
neighbouring clusters together, so mergers happen by drift on a time
`~ R / (gamma beta m)` that does not depend on `N`.

**Choosing parameters.** Three conditions have to hold together:

- The core, not the tail, must carry the fastest growing mode. The tail term
  is largest at the box scale, so `beta` has a narrow window: too large and the
  uniform state collapses into one box-sized cluster.
- The core-driven cluster spacing must be less than `R`, so that neighbours
  sit inside each other's tail.
- The clusters must be resolved: `sigma = sqrt(R/(gamma alpha)) >= 2 dx`, the
  core `a R >= 4 dx` and `R >= 8 dx` (the code warns for the last two).

The first two need `gamma` well above the core's own HK threshold
`~ 4 / (pi (a R)^3)`, while the third caps `gamma` from above, so in 2D the core cannot be smaller than about `1.5 sqrt(dx)`. On
`128^2` there is no usable window; `256^2` is the smallest grid that works.

**Example** (`exec/dean_kow/interaction/inputs_fv_pp_merge`): `256^2`,
`gamma = 7300`, `R = 0.45`, `a = 0.27`, `alpha = 1`, `beta = 0.02`, `D = 1`,
5 particles per cell. The fastest mode is (1,3) at growth rate 299, giving
about 10 clusters 0.32 apart, inside the tail. In one run (seed 11) there were
10 clusters at `t = 0.019`, 3 at 0.035 and a single cluster at 0.038. With
`beta = 0.03` the fastest mode is already the box mode (1,0).

| beta | fastest mode | clusters | tail pull e-folding time |
| --- | --- | --- | --- |
| 0.01 | (3,1) | ~10 | 0.06 |
| 0.02 | (3,1) | ~10 | 0.03 |
| 0.03 | (1,0) | 1 | - |

## Generalized Morse kernel

`amr.ip_type = morse` selects the sum of an attractive and a repulsive
exponential:

$$U(r) = -\epsilon_{a} e^{-(r - r_{e})/R_{a}} + \epsilon_{r} e^{-(r - r_{e})/R_{r}}
= -A e^{-r/R_{a}} + B e^{-r/R_{r}}$$

with `A = eps_att e^(re/R_att)` and `B = eps_rep e^(re/R_rep)`. The inputs are
`amr.ip_eps_att`, `ip_eps_rep`, `ip_R_att`, `ip_R_rep` and `ip_re` (default
0). The code requires `eps_att > eps_rep >= 0` and `R_att > R_rep > 0`, so the
attractive part is both stronger and longer ranged, and `U < 0` at large `r`.

**Transform.** The d-dimensional transform of `e^(-r/a)` is closed-form:

$$\hat{U}(k) = 2\pi \left[ -\frac{A R_{a}^{2}}{(1 + k^{2}R_{a}^{2})^{3/2}}
+ \frac{B R_{r}^{2}}{(1 + k^{2}R_{r}^{2})^{3/2}} \right] \; (2D), \qquad
8\pi \left[ -\frac{A R_{a}^{3}}{(1 + k^{2}R_{a}^{2})^{2}}
+ \frac{B R_{r}^{3}}{(1 + k^{2}R_{r}^{2})^{2}} \right] \; (3D)$$

so `Uhat(0) = 2 pi (B R_rep^2 - A R_att^2)` in 2D and
`8 pi (B R_rep^3 - A R_att^3)` in 3D.

**There is no finite-wavelength instability.** Write
`Uhat = g_r(k) [B - A (R_att/R_rep)^d h(k)]`, where `g_r` is the repulsive
transform per unit amplitude and `h(k) = g_a(k)/g_a(0) / (g_r(k)/g_r(0))`.
Both `g_r` and `h` decrease in `k`, because `R_att > R_rep`. Wherever the
bracket is negative, both terms of the derivative are positive, so `Uhat` is
increasing wherever it is negative. Two cases follow:

- `A R_att^d > B R_rep^d`: `Uhat(0) < 0` and the minimum is at `k = 0`. This is
  mean-field demixing, as for the attractive GEM kernel, with no selected
  wavelength.
- `A R_att^d <= B R_rep^d`: `Uhat > 0` for every `k`, and the uniform state is
  stable at any diffusivity.

With `re = 0` the input constraints force the first case. A positive `re`
multiplies `B/A` by `e^(re (1/R_rep - 1/R_att))` and can reach the second.

**Thresholds.** In an infinite domain the uniform state is unstable once
`D < D_crit = -phi0 Uhat(0)`. In the periodic box `k = 0` is the conserved
mass, so the lowest mode that can grow is `k1 = 2 pi / L`, and the box
threshold is

$$D < D_{\mathrm{crit}}^{\mathrm{box}} = -\phi_{0} \, \hat{U}(2\pi/L)$$

Because the minimum is at `k = 0`, this box threshold is always the smaller of
the two, and noticeably so when `R_att` is not small against `L`. Example: in
2D with `eps_att = 2, eps_rep = 1, R_att = 0.04, R_rep = 0.02, L = 1`, the
infinite-domain `D_crit = 0.01759`. The box value is `0.01589`, which the
discrete `Uhat` at 256^2 reproduces to 5e-6 relative. `PrintUhatMinMax`
prints both values.

**Range.** The tail decays like `e^(-r/R_att)` and never vanishes, so there is
no default `ip_range`.
- The truncation check compares `|U(ip_range)|` with `max |U|` on
  `[0, ip_range]`, because `U(0) = B - A` can be small or zero.
- The code prints a suggested cutoff,
  `re + R_att ln(1e5 eps_att / max|U|)`.
- The particle solver needs `ip_range < L/2`. The mesh kernel is wrapped at
  `L/2`, and the code warns when `|U(L/2)| / max|U| > 1e-4`.

In practice this means `R_att` must be well below `L / (2 ln(1e4 ...))`. With
`R_att = 0.1` in a unit box, the tail at `L/2` is still 1.3 percent of `max|U|`
and periodic images interact. `R_att = 0.04` or less keeps it negligible.

**Force at r = 0.** `dU/dr = eps_att/R_att e^(...) - eps_rep/R_rep e^(...)` is
generally nonzero at `r = 0`: the kernel has a cusp. So `(dU/dr)/r` diverges,
and the pair force is set to zero for `r < 1e-6` (`ip_morse_rmin` in
`Potentials.H`).

**Resolution.** Because of the cusp at the origin, `Uhat` decays only like
`k^-(d+1)`. Keep `R_rep/dx` at 4 or more; the code warns below that.

## Numerical limits

`EstTimeStep` in `AmrCoreAdv.cpp` sets `dt = cfl * 2 / lambda`, the
fraction `cfl` of the forward Euler limit for the largest eigenvalue `lambda`:

$$\lambda = D \sum_d \frac{4}{dx_d^{2}} + \max(\phi) \, S_{\mathrm{int}}
\; \left[ + \; 4 \sum_d \frac{\max|w_d|}{dx_d} \right]$$

- The first term is the explicit diffusion limit.
- `S_int` (`int_stiff`) is `max_k keff^2 max(Uhat(k) cellvol, 0)`, the interaction
  operator linearized about a uniform state, scaled by the largest density. It
  is computed from the discrete `Uhat` of whichever kernel is in use, but it
  only sees the positive part of `Uhat`: for an attractive kernel with
  `Uhat <= 0` at every k (attractive gema with `alpha <= 2`, some morse
  settings) it is zero. It also does not bound the drift in a clustered state.
  As clusters form, `max(phi)` grows like `gamma alpha m^2 / (2 pi R)` for a
  pp or hk cluster of mass m, and dt falls with it.
- The bracketed term is added only with `drift_flux_type = 1`. With
  `w = grad(C + V_ext)` the face drift velocity, the Scharfetter-Gummel update
  keeps `phi >= 0` when `dt sum_d (2D/dx_d^2 + 2|w_d|/dx_d) <= 1`; the Bernoulli
  weights on a cell's two faces add to at most `2 + |Pe_lo| + |Pe_hi|`, the worst
  case being flow out of both faces. With `cfl <= 1` the factor 4 guarantees it.

`max|w_d|` is measured over all faces in the previous step (`MeasureDrift`,
which forms `w` exactly as the flux kernels do, external potential included),
so dt lags the drift by one step. Before the first measurement (the first step,
or the first step after a restart) the interaction drift is bounded instead by
`max|dU/dr| * sum|phi| cellvol`.

Two further controls, both off by default:

| Input | Effect |
| --- | --- |
| `drift_cfl = c` | Also limit `dt * sum_d max|w_d|/dx_d <= c`, for either drift flux |
| `drift_diag_int = n` | Every n steps print `max|w_d|`, the drift CFL number `dt sum_d max|w_d|/dx_d` and the largest cell Peclet number `max|w_d| dx_d / D` |

The centered drift flux (`drift_flux_type = 0`) is only monotone for cell Peclet
number `<= 2`; the diagnostic flags larger values. In the 256^2 pp merger
(`inputs_fv_pp_merge`) it exceeds 2 at the core edge of the merged cluster, and
in the 128^2 hk cluster run (`inputs_fv_hk_cluster` at 128^2) it reaches about
2.5. Neither run is limited by the drift CFL number, which stays well below 1.

Resolution matters too, and only for the repulsive case. The selected
wavelength is `1.17 R`, which at `n_cell = 64` with `R = 0.1` is 7.3 cells: too
coarse for a cluster lattice. 128 cells per side puts `R` at 12.8 cells and the
wavelength at about 15. The kernel itself is fine at 64, since `Uhat` is
unchanged between 64 and 128; it is the emergent pattern that is not.

All interaction parameters are read from the inputs file by
`read_potential_params` in `Potentials.H`, and the mesh and particle solvers
share them:

| Parameter | Location |
| --- | --- |
| amr.ip_type, amr.ip_eps, amr.ip_R, amr.ip_alpha | inputs file (`Potentials.H`); gema, hk, pp |
| amr.ip_pp_alpha, ip_pp_beta, ip_pp_a | inputs file (`Potentials.H`); pp |
| amr.ip_eps_att, ip_eps_rep, ip_R_att, ip_R_rep, ip_re | inputs file (`Potentials.H`); morse |
| amr.ip_range (particle cutoff) | inputs file; defaults to ip_R for hk and pp, required for morse particles |
| diff_coeff, cfl | the inputs files |

For gema, hk and pp, the sign of the kernel is the sign of `amr.ip_eps`: positive
is repulsive, negative is attractive.

All figures here are linear, mean-field and deterministic, computed with
`dorand = 0`. With noise on, fluctuations smear each transition and produce
precursor structure somewhat above threshold, so the thresholds are where
growth changes sign, not where structure first becomes visible.
