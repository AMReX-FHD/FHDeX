# src_deankow/interaction

Dean–Kawasaki solver for diffusing particles with a pairwise interaction
potential and an optional external potential. The density φ is advanced either
as a fluctuating SPDE on the mesh or, with `amr.use_particles = 1`, with
particles on a refined level. Both use the same potential parameters, read from
the inputs file by `read_potential_params` in `Potentials.H`.

The mesh solver advances

  ∂φ/∂t = ∇·( D ∇φ + φ ∇(C + V_ext) ) + noise,  C = U ⋆ φ

where φ is normalized to integrate to 1 over the domain, `diff_coeff` = D,
U(r) is the interaction potential and V_ext the external potential.
`num_part` sets the number of particles N, which sets the noise amplitude and
the particle count; it doesn't change the mean-field strength. The convolution C
is computed by FFT on level 0.

Example inputs are in `exec/dean_kow/interaction`. `STABILITY.md` covers the
linear stability of each potential, choosing parameters, and the time step.

## Interaction potentials

The interaction is switched on with `amr.use_int_pot = 1` (the default) and
chosen with `amr.ip_type`. The force on particle i from particle j is
−∇U(|x_i − x_j|), weighted by 1/N.

For `gema`, `hk` and `pp` the sign of the potential is the sign of `amr.ip_eps`:
**positive is repulsive, negative is attractive.**

| `amr.ip_type` | U(r) | Range |
| --- | --- | --- |
| `gema` (default) | ε·exp(−(r/R)^α) | Unbounded (needs a cutoff for particles) |
| `hk` | ½·ε·R·(1 − (r/R)²) for r < R, 0 beyond | Compact, radius R |
| `morse` | −ε_att·exp(−(r − r_e)/R_att) + ε_rep·exp(−(r − r_e)/R_rep) | Unbounded |
| `pp` | −ε·R·w(r/R), with w as below | Compact, radius R |

### gema: generalized exponential model

U(r) = `ip_eps` · exp(−(r/`ip_R`)^`ip_alpha`)

A bounded, soft kernel. α ≤ 2 (α = 2 is the Gaussian core model) has a Fourier
transform that is positive everywhere. α > 2 has negative lobes, so a repulsive
kernel can form a cluster crystal (see `STABILITY.md`).

### hk: Hegselmann–Krause paraboloid

U(r) = ½ · `ip_eps` · `ip_R` · (1 − (r/`ip_R`)²) for r < `ip_R`, 0 for r ≥ `ip_R`

This is the prototype potential of Gerber, Gvalani, Hairer, Pavliotis and
Schlichting (arXiv:2510.17629), W = γℓ·½((x/ℓ)² − 1), with `ip_R` = ℓ and
`ip_eps` = −γ. Inside its range the pair force is a constant-stiffness spring,
(dU/dr)/r = −ε/R. In 2D, clusters form about 2.06·R apart, beyond the range of
the potential.

### morse: generalized Morse

U(r) = −`ip_eps_att`·exp(−(r − `ip_re`)/`ip_R_att`) + `ip_eps_rep`·exp(−(r − `ip_re`)/`ip_R_rep`)

A short-range repulsive core with a longer-range attractive tail, so there is
a potential well. The inputs must satisfy `ip_R_att` > `ip_R_rep` > 0 and
`ip_eps_att` > `ip_eps_rep` ≥ 0. The pair force is set to zero below
r = 10⁻⁶ to avoid round-off. At startup the code prints a suggested `ip_range`.

### pp: piecewise parabolic

U(r) = −`ip_eps` · `ip_R` · w(r/`ip_R`), with a = `ip_pp_a`:

- w(s) = α(s² − a²)/2 + β(a² − 1)/2 for s ≤ a
- w(s) = β(s² − 1)/2 for a < s ≤ 1
- w(s) = 0 for s > 1

with α = `ip_pp_alpha` (core stiffness) and β = `ip_pp_beta` (tail stiffness).
It is the piecewise parabolic potential of the same paper, again with
`ip_R` = ℓ and `ip_eps` = −γ. α = β = 1 is exactly `hk`, bit for bit.

With a strong core and a weak, longer tail, the core sets the cluster pattern
while the tail pulls neighbouring clusters together, so clusters merge by drift
instead of by slow random walk. `exec/dean_kow/interaction/inputs_fv_pp_merge`
is an example, and `STABILITY.md` explains how to choose the parameters.

### Interaction potential inputs

All use the `amr.` prefix.

| Input | Default | Used by | Meaning |
| --- | --- | --- | --- |
| `use_int_pot` | 1 | all | 1 = interaction on, 0 = off |
| `ip_type` | `gema` | all | `gema`, `hk`, `morse` or `pp` (case-insensitive) |
| `ip_eps` | 666 | gema, hk, pp | Strength. Positive is repulsive, negative attractive. |
| `ip_R` | 0.1 | gema, hk, pp | Length scale: kernel radius for hk and pp. Must be > 0. |
| `ip_alpha` | 3 | gema | Exponent α. Must be > 0. |
| `ip_eps_att` | 1 | morse | Attractive amplitude. Must be > `ip_eps_rep`. |
| `ip_eps_rep` | 0.5 | morse | Repulsive amplitude. Must be ≥ 0. |
| `ip_R_att` | 0.1 | morse | Attractive decay length. Must be > `ip_R_rep`. |
| `ip_R_rep` | 0.05 | morse | Repulsive decay length. Must be > 0. |
| `ip_re` | 0 | morse | Offset r_e. Must be ≥ 0. |
| `ip_pp_alpha` | 1 | pp | Core stiffness α. Must be > 0. |
| `ip_pp_beta` | 1 | pp | Tail stiffness β. Must be > 0. |
| `ip_pp_a` | 0.5 | pp | Core radius as a fraction of `ip_R`. Must be in (0, 1). |
| `ip_range` | 0 | particles | Particle pair cutoff; ≤ 0 means none. Defaults to `ip_R` for hk and pp. Required (> 0) for gema and morse when `use_particles = 1`. The mesh solver always uses the untruncated kernel. |

At startup the code prints the potential's linear stability and resolution
warnings. It stops if `ip_range` truncates the potential by more than 10⁻⁴ of
its maximum. The old names `interaction_strength`, `interaction_scale` and
`interaction_range` are rejected.

## External potential

Switched on with `amr.use_ext_pot = 1` (default 0). The drift it adds is

- x: dV/dx = 2(x − `ep_beta`)(x − `ep_alpha`)(2x − `ep_alpha` − `ep_beta`)/`ep_gamma`,
  the double well V = (x − α)²(x − β)²/γ;
- y: dV/dy = 2(y − ½)³/`ep_gamma`, the quartic well V = (y − ½)⁴/(2γ).

In 3D the z direction reuses the y form; this is a placeholder (see the
`TODO(3D)` in `mykernel.H`).

| Input | Default | Meaning |
| --- | --- | --- |
| `amr.use_ext_pot` | 0 | 1 = external potential on |
| `amr.ep_alpha` | 0.25 | First well position in x |
| `amr.ep_beta` | 0.75 | Second well position in x |
| `amr.ep_gamma` | 1.1e-3 | Well depth scale; smaller is deeper. Must be > 0. |

## Discretization options

These inputs have no prefix, like `diff_coeff`. All defaults reproduce the
original scheme bit for bit.

### Noise modification

| Input | Values | Default | Effect |
| --- | --- | --- | --- |
| `noise_avg_type` | 0 or 1 | 0 | How the noise amplitude at a face is formed from the two adjacent densities. **0:** average of √φ (original). **1:** √(arithmetic average of φ) × min(n_L, 1) × min(n_R, 1), with n = `num_part`·φ·cell volume the particle count in each neighbouring cell. The noise turns off smoothly next to cells holding less than one particle and vanishes next to empty ones (as reactDiff `avg_type = 1`). |

It uses `num_part` and the cell size, so nothing else needs setting.

### Scharfetter–Gummel drift flux

| Input | Values | Default | Effect |
| --- | --- | --- | --- |
| `drift_flux_type` | 0 or 1 | 0 | Discretization of the diffusion plus drift flux, where the drift w = ∇(C + V_ext) includes both the interaction and the external potential. **0:** centred, D(φ_R − φ_L)/dx + w(φ_L + φ_R)/2 (original). **1:** Scharfetter–Gummel, (D/dx)·[B(−Pe)·φ_R − B(Pe)·φ_L], with B(z) = z/(e^z − 1) and Pe = w·dx/D. It keeps φ ≥ 0 under the time step it imposes, and reduces to the centred flux at small Pe. It requires `diff_coeff` > 0. |

`drift_flux_type = 1` also tightens the time step automatically: `EstTimeStep`
adds 4·Σ_d max|w_d|/dx_d to λ, with w measured over the faces in the previous
step (a worst-case bound is used before the first measurement).

### Related diagnostics and limits

| Input | Default | Effect |
| --- | --- | --- |
| `neg_diag_int` | 0 (off) | Every n steps print min φ, the fraction of negative cells and the negative mass |
| `drift_diag_int` | 0 (off) | Every n steps print max\|w\| per direction, the drift CFL number dt·Σ_d max\|w_d\|/dx_d and the largest cell Péclet number max\|w_d\|·dx_d/D. Péclet numbers above 2 are flagged, because the centred flux is not monotone there. |
| `drift_cfl` | 0 (off) | If > 0, also limit dt·Σ_d max\|w_d\|/dx_d ≤ `drift_cfl`, for either drift flux |

Example, both fixes with diagnostics:

```
noise_avg_type  = 1
drift_flux_type = 1
neg_diag_int    = 1000
drift_diag_int  = 1000
```

Negative densities come mainly from the noise when there are few particles per
cell. Raising `num_part` reduces them more than either option does.

## Postprocessing

`python/stack_plotfiles_time.py` stacks a series of 2D plotfiles into a 3D
space–time plotfile. `python/cluster_diagnostics.py` finds and tracks
clusters, logs mergers, and computes radial distribution functions. Each has a
README in `python/`.
