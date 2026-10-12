# Interacting Brownian particles: formation, merging and splitting of clusters

This note describes the model solved in `src_deankow/interaction` and gives
estimates of four things, for each interaction potential in the code:
1. the linear stability of the uniform state and the time for clusters to
   form;
2. how many clusters form initially;
3. the time scales for clusters to merge;
4. the time scales for clusters to split.

The potentials are generalized exponential (`gema`), Hegselmann–Krause
(`hk`), generalized Morse (`morse`) and piecewise parabolic (`pp`). For `gema`
and `hk`, both attractive and repulsive cases are treated.

`README.md` lists the inputs, and `STABILITY.md` has more detail on the
linear stability and the time step. All numbers below are for 2D and the
unit periodic box. They come from the representative parameter sets in §2.3.

## 1. Model

### 1.1 Particles

N Brownian particles at positions Xᵢ in the periodic box interact through a
pair potential U, weighted by 1/N:

  dXᵢ = −(1/N) Σⱼ ∇U(|Xᵢ − Xⱼ|) dt + √(2D) dWᵢ

The mobility is 1, so D plays the role of temperature. A particle at x feels
the potential C(x) = (1/N) Σⱼ U(|x − Xⱼ|).

### 1.2 Mean field and Dean–Kawasaki SPDE

With the empirical density φ = (1/N) Σᵢ δ(x − Xᵢ), normalized to ∫φ = 1, so
that the mean is μ₀ = 1/A = 1,

  ∂ₜφ = ∇·( D∇φ + φ∇C ) + N^(−1/2) ∇·( √(2Dφ) ξ ),  C = U ⋆ φ

This is the McKean–Vlasov equation plus Dean–Kawasaki noise.
- **Free energy:** without noise it is the gradient flow of
  F[φ] = D ∫ φ ln φ + ½ ∫∫ U(x − y) φ(x) φ(y).
- **Steady states:** the stationary states are fixed points of
  φ = e^(−C/D)/Z.
- **Finite N:** N enters only through the noise. The mean-field strength is
  set by U, D and μ₀, not by N.

In the code, C is computed by FFT (`AdvancePhiAtLevel.cpp`), and the flux
D∇φ + φ∇C is discretized on faces (`mykernel.H`). The particle version
(`StochasticPC`) sums the pair forces within the cutoff `amr.ip_range`.

### 1.3 Potentials

| `ip_type` | U(r) | Range | U(0) | Attractive when |
| --- | --- | --- | --- | --- |
| `gema` | ε e^(−(r/R)^α) | unbounded, fast decay | ε | ε < 0 |
| `hk` | ½ε R (1 − r²/R²) for r < R | R | ½εR | ε < 0 |
| `morse` | −ε_att e^(−(r−r_e)/R_att) + ε_rep e^(−(r−r_e)/R_rep) | unbounded, exponential | −ε_att + ε_rep (at r_e = 0) | always: the code requires ε_att > ε_rep and R_att > R_rep |
| `pp` | −ε R w(r/R), with a parabolic core (stiffness α, radius aR) and tail (stiffness β) | R | ½εR[α a² + β(1 − a²)] | ε < 0 |

For `hk` and `pp`, ε = −γ maps onto the Hegselmann–Krause and
piecewise-parabolic potentials of Gerber et al. (arXiv:2510.17629). Their
γℓ·w(x/ℓ) is the code's U, with ℓ = R.

## 2. General results

### 2.1 Linear stability and the formation time

Linearizing about φ = μ₀ gives, for each Fourier mode k,

  σ(k) = −k² ( D + μ₀ Û(k) )

where Û(k) = 2π ∫ U(r) J₀(kr) r dr is the 2D Fourier transform. The uniform
state is unstable at k when μ₀Û(k) < −D. On the periodic box only the
lattice modes k = 2π(n_x, n_y) exist, and k = 0 (total mass) is excluded.
The code prints these quantities at startup:
- the continuum and discrete Û;
- the discrete minimum Λ = −μ₀ Û_min/D, with instability when Λ > 1;
- the fastest box mode and its growth rate.

There are two kinds of instability:
- **Attractive kernels** (Û(0) < 0): the band reaches down to the smallest box
  mode. Growth peaks where k²(−D − μ₀Û(k)) is largest, which is a finite k set
  by the kernel's shape. The system clusters and then coarsens.
- **Repulsive kernels** (Û(0) > 0): only kernels whose Û has a negative lobe
  can go unstable, at the lobe's k*, and only when D < D_crit = μ₀|Û_min|.
  The result is a cluster crystal, a lattice of clusters with spacing set by
  k* (Likos et al.).

**Formation time.** The noise seeds each mode with relative amplitude about
N^(−1/2). Nonlinear clusters appear once the fastest mode has grown to O(1):

  t_form ≈ ln(√N) / σ_max = (ln N)/(2σ_max)

This holds up to an O(1)/σ_max correction for the initial transient and the
saturation level. The time grows only logarithmically with N.

### 2.2 Number of clusters

The fastest-growing mode sets the initial spacing λ = 2π/|k_max|. Its
nonlinear pattern breaks into blobs on a lattice:

  n_square = (|k_max| L/2π)²,  n_hex = (√3/8π²) A |k_max|²  ≈ 0.87 n_square

For a hexagonal lattice whose reciprocal vectors have magnitude |k_max|, the
nearest-neighbour spacing is a = 4π/(√3|k_max|). In a small box the lattice
mode that wins is a specific (n_x, n_y), often a square or rectangular one,
so n_square is usually closer. With several modes of nearly equal σ, the
pattern is mixed and the count is uncertain by about ±20%.

The cluster mass is about m ≈ 1/n. The cluster width follows from the
self-consistent well C(r) ≈ m U(r) around the cluster:
- **Harmonic core** (`hk`, `pp`, `morse` with ε_att/R_att = ε_rep/R_rep):
  σ_c² = D/(m U″(0)).
- **Flat core** (`gema` with α > 2, where U ≈ U(0) + |ε|(r/R)^α):
  r_c ≈ R (D/(m|ε|))^(1/α).

### 2.3 Representative parameter sets

| Case | Parameters | D | N | Source |
| --- | --- | --- | --- | --- |
| gema attractive | ε = −666, R = 0.1, α = 3 | 0.5 | 81920 | `inputs_fv` with the sign flipped |
| gema repulsive | ε = +666, R = 0.1, α = 3 | 0.25 (D_crit = 0.55) | 81920 | `inputs_fv`, lower D |
| hk attractive | ε = −3000 (γ = 3000), R = 0.1 | 1 | 81920 | `inputs_fv_hk_cluster` (128²) |
| hk repulsive | ε = +2·10⁴, R = 0.1 | 0.5 (D_crit = 0.92) | 81920 | — |
| morse | ε_att = 200, ε_rep = 100, R_att = 0.1, R_rep = 0.05, r_e = 0 | 0.5 | 81920 | — |
| pp | ε = −7300, R = 0.45, a = 0.27, α = 1, β = 0.02 | 1 | 327680 | `inputs_fv_pp_merge` |

### 2.4 Merging mechanisms

Three mechanisms merge clusters, with very different time scales.

**(M1) Drift.** Clusters within each other's interaction range attract. Each
moves toward the other at speed m|U′(d)| (unit mobility), so a pair at
distance d₀ merges in about

  t_drift ≈ ∫ dd / (2m|U′(d)|)

taken from the core size up to d₀. This is independent of N. It is the fast
route when the initial spacing λ lies inside the range: `pp` with a weak
tail, `morse`, and `gema` if λ is under about 1.5R.

**(M2) Brownian coalescence.** Clusters beyond each other's range move as
Brownian particles of diffusivity D/(N m). They merge when random motion
closes the gap g (at most λ):

  t_brown ≈ g² N m / (8D)

This grows linearly with N, which is the coalescence time scale of
Gerber et al.

**(M3) Mass exchange (Ostwald-type ripening).** Particles escape over a
barrier about equal to the well depth m|U(0)|, at a per-particle rate of
roughly (D/σ_c²)·e^(−m|U(0)|/D). Smaller clusters have shallower wells, lose
mass faster, and eventually vanish:

  t_exch ≈ (σ_c²/D) · e^(m|U(0)|/D)

This is independent of N and present in the mean-field PDE; it is the
dynamical metastability of Gerber et al. Differences in mass enter through
the exponent, e^(Δm|U(0)|/D).

### 2.5 Splitting

- **Attractive kernels.** Splitting a cluster of mass m into two halves out of
  each other's range costs a total interaction energy N(m/2)²|U(0)|. The
  rate is about e^(−N m²|U(0)|/4D): exponentially small in N, and absent from
  the mean-field PDE, which is a gradient flow.
  - Genuine splitting is therefore not observed at these N.
  - What looks like splitting during formation is the break-up of stripes or
    elongated blobs into round clusters. That is part of the nonlinear stage
    of the instability, on a time of about 1/σ_max.
- **Repulsive (cluster-crystal) kernels.** The number of clusters is fixed by
  the lattice spacing 2π/k*, and the occupancy by N/n_clusters. Occupancies
  change by single-particle hops between neighbouring clusters, over the
  lattice barrier ΔC, at a rate of about (D/a²)e^(−ΔC/D). Splitting or merging
  a whole cluster means creating or annihilating a lattice site. That is a
  collective, activated rearrangement, much slower than hopping. It happens
  only if the initial lattice is badly mismatched to k*, for example when it
  is incommensurate with the box.

## 3. Results for each potential

### 3.1 Summary

| Case | σ_max | Fastest mode | t_form | Clusters (square / hex) | Spacing λ | m | Width | m\|U(0)\|/D |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gema attractive | 3633 | (3,2) | 0.0016 | 13 / 11 | 0.277 (2.8R) | 0.077 | r_c ≈ 0.021 | 102 |
| gema repulsive | 757 | (7,4) | 0.0075 | 65 / 56 | 0.124; hex a = 0.143 (1.43R) | 0.015 (≈1300 particles) | 0.008 | n/a |
| hk attractive | 265 | (3,1) | 0.021 | 10 / 9 | 0.316 (3.2R) | 0.10 | σ_c ≈ 0.018 | 15 |
| hk repulsive | 1750 | (10,3) | 0.0032 | 109 / 94 | 0.096; hex a = 0.111 (1.11R) | 0.009 (≈870 particles) | ring, radius 0.012 | n/a |
| morse | 283 | (1,1) | 0.020 | 2 | 0.707 (7.1R_att) | 0.5 | σ_c ≈ 0.007 | 100 |
| pp | 299 | (3,1) | 0.021 | 10 / 9 | 0.316 (0.70R, inside the tail) | 0.10 | σ_c ≈ 0.025 | 15 |

| Case | Drift merging (M1) | Brownian merging (M2) | Mass exchange (M3) | Splitting |
| --- | --- | --- | --- | --- |
| gema attractive | none: force at λ ≈ 7·10⁻⁶ | ≈ 120 | ∝ e^102: never | e^(−1.6·10⁵): never |
| gema repulsive | clusters repel | — | hops: ΔC/D = 17, τ_hop ≈ 3·10⁶ | lattice change: never |
| hk attractive | none: λ > R | ≈ 100 | ≈ 3·10⁻⁴ × e^15 ≈ 10³ | e^(−3·10⁴): never |
| hk repulsive | clusters repel | — | hops: ΔC/D = 12, τ_hop ≈ 4·10³ | lattice change: ≫ τ_hop |
| morse | ≈ 0.06 | ≈ 5·10³ | ∝ e^100: never | e^(−10⁶): never |
| pp | ≈ 0.015 | ≈ 400 | ≈ 6·10⁻⁴ × e^15 ≈ 2·10³ | e^(−1.2·10⁵): never |

For scale, in simulation steps: at 128², dt ≈ 2–4·10⁻⁶ for these runs, so
t = 100 is about 3·10⁷ steps. Only times up to about 0.1–1 are reachable
routinely.

### 3.2 GEM-α, attractive (ε < 0)

- **Stability.** Û(0) = 0.0284ε for α = 3, so Λ = 0.0284|ε|/D = 38 here. The
  band runs from the smallest box mode up to about 6·2π.
- **Formation.** It is very fast, t_form ≈ 0.0016. Growth peaks at
  kR ≈ 2.3 (λ ≈ 0.28), which does not depend on |ε| once Λ ≫ 1 (see
  `STABILITY.md`).
- **Number of clusters.** About 11–13 (λ ≈ 2.8R). This scales like 1/R² and
  depends only weakly on ε and D, through the band edge.
- **Cluster shape.** For α = 3 the core is flat (U″(0) = 0), so the clusters
  are compact blobs of radius r_c ≈ R(D/m|ε|)^(1/3) ≈ 0.02, with steep edges.
- **Merging.**
  - Drift: clusters 2.8R apart see a force of only m|U′(λ)| ≈ 7·10⁻⁶, because
    of the factor e^(−(2.8)³) ≈ 3·10⁻¹⁰. That is effectively none.
  - Mass exchange: the well depth m|ε|/D ≈ 100 rules it out.
  - What remains is Brownian coalescence, with t ≈ 120 at N = 81920, growing
    like N.

  So the initial cluster pattern is effectively frozen on simulation times.
  It coarsens to the single-cluster state (`STABILITY.md`, "ends as one
  blob") only over times ∝ N. Weaker attraction or higher D moves λ toward
  R, where drift merging becomes possible.
- **Splitting:** never.

### 3.3 GEM-α, repulsive (ε > 0)

- **Stability.** Û(0) > 0, so there is no long-wave instability. For α > 2
  the transform has a negative lobe: its minimum is −0.029·Û(0) at kR = 5.03.
  The uniform state is unstable for D < D_crit = μ₀|Û_min| = 0.55 at ε = 666.
- **Formation.** At D = 0.25, σ_max ≈ 760 at k ≈ 8·2π, so t_form ≈ 0.0075.
  Near D_crit the growth rate vanishes like (D_crit − D)·k*², and t_form
  diverges.
- **Number of clusters.** A hexagonal cluster crystal with spacing
  a ≈ 1.43R = 0.143 forms, about 56 clusters of about 1300 particles each in
  the unit box. The number is set by R, not by N or ε: n ≈ 0.22·(k*)²A with
  k*R ≈ 5. The occupancy is N/n.
- **Cluster shape.** Each cluster sits at a minimum of the lattice potential.
  The curvature there gives a width of about sqrt(D/κ) ≈ 0.008.
- **Merging and splitting.** Neighbouring clusters repel, so there is no drift
  or Brownian merging. Particles hop between neighbouring sites over a
  barrier ΔC ≈ 4.4 (ΔC/D ≈ 17), taking τ_hop ≈ (a²/D)e^(ΔC/D) ≈ 3·10⁶. The
  lattice is therefore static on any simulation time. Merging or splitting
  clusters, meaning a change in the number of lattice sites, does not occur
  unless the initial pattern is strongly mismatched to k*. Defects created
  during formation anneal only by hopping, so they persist.

### 3.4 HK, attractive (ε < 0)

- **Stability.** Û(k) = εR³ Ŝ(kR), with Ŝ(0) = π/4. Λ = 2.28 here (2.3× above
  onset). The onset is ε_crit ≈ −4D/(πR³) for small R.
- **Formation.** The fastest mode is (3,1), with σ ≈ 265 and t_form ≈ 0.021.
  As |ε| grows, the fastest kR moves toward the maximum of J₂, kR ≈ 3.05.
- **Number of clusters.** About 10 (9 hexagonal), with spacing λ ≈ 3.2R. In
  the strong-coupling limit the spacing tends to 2π/3.05·R ≈ 2.06R. Either
  way it lies **outside** the range R. The compact support is what freezes
  the pattern.
- **Cluster shape.** Gaussian clusters, with σ_c² = DR/(m|ε|), giving
  σ_c ≈ 0.018.
- **Merging.**
  - Drift: none, because neighbours are beyond R.
  - Brownian coalescence: t ≈ 100 at N = 81920 (∝ N).
  - Mass exchange: well depth m|U(0)|/D = m·γR/(2D) ≈ 15, so
    t_exch ≈ 3·10⁻⁴·e^15 ≈ 10³. That is Gerber et al.'s e^(γℓΔm)
    metastability, with D = 1.

  Both merging routes are out of reach of simulations. At fixed 5 particles
  per cell, the faster Brownian route needs γ well above onset for small λ,
  but the grid must also resolve σ_c. Both cannot hold on practical grids;
  see the discussion leading to `pp`.
- **Splitting:** never (e^(−3·10⁴)).

### 3.5 HK, repulsive (ε > 0)

- **Stability.** The transform's negative lobe, Ŝ_min = −0.046 at kR = 6.38,
  gives D_crit = μ₀ε R³|Ŝ_min| = 0.92 at ε = 2·10⁴, R = 0.1. The lobe is about
  twice as deep relative to Û(0) as GEM's (`STABILITY.md`).
- **Formation.** At D = 0.5, σ_max ≈ 1750 at k ≈ 10.4·2π, so t_form ≈ 0.003.
- **Number of clusters.** A hexagonal crystal with spacing a ≈ 1.11R, about 94
  clusters of about 870 particles.
- **Cluster shape.** The potential is concave inside the range: each parabolic
  piece has curvature −ε/R. So the lattice potential has its minimum off the
  site, at about 0.012. The clusters are small rings, confined by the kernel
  edges at the neighbouring sites. Resolving them needs dx well below 0.01.
- **Merging and splitting.** Hops cross a barrier ΔC/D ≈ 12, with
  τ_hop ≈ 4·10³. The crystal is static on simulation times, as for GEM; the
  cluster number changes only through collective lattice rearrangements.

### 3.6 Morse

- **Stability.** With the code's constraints (ε_att > ε_rep, R_att > R_rep),
  Û is increasing wherever it is negative (`STABILITY.md`). So either
  Û(0) < 0, giving long-wave demixing, or Û > 0 everywhere, which is stable.
  There is no cluster crystal.
  - At the code defaults (ε_att = 1, ε_rep = 0.5) Û(0) = −0.055, which is
    stable for any D ≳ 0.06.
  - With ε_att = 200, ε_rep = 100: Û(0) = −10.5 and Λ ≈ 21.
- **Formation.** Since Û(k) ≈ Û(0) + c k², growth peaks at a finite
  k_max² ≈ (|μ₀Û(0)| − D)/(2c), where c reflects the range R_att. Here the
  fastest mode is (1,1), with σ ≈ 283 and t_form ≈ 0.02.
- **Number of clusters.** About 2 here. More clusters need a larger |Û(0)|/D
  or a shorter R_att. Near threshold only the box mode grows: one band or
  blob.
- **Cluster shape.** With ε_att/R_att = ε_rep/R_rep (200/0.1 = 100/0.05) the
  cusp cancels and the core is harmonic, U″(0) = −ε_att/R_att² + ε_rep/R_rep².
  That gives σ_c ≈ 0.007. If ε_att/R_att > ε_rep/R_rep the core has a cusp,
  and the clusters are sharply peaked with radius D/(m|U′(0)|).
- **Merging.** The exponential attractive tail reaches between clusters, so
  drift merging takes ≈ 0.06. Morse coarsens quickly by drift, toward a
  single cluster. Brownian merging (≈ 5·10³) and mass exchange (e^100) are
  irrelevant.
- **Splitting:** never.

### 3.7 PP (piecewise parabolic)

- **Stability.** Û is the sum of two HK transforms, a tail of radius R and a
  core of radius aR. The core sets the fastest mode, provided the tail is weak
  (β ≲ 0.02 here; β ≥ 0.03 lets the box mode win).
- **Formation.** The fastest mode is (3,1), with σ ≈ 299 and t_form ≈ 0.021.
  **Observed:** clusters appear from t ≈ 0.014 and number 10 by t ≈ 0.018
  (`inputs_fv_pp_merge`, seed 11).
- **Number of clusters.** About 10 predicted, and 10 observed. The spacing,
  0.32 = 0.70R, lies **inside** the tail.
- **Cluster shape.** Set by the core stiffness: σ_c² = DR/(m|ε|α), giving
  σ_c ≈ 0.025.
- **Merging.** The tail pulls neighbours together: m|U′(λ)| ≈ 10, so drift
  merging takes ≈ 0.015, independent of N. **Observed:** 9 mergers, from 10
  clusters to 1, between t ≈ 0.027 and 0.036. Brownian merging (≈ 400) and
  mass exchange (≈ 2·10³) are irrelevant. This is why `pp` was added:
  merging by drift on simulation times at SPDE particle counts.
- **Splitting:** never (e^(−10⁵)).

## 4. How the time scales depend on the parameters

| Quantity | Scaling |
| --- | --- |
| t_form | ln N / (2σ_max). σ_max grows like (μ₀\|Û\| − D)k², so t_form shrinks with coupling and depends only logarithmically on N |
| Number of clusters | ∝ k_max² A: about A/R² times a kernel-dependent constant. Nearly independent of N, and of ε far above onset |
| Drift merging | ∝ 1/(m\|U′(λ)\|); independent of N; requires λ inside the range |
| Brownian merging | ∝ N m λ²/D; linear in N |
| Mass exchange | ∝ e^(m\|U(0)\|/D); independent of N; exponential in the well depth |
| Splitting (attractive) | ∝ e^(N m²\|U(0)\|/4D): never at these N |
| Hopping (repulsive crystals) | ∝ (a²/D) e^(ΔC/D) |

**How the estimates were made.**
- **Growth rates:** σ(k) is evaluated from the numerical Hankel transform of U
  on the box lattice modes. These are continuum Û values; the code's startup
  printout gives the discrete equivalents, which agree to a few percent at
  the resolutions used.
- **Mass, width and merging:** these use the leading estimates above,
  m ≈ 1/n and a harmonic or flat well. The mass-exchange and hopping times
  are Arrhenius estimates with an O(1) prefactor σ_c²/D or a²/D. Treat them
  as order-of-magnitude values; the exponents are the reliable part.
- **Lattice barriers** for the repulsive crystals are from direct sums over a
  hexagonal lattice of point clusters.

## References

- D. S. Dean, *J. Phys. A* 29, L613 (1996).
- C. N. Likos, A. Lang, M. Watzlawek and H. Löwen, *Phys. Rev. E* 63, 031206
  (2001): criterion for cluster crystals of bounded repulsive potentials.
- N. Gerber, R. S. Gvalani, M. Hairer, G. A. Pavliotis and A. Schlichting,
  *Formation of clusters and coarsening in weakly interacting diffusions*,
  arXiv:2510.17629: the hk and pp potentials, coalescence ∝ N and mass
  exchange ∝ e^(γℓΔm).
