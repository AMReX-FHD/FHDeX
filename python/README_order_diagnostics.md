# Order diagnostics in cluster_diagnostics.py and cluster_diagnostics_3d.py

This note explains the orientational-order diagnostics written by
`cluster_diagnostics.py` (2D) and `cluster_diagnostics_3d.py` (3D): what each
quantity is and how to read it. The LaTeX version is
`order_diagnostics.tex`.

## Why g(r) is not enough

The radial distribution function g(r) counts pairs of clusters at distance r,
averaged over every direction. It records how far apart clusters are, but not
in which directions their neighbours lie. A crystal and a dense amorphous
packing with the same spacing therefore have g(r) peaks in the same places;
the crystal's peaks are only somewhat sharper and its higher peaks a little
more structured. That is the case in `exec/dean_kow/interaction`:

- the 2D run `_new_repul_2` looks like a hexagonal lattice;
- the 3D run `_3d_noise` looks amorphous;
- the cluster g(r) of both has the same nearest-neighbour peak.

Telling these states apart needs angular information. The tools measure it in
two independent ways:

1. **Bond-orientational order in real space:** the directions from each
   cluster centre to its neighbours, compared within each neighbourhood
   (local order) and across the box (global order).
2. **Bragg spots in k space:** whether the main peak of the structure factor
   S(**k**) sits on a few wave vectors (spots) or is spread around the ring.

Both are written for every frame to `order.txt`. The per-cluster bond order is
added as extra columns of `clusters.txt`.

## 1. Neighbours and bonds

Order is measured on the cluster centres, the circular-mean centres of mass in
`clusters.txt`. For each centre **x**ᵢ:

- **Bond vectors:** **r**ᵢⱼ = **x**ⱼ − **x**ᵢ, using the minimum image in each
  periodic direction.
- **Neighbours, by default:** the `--nn` nearest centres (6 in 2D, 12 in 3D).
- **Neighbours with `--nn-cut R`:** every centre with |**r**ᵢⱼ| < R. A good
  choice of R is the first minimum of the cluster g(r) (`rdf_clusters.txt`),
  after its first peak. The number of neighbours nᵢ then varies from cluster
  to cluster. A cluster with no neighbours gets NaN and is left out of the
  averages.

The neighbour count is a choice, and it should match the coordination of the
packing you are testing for:

| Packing | Natural neighbour count |
| --- | --- |
| 2D hexagonal | 6 |
| 3D fcc or hcp | 12 |
| 3D bcc | 8 (or 14 including the second shell) |
| 3D simple cubic | 6 |

Bonds are not symmetrized: j can be among i's nearest neighbours without i
being among j's. With a fixed `--nn` this is the usual convention and has
little effect.

## 2. Local and global bond order in 2D: ψ₆

Let θᵢⱼ be the angle of bond **r**ᵢⱼ to the x axis. Each cluster gets the
complex number

  ψ₆,ᵢ = (1/nᵢ) Σⱼ exp(6 i θᵢⱼ).

The factor 6 makes the six bonds of a perfect hexagon all contribute the same
phase: rotating a bond by 60° changes 6θ by 360°. So:

- **The modulus |ψ₆,ᵢ|** says how hexagonal the neighbourhood is. It is 1 for
  a perfect hexagon and small when the bond angles are irregular.
- **The phase arg ψ₆,ᵢ / 6** is the orientation of that hexagon.

Two averages over the N clusters of a frame are reported:

| Quantity | Definition | Measures |
| --- | --- | --- |
| `psi6_local` | (1/N) Σᵢ \|ψ₆,ᵢ\| | how hexagonal each neighbourhood is, whatever its orientation |
| `psi6_global` | \|(1/N) Σᵢ ψ₆,ᵢ\| | whether the hexagons share one orientation across the box |

Taking the modulus before averaging ignores orientation. Averaging the complex
numbers first lets differently oriented hexagons cancel. Since the mean of a
set of numbers is never larger in modulus than the mean of their moduli,
`psi6_global` ≤ `psi6_local` always holds.

`clusters.txt` gets a `psi6` column: |ψ₆,ᵢ| of each cluster.

### Reference values (6 nearest neighbours)

| Configuration | psi6_local | psi6_global |
| --- | --- | --- |
| perfect hexagonal lattice | 1 | 1 |
| hexagonal lattice, random displacements of 5% of the spacing | 0.92 | 0.92 |
| random points (Poisson) | 0.36 | about 0.4/√N (0.054 for N = 80, 0.018 for N = 500) |

The random value 0.36 comes from averaging six random phases: |mean| ≈ √(π/24)
≈ 0.36. The global floor about 0.4/√N is the mean of N random unit-scale
complex numbers. It is not zero in a finite system, so a global value must be
compared with 0.4/√N, not with 0.

## 3. Local and global bond order in 3D: Steinhardt q₄, q₆

In 3D a bond has a direction r̂ᵢⱼ (polar angle θ, azimuth φ) rather than a
single angle. The angular pattern of each neighbourhood is expanded in
spherical harmonics:

  q_lm,ᵢ = (1/nᵢ) Σⱼ Y_lm(r̂ᵢⱼ),  m = −l … l,

  q_l,ᵢ = √( 4π/(2l+1) · Σ_m |q_lm,ᵢ|² ).

These are the Steinhardt–Nelson–Ronchetti bond-order parameters:

- The set q_lm,ᵢ (2l+1 complex numbers) is the 3D analogue of ψ₆,ᵢ: it encodes
  both the shape and the orientation of the neighbourhood.
- The sum over m makes q_l,ᵢ rotationally invariant: it describes the shape of
  the neighbourhood regardless of its orientation, as |ψ₆| does in 2D.
- The normalization gives q_l = 1 when all bonds point in the same direction.
- l = 4 and l = 6 are the standard pair: cubic and icosahedral symmetries
  have strong l = 4 and l = 6 components, and odd l vanish for
  centrosymmetric neighbourhoods.
- Y_lm is computed from Legendre polynomials in `numpy.polynomial`, without
  scipy.

The local and global averages follow the 2D logic:

| Quantity | Definition | Measures |
| --- | --- | --- |
| `q4`, `q6` | (1/N) Σᵢ q_l,ᵢ | how regular each neighbourhood is (rotationally invariant) |
| `Q4_global`, `Q6_global` | √(4π/(2l+1) Σ_m \|(1/N) Σᵢ q_lm,ᵢ\|²) | whether the neighbourhoods share one orientation |

In the global value the q_lm are averaged over clusters before the norm, so
differently oriented neighbourhoods cancel. It equals the local value for a
single perfect crystal of any orientation, and is near zero when orientations
are random.

`clusters.txt` gets `q4` and `q6` columns: q₄,ᵢ and q₆,ᵢ of each cluster.

### Crystal packings

The reference lattices are the standard ways of packing equal spheres (here,
cluster centres) in 3D:

| Name | Structure | Nearest neighbours | Packing fraction |
| --- | --- | --- | --- |
| **fcc**, face-centred cubic | a cube with a point at each corner and at the centre of each face (4 points per cubic cell) | 12 | 0.74, the densest possible |
| **hcp**, hexagonal close-packed | hexagonal layers stacked ABAB… (see below) | 12 | 0.74 |
| **bcc**, body-centred cubic | a cube with a point at each corner and one at its centre (2 points per cell) | 8, plus 6 others only 15% farther | 0.68 |
| **sc**, simple cubic | points only at the cube corners (1 point per cell) | 6 | 0.52 |

**fcc and hcp.** Both are stacks of close-packed layers, each identical to
the 2D hexagonal lattice. Each layer sits in the hollows of the one below,
and there are two sets of hollows, so three layer positions A, B and C are
possible. ABAB… stacking gives hcp; ABCABC… gives fcc. Random stackings,
mixtures of the two, are common in colloidal crystals. In both lattices every
point has 12 neighbours at the same distance, so g(r) barely tells them apart.
Their angular arrangements differ, and so do q₄ and q₆.

**bcc** is less dense but more open: each point sits at the centre of a cube
of 8 neighbours, with 6 more at 2/√3 ≈ 1.15 times that distance. Its
reference values are for the 8 nearest neighbours. Using 12 mixes the first
and second shells, so use `--nn 8` (or `--nn 14`, or `--nn-cut` between the
shells) when testing for bcc.

**sc** is rarely stable for equal spheres. It is listed because its q₄ and q₆
are very distinct.

Which lattice forms depends on the interaction:
- **Hard or short-range repulsion:** typically fcc or hcp, the close packings.
- **Soft, long-range repulsion:** often bcc, especially near melting.
- **Soft, penetrable potentials** (for example the gema and HK potentials of
  `src_deankow/interaction`): cluster crystals can be fcc or bcc depending on
  density and noise.

A **polycrystal** is many crystalline grains with different orientations.
An **amorphous** (glassy or liquid-like) packing has no lattice at all, only
short-range order.

### Reference values

Perfect lattices, checked with the code; the values are independent of the
lattice's orientation, and global = local:

| Packing (neighbours) | q₄ | q₆ |
| --- | --- | --- |
| fcc (12) | 0.191 | 0.575 |
| hcp (12) | 0.097 | 0.485 |
| bcc (8) | 0.509 | 0.629 |
| simple cubic (6) | 0.764 | 0.354 |
| icosahedral cluster (12), for comparison | 0 | 0.663 |
| random points (12) | 0.28 | 0.28; global Q₆ about 0.4/√N (0.017 for N = 500) |

Thermal displacements lower all crystal values and broaden their
distributions. A dense liquid or glass typically has q₆ of about 0.35–0.45,
above the random-point value because neighbours keep a minimum distance but
below the crystals.

## 4. Bragg spots or a ring: S(k) on the peak shell

For each frame, the tool computes the full structure factor
S(**k**) = |n̂(**k**)|²/N_tot on the FFT grid, before any shell average, as in
`sf.txt`. It then takes the shell (ring in 2D, spherical shell in 3D) with the
largest mean S, which is normally the main peak at 2π/(cluster spacing). On
that shell, with n wave vectors and values S₁ … Sₙ:

  PR = (Σ Sₐ)² / (n · Σ Sₐ²),

  top_share = (sum of the `--top-modes` largest Sₐ) / Σ Sₐ.

| Column | Meaning |
| --- | --- |
| `kpeak` | mean \|**k**\| of the peak shell, about 2π/(nearest-neighbour spacing) times a lattice factor |
| `nmodes` | n, the number of FFT wave vectors on the shell |
| `PR` | participation ratio: the fraction of the shell's modes that carry the signal |
| `top_share` | fraction of the shell's S in its largest modes (default 6 in 2D, 12 in 3D) |

PR is the inverse-participation idea: if S were equal on k of the n modes and
zero elsewhere, PR = k/n exactly. It assumes no particular lattice.

- **Crystal:** a few Bragg spots carry the shell, so PR ≈ (number of spots)/n
  and top_share ≈ 1 when `--top-modes` is at least the number of spots. Spots
  come in ±**k** pairs, so a 2D hexagonal lattice has 6 on its first ring, an
  fcc lattice 8 on its first shell, and bcc 12.
- **Isotropic ring (amorphous, liquid):** with random phases each Sₐ is
  approximately exponentially distributed about the ring average, which gives
  PR ≈ 0.5 (because ⟨S²⟩ = 2⟨S⟩²) and top_share about `top_modes`/n, somewhat
  more due to fluctuations.
- **Uniform field** (no structure, e.g. t = 0): PR and top_share are NaN.

Reference values, measured with the code:

| Field | n | PR | top_share |
| --- | --- | --- | --- |
| 2D hexagonal lattice of Gaussian blobs | 60 | 0.100 = 6/60 | 1.00 |
| 2D random blobs | 12 | 0.61 | — |
| 3D fcc lattice of blobs, 48³ | 602 | 0.0133 = 8/602 | 1.00 |
| 3D random blobs | 98 | 0.55 | 0.36 |

This measure is independent of the cluster detection (`--threshold`,
`--mmin`) and of the neighbour rule, so it is a cross-check on the bond order.

### Seeing the spots

`--sf-full` writes the full S(**k**), averaged over the frames in
`--t-range`, as the plotfile `sf_full_avg`, with k = 0 at the centre and
components `S` and `log10S`. `--sf-full-per-frame` writes one per frame.
Open them in VisIt, ParaView or amrvis.

- **2D crystal:** distinct bright spots. A hexagonal lattice has six on the
  first ring, 60° apart, and more at √3 and 2 times that radius.
- **Polycrystal:** several sets of spots at different angles on the same ring.
- **Amorphous state:** a ring of uniform, speckled intensity.
- **3D:** view slices through k = 0 (planes kx = 0, ky = 0, kz = 0), or an
  isosurface of `log10S` just below the peak value. Crystal spots are about
  nmodes/nspots times the shell mean, which is 100 or more. Random speckle
  on a ring of n modes reaches only about ln n times the mean, about 7 for
  1000 modes.

## 5. The columns of order.txt

There is one row per plotfile, after a header line that gives the neighbour
rule and `--top-modes`.

### 2D (`cluster_diagnostics.py`)

| # | Column | Meaning |
| --- | --- | --- |
| 1 | `t` | plotfile time |
| 2 | `step` | plotfile step |
| 3 | `nclusters` | number of clusters in the frame (as in `summary.txt`) |
| 4 | `psi6_local` | mean over clusters of \|ψ₆,ᵢ\|: how hexagonal each neighbourhood is, whatever its orientation. 1 is a perfect hexagon; random points give about 0.36 |
| 5 | `psi6_global` | \|mean over clusters of ψ₆,ᵢ\|: whether the hexagons share one orientation across the box. About equal to `psi6_local` for a single crystal; about 0.4/√N with no common orientation |
| 6 | `kpeak` | mean \|k\| of the S(k) ring with the largest mean S in this frame, normally the main peak at about 2π/spacing |
| 7 | `nmodes` | number of FFT wave vectors on that ring |
| 8 | `PR` | participation ratio (ΣS)²/(nmodes·ΣS²) on that ring, roughly the fraction of modes carrying the signal: about (number of spots)/nmodes for Bragg spots, about 0.5 for a uniform ring |
| 9 | `top_share` | fraction of the ring's total S in its `--top-modes` largest modes (default 6): about 1 for a crystal, small for a ring |

Example, the last frame of `_new_repul_2`:

```
0.08773803711 400000 53 0.987506 0.986861 50.306732 48 0.125774 0.996647
```

This is t = 0.0877, step 400000, with 53 clusters:
- Local and global ψ₆ are both 0.99: one hexagonal crystal with a single
  orientation.
- The S(k) peak ring is at |k| ≈ 50.3 and has 48 modes.
- PR = 0.126 ≈ 6/48, and the 6 largest modes hold 99.7% of the ring: six
  Bragg spots.

The t = 0 row reads `0 0 0 nan nan 7.58 8 nan nan`. The field is still
uniform, so the order values are NaN, and `kpeak` and `nmodes` only describe
the first ring.

### 3D (`cluster_diagnostics_3d.py`)

| # | Column | Meaning |
| --- | --- | --- |
| 1–3 | `t`, `step`, `nclusters` | as in 2D |
| 4 | `q4` | mean over clusters of q₄,ᵢ (local; fcc 0.19, bcc 0.51 with 8 nn, random 0.28) |
| 5 | `q6` | mean over clusters of q₆,ᵢ (local; fcc 0.575, hcp 0.485, random 0.28) |
| 6 | `Q4_global` | Q₄ from q_4m averaged over all clusters (one orientation across the box) |
| 7 | `Q6_global` | Q₆ from q_6m averaged over all clusters: about `q6` for a single crystal, about 0.4/√N with no common orientation |
| 8–11 | `kpeak`, `nmodes`, `PR`, `top_share` | as in 2D; `--top-modes` defaults to 12 |

Example, the last frame of `_3d_noise`:

```
0.03797721181 60000 476 0.126799 0.424355 0.016313 0.073719 56.909256 1142 0.318627 0.096069
```

- 476 clusters with local q₆ = 0.42 and q₄ = 0.13: between random and the
  crystals, and matching no crystal.
- Global Q₆ = 0.074, about six times below local: no common orientation.
- The peak ring at |k| ≈ 56.9 has 1142 modes. PR = 0.32, and the 12 largest
  modes hold only 9.6% of it: a ring, not spots. This is an amorphous
  packing.

## 6. How to read them together

### Decision table

| Local order | Global order | S(k) peak shell | Interpretation |
| --- | --- | --- | --- |
| near the crystal value | ≈ local | PR ≈ spots/n, top_share ≈ 1 | **single crystal** spanning the box |
| near the crystal value | well below local, well above 0.4/√N | PR intermediate, several groups of spots | **polycrystal** of a few grains; global ≈ local/√(number of grains) |
| near the crystal value | ≈ 0.4/√N | ring, PR ≈ 0.3–0.5 | **fine polycrystal** of many small grains, hard to distinguish from amorphous; look at per-cluster histograms |
| between random and crystal | ≈ 0.4/√N | ring | **amorphous / dense liquid**: neighbours are kept at a distance but not on a lattice |
| ≈ random value | ≈ 0.4/√N | ring, or no clear peak | **disordered, gas-like** arrangement of clusters |

In 2D, dense liquids and hexatic phases can also have fairly high local ψ₆,
so local ψ₆ alone does not prove a crystal. The global value and the spots
decide.

### Time series

Plot the `order.txt` columns against t:

- **Crystallization:** local order rises first, then global order, while PR
  drops towards spots/n and top_share rises towards 1.
- **Coarsening of a polycrystal:** local order is flat, global order rises
  slowly as grains merge.
- **Amorphous or glassy state:** all quantities stay roughly constant. A slow
  drift may mean very slow ordering; in that case a longer run is needed before
  concluding.

### Per-cluster distributions

The `clusters.txt` columns (`psi6`, or `q4 q6`) give the distribution over
clusters in a frame:

- **One narrow peak at a crystal value:** crystal.
- **One broad peak below the crystal values:** amorphous.
- **Two peaks:** crystalline grains coexisting with disordered regions; the
  weight of the upper peak estimates the crystalline fraction.
- **3D, q₄ against q₆ scatter plot:** points cluster near (0.19, 0.575) for fcc,
  (0.10, 0.485) for hcp and (0.51, 0.63) for bcc (8 neighbours). This helps
  identify which crystal forms.

## 7. Worked example: `_new_repul_2` (2D) and `_3d_noise` (3D)

| | 2D `_new_repul_2` | 3D `_3d_noise` |
| --- | --- | --- |
| Frames | t = 0.022 to 0.088 | t = 0.032 to 0.038 |
| Clusters N | 41–53 | 465–480 |
| Local order | ψ₆ 0.63 → 0.89 → 0.91 → 0.988 | q₄ 0.13, q₆ 0.41–0.42 |
| Global order | ψ₆ 0.55 → 0.86 → 0.89 → 0.987 | Q₄ 0.012–0.016, Q₆ 0.058–0.074 |
| Global noise floor 0.4/√N | 0.055 | 0.018 |
| S(k) peak shell | 48 modes; PR 0.23 → 0.126; top_share (6) 0.66 → 0.997 | 1142 modes; PR 0.32; top_share (12) 0.08–0.10 |

**2D: a single hexagonal crystal, formed during the run.**
- Local and global ψ₆ agree and approach 1.
- Six Bragg spots carry essentially the whole ring: PR 0.126 against
  6/48 = 0.125.
- Early on the ordering is incomplete: ψ₆ 0.63 and PR 0.23. The lattice
  anneals by t ≈ 0.09.

**3D: amorphous.**
- Local q₆ of 0.42 is above the random value 0.28 but below hcp (0.485) and
  fcc (0.575), and q₄ = 0.13 matches no crystal. This is the signature of a
  dense, liquid-like packing.
- Global Q₆ of 0.07 is about 6 times smaller than the local value: the
  neighbourhoods have no common orientation. It is only 4 times the
  random-orientation floor, so at most there are very many tiny grains.
- The S(k) peak is a ring: the 12 largest of 1142 modes hold under 10% of it,
  and PR is 0.32, far above the ≈ 0.01 that Bragg spots would give.
- In `sf_full_avg` over t = 0.037–0.038, the brightest modes on the peak
  shell are about 14 times the shell mean (|k| ≈ 53–55). That is above
  pure speckle (about 7) but far below Bragg spots (about 100). The ring
  is weakly lumpy: some local ordering, no lattice.
- All values drift only slowly over the window, so slow ordering at later
  times cannot be excluded.

## 8. Caveats

- **Few clusters:** with N clusters the global floor is about 0.4/√N. With 10–20
  clusters, global values of 0.1–0.15 can arise by chance.
- **Neighbour rule:** use `--nn` matching the expected coordination. For a
  suspected bcc packing use `--nn 8`; for loose or polydisperse packings use
  `--nn-cut` at the first minimum of the cluster g(r). Different rules shift
  the reference values; compare against values computed with the same rule.
- **Cluster detection:** stray small clusters inside a lattice, or a merger in
  progress, lower the local order. Check `--mmin` and `--threshold` if local
  order is unexpectedly low. The S(k) measure does not depend on this.
- **Periodic box:** a crystal fits the box without strain only for certain
  sizes and orientations. A crystal that does not fit has defects or a twist,
  which lowers global order and smears the spots, so PR is moderate. The
  clean spots of `_new_repul_2` show that its lattice fits the box.
- **Peak shell:** the shell of largest mean S is used. If a low-k shell is
  larger, for example during coarsening, `kpeak` jumps, and PR then describes
  that shell instead. Check `kpeak` against 2π/(spacing).
- **Not implemented:** classifying individual clusters as crystalline with the
  ten Wolde criterion (normalized q₆,ᵢ·q₆,ⱼ > 0.5 for most neighbours), which
  would give crystalline fractions and grain sizes directly.

## References

- P. J. Steinhardt, D. R. Nelson and M. Ronchetti, Phys. Rev. B 28, 784 (1983):
  bond-orientational order, q₄ and q₆.
- D. R. Nelson and B. I. Halperin, Phys. Rev. B 19, 2457 (1979): ψ₆ and 2D
  melting.
- P. R. ten Wolde, M. J. Ruiz-Montero and D. Frenkel, J. Chem. Phys. 104, 9932
  (1996): per-particle crystallinity from q₆ correlations.
- W. Lechner and C. Dellago, J. Chem. Phys. 129, 114707 (2008): averaged
  q₄, q₆ for better separation of crystal types.
