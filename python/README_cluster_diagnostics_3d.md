# cluster_diagnostics_3d.py

The 3D analog of `cluster_diagnostics.py`. It computes cluster statistics and
radial distribution functions for a sequence of 3D AMReX plotfiles in which x,
y and z are all spatial. It was written for the 3D runs of
`src_deankow/interaction` (for example `exec/dean_kow/interaction/inputs_fv_3d`).

For each plotfile the script:
1. finds the clusters in the density field and measures each one's mass,
   centre, size and peak density;
2. matches them to the clusters in the previous plotfile, giving each cluster a
   stable id and logging mergers, splits, new clusters and clusters that
   dissolve;
3. accumulates two radial distribution functions: g(r) of the density field and
   g(r) of the cluster centres, and the structure factor S(k);
4. measures the orientational order of the cluster centres (Steinhardt q4, q6)
   and whether S(k) has Bragg spots or a ring, to tell a cluster crystal from
   an amorphous packing.

The results are plain whitespace-separated text tables with `#` header lines.
The script needs only Python 3 and numpy, and reuses the plotfile reader in
`stack_plotfiles_time.py`, which must be in the same directory.

The definitions, options and outputs match the 2D tool
(`README_cluster_diagnostics.md`). The differences are listed below.

## Differences from the 2D tool

| | 2D (`cluster_diagnostics.py`) | 3D (`cluster_diagnostics_3d.py`) |
| --- | --- | --- |
| Connectivity | 4 face neighbours | 6 face neighbours |
| Cluster centre | (xc, yc) | (xc, yc, zc); `clusters.txt` gains a `zc` column |
| Bond order | ψ₆ (6 nn); `clusters.txt` column `psi6` | Steinhardt q4, q6 (12 nn); `clusters.txt` columns `q4 q6` |
| Cell measure | ΔV = Δx Δy | ΔV = Δx Δy Δz |
| Cluster g(r) shell | π((r+dr)² − r²)/A | (4π/3)((r+dr)³ − r³)/V |
| Field g(r) | 2D FFT autocorrelation | 3D FFT autocorrelation (`rfftn`) |
| Labelling | per-cell flood fill | vectorized label propagation (fast at 128³) |
| Per-cluster sums | per-cluster loops | `np.bincount` / `np.maximum.at` over labels |
| `--mmin` default | 0.01 | 10⁻⁴; and a warning when most mass above the threshold is in regions lighter than `--mmin` |

At 128³ (2·10⁶ cells) a frame takes about 2 s and about 300 MB of memory.

## Usage

```
python3 cluster_diagnostics_3d.py plt_fv3d_* --num-part 2e7 -o diag3d
python3 cluster_diagnostics_3d.py plt_fv3d_* --num-part 2e7 -o diag3d --every 10
```

Plotfiles are sorted by the step number in their headers, not by file name.
All of them must have the same domain and cell size. A 2D plotfile is rejected
with a pointer to `cluster_diagnostics.py`.

## Assumptions and limitations

- **3D, periodic in all directions, level 0 only.** Finer levels are ignored,
  with a warning.
- **Density field only.** Both g(r) are computed from the density variable
  (`phi0` by default).
- **Touching clusters count as one.** Two clusters become one as soon as the
  density between them exceeds the threshold; this is the definition of a
  merger. In 3D the region above the threshold extends further from a
  cluster's centre than in 2D (for a Gaussian, out to where the density
  falls to c·φ̄), so merges are registered at larger separations.
- **Total mass includes negative densities,** as in 2D. The mass fraction in
  clusters can slightly exceed 1.
- **Choosing `--mmin`.** It is a fraction of the total mass. The default,
  10⁻⁴, comes from the 128³ interaction run (`_3d_noise/plt_fv3d_060000`).
  There, the regions above the 5× threshold split into:
  - noise: 1–5 cells, mass ≤ 3·10⁻⁵;
  - no regions at all between 3·10⁻⁵ and 10⁻⁴;
  - small, young or evaporating clusters: 10⁻⁴–10⁻³;
  - about 470 clusters of 10⁻³–3·10⁻³, holding 91% of the mass.

  10⁻⁴ keeps every cluster and drops only noise. The same gap appears at
  thresholds of 2, 10 and 20× the mean. For other runs, a rule of thumb is
  `--mmin` ≈ 0.05/(number of clusters). If most of the mass above the
  threshold is in regions lighter than `--mmin`, the script warns once.
  The 2D tool keeps 0.01, which suits runs of about 10 clusters.

## Definitions

Let the domain be the periodic box [x_lo, x_lo + L_x) × [y_lo, y_lo + L_y) ×
[z_lo, z_lo + L_z), with cells of volume ΔV = Δx Δy Δz and centres xᵢ, and let
φᵢ be the density in cell i. The minimum-image displacement in direction d is
mi_d(a) = ((a + L_d/2) mod L_d) − L_d/2.

### Cluster

With φ̄ = Σφᵢ ΔV / (L_x L_y L_z) and threshold factor c (`--threshold`), a
cluster is a set of cells connected through shared faces, including across the
periodic boundaries, in which φᵢ > c φ̄. It is kept only if its mass
M = Σ φᵢ ΔV is at least `--mmin` times the total mass.

**Labelling.** Every cell above the threshold starts with a unique label. Each
cell then repeatedly takes the largest label among itself and its six
neighbours (with periodic wrap), restricted to cells above the threshold, until
nothing changes. This takes about as many sweeps as the largest cluster's
diameter in cells. The surviving labels are then renumbered 1..n.

For each cluster C, with weights wᵢ = φᵢ ΔV:

| Quantity | Definition |
| --- | --- |
| mass M | Σ_{i∈C} wᵢ |
| centre (xc, yc, zc) | circular mean in each direction, x_c = (L_x/2π)·atan2(Σ wᵢ sin(2πxᵢ/L_x), Σ wᵢ cos(2πxᵢ/L_x)) mod L_x; correct for clusters straddling a boundary |
| Rg | radius of gyration, √(Σ wᵢ \|dᵢ\|² / M), with dᵢ the minimum-image displacement from the centre |
| peak | largest φ in the cluster |
| ncells | number of cells |
| wide | 1 if some cell is ≥ 0.45 L_d from the centre in any direction d (minimum image ambiguous; mainly for a single box-filling cluster) |

### Tracking between frames

These are the same rules as in 2D, with 3D minimum-image distances. Let
O_ba = Σ_{i ∈ b∩a} φᵢ^old ΔV.
1. **Matching.** Old cluster b goes to argmax_a O_ba if that exceeds 0.1 M_b.
   Otherwise it goes to the new cluster with the nearest centre, provided that
   is within `--search-radius`.
2. **Merger.** A new cluster with two or more old clusters assigned, provided
   it holds at least 0.7 of their summed mass. Old clusters that would break
   this are dissolved instead, starting with the farthest.
3. **Split.** An old cluster with O_ba > 0.1 M_b for two or more new clusters.
4. **Formed and dissolved.** Unassigned new and old clusters that are not part
   of a split.
5. **Ids.** A new cluster takes the id of its most massive parent. New clusters
   and split-off pieces get fresh ids.

Tracking needs plotfiles close enough in time that clusters move less than
about their own size between frames.

### Density-field g(r)

With cell counts nᵢ = φᵢ N ΔV, N = `--num-part`, M the number of cells and
N_tot = Σ nᵢ:

  G(Δ) = Σᵢ nᵢ n_{i+Δ} − δ_{Δ,0} Σᵢ nᵢ,  g(Δ) = M·G(Δ) / (N_tot(N_tot − 1))

G is a periodic 3D autocorrelation, computed by FFT; a Poisson field gives
g = 1 at every Δ. Each displacement is wrapped componentwise into [−M_d/2, M_d/2),
so r = |(Δx Δx, Δy Δy, Δz Δz)| is the minimum-image distance. g(r) in a bin
is the average of g(Δ) over the lattice displacements in it, which normalizes by
the exact count per shell. The result is averaged over the frames in
`--t-range`.

### Cluster-centre g(r)

For every frame in `--t-range` with n_c ≥ 2 clusters, minimum-image distances
between all pairs of centres are histogrammed, counting each pair once. With
V = L_x L_y L_z,

  g(r_k) = Σ_frames pairs in [r_k, r_k + Δr) / Σ_frames ½ n_c(n_c − 1) · (4π/3)((r_k + Δr)³ − r_k³)/V

The shell volume is exact for r ≤ L/2, because every such sphere lies inside
the minimum-image cube. Pair and expected counts are summed over frames before
dividing.

### Structure factor S(k) (as in 2D, with 3D wave vectors)

From the same cell counts as the field g(r) (nᵢ = φᵢ N ΔV, N_tot = Σ nᵢ):

  S(k) = |n̂(k)|² / N_tot,  n̂(k) = Σᵢ nᵢ e^(−ik·xᵢ),  k ≠ 0

A Poisson field gives S = 1, and S(k) = 1 + (N/V)∫(g − 1)e^(ik·r) dr. The wave
vectors come from `fftfreq`: index i > N_d/2 is the negative wavenumber
i − N_d, so k_d = 2π·(index)/L_d. S is averaged over shells: shell j holds the
modes with |k| in [(j − ½)Δk, (j + ½)Δk), with Δk = `--dk` (default
2π/max L_d, the smallest box wavenumber). Each shell is divided by its exact
number of modes. Only complete shells, up to the smallest Nyquist wavenumber
π/Δx_d, are written. The result is averaged over the frames in `--t-range`.
Statistics are limited by the number of plotfiles. For many samples, use the
in-code structure factor of `src_deankow/interaction` (`struct_fact_int`),
which writes the same shell average.

### Orientational order: crystal or amorphous

A detailed description and interpretation guide is in
`README_order_diagnostics.md` (LaTeX: `order_diagnostics.tex`).

g(r) averages over directions, so a cluster crystal and an amorphous packing
with the same spacing can have the same g(r) peaks. Two measures tell them
apart. Both are written to `order.txt` for every frame.

**Steinhardt bond order.** Each cluster centre i is joined to its `--nn`
nearest neighbours (default 12, minimum image), or to all neighbours closer
than `--nn-cut`. For l = 4 and 6:

  q_lm,ᵢ = (1/nᵢ) Σⱼ Y_lm(r̂ᵢⱼ),  q_l,ᵢ = √(4π/(2l+1) Σ_m |q_lm,ᵢ|²)

The local order is the mean of q_l,ᵢ. The global Q_l uses q_lm averaged over
all clusters before taking the norm. It equals the local value in a single
crystal of any orientation, and is about 1/√N for an amorphous packing or a
polycrystal of many grains. Y_lm comes from numpy's Legendre polynomials (no
scipy).

The packings are fcc (face-centred cubic), hcp (hexagonal close-packed), bcc
(body-centred cubic) and sc (simple cubic). They are defined and compared in
`README_order_diagnostics.md`, section 3.

| Packing (neighbours) | q4 | q6 |
| --- | --- | --- |
| fcc (12) | 0.191 | 0.575 |
| hcp (12) | 0.097 | 0.485 |
| bcc (8) | 0.509 | 0.629 |
| simple cubic (6) | 0.764 | 0.354 |
| random points (12) | 0.28 | 0.28; global Q6 about 0.4/√N |

**Bragg spots or a ring.** This is the same as in 2D:
PR = (Σ S)²/(n Σ S²) on the S(k) shell with the largest mean S. It is about
(number of spots)/n for a crystal, for example 8/602 = 0.0133 for an fcc
lattice of blobs, and about 0.5 for a ring. `top_share` is the fraction of the
shell in its `--top-modes` largest modes (default 12).

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `plotfiles` | (required) | 3D plotfile directories |
| `--num-part N` | (required) | `num_part` of the run; converts φ to particle counts for the field g(r) |
| `--var NAME` | `phi0` | plotfile variable used as the density |
| `--threshold C` | 5 | cluster cells have φ > C φ̄ |
| `--mmin F` | 10⁻⁴ | smallest cluster kept, as a fraction of the total mass; about 0.05/(number of clusters) |
| `--every K` | 1 | use every K-th plotfile, after sorting by step |
| `--rmax R` | L/2 | range of both g(r); values above L/2 (shortest side) are reduced to L/2 with a warning |
| `--dr D` | Δx | bin width of the field g(r) |
| `--dr-clusters D` | L/50 | bin width of the cluster-centre g(r) |
| `--search-radius R` | L/4 | largest distance a centre may move between frames and still be matched when its old cells overlap no new cluster |
| `--t-range T0 T1` | all frames | time window for both g(r) averages |
| `--rdf-per-frame` | off | also write the field g(r) of each frame in the window |
| `--dk D` | 2π/max L | shell width of the structure factor S(k) |
| `--sf-per-frame` | off | also write the shell-averaged S(k) of each frame in the window |
| `--nn K` | 12 | bond order: number of nearest neighbours of each cluster |
| `--nn-cut R` | off | bond order: use all neighbours closer than R instead, e.g. the first minimum of the cluster g(r) |
| `--top-modes K` | 12 | `order.txt` top_share: number of largest modes on the S(k) peak shell |
| `--sf-full` | off | Write the full S(**k**) on the k grid, averaged over the frames in `--t-range`, as the plotfile `sf_full_avg` |
| `--sf-full-per-frame` | off | Also write the full S(**k**) of each frame in the window as `sf_full_<step>` |
| `-o DIR`, `--output DIR` | `cluster_diag_3d` | output directory, created if needed; files are overwritten |

## Output files

These are the same files as the 2D tool. Times are plotfile times, lengths are
in the plotfile's units, and masses are Σφ ΔV, a fraction of the total.

- **`summary.txt`:** `t step nclusters mass_fraction largest_mass merges splits forms dissolves`
- **`clusters.txt`:** `t step id mass xc yc zc Rg peak ncells wide q4 q6`, one
  line per cluster per frame, sorted by id.
- **`events.txt`:** `t step type old_ids old_masses new_ids new_masses`, one line
  per merger or split.
- **`rdf_field.txt`:** a header giving the number of frames, the time window,
  `num_part`, Δr and r_max, then `r g ndisplacements` and, with
  `--rdf-per-frame`, one `g(t=...)` column per frame.
- **`rdf_clusters.txt`:** a header, then `r_lo r_hi g pairs expected`.
- **`sf.txt`:** a header, then `k S nmodes` (and `S(t=...)` per frame with
  `--sf-per-frame`): the shell-averaged structure factor.
- **`order.txt`:** one line per frame,
  `t step nclusters q4 q6 Q4_global Q6_global kpeak nmodes PR top_share`. The
  columns q4 and q6 are local means, Q4_global and Q6_global are global, and
  the S(k) columns are as in 2D.
- **`sf_full_avg`, `sf_full_<step>`** (with `--sf-full` or
  `--sf-full-per-frame`): the full S(**k**) as a plotfile on the 3D k grid,
  with k = 0 at the centre and components `S` and `log10S`. It is 33 MB at
  128³; see the 2D README for the layout.

See `README_cluster_diagnostics.md` for the meaning of each column.

## Checks

| Test | Result |
| --- | --- |
| Poisson field, 64³, N = 4·10⁶ | field g = 1.00000, rms deviation 1·10⁻⁵ |
| 4×4×4 lattice of Gaussian blobs, including sites on the boundaries | 64 clusters at the sites (xc = 1 ≡ 0). Cluster g(r) nonzero only at L/4, √2L/4, √3L/4 and L/2, with pair counts 192, 384, 256 and 96 |
| Lattice shifted by half a cell and by half the box | identical g(r) |
| `--rmax 0.7` | clipped to 0.5 with a warning |
| Two blobs merging across the z boundary, and the reverse | one merge and one split event, with masses 0.594 + 0.394 → 0.988 |
| 2D pattern extruded in z, against the 2D tool on the 2D pattern | same 4 clusters, masses and centres |
| 3D interaction run (`_3d_noise/plt_fv3d_*`, 128³), default `--mmin` | about 476 clusters (92% of the mass), mostly of mass 0.001–0.003; cluster g(r) peak at 0.14–0.16; about 2 s per frame |
| Bond order of sc, bcc, fcc and hcp point lattices; random points | the q4 and q6 values in the table above, global = local, invariant under rotation; random: q6 0.283, Q6 0.017 (N = 500) |
| fcc lattice of blobs, 48³; random blobs | PR = 0.0133 = 8/602, top_share 1.0; random: PR 0.55 |
| `_3d_noise`, t = 0.032–0.038 | q4 0.13, q6 0.41–0.42, Q4 0.012–0.016, Q6 0.058–0.074; PR 0.32 over 1142 modes, top 12 modes 8–10%. There is local packing order above random, but no long-range orientation and no Bragg spots: amorphous, unlike the 2D `_new_repul_2` crystal (ψ₆ 0.99 local and global, 6 spots) |
