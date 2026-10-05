# cluster_diagnostics.py

Cluster statistics and radial distribution functions for a sequence of 2D AMReX
plotfiles, written for the Dean–Kawasaki interaction runs in
`exec/dean_kow/interaction` (for example `inputs_fv_pp_merge` and
`inputs_fv_hk_cluster`).

For each plotfile the script:

1. finds the clusters in the density field and measures each one's mass,
   centre, size and peak density;
2. matches them to the clusters in the previous plotfile, giving each cluster a
   stable id and logging mergers, splits, new clusters and clusters that
   dissolve;
3. accumulates two radial distribution functions: g(r) of the density field and
   g(r) of the cluster centres.

The results are plain whitespace-separated text tables with `#` header lines,
ready for gnuplot, numpy or a spreadsheet.

It needs only Python 3 and numpy, and reuses the plotfile reader in
`stack_plotfiles_time.py`, which must be in the same directory.

## Usage

```
python3 cluster_diagnostics.py plt_pp_* --num-part 327680 -o diag
python3 cluster_diagnostics.py plt_pp_* --num-part 327680 -o diag --t-range 0.02 0.03 --rmax 0.4
```

Plotfiles are sorted by the step number in their headers, not by file name.
All of them must have the same domain and cell size.

## Assumptions and limitations

- **2D, periodic domain.** Clusters may wrap across the boundaries, and all
  distances use the minimum-image convention. The plotfile header doesn't record
  periodicity, so the script assumes it.
- **Level 0 only.** Finer levels are ignored, with a warning.
- **Density field only.** The plotfiles contain no particle positions, so both
  g(r) are computed from the density variable (`phi0` by default).
- **Touching clusters count as one.** Two clusters are one cluster as soon as
  the density between them exceeds the threshold. This is the definition of a
  merger used here. Clusters are not split at density saddles.
- **Total mass includes negative densities.** At low particle counts per cell,
  φ can be slightly negative between clusters. The total mass, used for
  `--mmin` and the mass fraction, is the plain sum Σφ·dV, so the mass fraction
  in clusters can slightly exceed 1.

## Definitions

### Cluster

A cluster is a set of cells connected through shared faces, including across
the periodic boundaries, where

  φ > threshold × φ̄,  with φ̄ = Σφ·dV / (domain area).

It is kept only if its mass Σφ·dV is at least `mmin` times the total mass.

For each cluster:

| Quantity | Definition |
| --- | --- |
| mass | Σ φ·dV over its cells |
| centre (xc, yc) | Circular mean in each direction: x_c = L/(2π)·atan2(Σ w sin(2πx/L), Σ w cos(2πx/L)), with w = φ·dV. Correct for clusters straddling the boundary. |
| Rg | Radius of gyration, √(Σ w·\|d\|² / Σ w), with d the minimum-image displacement of each cell from the centre |
| peak | Largest φ in the cluster |
| ncells | Number of cells |
| wide | 1 if some cell is at least 0.45·L from the centre in either direction. The minimum image, and so the centre and Rg, are then ambiguous; this mainly happens for a single cluster that fills the box. |

### Tracking between frames

For each pair of consecutive plotfiles:

1. **Matching.** Each old cluster is assigned to the new cluster that holds the
   largest part of its old mass, measured at the old time over the cells that
   are now in the new cluster, provided that part exceeds 10% of the old
   cluster's mass. If no new cluster does, the old cluster is assigned to the
   nearest new centre within `--search-radius`. The fallback catches the end of
   a merger, when a small cluster falls into a large one within one frame and
   its old cells no longer overlap anything.
2. **Merger.** A new cluster with two or more old clusters assigned is a merger,
   provided it holds at least 70% of their summed mass. Old clusters that would
   break that condition are treated as dissolved, starting with the farthest.
3. **Split.** An old cluster with more than 10% of its mass in each of two or
   more new clusters is a split.
4. **Formed / dissolved.** A new cluster with no old cluster assigned, and not
   split off one, is formed. An old cluster assigned to nothing, and not split,
   is dissolved.
5. **Ids.** A new cluster takes the id of its most massive assigned old cluster.
   Formed clusters and pieces split off get new ids, so ids are never reused.

The 10% and 70% thresholds are fixed in the code (`track(frac=0.1,
mass_keep=0.7)`). Tracking needs plotfiles close enough in time that clusters
move less than about their own size between frames, apart from the final plunge
of a merger. The cadence in `inputs_fv_pp_merge` (`plot_int = 1000`) is
sufficient.

### Density-field g(r)

Built from the cell particle counts n_i = φ_i · num_part · dV:

1. **Autocorrelation:** G(Δ) = Σ_i n_i n_{i+Δ} for every lattice displacement Δ,
   computed by FFT. It is periodic by construction.
2. **Self pairs:** Σ_i n_i is subtracted at Δ = 0.
3. **Normalization:** g(Δ) = G(Δ)·M / (N(N−1)), with M the number of cells and N =
   Σ n_i. A Poisson (ideal gas) field gives g = 1 at every Δ, including Δ = 0.
4. **Binning:** each Δ is wrapped componentwise into [−M_d/2, M_d/2), so r = |Δ·dx|
   is the minimum-image distance. g(r) in a bin is the average of g(Δ) over the
   displacement vectors that fall in it, which normalizes by the exact number of
   lattice vectors per shell rather than 2πr·dr.
5. **Averaging:** the result is averaged over the frames in `--t-range`.

`num_part` enters only through the self-pair correction, which is of relative
size 1/N.

### Cluster-centre g(r)

For every frame in `--t-range` with at least two clusters, the minimum-image
distances between all pairs of cluster centres are histogrammed, counting each
pair once. The result is

  g(r) = Σ_frames pairs in [r, r+dr) / Σ_frames ½ n_c(n_c−1) · π((r+dr)² − r²) / A

with n_c the number of clusters in the frame and A the domain area. The annulus
area is exact up to L/2, because every annulus up to that radius lies inside the
minimum-image cell. Pair counts and expected counts are summed over frames
before dividing, because individual frames often have only about ten clusters.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `plotfiles` | (required) | 2D plotfile directories, e.g. `plt_pp_*` |
| `--num-part N` | (required) | `num_part` of the run. Converts φ to particle counts for the field g(r). |
| `--var NAME` | `phi0` | Plotfile variable used as the density |
| `--threshold C` | 5 | Cluster cells have φ > C × mean φ |
| `--mmin F` | 0.01 | Smallest cluster kept, as a fraction of the total mass |
| `--every K` | 1 | Use every K-th plotfile, after sorting by step |
| `--rmax R` | L/2 | Range of both g(r). Values above L/2 (the shortest side) are reduced to L/2 with a warning. |
| `--dr D` | dx | Bin width of the field g(r). With D ≤ dx the first bin holds only Δ = 0. |
| `--dr-clusters D` | L/50 | Bin width of the cluster-centre g(r) |
| `--search-radius R` | L/4 | Furthest a cluster centre may move between frames and still be matched when its old cells overlap no new cluster |
| `--t-range T0 T1` | all frames | Time window for both g(r) averages. Cluster statistics and tracking always use every frame. |
| `--rdf-per-frame` | off | Also write the field g(r) of each frame in the window |
| `-o DIR`, `--output DIR` | `cluster_diag` | Output directory, created if needed. Existing files are overwritten. |

The script also prints one line per frame with the cluster count, the mass
fraction in clusters and the events since the previous frame.

## Output files

Times are the plotfile times. Lengths and positions are in the plotfile's
physical units. Masses are Σφ·dV, which is a fraction of the total because the
code normalizes φ to integrate to 1.

### summary.txt: one line per plotfile

| Column | Meaning |
| --- | --- |
| t | Time |
| step | Step number |
| nclusters | Number of clusters |
| mass_fraction | Total mass in clusters / total mass |
| largest_mass | Mass of the largest cluster (0 if none) |
| merges | Merger events since the previous plotfile |
| splits | Split events since the previous plotfile |
| forms | Clusters formed since the previous plotfile |
| dissolves | Clusters dissolved since the previous plotfile |

### clusters.txt: one line per cluster per plotfile, sorted by id

| Column | Meaning |
| --- | --- |
| t, step | Time and step |
| id | Cluster id, stable across frames |
| mass | Cluster mass |
| xc, yc | Centre of mass, in [prob_lo, prob_hi) |
| Rg | Radius of gyration |
| peak | Largest φ in the cluster |
| ncells | Number of cells |
| wide | 1 if the cluster spans about half the box or more, so xc, yc and Rg are unreliable |

To follow one cluster, filter by id, e.g. `awk '$3==4' clusters.txt`.

### events.txt: one line per merger or split

| Column | Meaning |
| --- | --- |
| t, step | Time and step of the first plotfile in which the event shows |
| type | `merge` or `split` |
| old_ids | Comma-separated ids before the event: the merging clusters, or the cluster that split |
| old_masses | Their masses at the previous plotfile |
| new_ids | Id(s) after the event: the merged cluster, or the pieces of the split |
| new_masses | Their masses |

Formed and dissolved clusters appear only as counts in `summary.txt`.

### rdf_field.txt: density-field g(r)

The first header line gives the number of frames, the time window,
`num_part`, `dr` and `rmax`. Then one line per non-empty bin:

| Column | Meaning |
| --- | --- |
| r | Mean distance of the lattice displacements in the bin |
| g | g(r) averaged over the frames in the window |
| ndisplacements | Number of lattice displacements in the bin |
| g(t=…) | With `--rdf-per-frame`, one column per frame |

The first bin, at r = 0, is the same-cell value. In clustered states it is
large, because of the cluster's own density.

### rdf_clusters.txt: cluster-centre g(r)

The first header line gives the number of frames, the time window, `dr` and
`rmax`. Then one line per bin:

| Column | Meaning |
| --- | --- |
| r_lo, r_hi | Bin edges |
| g | Σ pairs / Σ expected (0 where nothing is expected) |
| pairs | Cluster pairs in the bin, summed over frames |
| expected | Ideal-gas expectation for the bin, summed over frames |

`pairs` shows how much each bin's g rests on. With a few clusters per frame,
bins with only a handful of pairs are noisy.

## Example

The `inputs_fv_pp_merge` run (256², `num_part = 327680`, seed 11) gives:

- clusters from t ≈ 0.014, 10 of them by t ≈ 0.018;
- 9 mergers between t ≈ 0.027 and 0.036, ending in one cluster with 96–99% of
  the mass;
- with `--t-range 0.02 0.03`, a cluster-centre g(r) with almost no pairs below
  r ≈ 0.2 and a peak over 0.24–0.36 (the expected spacing is about 0.32);
- a field g(r) of about 14 at r = 0, a dip to about 0.2 near r ≈ 0.13 and a
  second peak of about 1.5 at r ≈ 0.3.

The cluster tracks can be viewed in space–time by stacking the same plotfiles
with `stack_plotfiles_time.py`.
