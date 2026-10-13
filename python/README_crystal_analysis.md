# Cluster-crystal analysis tools

These scripts go with `cluster_crystals.tex`, which describes what the 2D and 3D
interaction runs of `exec/dean_kow/interaction` show near the instability
threshold, and how to interpret it. Three scripts work on plotfiles, and two are
theory tools that need only the potential's parameters.

| script | input | what it does |
| --- | --- | --- |
| `sk_growth.py` | 2D or 3D plotfiles | S(k) on the critical shell against time; onset times and fitted growth rate against linear theory |
| `sf_spots.py` | 2D or 3D plotfiles | time-averaged full S(k): shell profile, Bragg-spot statistics, brightest spots and the angles between them |
| `order_reference.py` | none | reference Steinhardt q4, q6, Q4, Q6 for noisy and polycrystalline bcc, fcc and sc lattices |
| `amplitude_eqs.py` | potential | weakly nonlinear coefficients γ, g, h (or γ_eff, G) and the saddle-node D_SN for each lattice |
| `lattice_continuation.py` | potential | follows a hexagonal, stripe, bcc or fcc steady state of the deterministic equation in D, with its free energy |

`ip_kernels.py` is a shared module, not a script. It holds the potentials
(gema, hk, morse, pp, as in `src_deankow/interaction/Potentials.H`), their 2D and
3D Fourier transforms, and a reader for `amr.ip_*` in an inputs file.

The scripts need only numpy. The plotfile tools use `stack_plotfiles_time.py`,
and `order_reference.py` uses `cluster_diagnostics_3d.py`, from this directory.
Every script has `--help`. Run the examples from `exec/dean_kow/interaction`:

```
P=../../../python
```

## sk_growth.py: linear or nonlinear onset

```
python3 $P/sk_growth.py _new_2d_10_vvclose/plt_fv_* --num-part 327680 \
    --D 0.545 --sigma 5.61 -o vvclose_sk.txt
```

- **Normalization.** S(k) = |Σ n e^{−ik·x}|²/N with n = Nφ/Σφ, so a Poisson
  field has S = 1. `--num-part` must be the run's `num_part`.
- **Shell.** By default it is the |m| of the largest S in the last frame
  (k = 2πm on the unit box), of width `--width`. Use `--shell 8` to fix it.
- **Linear theory.** With `--D` and `--sigma` (the "max growth rate sigma" of
  the startup output), it reports the time S first exceeds each `--thresholds`
  value, against noise-driven linear growth dS/dt = 2σS + 2k²D from S = 0. It
  also fits the growth rate over `--fit-range` and gives it as a ratio to σ.
- **Reading it:**
  - **Ratio about 1:** linear growth.
  - **Ratio well above 1:** the quadratic (triad) term has taken over, as it
    does for cluster crystals near D_crit.
  - **Linearly stable run (σ < 0):** the linear prediction is "never". S
    should sit near 1/(1 − Λ) unless clusters nucleate.
- **Output.** `-o` writes the per-frame table: t, step, shell mean and max S,
  overall max S and its |m|, and φ_max/mean.

## sf_spots.py: crystal or ring, and which lattice

```
python3 $P/sf_spots.py near_3d/plt_fv3d_* --every 10
python3 $P/sf_spots.py _new_2d_10_vvclose/plt_fv_12* --every 4
```

It averages |FFT(φ/φ̄ − 1)|² over the frames (`--every`, `--last`) and prints:

- **Shell profile:** mean and max S per |m| bin, up to `--mmax`.
- **Peak-shell statistics:** the modes in the shell, max/mean, the share of
  the shell's S in the brightest 12, 24 and 48 modes, and the participation
  ratio PR = (ΣS)²/(nΣS²).
- **Brightest ± pairs:** each one's integer mode (x, y[, z]), its S relative
  to the brightest, and the angles to the brightest six, folded to ≤ 90°.

How to read it:

| quantity | isotropic ring (liquid, amorphous) | single crystal |
| --- | --- | --- |
| max/mean on the peak shell | ≈ ln n | ≈ n/(number of spots) |
| PR | ≈ 0.3–0.6 | ≈ (number of spots)/n |
| top-12 share | ≈ 12/n | ≈ 1 |

The spot angles identify the lattice:

| lattice | first-shell spots | angles between spots | next shells (in units of the first) |
| --- | --- | --- | --- |
| 2D hexagonal | 6 | 60° | √3, 2 |
| bcc | 12 ⟨110⟩ | 60°, 90° | √2, √3, 2 |
| fcc | 8 ⟨111⟩ | 70.5° | 1.155, 1.633 |

Several spot sets at fixed rotations to each other mean grains or twins. A
bcc Σ3 twin is a 180° rotation about ⟨111⟩, which maps (6,6,0) to (2,2,8).
`--save S.npy` keeps the averaged S for further analysis.

## order_reference.py: what q4, q6, Q4, Q6 should be

```
python3 $P/order_reference.py --lattice bcc --cells 6 --nn-cut 0.2
python3 $P/order_reference.py --lattice bcc fcc --cells 6 --nn 12
```

This evaluates noisy lattices (`--noise`, as a fraction of the nearest-neighbour
distance), Voronoi polycrystals (`--grains`) and random points. It uses the same
`steinhardt()` as `cluster_diagnostics_3d.py`, so the columns compare directly
with `order.txt`.

- **Match the neighbour rule.** Use the same rule as the diagnostics run.
  For bcc, `--nn-cut` between the second shell (a) and the third (√2 a) is
  best, because 12 neighbours cuts through its 8 + 6 first two shells.
- **Match the density.** Choose `--cells` so the number of points matches the
  cluster count. For example, 6 gives 432 bcc points, the lattice that fits
  the 3D runs.
- **Reading it.** Thermal noise lowers q6 and Q6 together. Misoriented grains
  lower Q6 much more than q6. A single noisy crystal has Q6/q6 ≈ 0.9, two
  random grains about 0.45, and four about 0.3. fcc keeps a large global Q4;
  bcc's is small.

## amplitude_eqs.py: how far above D_crit a lattice survives

```
python3 $P/amplitude_eqs.py --dim 2 --inputs inputs_fv --k 8 --dcrit 0.5472
python3 $P/amplitude_eqs.py --dim 3 --inputs inputs_fv_3d
```

It reduces the equation, for equal amplitudes on a lattice's critical shell, to
dA/dt = σA + γ_eff A² − G A³, and prints γ_eff, G, the amplitude at D_crit, and
D_SN = D_crit + γ_eff²/(4Gk_c²).

- **Lattices.**
  - 2D: hex, stripes, squares.
  - 3D: bcc, rods, fcc, sc, lamellae.
- **Potential.** Set it with `--inputs` or the `--ip-type`/`--eps`/... options.
- **k_c and D_crit.** By default they are continuum values. `--k 8 --dcrit 0.5472`
  uses the box mode (8,0) and the code's discrete D_crit.
- **Reading it.**
  - γ_eff = 0: no triads on the shell. The lattice is supercritical if G > 0.
  - G < 0: a slaved mode is nearly unstable, and one-shell theory fails.
  - For 2D hexagons it also prints g and h; g = h = γ/2 gives
    D_SN = (7/6) D_crit.

## lattice_continuation.py: the full deterministic branch

```
python3 $P/lattice_continuation.py --lattice hex --inputs inputs_fv --k 8 \
    --D 0.50 0.55 0.60 0.62 0.64 0.645 0.6475 0.65
python3 $P/lattice_continuation.py --lattice bcc --inputs inputs_fv_3d \
    --D 0.30 0.36 0.40 0.42 0.43 0.44 0.45 --tmax 0.06
python3 $P/lattice_continuation.py --lattice fcc --inputs inputs_fv_3d --scale 1.05 \
    --D 0.30 0.36 0.39 0.42 --tmax 0.06
```

It solves the noise-free equation pseudo-spectrally on one periodic cell of the
lattice, starting from Gaussian clusters. It then runs `--tmax` at each D in
turn, starting from the previous state.

- **Columns.** For each D it prints φ_max, φ_min, the first-harmonic
  amplitude A, the free energy relative to the uniform state (per unit volume,
  μ₀ = 1), and the residual max|∂φ/∂t|.
- **Reading it.**
  - **Steady:** a residual below about 1e-6.
  - **Near or past the saddle-node:** a growing residual, then decay to
    uniform (the scan stops).
  - **D_F:** the free-energy column changes sign there, the D up to which the
    lattice is the global minimum.
- **Options.** `--scale` changes the lattice spacing; fcc prefers about 1.05
  for the 3D gema kernel. `--n` sets the grid.
- **Cost.** 2D cells take seconds per D. 3D cells take 5–20 minutes per D
  for `--tmax 0.06` at the default `--dt 2e-6`.
