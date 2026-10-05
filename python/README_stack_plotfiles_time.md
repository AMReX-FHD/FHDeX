# stack_plotfiles_time.py

Stacks a sequence of 2D AMReX plotfiles into one 3D plotfile, with time as the
third direction (z). The result is a standard single-level AMReX plotfile, so
amrvis (3D build), AMReXplorer, VisIt, ParaView and the AMReX
`Tools/Plotfile` utilities (`fboxinfo`, `fextrema`, `fsnapshot`, `fcompare`,
…) can all read it.

In the space–time volume, a cluster that sits still is a tube along z, a moving
cluster is a tilted or curved tube, and a merger is two tubes joining.

It needs only Python 3 and numpy.

## Usage

```
python3 stack_plotfiles_time.py plt_hk_* -o spacetime_plt --tscale 10
python3 stack_plotfiles_time.py plt_pp_* -o spacetime_pp --vars phi0 --tscale 25 --selftest
```

Then open `spacetime_plt` like any 3D plotfile, e.g.
`amrvis3d...ex spacetime_plt`, or take a space–time slice with
`fsnapshot.gnu.ex -v phi0 -n 1 spacetime_plt`.

## What it does

1. Reads the `Header` of every input plotfile and sorts the files by step
   number, not by name.
2. Checks that all files are 2D and have the same domain, cell size,
   prob_lo/prob_hi, coordinate system and variable list. It stops with an
   error on a mismatch.
3. Reads the level-0 data of each file. Slice k of the output (z cell index k)
   is the k-th file, copied cell for cell; nothing is interpolated or resampled.
4. Writes the 3D plotfile one chunk of z slices at a time, so memory use
   doesn't grow with the number of plotfiles.
5. Writes `times.txt` inside the output plotfile with the actual step and time
   of every slice.

## The time axis

A plotfile needs uniform cell spacing in every direction, but output times are
often slightly uneven: here dt changes with the density, so `plot_int` steps
don't span equal times. So:

- **Default:** the z cell size is the mean interval between outputs,
  dz = tscale·(t_last − t_first)/(nt − 1), and slice k is centred at
  z = tscale·t_first + k·dz. The z range is [tscale·t_first − dz/2,
  tscale·t_first + (nt − ½)·dz].
- **`--tscale s`:** stretches z by s. Typically x and y span [0, 1] while the
  run covers a much shorter time (about 0.05 for the clustering runs), and
  viewers draw axes in true proportion. `--tscale` gives the volume a usable
  aspect ratio; for example, s = 20 makes t ∈ [0, 0.05] span z ∈ [0, 1]. Then
  z = s·t, not t.
- **`--index`:** ignores time, setting z = slice index and dz = 1, so the z range
  is [0, nt]. Use it when output times are very uneven or not increasing.
- **Check:** if any output time differs from the uniform spacing by more than
  1% of the mean interval, the script prints a warning. `times.txt` always
  gives the true time of each slice.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `plotfiles` | (required) | 2D plotfile directories, e.g. `plt_hk_*`. At least two. |
| `-o NAME`, `--output NAME` | `spacetime_plt` | Output plotfile directory. It must not already exist. |
| `--vars V1 V2 …` | all | Variables to keep, e.g. `--vars phi0` |
| `--every K` | 1 | Use every K-th plotfile, after sorting by step |
| `--tscale S` | 1 | Scale factor for the z (time) axis; must be > 0 |
| `--index` | off | Use the slice index as z (dz = 1) instead of time |
| `--kchunk K` | 32 | Number of z slices per output box (FAB). Larger means fewer, bigger boxes and more memory while writing. |
| `--selftest` | off | After writing, read the output back and check that every z slice equals its input plotfile bit for bit |

## Output

The output directory is a single-level 3D plotfile:

| File | Contents |
| --- | --- |
| `Header` | `HyperCLaw-V1.1` header. Variables are those selected, spacedim 3, time = the last input's time, step = the last input's step, prob_lo/prob_hi = the input's x, y extent plus the z range above, domain nx × ny × nt, cell size dx, dy, dz. |
| `Level_0/Cell_H` | VisMF header: one box per chunk of `--kchunk` slices, `((0,0,k0) (nx-1,ny-1,k1) (0,0,0))`, with each box's min and max per variable |
| `Level_0/Cell_D_00000`, … | One file per box: a FAB header line, then the data as little-endian 64-bit reals, x fastest, one variable after another |
| `times.txt` | One line per slice (see below) |

`times.txt` columns:

| Column | Meaning |
| --- | --- |
| k | z cell index of the slice |
| step | Step number of the input plotfile |
| time | Time of the input plotfile |
| z_center | z coordinate of the slice's cell centre in the output |

Data is always written as 64-bit reals, whatever precision the inputs used.

## Assumptions and limitations

- **2D inputs, level 0 only.** If a file has finer levels, a warning is printed
  and they are ignored. Level 0 must cover the whole domain.
- **Uniform z.** The stacked file is evenly spaced in z, while real output
  times may not be. Read true times from `times.txt`.
- **No time interpolation.** Each slice is one plotfile. Fewer plotfiles give a
  coarser space–time volume; use a smaller `plot_int` in the run for finer
  time resolution.
- **Large runs.** The output size is about nx·ny·nt·(number of variables)·8
  bytes. Use `--vars` and `--every` to reduce it.

## Verifying a stacked file

```
python3 stack_plotfiles_time.py plt_* -o st --selftest     # bit-for-bit slice check
fboxinfo.gnu.ex st                                         # domain nx x ny x nt, 100% coverage
fextrema.gnu.ex st                                         # min/max = min/max over the inputs
```

For example, 11 plotfiles of a 128² `inputs_fv_hk_cluster` run (16 boxes each)
stack into a 128 × 128 × 11 file that passes `--selftest`. Its `fextrema` range
equals the range over the inputs.

## See also

`cluster_diagnostics.py` (README: `README_cluster_diagnostics.md`) computes
cluster counts, tracks, mergers and radial distribution functions from the same
plotfiles.
