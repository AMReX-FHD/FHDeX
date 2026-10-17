#!/usr/bin/env python3
"""Histogram of cluster masses from the clusters.txt of cluster_diagnostics.py
(2D) or cluster_diagnostics_3d.py (3D).

All clusters of all frames in the chosen range are pooled. Masses are the
"mass" column: the fraction of the total mass in the cluster. The output file
has one line per bin: the bin centre and the fraction of hits in the bin, so the
fractions sum to 1.

Options:
  --t-range T0 T1 or --step-range S0 S1   frames to use (default: all)
  --bin W         bin width in mass (default: (max - min)/50)
  --mmin, --mmax  histogram range (default: the range of the data)
  --no-wide       drop clusters flagged wide (spanning most of the box)

Example:
    python3 cluster_mass_hist.py near_3d/diag2/clusters.txt --bin 2e-5 -o mass_hist.txt
"""

import argparse

import numpy as np


def read_clusters(path):
    """Column name -> array for a clusters.txt file."""
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                names = line[1:].split()
                break
    data = np.loadtxt(path, ndmin=2)
    return {n: data[:, i] for i, n in enumerate(names)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("clusters", help="clusters.txt")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--t-range", type=float, nargs=2, metavar=("T0", "T1"))
    g.add_argument("--step-range", type=float, nargs=2, metavar=("S0", "S1"))
    ap.add_argument("--bin", type=float, help="bin width (default: (max - min)/50)")
    ap.add_argument("--mmin", type=float, help="lower end of the histogram")
    ap.add_argument("--mmax", type=float, help="upper end of the histogram")
    ap.add_argument("--no-wide", action="store_true", help="drop clusters flagged wide")
    ap.add_argument("-o", "--output", default="mass_hist.txt")
    args = ap.parse_args()

    c = read_clusters(args.clusters)
    keep = np.ones(len(c["mass"]), dtype=bool)
    if args.t_range:
        keep &= (c["t"] >= args.t_range[0]) & (c["t"] <= args.t_range[1])
    if args.step_range:
        keep &= (c["step"] >= args.step_range[0]) & (c["step"] <= args.step_range[1])
    if args.no_wide and "wide" in c:
        keep &= c["wide"] == 0
    m = c["mass"][keep]
    if m.size == 0:
        raise SystemExit("no clusters in the selected range")
    nframes = np.unique(c["step"][keep]).size

    lo = args.mmin if args.mmin is not None else m.min()
    hi = args.mmax if args.mmax is not None else m.max()
    w = args.bin if args.bin else (hi - lo)/50 if hi > lo else max(lo, 1e-12)*0.1
    nb = max(1, int(np.ceil((hi - lo)/w)))
    edges = lo + w*np.arange(nb + 1)
    if edges[-1] <= hi:              # make the top value fall inside the last bin
        edges = np.append(edges, edges[-1] + w)
    counts, edges = np.histogram(m, bins=edges)
    frac = counts/counts.sum()
    centres = 0.5*(edges[:-1] + edges[1:])

    with open(args.output, "w") as f:
        f.write("# cluster mass histogram from %s: %d clusters in %d frames, t = %.6g to %.6g\n"
                % (args.clusters, m.size, nframes, c["t"][keep].min(), c["t"][keep].max()))
        f.write("# bin width %.6g; mean mass %.6g, median %.6g; %d hits outside the range\n"
                % (w, m.mean(), np.median(m), m.size - counts.sum()))
        f.write("# bin_centre fraction\n")
        for x, p in zip(centres, frac):
            f.write("%.8g %.8g\n" % (x, p))
    print("%d clusters in %d frames, %d bins of width %.4g -> %s"
          % (m.size, nframes, len(centres), w, args.output))


if __name__ == "__main__":
    main()
