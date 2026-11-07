#!/usr/bin/env python3
"""Zone statistics for the multilevel Dean-Kawasaki hybrid (src_deankow/multilevel).

For each plotfile the level-0 data (which holds the averaged-down particle
histogram in the particle region) is split into three zones taken from the
plotfile itself: the level-1 grids (particles), the one-cell ring around them
(ghost cells, periodic) and the rest (SPDE). Per zone it reports the
time-averaged particles per cell and the cell variance over the Poisson value
(var(n)/mean(n)), with the standard error over snapshots.

With --pot it also prints the y-averaged x-profile against the Boltzmann
distribution of the potential V = (x-a)^2 (x-b)^2 + ky (y-1/2)^4 (ky from
--pot-ky, default 0.5), rho ~ exp(-2 V/gamma), normalised to the measured total.

Usage:
  multilevel_hybrid_stats.py 'plt_eq_*' [--skip 4]
  multilevel_hybrid_stats.py 'plt_pot_*' --pot 0.25 0.75 2.2e-3 --pot-ky 0 --skip 4
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stack_plotfiles_time import read_header, read_level0, read_cell_h  # noqa: E402


def level1_mask(plt, h, shape):
    """Boolean mask of the level-1 cells on the level-0 index space (ratio 1 in the
    code; the plotfile writer refines level 1 by 2 for display)."""
    mask = np.zeros(shape, dtype=bool)
    d = os.path.join(plt, "Level_1")
    if h["finest_level"] < 1 or not os.path.isdir(d):
        return mask
    boxes, _ = read_cell_h(d)
    for lo, hi in boxes:
        i0, j0 = lo[0] // 2, lo[1] // 2
        i1, j1 = (hi[0] + 1) // 2 - 1, (hi[1] + 1) // 2 - 1
        mask[j0:j1 + 1, i0:i1 + 1] = True
    return mask


def ring_mask(region):
    """Cells within one cell (including diagonals, periodic) of the region, minus it."""
    g = region.copy()
    for dj in (-1, 0, 1):
        for di in (-1, 0, 1):
            g |= np.roll(np.roll(region, dj, axis=0), di, axis=1)
    return g & ~region


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("glob", help="plotfile glob, e.g. 'plt_eq_*'")
    ap.add_argument("--skip", type=int, default=0, help="snapshots to skip at the start (burn-in)")
    ap.add_argument("--pot", nargs=3, type=float, metavar=("ALPHA", "BETA", "GAMMA"),
                    help="compare the x-profile with exp(-2V/gamma)")
    ap.add_argument("--pot-ky", type=float, default=0.5, help="coefficient of (y-1/2)^4 in V")
    args = ap.parse_args()

    files = sorted(f for f in glob.glob(args.glob) if ".old." not in f and os.path.isdir(f))
    files = files[args.skip:]
    if not files:
        sys.exit("no plotfiles match " + args.glob)

    zones = ["particles", "ring", "SPDE"]
    means = {z: [] for z in zones}
    vp = {z: [] for z in zones}
    ncells = {z: [] for z in zones}
    profile = []
    total = []
    for f in files:
        h = read_header(f)
        d = read_level0(f, h)
        phi = d[h["varnames"].index("phi")]
        dx, dy = h["dx"][0], h["dx"][1]
        n = phi*dx*dy                     # particles per cell
        region = level1_mask(f, h, n.shape)
        ring = ring_mask(region)
        rest = ~(region | ring)
        for z, m in zip(zones, (region, ring, rest)):
            if m.sum() == 0:
                continue
            v = n[m]
            means[z].append(v.mean())
            vp[z].append(v.var(ddof=1)/v.mean() if v.mean() > 0 else np.nan)
            ncells[z].append(m.sum())
        profile.append(n.mean(axis=0))
        total.append(n.sum())

    print("# %d snapshots from %s (t = %.4g .. %.4g); total particles %.1f +- %.1f"
          % (len(files), files[0], read_header(files[0])["time"], read_header(files[-1])["time"],
             np.mean(total), np.std(total)))
    print("# zone        cells   <n>/cell  +- s.e.   var/Poisson  +- s.e.")
    for z in zones:
        if not means[z]:
            print("  %-10s  (none)" % z)
            continue
        m = np.array(means[z])
        p = np.array(vp[z])
        k = len(m)
        print("  %-10s %7.1f   %8.4f  %7.4f    %7.4f    %7.4f"
              % (z, np.mean(ncells[z]), m.mean(), m.std(ddof=1)/np.sqrt(k) if k > 1 else 0.,
                 np.nanmean(p), np.nanstd(p, ddof=1)/np.sqrt(k) if k > 1 else 0.))

    if args.pot:
        a, b, g = args.pot
        h = read_header(files[0])
        nx = len(profile[0])
        ny = (h["domain"][1][1] - h["domain"][0][1] + 1)
        dx, dy = h["dx"][0], h["dx"][1]
        x = h["prob_lo"][0] + (np.arange(nx) + 0.5)*dx
        y = h["prob_lo"][1] + (np.arange(ny) + 0.5)*dy
        V = (x[None, :] - a)**2*(x[None, :] - b)**2 + args.pot_ky*(y[:, None] - 0.5)**4
        w = np.exp(-2.*V/g)
        prof = np.mean(profile, axis=0)
        err = np.std(profile, axis=0, ddof=1)/np.sqrt(len(profile))
        exact = w.mean(axis=0)*np.mean(total)/w.sum()   # y-averaged, same total
        print("# x-profile (y-averaged particles per cell), time-averaged, vs exp(-2V/gamma)")
        print("#   i      x       <n>      s.e.     exact    ratio")
        for i in range(nx):
            print("  %3d  %7.4f  %8.3f  %7.3f  %8.3f  %7.4f" % (i, x[i], prof[i], err[i], exact[i],
                                                             prof[i]/exact[i] if exact[i] > 0 else np.nan))
        low = exact < 10.
        print("# cells with exact < 10/cell: measured %.1f, exact %.1f (ratio %.4f); rest %.1f vs %.1f"
              % (prof[low].sum()*ny, exact[low].sum()*ny, prof[low].sum()/exact[low].sum(),
                 prof[~low].sum()*ny, exact[~low].sum()*ny))


if __name__ == "__main__":
    main()
