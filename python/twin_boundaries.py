#!/usr/bin/env python3
"""Grains of a twinned bcc cluster crystal and the boundaries between them.

Reads clusters.txt files from cluster_diagnostics_3d.py (unit periodic box).
Two grain orientations are distinguished by Bragg vectors (integer box modes
m, k = 2 pi m) that only one grain has. The defaults are the box-aligned bcc
("parent", (6,6,0), (6,0,6), (0,6,6)) and its Sigma3 twin, the 180 degree
rotation about [111] ((2,2,8), (2,8,2), (8,2,2)); the vectors perpendicular to
[111] are shared and are not used.

For each cluster i and grain g the local coherence
  c_g,i = < cos(G . (x_j - x_i)) >  over G of grain g and neighbours j with |x_j - x_i| < --rloc
is about 1 inside grain g and about 0 elsewhere. A cluster is assigned to the
twin if c_twin - c_parent > --margin, to the parent if below -margin, and
otherwise to a boundary. Per frame the script prints:

  * the fractions of parent, twin and boundary clusters;
  * the boundary positions along --axis: zero crossings of the
    Gaussian-smoothed (width --width) profile of c_twin - c_parent, labelled
    P->T or T->P in the direction of increasing coordinate;
  * with --map n, the parent thickness (fraction of the axis where the
    smoothed profile is negative) on an n x n grid over the other two axes;
    a zero means the parent grain is perforated there.

Example:
    python3 twin_boundaries.py near_3d/diag2/clusters.txt --map 6 --every 5
"""

import argparse

import numpy as np


def coherence(xs, gs, rloc):
    dx = xs[None] - xs[:, None]
    dx -= np.round(dx)
    r = np.linalg.norm(dx, axis=2)
    np.fill_diagonal(r, np.inf)
    w = r < rloc
    return np.mean([(np.cos(2*np.pi*dx @ g)*w).sum(axis=1)/w.sum(axis=1) for g in gs], axis=0)


def profile(x, f, nb, width):
    xb = (np.arange(nb) + 0.5)/nb
    dx = xb[:, None] - x[None]
    dx -= np.round(dx)
    k = np.exp(-0.5*(dx/width)**2)
    return xb, (k*f).sum(axis=1)/k.sum(axis=1)


def crossings(xb, p):
    nb = len(xb)
    out = []
    for i in range(nb):
        a, b = p[i], p[(i + 1) % nb]
        if a*b < 0:
            out.append(((xb[i] + a/(a - b)/nb) % 1.0, "P->T" if a < 0 else "T->P"))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("clusters", nargs="+", help="clusters.txt files")
    ap.add_argument("--parent", type=float, nargs="+", default=[6, 6, 0, 6, 0, 6, 0, 6, 6],
                    help="parent-only Bragg modes, flattened triples")
    ap.add_argument("--twin", type=float, nargs="+", default=[2, 2, 8, 2, 8, 2, 8, 2, 2],
                    help="twin-only Bragg modes, flattened triples")
    ap.add_argument("--axis", type=int, default=0, help="axis normal to the boundaries (0 = x)")
    ap.add_argument("--rloc", type=float, default=0.26, help="neighbourhood radius for the coherence")
    ap.add_argument("--margin", type=float, default=0.15)
    ap.add_argument("--width", type=float, default=0.06, help="smoothing width of the profile")
    ap.add_argument("--map", type=int, default=0, help="parent-thickness map on an n x n grid")
    ap.add_argument("--every", type=int, default=1)
    args = ap.parse_args()

    gp = np.array(args.parent).reshape(-1, 3)
    gt = np.array(args.twin).reshape(-1, 3)
    ax = args.axis
    others = [a for a in range(3) if a != ax]
    print("# t step N parent twin boundary | boundary positions along axis", ax)
    for fn in args.clusters:
        d0 = np.loadtxt(fn)
        steps = np.unique(d0[:, 1])[::args.every]
        for s in steps:
            d = d0[d0[:, 1] == s]
            xs = d[:, 4:7]
            f = coherence(xs, gt, args.rloc) - coherence(xs, gp, args.rloc)
            xb, p = profile(xs[:, ax], f, 48, args.width)
            c = crossings(xb, p)
            print(f"{d[0, 0]:.4f} {int(s)} {len(xs)} {np.mean(f < -args.margin):.2f} "
                  f"{np.mean(f > args.margin):.2f} {np.mean(abs(f) <= args.margin):.2f} | "
                  + "  ".join(f"{x:.3f} {k}" for x, k in c))
            if args.map:
                n = args.map
                half = max(0.17, 1.0/n)
                th = np.empty((n, n))
                for i in range(n):
                    for j in range(n):
                        du = xs[:, others[0]] - (i + 0.5)/n
                        dv = xs[:, others[1]] - (j + 0.5)/n
                        du -= np.round(du)
                        dv -= np.round(dv)
                        m = (abs(du) < half) & (abs(dv) < half)
                        _, pm = profile(xs[m, ax], f[m], 60, 0.04)
                        th[i, j] = np.mean(pm < 0)
                print(f"  parent thickness on {n}x{n} grid (rows axis {others[0]}, columns axis {others[1]}): "
                      f"mean {th.mean():.3f}, min {th.min():.3f}")
                for row in th:
                    print("   ", " ".join(f"{v:.2f}" for v in row))


if __name__ == "__main__":
    main()
