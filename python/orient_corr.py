#!/usr/bin/env python3
"""Bond-orientational correlation C6(r) of 3D cluster centres.

Reads clusters.txt from cluster_diagnostics_3d.py (columns t step id mass xc
yc zc ...; unit periodic box) and, for one frame (default the last), computes
for each cluster q_6m,i = (1/n) sum_j Y_6m(bond ij) over its --nn nearest
neighbours, then

  C6(r) = < Re sum_m q_6m,i conj(q_6m,j) / (|q_6,i| |q_6,j|) >   over pairs at r.

C6 is 1 at all r in a single crystal and decays to near 0 within about one
neighbour shell in an amorphous packing; in a polycrystal it decays over the
grain size. The script also prints the grain-count estimate
(q6 / Q6)^2, the number of randomly oriented grains that would give the
observed global Q6.

Example:
    python3 orient_corr.py _3d_.29/diag/clusters.txt justabove_3d/diag/clusters.txt
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_diagnostics_3d import neighbours, ylm  # noqa: E402


def q6m(xs, nn):
    """Normalized per-cluster q_6m (rows m = 0..6), the weights of m and -m,
    and the local q6 of each cluster."""
    nc = len(xs)
    i, b = neighbours(xs.copy(), [1.0, 1.0, 1.0], nn)
    r = np.linalg.norm(b, axis=1)
    y = ylm(6, b[:, 2]/r, np.arctan2(b[:, 1], b[:, 0]))
    q = np.array([np.bincount(i, weights=ym.real, minlength=nc)
                  + 1j*np.bincount(i, weights=ym.imag, minlength=nc) for ym in y])/nn
    w = np.full(7, 2.0)
    w[0] = 1.0
    norm = np.sqrt((w[:, None]*np.abs(q)**2).sum(axis=0))
    q6 = np.sqrt(4*np.pi/13)*norm
    Q6 = np.sqrt(4*np.pi/13*(w*np.abs(q.mean(axis=1))**2).sum())
    return q/norm, w, q6, Q6


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("clusters", nargs="+", help="clusters.txt files")
    ap.add_argument("--step", type=float, help="frame to use (default: last in each file)")
    ap.add_argument("--nn", type=int, default=14, help="nearest neighbours (14 = bcc 8 + 6)")
    ap.add_argument("--dr", type=float, default=0.04)
    ap.add_argument("--rmax", type=float, default=0.5)
    args = ap.parse_args()

    bins = np.arange(0, args.rmax + 1e-12, args.dr)
    print("# r bin centres:", " ".join(f"{0.5*(a + b):5.2f}" for a, b in zip(bins[:-1], bins[1:])))
    for fn in args.clusters:
        d = np.loadtxt(fn)
        step = args.step if args.step is not None else d[:, 1].max()
        d = d[d[:, 1] == step]
        xs = d[:, 4:7]
        qn, w, q6, Q6 = q6m(xs, args.nn)
        C = np.real(np.einsum("m,mi,mj->ij", w, qn, np.conj(qn)))
        dx = xs[None] - xs[:, None]
        dx -= np.round(dx)
        r = np.linalg.norm(dx, axis=2)
        np.fill_diagonal(r, np.inf)
        prof = []
        for a, b in zip(bins[:-1], bins[1:]):
            m = (r >= a) & (r < b)
            prof.append(f"{C[m].mean():5.2f}" if m.any() else "    -")
        print(f"{fn}: t = {d[0, 0]:.4f}, step {int(step)}, {len(xs)} clusters, "
              f"q6 = {q6.mean():.3f}, Q6 = {Q6:.3f}, grains ~ (q6/Q6)^2 = {(q6.mean()/Q6)**2:.1f}")
        print("  C6(r):", " ".join(prof))


if __name__ == "__main__":
    main()
