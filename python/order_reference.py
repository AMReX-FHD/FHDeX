#!/usr/bin/env python3
"""Reference Steinhardt order values for noisy and polycrystalline 3D lattices.

Builds point sets in the unit periodic box and evaluates them with the same
steinhardt() as cluster_diagnostics_3d.py, so the numbers can be compared
directly with its order.txt (q4, q6 local; Q4, Q6 global):

  * perfect lattices (bcc, fcc, sc) of --cells conventional cells per side,
    with Gaussian displacements of each --noise value (a fraction of the
    nearest-neighbour distance), mimicking thermal motion of clusters;
  * polycrystals: --grains Voronoi grains, each a randomly rotated copy of
    the lattice at the same density, with --poly-noise displacement;
  * uniformly random points of the same number.

Thermal noise lowers q6 and Q6 together; misoriented grains lower the
global Q6 much more than the local q6. So Q6/q6 is near 0.9 for a single
noisy crystal and falls toward 0.4/(q6 sqrt(N)) for many small grains. Global
Q4 is large for fcc (0.19) and small for bcc (0.04) even in a single crystal.

Example (the bcc lattice commensurate with the 3D runs, 6^3 cells = 432
clusters, neighbours within 0.2 as in --nn-cut 0.2):
    python3 order_reference.py --lattice bcc --cells 6 --nn-cut 0.2
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_diagnostics_3d import steinhardt  # noqa: E402

BASIS = {"bcc": [[0, 0, 0], [.5, .5, .5]],
         "fcc": [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]],
         "sc": [[0, 0, 0]]}
NN = {"bcc": np.sqrt(3)/2, "fcc": 1/np.sqrt(2), "sc": 1.0}   # nn distance / cell size


def lattice(kind, n):
    g = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)], float)
    return ((g[:, None, :] + np.array(BASIS[kind])[None]).reshape(-1, 3))/n


def random_rotation(rng):
    q = rng.normal(size=4)
    a, b, c, d = q/np.linalg.norm(q)
    return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                     [2*(b*c + a*d), a*a - b*b + c*c - d*d, 2*(c*d - a*b)],
                     [2*(b*d - a*c), 2*(c*d + a*b), a*a - b*b - c*c + d*d]])


def polycrystal(kind, n, grains, rng):
    """Voronoi grains (periodic) of randomly rotated copies of the lattice."""
    m = int(np.ceil(n*np.sqrt(3))) + 1
    big = (lattice(kind, m) - 0.5)*m/n               # same spacing, covers the rotated cube
    seeds = rng.random((grains, 3))
    out = []
    for g in range(grains):
        p = big @ random_rotation(rng).T + 0.5
        p = p[(p >= 0).all(1) & (p < 1).all(1)]
        d = p[:, None, :] - seeds[None]
        d -= np.round(d)
        out.append(p[np.argmin((d*d).sum(2), axis=1) == g])
    return np.vstack(out)


def order(xs, nn, rcut):
    cl = [{"x": x, "y": y, "z": z} for x, y, z in xs]
    o = steinhardt(cl, np.ones(3), nn, rcut)
    return np.nanmean(o[4][0]), np.nanmean(o[6][0]), o[4][1], o[6][1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lattice", nargs="+", default=["bcc", "fcc"], choices=list(BASIS))
    ap.add_argument("--cells", type=int, default=6, help="conventional cells per box side")
    ap.add_argument("--noise", type=float, nargs="+", default=[0, 0.05, 0.10, 0.15],
                    help="displacement rms, as a fraction of the nearest-neighbour distance")
    ap.add_argument("--grains", type=int, nargs="*", default=[1, 2, 4, 8, 16])
    ap.add_argument("--poly-noise", type=float, default=0.10)
    ap.add_argument("--nn", type=int, default=12, help="nearest neighbours (as in cluster_diagnostics_3d)")
    ap.add_argument("--nn-cut", type=float, help="use all neighbours closer than this instead")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    nb = "neighbours closer than %g" % args.nn_cut if args.nn_cut else "%d nearest neighbours" % args.nn
    print("# %s; columns: q4 q6 Q4_global Q6_global" % nb)
    for kind in args.lattice:
        a_nn = NN[kind]/args.cells
        base = lattice(kind, args.cells)
        print("\n%s, %d^3 cells, %d points, nn distance %.4f" % (kind, args.cells, len(base), a_nn))
        for s in args.noise:
            xs = (base + rng.normal(scale=s*a_nn, size=base.shape)) % 1.0
            print("  noise %.2f           %.3f %.3f %.4f %.3f" % (s, *order(xs, args.nn, args.nn_cut)))
        for g in args.grains:
            xs = polycrystal(kind, args.cells, g, rng)
            xs = (xs + rng.normal(scale=args.poly_noise*a_nn, size=xs.shape)) % 1.0
            print("  %2d grains (N=%4d)   %.3f %.3f %.4f %.3f" % (g, len(xs), *order(xs, args.nn, args.nn_cut)))
    xs = rng.random((len(lattice(args.lattice[0], args.cells)), 3))
    print("\nrandom points (N=%d)   %.3f %.3f %.4f %.3f" % (len(xs), *order(xs, args.nn, args.nn_cut)))


if __name__ == "__main__":
    main()
