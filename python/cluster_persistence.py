#!/usr/bin/env python3
"""Persistence and direction of cluster motion over the last frames of a 3D run.

Reads clusters.txt from cluster_diagnostics_3d.py (columns t step id mass xc
yc zc ...; unit periodic box) and follows every cluster present (by id) in
each of the last --intervals + 1 frames. For each it computes

  straightness = |net displacement| / path length,

about 1 for steady drift and about 1/sqrt(K) for a random walk of K steps,
the net speed, and the cosine between the net displacement and the pull of
all other clusters, sum_j m_j r_ij (r_ij the minimum-image vector to j, in
the first of those frames): the direction a long-range attraction would
move it. It prints the median straightness and net speed (also in cells of
--n-cell) and the mean cosine and the fraction of clusters with a positive
one.

Example:
    python3 cluster_persistence.py _ppc1_driftmerge_3d_B/diag/clusters.txt --intervals 6
"""

import argparse

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("clusters", help="clusters.txt from cluster_diagnostics_3d.py")
    ap.add_argument("--intervals", type=int, default=6,
                    help="number K of last frame intervals to follow (default 6)")
    ap.add_argument("--n-cell", type=int, default=128,
                    help="cells per side, to express the net displacement in cells (default 128)")
    a = ap.parse_args()

    rows = np.array([[float(v) for v in line.split()[:7]]
                     for line in open(a.clusters) if not line.startswith("#")])
    ts = np.unique(rows[:, 0])[-a.intervals-1:]
    frames = [{int(r[2]): r for r in rows[rows[:, 0] == t]} for t in ts]
    ids = set.intersection(*[set(f) for f in frames])

    st, cs, sp = [], [], []
    for k in ids:
        x = np.array([f[k][4:7] for f in frames])
        steps = np.diff(x, axis=0)
        steps -= np.round(steps)
        net = steps.sum(0)
        path_len = np.linalg.norm(steps, axis=1).sum()
        if path_len == 0:
            continue
        st.append(np.linalg.norm(net)/path_len)
        sp.append(np.linalg.norm(net)/(ts[-1]-ts[0]))
        pull = np.zeros(3)
        for j, r in frames[0].items():
            if j != k:
                d = r[4:7] - frames[0][k][4:7]
                d -= np.round(d)
                pull += r[3]*d
        cs.append(net @ pull/np.linalg.norm(net)/np.linalg.norm(pull))
    st, cs, sp = map(np.array, (st, cs, sp))
    span = ts[-1] - ts[0]
    print(f"t {ts[0]:.5f}->{ts[-1]:.5f} ({len(ts)-1} intervals), {len(st)} clusters: "
          f"straightness median {np.median(st):.2f} (random walk ~{1/np.sqrt(len(ts)-1):.2f}); "
          f"net speed median {np.median(sp):.2f} = {np.median(sp)*span*a.n_cell:.2f} cells; "
          f"cos(net, pull) mean {cs.mean():+.2f}, fraction > 0 {np.mean(cs > 0):.2f}")


if __name__ == "__main__":
    main()
