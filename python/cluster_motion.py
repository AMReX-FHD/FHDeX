#!/usr/bin/env python3
"""Cluster masses and motion relative to Brownian drift, for 3D cluster runs.

Reads clusters.txt from cluster_diagnostics_3d.py (columns t step id mass xc
yc zc ...; unit periodic box). For the last frame it prints the cluster
count, the largest masses, the median mass and the share of the total
cluster mass in the largest one and the largest five.

Then, for each of the last --pairs pairs of consecutive frames, it compares
the displacement |dx| of every cluster tracked in both (by id, minimum
image) with its Brownian rms displacement sqrt(6 D dt / n_p), where
n_p = mass * num_part is the number of particles in the cluster, and reports
the median and maximum ratio. A ratio of order 1 means the clusters only
diffuse; a ratio much larger than 1 means they drift. It also reports how
the displacement points relative to the largest cluster: the mean cosine of
the angle between them and the fraction of clusters moving toward it.

Example:
    python3 cluster_motion.py _ppc1_collapse_3d_B/diag/clusters.txt \\
        --num-part 83886080 --diff-coeff 1.0
"""

import argparse

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("clusters", help="clusters.txt from cluster_diagnostics_3d.py")
    ap.add_argument("--num-part", type=float, required=True,
                    help="total number of particles in the run")
    ap.add_argument("--diff-coeff", type=float, required=True,
                    help="particle diffusion coefficient D")
    ap.add_argument("--pairs", type=int, default=6,
                    help="number of last frame pairs to report (default 6)")
    a = ap.parse_args()

    rows = np.array([[float(v) for v in line.split()[:7]]
                     for line in open(a.clusters) if not line.startswith("#")])
    ts = np.unique(rows[:, 0])

    def frame(t):
        return {int(r[2]): r for r in rows[rows[:, 0] == t]}

    last = frame(ts[-1])
    m = np.array(sorted([r[3] for r in last.values()], reverse=True))
    print(f"t = {ts[-1]:.5f}: N = {len(m)}, largest masses {np.round(m[:6], 4)}, "
          f"median {np.median(m):.4f}, top-1 share {m[0]/m.sum():.2f}, "
          f"top-5 share {m[:5].sum()/m.sum():.2f}")

    for i in range(max(1, len(ts) - a.pairs), len(ts)):
        f0, f1 = frame(ts[i-1]), frame(ts[i])
        dt = ts[i] - ts[i-1]
        big = max(f1.values(), key=lambda r: r[3])
        ratio, cosang = [], []
        for k in set(f0) & set(f1):
            if k == int(big[2]):
                continue
            d = f1[k][4:7] - f0[k][4:7]
            d -= np.round(d)
            rms = np.sqrt(6*a.diff_coeff*dt/(f0[k][3]*a.num_part))
            ratio.append(np.linalg.norm(d)/rms)
            to_big = big[4:7] - f0[k][4:7]
            to_big -= np.round(to_big)
            if np.linalg.norm(d) > 0:
                cosang.append(d @ to_big/np.linalg.norm(d)/np.linalg.norm(to_big))
        if not ratio:
            continue
        ratio = np.array(ratio)
        cosang = np.array(cosang)
        print(f"  {ts[i-1]:.5f}->{ts[i]:.5f}: {len(ratio)} tracked, |dx|/Brownian rms median "
              f"{np.median(ratio):6.1f} max {ratio.max():7.1f}; cos(displacement, direction to "
              f"largest) mean {cosang.mean():+.2f}, fraction moving toward it {np.mean(cosang > 0):.2f}")


if __name__ == "__main__":
    main()
