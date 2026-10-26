#!/usr/bin/env python3
"""Drift-merge or collective collapse? Merger generations and the coherence of
cluster motion, from the clusters.txt and events.txt of cluster_diagnostics.py
(2D) or cluster_diagnostics_3d.py (3D).

In drift-merge each cluster is pulled only by its nearest neighbours: mergers
come in generations, each slower than the last, and stop once neighbours are
farther apart than the tail range. In a collapse the cluster gas is unstable
at the box scale: all clusters move together up the large-scale density
gradient and the mergers cascade through the contracting cloud. Measured here:

Generations (generations.txt, one line per merger). A cluster present in the
first frame, or formed later, has generation 0; a merger product has 1 + the
largest generation of its parents; split pieces keep the parent's generation.
For each merger: the parents' masses and generations, and their separation
when the younger of the two first appeared ("birth separation"), in units of
nn0, the mean nearest-neighbour distance of the first frame. Per generation
the summary gives the count, the first, median and last merger time and the
median birth separation.

Motion (modes.txt, per frame). Cluster velocities v_i are centre displacements
over --lag frames (minimum image) of ids that take part in no event in
between, divided by the elapsed time.
  nn         mean nearest-neighbour distance
  A1         box-mode amplitude of the cluster mass, max over the lowest
             wavevectors (the axes and, in 2D, the diagonals) of
             |sum_j m_j exp(i k.x_j)|/sum_j m_j. About N^-1/2 for clusters
             spread over the box, 1 for all mass in one place.
  cos_nn     mean cosine between v_i and the direction to the nearest neighbour
  cos_box    mean cosine between v_i and the gradient, at x_i, of the cluster
             mass density smoothed to those lowest modes (up the box-scale
             density gradient)
  C1, C2, C3 mean cosine between v_i and v_j over pairs at separations
             0.5-1.5, 1.5-2.5, 2.5-3.5 nn0
Drift-merge: cos_nn near 1 for the pairs that are closing, cos_box about 0,
C2 and C3 about 0 (pairs approaching each other make C1 negative). Collapse:
cos_box well above 0 and C2, C3 non-zero (a converging flow makes clusters
across the cloud anti-correlated), with A1 growing.

Flicker events (a merger undone by the opposite split within --flicker-steps)
are removed as in cluster_tracks.py.

Example:
    python3 merge_modes.py diag/clusters.txt --flicker-steps 2000 -o modes
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_tracks import filter_flicker, read_clusters, read_events  # noqa: E402


def min_image(d, L):
    return d - L*np.round(d/L)


def low_modes(L, dim):
    """Lowest wavevectors: the axes and, in 2D, the two diagonals."""
    ks = []
    for d in range(dim):
        k = np.zeros(dim)
        k[d] = 2*np.pi/L[d]
        ks.append(k)
    if dim == 2:
        ks += [np.array([2*np.pi/L[0], 2*np.pi/L[1]]), np.array([2*np.pi/L[0], -2*np.pi/L[1]])]
    return np.array(ks)


def frame_stats(x, m, v, L, ks, nn0):
    """nn, A1, cos_nn, cos_box, C1..C3 for one frame (v: velocities, nan where unknown)."""
    n = len(x)
    out = dict(nn=np.nan, A1=np.nan, cos_nn=np.nan, cos_box=np.nan, C1=np.nan, C2=np.nan, C3=np.nan)
    if n < 2:
        return out
    d = min_image(x[None, :, :] - x[:, None, :], L)      # d[i, j] = x_j - x_i
    r = np.linalg.norm(d, axis=2)
    np.fill_diagonal(r, np.inf)
    jn = r.argmin(axis=1)
    out["nn"] = r[np.arange(n), jn].mean()

    ph = np.exp(1j*(x @ ks.T))                          # (n, nk)
    rho = (m[:, None]*ph).sum(axis=0)/m.sum()           # (nk,)
    out["A1"] = np.abs(rho).max()

    ok = np.isfinite(v).all(axis=1)
    sp = np.linalg.norm(v, axis=1)
    ok &= sp > 0
    if ok.sum() < 2:
        return out
    u = v[ok]/sp[ok, None]
    dn = d[np.arange(n), jn][ok]
    out["cos_nn"] = np.mean((u*dn).sum(axis=1)/np.linalg.norm(dn, axis=1))
    # gradient of sum_k rho_k exp(-i k.x) + c.c. (the smoothed density, up to a constant)
    g = np.real((1j*np.conj(rho)[None, :, None]*ks[None, :, :]*np.conj(ph)[:, :, None]).sum(axis=1))[ok]
    gn = np.linalg.norm(g, axis=1)
    good = gn > 0
    if good.any():
        out["cos_box"] = np.mean((u[good]*g[good]).sum(axis=1)/gn[good])
    ro = r[np.ix_(ok, ok)]
    uu = u @ u.T
    for b, (lo, hi) in enumerate(((0.5, 1.5), (1.5, 2.5), (2.5, 3.5))):
        sel = (ro >= lo*nn0) & (ro < hi*nn0)
        if sel.any():
            out["C%d" % (b + 1)] = uu[sel].mean()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("clusters", help="clusters.txt")
    ap.add_argument("--events", help="events.txt (default: next to clusters.txt)")
    ap.add_argument("--lengths", type=float, nargs="+", help="box lengths (default: unit box)")
    ap.add_argument("--lag", type=int, default=1, help="frames between the two positions of a velocity")
    ap.add_argument("--flicker-steps", type=float, default=0,
                    help="cancel merge/split pairs of the same clusters within this many steps (0 = off)")
    ap.add_argument("--range", type=float, help="tail range R; the summary reports final nn/R")
    ap.add_argument("-o", "--output", default="modes", help="output directory")
    args = ap.parse_args()

    c = read_clusters(args.clusters)
    axes = [a for a in ("xc", "yc", "zc") if a in c]
    dim = len(axes)
    L = np.array(args.lengths if args.lengths else [1.0]*dim, float)
    if len(L) != dim:
        raise SystemExit("--lengths needs %d values" % dim)
    events = read_events(args.events or os.path.join(os.path.dirname(args.clusters) or ".", "events.txt"))
    ev, held, transient = filter_flicker(events, args.flicker_steps)
    rename = {s2: (b, cc) for s1, s2, a, b, cc in held}   # after a cancelled merge, B continues as C

    pos = np.stack([c[a] for a in axes], axis=1)
    frames = []
    for s in np.unique(c["step"]):
        sel = c["step"] == s
        frames.append((int(s), c["t"][sel][0], c["id"][sel].astype(int), pos[sel], c["mass"][sel]))
    fx = [{i: p for i, p in zip(ids, x)} for _, _, ids, x, _ in frames]
    fm = [{i: mm for i, mm in zip(ids, m)} for _, _, ids, _, m in frames]
    steps = [f[0] for f in frames]

    # generations and birth separations
    gen, born = {}, {}                                   # id -> generation, frame index of first appearance
    for i in frames[0][2]:
        gen[i], born[i] = 0, 0
    merges = []
    touched = {}                                         # frame index -> ids in an event there
    for k in range(1, len(frames)):
        s = steps[k]
        if s in rename:
            b, cc = rename[s]
            gen[cc], born[cc] = gen.get(b, 0), born.get(b, k)
        for kind, old, new in ev.get(s, []):
            touched.setdefault(k, set()).update(old + new)
            if kind == "merge" and len(old) >= 2:
                g = 1 + max(gen.get(o, 0) for o in old)
                kb = max(born.get(o, k - 1) for o in old)
                xb = [fx[kb].get(o) for o in old]
                sep = np.nan
                if len(old) == 2 and all(p is not None for p in xb):
                    sep = np.linalg.norm(min_image(xb[1] - xb[0], L))
                merges.append((frames[k][1], s, g, [fm[k - 1].get(o, np.nan) for o in old],
                               [gen.get(o, 0) for o in old], sep, frames[k][1] - frames[kb][1]))
                gen[new[0]], born[new[0]] = g, k
            elif kind == "split":
                g = gen.get(old[0], 0)
                for n in new:
                    gen[n] = g
                    if n != old[0]:
                        born[n] = k
        for i in frames[k][2]:
            if i not in gen:
                gen[i], born[i] = 0, k

    nn0 = frame_stats(frames[0][3], frames[0][4], np.full_like(frames[0][3], np.nan), L,
                      low_modes(L, dim), 1.0)["nn"]
    os.makedirs(args.output, exist_ok=True)
    with open(os.path.join(args.output, "generations.txt"), "w") as f:
        f.write("# nn0 = %.6g\n# t step generation masses parent_generations birth_sep birth_sep/nn0 time_since_birth\n"
                % nn0)
        for t, s, g, ms, gs, sep, dtb in merges:
            f.write("%.8g %d %d %s %s %.5g %.4g %.5g\n" % (t, s, g, ",".join("%.4g" % v for v in ms),
                                                          ",".join(str(v) for v in gs), sep, sep/nn0, dtb))

    # motion statistics per frame
    ks = low_modes(L, dim)
    rows = []
    for k in range(len(frames)):
        s, t, ids, x, m = frames[k]
        v = np.full_like(x, np.nan)
        k2 = k + args.lag
        if k2 < len(frames):
            bad = set().union(*[touched.get(j, set()) for j in range(k + 1, k2 + 1)])
            dt = frames[k2][1] - t
            for n, i in enumerate(ids):
                if i in fx[k2] and i not in bad and dt > 0:
                    v[n] = min_image(fx[k2][i] - x[n], L)/dt
        st = frame_stats(x, m, v, L, ks, nn0)
        rows.append((t, s, len(ids), m.sum(), st))
    with open(os.path.join(args.output, "modes.txt"), "w") as f:
        f.write("# t step nclusters mass nn A1 cos_nn cos_box C1 C2 C3\n")
        for t, s, n, mt, st in rows:
            f.write("%.8g %d %d %.5g %.5g %.5g %.4g %.4g %.4g %.4g %.4g\n" % (
                t, s, n, mt, st["nn"], st["A1"], st["cos_nn"], st["cos_box"], st["C1"], st["C2"], st["C3"]))

    # summary
    lines = ["%d frames, t = %.5g to %.5g; clusters %d -> %d; nn0 = %.4g" % (
        len(frames), frames[0][1], frames[-1][1], len(frames[0][2]), len(frames[-1][2]), nn0)]
    if args.range:
        lines.append("final nn = %.4g = %.2f R" % (rows[-1][4]["nn"], rows[-1][4]["nn"]/args.range))
    lines.append("generation  mergers  t_first  t_median  t_last  median birth_sep/nn0")
    for g in sorted({mg[2] for mg in merges}):
        ts = np.array([mg[0] for mg in merges if mg[2] == g])
        sp = np.array([mg[5] for mg in merges if mg[2] == g])
        lines.append("%10d %8d %8.5g %9.5g %7.5g %10.3g" % (g, len(ts), ts.min(), np.median(ts), ts.max(),
                                                          np.nanmedian(sp)/nn0 if np.isfinite(sp).any() else np.nan))
    # motion averages over the quarters of the run
    ts = np.array([r[0] for r in rows])
    q = np.linspace(ts[0], ts[-1], 5)
    lines.append("time window          clusters  A1     cos_nn  cos_box  C1     C2     C3")
    for a, b in zip(q[:-1], q[1:]):
        sel = [r for r in rows if a <= r[0] <= b]
        avg = lambda key: np.nanmean([r[4][key] for r in sel]) if sel else np.nan
        lines.append("%8.5g - %-8.5g %8.1f  %.3f  %6.3f  %6.3f  %6.3f %6.3f %6.3f" % (
            a, b, np.mean([r[2] for r in sel]), avg("A1"), avg("cos_nn"), avg("cos_box"),
            avg("C1"), avg("C2"), avg("C3")))
    text = "\n".join(lines) + "\n"
    open(os.path.join(args.output, "summary.txt"), "w").write(text)
    print(text, end="")


if __name__ == "__main__":
    main()
