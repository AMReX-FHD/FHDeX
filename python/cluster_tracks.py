#!/usr/bin/env python3
"""Distance each cluster has moved from where it started, from the
clusters.txt and events.txt of cluster_diagnostics.py (2D) or
cluster_diagnostics_3d.py (3D), and the average over all traces.

The diagnostics carry a cluster's id from frame to frame. A trace follows one
id and accumulates its displacement frame by frame (minimum image, so moving
across a periodic boundary is fine as long as a cluster moves less than half
the box between frames). The distance written is |x(t) - x(origin)| of the
unwrapped path.

  * Start: every cluster of the first frame starts a trace. A cluster that
    appears later without a split (a new cluster forming) starts its own
    trace unless --initial-only is given.
  * Split: the parent's trace continues on the child that kept its id. Each
    other child starts a branch trace that inherits the parent's origin and
    path; the branch begins with the parent's last point, so it joins on.
  * Merger: one of the traces on the merging ids continues on the merged id,
    and the others stop. --merge-rule chooses how:
      farther  (default) the trace that has moved farthest continues and moves
               to the merged centre;
      heavier  the trace on the heaviest parent continues and moves to the
               merged centre, a jump of only (other mass/total) x separation;
      nojump   the trace that has moved farthest continues, but takes no jump:
               from then on it adds only the merged cluster's own motion.
    The merged centre is the mass-weighted mean of the parents, so with
    farther or heavier the continuing trace jumps by up to their separation.
  * A cluster that disappears without an event (dissolves) ends its trace.

Flicker (--flicker-steps S). Two clusters joined by a bridge at the detection
threshold are seen as one cluster in some frames and two in others. A
two-way merge and the opposite split of the same pair, within S steps of each
other, cancel:
  * merge (A + B -> A) then split (A -> A + C): B and C are the same cluster.
    While joined, the traces of A and B hold their last positions (the merged
    centre is a mass-weighted mean, not a cluster position). Afterwards the
    trace of B continues on C.
  * split (A -> A + C) then merge (A + C -> A): the piece C is transient; no
    branch is started and the trace stays on A.
Pairs are formed in time order, each event with the next opposite event on the
same kept id. S of about 2000 steps removes the flicker in the 3D gema runs,
where genuine mergers persist for 6000 steps or more.

Output:
  -o (default tracks.txt): for each trace a comment line (start id, parent
     trace, how it ended), then lines "t distance id", and a blank line
     between traces.
  --mean (default: -o with _mean before the extension): per frame
     "t mean_dist median_dist rms_dist mean_dx mean_dy [mean_dz] ntraces" over
     the traces alive in that frame. mean_dist is the average drift distance,
     rms_dist its root mean square (the mean squared displacement is
     rms_dist^2), and mean_d* the average displacement vector, which is about
     0 unless the clusters drift collectively. A merger moves the continuing
     trace to the mass-weighted merged centre, a jump of up to the parents
     separation; the median is insensitive to these few jumps.

Example:
    python3 cluster_tracks.py _new_3d_noise/diag3d2/clusters.txt --flicker-steps 2000 -o tracks.txt
"""

import argparse
import os

import numpy as np


def read_clusters(path):
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                names = line[1:].split()
                break
    data = np.loadtxt(path, ndmin=2)
    return {n: data[:, i] for i, n in enumerate(names)}


def read_events(path):
    """List of (step, type, old_ids, new_ids) in time order."""
    ev = []
    if not os.path.exists(path):
        return ev
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            w = line.split()
            ids = lambda s: [int(v) for v in s.split(",")]
            ev.append((int(w[1]), w[2], ids(w[3]), ids(w[5])))
    return sorted(ev, key=lambda e: e[0])


def filter_flicker(events, smax):
    """Cancel two-way merge/split pairs of the same kept id within smax steps.
    Returns the remaining events (step -> list), the held pairs
    (list of (s1, s2, A, B, C) for merge-then-split) and the transient pieces
    (id -> (s1, s2)) of split-then-merge."""
    keep = list(range(len(events)))
    held, transient = [], {}
    pending = {}                                     # kept id -> index of last unpaired event
    for i, (s, kind, old, new) in enumerate(events):
        two = (kind == "merge" and len(old) == 2) or (kind == "split" and len(new) == 2)
        if not two or smax <= 0:
            continue
        a = new[0] if kind == "merge" else old[0]
        j = pending.get(a)
        if j is not None:
            sj, kj, oj, nj = events[j]
            if kj != kind and s - sj <= smax:
                if kj == "merge":                    # merge A+B, then split A -> A + C
                    b = [x for x in oj if x != a]
                    c = [x for x in new if x != a]
                    if b and c:
                        held.append((sj, s, a, b[0], c[0]))
                        keep.remove(i)
                        keep.remove(j)
                        del pending[a]
                        continue
                else:                                # split A -> A + C, then merge A + C
                    c = [x for x in nj if x != a]
                    if c and c[0] in old:
                        transient[c[0]] = (sj, s)
                        keep.remove(i)
                        keep.remove(j)
                        del pending[a]
                        continue
        pending[a] = i
    out = {}
    for i in keep:
        s, kind, old, new = events[i]
        out.setdefault(s, []).append((kind, old, new))
    return out, held, transient


class Trace:
    def __init__(self, k, cid, x, t, parent=None):
        self.k, self.cid, self.start_id, self.parent = k, cid, cid, parent
        self.x_last = x.copy()                      # last raw (wrapped) position
        if parent is None:
            self.x0, self.xu = x.copy(), x.copy()
            self.pts = [(t, np.zeros_like(x), cid)]
        else:                                       # inherit the parent's state before this frame
            self.x0 = parent.x0.copy()
            self.xu = parent.xu_prev.copy()
            self.pts = [parent.pt_prev]
        self.end = "end of data"
        self.hold_until = None

    def dist(self):
        return float(np.linalg.norm(self.xu - self.x0))

    def move(self, x, L, t):
        d = x - self.x_last
        d -= L*np.round(d/L)
        self.xu = self.xu + d
        self.x_last = x.copy()
        self.pts.append((t, self.xu - self.x0, self.cid))

    def hold(self, t):
        self.pts.append((t, self.xu - self.x0, self.cid))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("clusters", help="clusters.txt")
    ap.add_argument("--events", help="events.txt (default: next to clusters.txt)")
    ap.add_argument("--lengths", type=float, nargs="+", help="box lengths (default: unit box)")
    ap.add_argument("--initial-only", action="store_true",
                    help="only follow clusters present in the first frame (and their split branches)")
    ap.add_argument("--merge-rule", choices=["farther", "heavier", "nojump"], default="farther",
                    help="which trace continues at a merger, and whether it jumps (default farther)")
    ap.add_argument("--flicker-steps", type=float, default=0,
                    help="cancel merge/split pairs of the same clusters within this many steps (0 = off)")
    ap.add_argument("-o", "--output", default="tracks.txt")
    ap.add_argument("--mean", help="per-frame averages over the traces (default: <output>_mean)")
    args = ap.parse_args()

    c = read_clusters(args.clusters)
    axes = [a for a in ("xc", "yc", "zc") if a in c]
    dim = len(axes)
    L = np.array(args.lengths if args.lengths else [1.0]*dim, float)
    if len(L) != dim:
        raise SystemExit("--lengths needs %d values" % dim)
    events = read_events(args.events or os.path.join(os.path.dirname(args.clusters) or ".", "events.txt"))
    ev, held, transient = filter_flicker(events, args.flicker_steps)
    hold_start = {}                                  # step -> [(s2, A, B, C)]
    for s1, s2, a, b, cc in held:
        hold_start.setdefault(s1, []).append((s2, a, b, cc))

    steps = np.unique(c["step"])
    pos = np.stack([c[a] for a in axes], axis=1)
    frames = []
    for s in steps:
        sel = c["step"] == s
        frames.append((int(s), c["t"][sel][0], {int(i): p for i, p in zip(c["id"][sel], pos[sel])},
                       {int(i): m for i, m in zip(c["id"][sel], c["mass"][sel])}))

    traces, on = [], {}                              # on: cluster id -> live trace
    renew = {}                                       # step -> [(trace, new id)] at the end of a hold
    s0, t0, f0, mass_prev = frames[0]
    for cid, x in f0.items():
        on[cid] = Trace(len(traces), cid, x, t0)
        traces.append(on[cid])

    for s, t, fr, mass in frames[1:]:
        merged_into, split_from = {}, {}
        for kind, old, new in ev.get(s, []):
            if kind == "merge":
                merged_into[new[0]] = old
            elif kind == "split":
                for n in new:
                    if n != old[0]:
                        split_from[n] = old[0]
        for tr in on.values():                       # state before this frame, for split branches
            tr.xu_prev, tr.x_prev, tr.pt_prev = tr.xu.copy(), tr.x_last.copy(), tr.pts[-1]

        # flicker: start holds (merge cancelled) and end them (split cancelled)
        held_ids = set()
        for s2, a, b, cc in hold_start.get(s, []):
            for cid in (a, b):
                if cid in on:
                    on[cid].hold_until = s2
            if b in on:                              # B continues as C after the hold
                renew.setdefault(s2, []).append((on[b], cc))
        for tr, cc in renew.pop(s, []):
            if tr.cid in on and on[tr.cid] is tr:
                del on[tr.cid]
            tr.cid = cc
            on[cc] = tr
        for tr in on.values():
            if tr.hold_until is not None and s >= tr.hold_until:
                tr.hold_until = None
            if tr.hold_until is not None:
                held_ids.add(tr.cid)

        # mergers: one trace continues on the merged id (--merge-rule)
        for nid, olds in merged_into.items():
            live = [on[o] for o in olds if o in on]
            if not live:
                continue
            if args.merge_rule == "heavier":
                best = max(live, key=lambda tr: mass_prev.get(tr.cid, 0.0))
            else:
                best = max(live, key=lambda tr: tr.dist())
            if args.merge_rule == "nojump" and nid in fr:
                best.x_last = fr[nid].copy()         # no displacement at the merger itself
            for tr in live:
                del on[tr.cid]
                if tr is not best:
                    tr.end = "merged into id %d at t = %.6g" % (nid, t)
            best.cid = nid
            on[nid] = best

        # clusters that vanished without an event
        for cid in [cid for cid in on if cid not in fr and cid not in held_ids]:
            on.pop(cid).end = "id %d vanished at t = %.6g" % (cid, t)

        for cid, tr in on.items():
            if cid in held_ids:
                tr.hold(t)
            else:
                tr.move(fr[cid], L, t)

        # split branches and newly formed clusters
        for cid, x in fr.items():
            if cid in on or (cid in transient and transient[cid][0] <= s < transient[cid][1]):
                continue
            parent = on.get(split_from.get(cid))
            if parent is not None:
                tr = Trace(len(traces), cid, parent.x_prev, t, parent)
                tr.move(x, L, t)
            elif args.initial_only or cid in split_from:
                continue
            else:
                tr = Trace(len(traces), cid, x, t)
            traces.append(tr)
            on[cid] = tr
        mass_prev = mass

    with open(args.output, "w") as f:
        f.write("# cluster traces from %s: %dD, box %s, flicker window %g steps, merge rule %s;"
                " columns t distance id\n"
                % (args.clusters, dim, " ".join("%g" % v for v in L), args.flicker_steps, args.merge_rule))
        for tr in traces:
            par = "" if tr.parent is None else ", branch of trace %d" % tr.parent.k
            f.write("# trace %d: start id %d%s, ended by %s\n" % (tr.k, tr.start_id, par, tr.end))
            for tt, dx, cid in tr.pts:
                f.write("%.10g %.8g %d\n" % (tt, np.linalg.norm(dx), cid))
            f.write("\n")

    # per-frame averages over the traces alive in each frame
    by_t = {}
    for tr in traces:
        for tt, dx, _ in tr.pts:
            by_t.setdefault(tt, []).append(dx)
    mean_out = args.mean or "%s_mean%s" % os.path.splitext(args.output)
    with open(mean_out, "w") as f:
        f.write("# averages over the traces alive in each frame, from %s\n" % args.clusters)
        f.write("# t mean_dist median_dist rms_dist " + " ".join("mean_d" + a[0] for a in axes) + " ntraces\n")
        for tt in sorted(by_t):
            d = np.array(by_t[tt])
            r = np.linalg.norm(d, axis=1)
            f.write("%.10g %.8g %.8g %.8g %s %d\n" % (tt, r.mean(), np.median(r), np.sqrt((r*r).mean()),
                                                 " ".join("%.6g" % v for v in d.mean(axis=0)), len(d)))

    ends = {}
    for tr in traces:
        key = tr.end.split()[0]
        ends[key] = ends.get(key, 0) + 1
    print("%d frames, %d traces (%s), %d flicker pairs cancelled -> %s, %s"
          % (len(frames), len(traces), ", ".join("%s %d" % kv for kv in sorted(ends.items())),
             len(held) + len(transient), args.output, mean_out))


if __name__ == "__main__":
    main()
