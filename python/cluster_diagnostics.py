#!/usr/bin/env python3
"""Cluster diagnostics and radial distribution functions for 2D AMReX plotfiles.

For each plotfile (level 0, periodic domain assumed) this finds clusters,
tracks them from frame to frame, and accumulates two radial distribution
functions, writing plain text tables to the output directory:

  summary.txt        per frame: cluster count, mass fraction in clusters,
                     largest mass, and merge/split/form/dissolve events
  clusters.txt       per frame and cluster: id, mass, centre, Rg, peak, cells
  events.txt         one line per merger or split
  rdf_field.txt      pair correlation g(r) of the density field
  rdf_clusters.txt   g(r) of the cluster centres of mass

A cluster is a face-connected region (wrapping across the periodic
boundaries) where phi > threshold * mean(phi), kept if its mass is at least
mmin times the total mass. Cluster centres use the circular mean in each
direction, so clusters that straddle the boundary are handled.

The field g(r) uses cell counts n = phi * num_part * dV. Their circular
autocorrelation (by FFT) minus the self pairs, normalized so that a Poisson
field gives g = 1, is binned over minimum-image displacements. Each bin is
normalized by the exact number of lattice displacements in it. The cluster
g(r) bins minimum-image distances between centres and is normalized by the
ideal-gas pair count in each annulus, accumulated over all frames. Both are
limited to r <= L/2, where every annulus lies inside the minimum-image cell.

Example:
    python3 cluster_diagnostics.py plt_pp_* --num-part 327680 -o diag
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stack_plotfiles_time import read_header, read_level0  # noqa: E402


def label_periodic(mask):
    """Face-connected components of a 2D mask with periodic wrap.

    Returns (labels, n) with labels 1..n and 0 outside the mask."""
    ny, nx = mask.shape
    lab = np.zeros(mask.shape, dtype=np.int64)
    n = 0
    for j0, i0 in zip(*np.nonzero(mask)):
        if lab[j0, i0]:
            continue
        n += 1
        lab[j0, i0] = n
        stack = [(j0, i0)]
        while stack:
            j, i = stack.pop()
            for jj, ii in ((j + 1) % ny, i), ((j - 1) % ny, i), (j, (i + 1) % nx), (j, (i - 1) % nx):
                if mask[jj, ii] and not lab[jj, ii]:
                    lab[jj, ii] = n
                    stack.append((jj, ii))
    return lab, n


def min_image(d, length):
    """Wrap displacements into [-L/2, L/2)."""
    return (d + 0.5 * length) % length - 0.5 * length


def find_clusters(phi, x, y, dv, lengths, threshold, mmin):
    """Label clusters and compute their properties.

    Returns (labels, clusters, mass fraction in clusters); labels are 1..n in
    the order of the clusters list."""
    total = phi.sum() * dv
    mean = total / (lengths[0] * lengths[1])
    lab, n = label_periodic(phi > threshold * mean)
    mass = np.bincount(lab.ravel(), weights=phi.ravel(), minlength=n + 1) * dv
    keep = [c for c in range(1, n + 1) if mass[c] >= mmin * total]
    relabel = np.zeros(n + 1, dtype=np.int64)
    relabel[keep] = np.arange(1, len(keep) + 1)
    lab = relabel[lab]

    clusters = []
    for c in range(1, len(keep) + 1):
        sel = lab == c
        w = phi[sel] * dv
        m = w.sum()
        cen, d2, wide = [], 0.0, False
        for coord, length in ((x[sel], lengths[0]), (y[sel], lengths[1])):
            th = 2.0 * np.pi * coord / length
            xc = (length / (2.0 * np.pi) * np.arctan2((w * np.sin(th)).sum(),
                                                      (w * np.cos(th)).sum())) % length
            d = min_image(coord - xc, length)
            d2 += (w * d * d).sum() / m
            wide = wide or np.abs(d).max() >= 0.45 * length
            cen.append(xc)
        clusters.append({"mass": m, "x": cen[0], "y": cen[1],
                         "rg": np.sqrt(d2), "peak": phi[sel].max(),
                         "ncells": int(sel.sum()), "wide": wide})
    return lab, clusters, sum(c["mass"] for c in clusters) / total


def track(lab_old, old, phi_old, lab_new, new, dv, lengths, search_radius, next_id,
          frac=0.1, mass_keep=0.7):
    """Match clusters between frames.

    Each old cluster is assigned to the new cluster holding the largest part of
    its old mass (if that part exceeds frac of it), or else to the nearest new
    centre within search_radius. A new cluster with two or more old clusters
    assigned is a merger, provided it holds at least mass_keep of their summed
    mass; old clusters that would break that are dissolved instead, smallest
    and farthest first. The overlap fallback catches the last stage of a
    merger, when a small cluster falls into a large one within one frame and
    its old cells no longer overlap anything. An old cluster whose mass is
    split significantly over two or more new clusters is a split.

    Sets "id" on each new cluster and returns (events, next_id, counts) where
    counts holds the number of merges, splits, forms and dissolves."""
    counts = {"merge": 0, "split": 0, "form": 0, "dissolve": 0}
    events = []
    if lab_old is None:
        for c in new:
            c["id"] = next_id
            next_id += 1
        return events, next_id, counts

    no, nn = len(old), len(new)
    # ov[b, a]: mass of old cluster b (at the old time) in cells now in new cluster a
    both = (lab_old > 0) & (lab_new > 0)
    key = lab_old[both] * (nn + 1) + lab_new[both]
    ov = np.bincount(key, weights=phi_old[both] * dv,
                     minlength=(no + 1) * (nn + 1)).reshape(no + 1, nn + 1)[1:, 1:]

    def dist(c1, c2):
        dx = min_image(c1["x"] - c2["x"], lengths[0])
        dy = min_image(c1["y"] - c2["y"], lengths[1])
        return np.hypot(dx, dy)

    assign = {}
    for b in range(no):
        a = int(np.argmax(ov[b])) if nn else -1
        if nn and ov[b, a] > frac * old[b]["mass"]:
            assign[b] = (a, 0.0)
        elif nn:
            d = [dist(old[b], c) for c in new]
            a = int(np.argmin(d))
            if d[a] <= search_radius:
                assign[b] = (a, d[a])

    groups = {}
    for b, (a, d) in assign.items():
        groups.setdefault(a, []).append((b, d))
    for a, members in groups.items():
        # keep the overlapping and nearest, most massive old clusters while the
        # new cluster can account for their mass
        members.sort(key=lambda bd: (bd[1], -old[bd[0]]["mass"]))
        kept, msum = [], 0.0
        for b, d in members:
            if not kept or new[a]["mass"] >= mass_keep * (msum + old[b]["mass"]):
                kept.append(b)
                msum += old[b]["mass"]
            else:
                del assign[b]
        groups[a] = kept

    # pieces split off an old cluster, other than the one it was assigned to
    split_children = {}
    for b in range(no):
        children = np.nonzero(ov[b] > frac * old[b]["mass"])[0]
        if len(children) > 1:
            split_children[b] = [int(a) for a in children]
    split_off = {a for b, ch in split_children.items() for a in ch if a != assign.get(b, (None,))[0]}

    for a in range(nn):
        parents = groups.get(a, [])
        if not parents:
            new[a]["id"] = next_id
            next_id += 1
            if a not in split_off:
                counts["form"] += 1
            continue
        # the new cluster carries the id of its largest parent
        new[a]["id"] = old[max(parents, key=lambda b: old[b]["mass"])]["id"]
        if len(parents) > 1:
            counts["merge"] += 1
            events.append(("merge", [old[b]["id"] for b in parents],
                           [old[b]["mass"] for b in parents], new[a]["id"], new[a]["mass"]))

    for b in range(no):
        if b in split_children:
            counts["split"] += 1
            events.append(("split", [old[b]["id"]], [old[b]["mass"]],
                           [new[a]["id"] for a in split_children[b]],
                           [new[a]["mass"] for a in split_children[b]]))
        elif b not in assign:
            counts["dissolve"] += 1
    return events, next_id, counts


class FieldRDF:
    """Pair correlation g(r) of cell counts, binned over minimum-image displacements."""

    def __init__(self, shape, dx, dr, rmax):
        ny, nx = shape
        ix = np.fft.fftfreq(nx) * nx   # displacement indices wrapped into [-n/2, n/2)
        iy = np.fft.fftfreq(ny) * ny
        r = np.sqrt((iy[:, None] * dx[1]) ** 2 + (ix[None, :] * dx[0]) ** 2)
        self.inside = r < rmax
        self.bins = np.floor(r[self.inside] / dr).astype(np.int64)
        nb = int(np.ceil(rmax / dr))
        self.count = np.bincount(self.bins, minlength=nb)
        self.r = np.bincount(self.bins, weights=r[self.inside], minlength=nb) / np.maximum(self.count, 1)
        self.shape = shape
        self.frames = []

    def add(self, n):
        ntot = n.sum()
        fk = np.fft.rfft2(n)
        corr = np.fft.irfft2(fk * np.conj(fk), s=self.shape)
        corr[0, 0] -= ntot   # remove self pairs
        g = corr * n.size / (ntot * (ntot - 1.0))
        gb = np.bincount(self.bins, weights=g[self.inside], minlength=len(self.count))
        self.frames.append(gb / np.maximum(self.count, 1))
        return self.frames[-1]


class ClusterRDF:
    """g(r) of cluster centres, accumulated over frames."""

    def __init__(self, lengths, dr, rmax):
        self.lengths = lengths
        self.edges = np.arange(0.0, rmax + 0.5 * dr, dr)
        self.edges = self.edges[self.edges <= rmax + 1e-12]
        if self.edges[-1] < rmax - 1e-12:
            self.edges = np.append(self.edges, rmax)
        self.pairs = np.zeros(len(self.edges) - 1)
        self.expected = np.zeros(len(self.edges) - 1)
        self.area = lengths[0] * lengths[1]
        self.shell = np.pi * (self.edges[1:] ** 2 - self.edges[:-1] ** 2)

    def add(self, clusters):
        nc = len(clusters)
        if nc < 2:
            return
        xs = np.array([[c["x"], c["y"]] for c in clusters])
        iu, ju = np.triu_indices(nc, k=1)
        d = xs[iu] - xs[ju]
        for k in range(2):
            d[:, k] = min_image(d[:, k], self.lengths[k])
        r = np.sqrt((d ** 2).sum(axis=1))
        self.pairs += np.histogram(r, bins=self.edges)[0]
        self.expected += 0.5 * nc * (nc - 1) * self.shell / self.area


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfiles", nargs="+", help="2D plotfile directories")
    ap.add_argument("--num-part", type=float, required=True,
                    help="num_part of the run (sets the cell counts for the field g(r))")
    ap.add_argument("--var", default="phi0", help="density variable (default phi0)")
    ap.add_argument("--threshold", type=float, default=5.0,
                    help="cluster cells have phi > threshold * mean(phi) (default 5)")
    ap.add_argument("--mmin", type=float, default=0.01,
                    help="minimum cluster mass as a fraction of the total (default 0.01)")
    ap.add_argument("--every", type=int, default=1, help="use every n-th plotfile")
    ap.add_argument("--rmax", type=float, help="RDF range (default and maximum: L/2)")
    ap.add_argument("--dr", type=float, help="field RDF bin width (default: dx)")
    ap.add_argument("--dr-clusters", type=float, help="cluster RDF bin width (default: L/50)")
    ap.add_argument("--search-radius", type=float,
                    help="largest distance a cluster centre may jump between frames when its "
                         "old cells no longer overlap any cluster (default: L/4)")
    ap.add_argument("--t-range", type=float, nargs=2, metavar=("T0", "T1"),
                    help="time window for the RDF averages (default: all frames)")
    ap.add_argument("--rdf-per-frame", action="store_true",
                    help="also write the field g(r) of every frame in the window")
    ap.add_argument("-o", "--output", default="cluster_diag", help="output directory")
    args = ap.parse_args()
    if args.every < 1 or args.threshold <= 0.0 or args.mmin < 0.0:
        ap.error("--every must be >= 1, --threshold > 0 and --mmin >= 0")

    files = sorted(((read_header(p), p) for p in args.plotfiles), key=lambda hp: hp[0]["step"])
    files = files[::args.every]
    h0 = files[0][0]
    if h0["dim"] != 2:
        sys.exit("only 2D plotfiles are supported")
    if args.var not in h0["varnames"]:
        sys.exit("unknown variable %s; available: %s" % (args.var, " ".join(h0["varnames"])))
    ivar = h0["varnames"].index(args.var)
    dlo, dhi = h0["domain"]
    shape = (dhi[1] - dlo[1] + 1, dhi[0] - dlo[0] + 1)
    dx = h0["dx"][:2]
    plo = h0["prob_lo"][:2]
    lengths = [h0["prob_hi"][d] - plo[d] for d in range(2)]
    dv = dx[0] * dx[1]
    x = plo[0] + (np.arange(shape[1]) + 0.5) * dx[0]
    y = plo[1] + (np.arange(shape[0]) + 0.5) * dx[1]
    xx, yy = np.meshgrid(x, y)

    rcap = 0.5 * min(lengths)
    rmax = rcap if args.rmax is None else args.rmax
    if rmax > rcap:
        print("warning: --rmax %g exceeds L/2 = %g; using L/2" % (rmax, rcap), file=sys.stderr)
        rmax = rcap
    if rmax <= 0.0:
        ap.error("--rmax must be positive")
    dr = args.dr if args.dr else min(dx)
    drc = args.dr_clusters if args.dr_clusters else min(lengths) / 50.0
    frdf = FieldRDF(shape, dx, dr, rmax)
    crdf = ClusterRDF(lengths, drc, rmax)
    t0, t1 = args.t_range if args.t_range else (-np.inf, np.inf)
    search_radius = args.search_radius if args.search_radius else 0.25 * min(lengths)

    os.makedirs(args.output, exist_ok=True)
    fsum = open(os.path.join(args.output, "summary.txt"), "w")
    fcl = open(os.path.join(args.output, "clusters.txt"), "w")
    fev = open(os.path.join(args.output, "events.txt"), "w")
    fsum.write("# t step nclusters mass_fraction largest_mass merges splits forms dissolves\n")
    fcl.write("# t step id mass xc yc Rg peak ncells wide\n")
    fev.write("# t step type old_ids old_masses new_ids new_masses\n")

    lab_old, old, phi_old, next_id, times = None, [], None, 1, []
    for h, p in files:
        if h["domain"] != h0["domain"] or h["dx"] != h0["dx"]:
            sys.exit("%s: domain or dx differs from %s" % (p, files[0][1]))
        if h["finest_level"] > 0:
            print("warning: %s has finer levels; only level 0 is used" % p, file=sys.stderr)
        phi = read_level0(p, h)[ivar]
        lab, clusters, frac_in = find_clusters(phi, xx, yy, dv, lengths, args.threshold, args.mmin)
        events, next_id, counts = track(lab_old, old, phi_old, lab, clusters, dv, lengths,
                                        search_radius, next_id)
        lab_old, old, phi_old = lab, clusters, phi

        t = h["time"]
        largest = max((c["mass"] for c in clusters), default=0.0)
        fsum.write("%.10g %d %d %.6f %.6f %d %d %d %d\n"
                   % (t, h["step"], len(clusters), frac_in, largest, counts["merge"],
                      counts["split"], counts["form"], counts["dissolve"]))
        for c in sorted(clusters, key=lambda c: c["id"]):
            fcl.write("%.10g %d %d %.6g %.8f %.8f %.6g %.6g %d %d\n"
                      % (t, h["step"], c["id"], c["mass"], c["x"], c["y"], c["rg"],
                         c["peak"], c["ncells"], int(c["wide"])))
        for kind, oid, om, nid, nm in events:
            fmt = lambda v: ",".join(str(e) for e in v) if isinstance(v, list) else str(v)
            fmtm = lambda v: ",".join("%.4g" % e for e in v) if isinstance(v, list) else "%.4g" % v
            fev.write("%.10g %d %s %s %s %s %s\n"
                      % (t, h["step"], kind, fmt(oid), fmtm(om), fmt(nid), fmtm(nm)))

        if t0 <= t <= t1:
            frdf.add(phi * args.num_part * dv)
            crdf.add(clusters)
            times.append(t)
        print("t = %.5g step %d: %d clusters, %.3f of the mass in clusters, %s"
              % (t, h["step"], len(clusters), frac_in,
                 ", ".join("%d %s" % (v, k) for k, v in counts.items() if v) or "no events"))
    fsum.close()
    fcl.close()
    fev.close()

    if not times:
        sys.exit("no frames in --t-range")
    gf = np.mean(frdf.frames, axis=0)
    with open(os.path.join(args.output, "rdf_field.txt"), "w") as f:
        f.write("# field g(r) over %d frames, t in [%g, %g], num_part = %g, dr = %g, rmax = %g\n"
                % (len(times), times[0], times[-1], args.num_part, dr, rmax))
        f.write("# r(mean over bin) g ndisplacements" +
                ("".join(" g(t=%.6g)" % t for t in times) if args.rdf_per_frame else "") + "\n")
        for b in range(len(gf)):
            if frdf.count[b] == 0:
                continue
            f.write("%.8g %.8g %d" % (frdf.r[b], gf[b], frdf.count[b]))
            if args.rdf_per_frame:
                f.write("".join(" %.8g" % g[b] for g in frdf.frames))
            f.write("\n")
    with open(os.path.join(args.output, "rdf_clusters.txt"), "w") as f:
        f.write("# cluster-centre g(r) over %d frames, t in [%g, %g], dr = %g, rmax = %g\n"
                % (len(times), times[0], times[-1], drc, rmax))
        f.write("# r_lo r_hi g pairs expected\n")
        for b in range(len(crdf.pairs)):
            g = crdf.pairs[b] / crdf.expected[b] if crdf.expected[b] > 0 else 0.0
            f.write("%.8g %.8g %.8g %d %.6g\n" % (crdf.edges[b], crdf.edges[b + 1], g,
                                                 crdf.pairs[b], crdf.expected[b]))
    print("wrote %s/{summary,clusters,events,rdf_field,rdf_clusters}.txt" % args.output)


if __name__ == "__main__":
    main()
