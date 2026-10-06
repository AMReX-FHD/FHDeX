#!/usr/bin/env python3
"""Cluster diagnostics and radial distribution functions for 3D AMReX plotfiles.

The 3D (x, y, z all spatial) analog of cluster_diagnostics.py, for example for
the 3D runs of src_deankow/interaction. For each plotfile (level 0, periodic
domain assumed) it finds clusters, tracks them from frame to frame, and
accumulates two radial distribution functions, writing plain text tables:

  summary.txt        per frame: cluster count, mass fraction in clusters,
                     largest mass, and merge/split/form/dissolve events
  clusters.txt       per frame and cluster: id, mass, centre (xc, yc, zc),
                     Rg, peak, cells, wide flag
  events.txt         one line per merger or split
  rdf_field.txt      pair correlation g(r) of the density field
  rdf_clusters.txt   g(r) of the cluster centres of mass
  sf.txt             shell-averaged structure factor S(k) of the density field
  order.txt          per frame: Steinhardt bond order q4, q6 of the cluster
                     centres (local and global) and Bragg-spot statistics of S(k)

A cluster is a face-connected (6-neighbour) region, wrapping across the
periodic boundaries, where phi > threshold * mean(phi), kept if its mass is at
least mmin times the total mass (default 1e-4, smaller than in 2D because 3D
boxes hold many more clusters). Cluster centres use the circular mean in each
direction. Labelling and the per-cluster sums are vectorized, so 128^3
plotfiles take seconds.

The field g(r) uses cell counts n = phi * num_part * dV: their circular
autocorrelation (by FFT) minus the self pairs, normalized so that a Poisson
field gives g = 1, binned over minimum-image displacements with each bin
normalized by its exact number of lattice displacements. The cluster g(r)
bins minimum-image distances between centres, normalized by the ideal-gas
pair count in each spherical shell, accumulated over frames. Both are limited
to r <= L/2, where every shell lies inside the minimum-image cell.

A regular cluster lattice and an amorphous packing can have the same g(r),
which only measures distances. The Steinhardt bond order
q_l,i = sqrt(4 pi/(2l+1) sum_m |<Y_lm(bond)>|^2) over the nearest neighbours of
each centre measures local packing (fcc: q4 0.19, q6 0.57; random points
about 0.28 for q6); Q_l, from <Y_lm> averaged over all clusters, measures
long-range orientational order and is about the local value in a single
crystal but near 0 in an amorphous state. A crystal also shows a few Bragg
spots on the main S(k) shell instead of a ring.

Example:
    python3 cluster_diagnostics_3d.py plt_fv3d_* --num-part 81920 -o diag3d
"""

import argparse
import os
import sys
from math import factorial

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stack_plotfiles_time import read_header, read_level0, write_plotfile  # noqa: E402


def min_image(d, length):
    """Wrap displacements into [-L/2, L/2)."""
    return (d + 0.5*length) % length - 0.5*length


def label_periodic_3d(mask):
    """Face-connected (6-neighbour) components of a 3D mask with periodic wrap.

    Max-label propagation over the six neighbours until nothing changes; the
    number of sweeps is about the largest cluster's diameter in cells.
    Returns (labels, n) with labels 1..n and 0 outside the mask."""
    lab = np.where(mask, np.arange(1, mask.size + 1).reshape(mask.shape), 0)
    while True:
        new = lab
        for axis in range(3):
            for shift in (1, -1):
                new = np.maximum(new, np.roll(lab, shift, axis=axis))
        new = np.where(mask, new, 0)
        if np.array_equal(new, lab):
            break
        lab = new
    uniq, inv = np.unique(lab, return_inverse=True)
    inv = inv.reshape(lab.shape)
    # uniq[0] is 0 when some cell is outside the mask
    if uniq.size and uniq[0] == 0:
        return inv, uniq.size - 1
    return inv + 1, uniq.size


def find_clusters(phi, coords, dv, lengths, threshold, mmin):
    """Label clusters and compute their properties (vectorized over cells).

    coords are the cell-centre coordinates (x, y, z), each broadcastable to
    phi's [nz, ny, nx] shape. Returns (labels, clusters, mass fraction in
    clusters, mass fraction above the threshold in regions dropped by mmin)."""
    total = phi.sum()*dv
    mean = total/np.prod(lengths)
    lab, n = label_periodic_3d(phi > threshold*mean)
    w = (phi*dv).ravel()
    lr = lab.ravel()
    mass = np.bincount(lr, weights=w, minlength=n + 1)
    keep = np.nonzero(mass[1:] >= mmin*total)[0] + 1
    dropped = (mass[1:].sum() - mass[keep].sum())/total
    relabel = np.zeros(n + 1, dtype=np.int64)
    relabel[keep] = np.arange(1, len(keep) + 1)
    lab = relabel[lab]
    lr = lab.ravel()
    nk = len(keep)
    if nk == 0:
        return lab, [], 0.0, dropped

    m = np.bincount(lr, weights=w, minlength=nk + 1)[1:]
    sel = lr > 0
    ls = lr[sel] - 1
    ws = w[sel]
    centres, d2 = [], np.zeros(nk)
    dmax = np.zeros(nk)
    wide = np.zeros(nk, dtype=bool)
    for k in range(3):
        xk = np.broadcast_to(coords[k], phi.shape).ravel()[sel]
        th = 2*np.pi*xk/lengths[k]
        s = np.bincount(ls, weights=ws*np.sin(th), minlength=nk)
        c = np.bincount(ls, weights=ws*np.cos(th), minlength=nk)
        xc = (lengths[k]/(2*np.pi)*np.arctan2(s, c)) % lengths[k]
        d = min_image(xk - xc[ls], lengths[k])
        d2 += np.bincount(ls, weights=ws*d*d, minlength=nk)
        dmax[:] = 0.0
        np.maximum.at(dmax, ls, np.abs(d))
        wide |= dmax >= 0.45*lengths[k]
        centres.append(xc)
    peak = np.full(nk, -np.inf)
    np.maximum.at(peak, ls, phi.ravel()[sel])
    ncells = np.bincount(ls, minlength=nk)
    clusters = [{"mass": m[i], "x": centres[0][i], "y": centres[1][i], "z": centres[2][i],
                 "rg": np.sqrt(d2[i]/m[i]), "peak": peak[i], "ncells": int(ncells[i]),
                 "wide": bool(wide[i])} for i in range(nk)]
    return lab, clusters, m.sum()/total, dropped


def track(lab_old, old, phi_old, lab_new, new, dv, lengths, search_radius, next_id,
          frac=0.1, mass_keep=0.7):
    """Match clusters between frames; the rules of cluster_diagnostics.track,
    with 3D minimum-image distances.

    Each old cluster is assigned to the new cluster holding the largest part
    of its old mass (if that part exceeds frac of it), or else to the nearest
    new centre within search_radius. A new cluster with two or more old
    clusters assigned is a merger if it holds at least mass_keep of their
    summed mass (old clusters that would break that are dissolved instead).
    An old cluster whose mass is split significantly over two or more new
    clusters is a split. Sets "id" on each new cluster and returns
    (events, next_id, counts)."""
    counts = {"merge": 0, "split": 0, "form": 0, "dissolve": 0}
    events = []
    if lab_old is None:
        for c in new:
            c["id"] = next_id
            next_id += 1
        return events, next_id, counts

    no, nn = len(old), len(new)
    both = (lab_old > 0) & (lab_new > 0)
    key = lab_old[both]*(nn + 1) + lab_new[both]
    ov = np.bincount(key, weights=phi_old[both]*dv,
                     minlength=(no + 1)*(nn + 1)).reshape(no + 1, nn + 1)[1:, 1:]

    def dist(c1, c2):
        return np.sqrt(sum(min_image(c1[q] - c2[q], lengths[k])**2
                           for k, q in enumerate(("x", "y", "z"))))

    assign = {}
    for b in range(no):
        if nn == 0:
            break
        a = int(np.argmax(ov[b]))
        if ov[b, a] > frac*old[b]["mass"]:
            assign[b] = (a, 0.0)
        else:
            d = [dist(old[b], c) for c in new]
            a = int(np.argmin(d))
            if d[a] <= search_radius:
                assign[b] = (a, d[a])

    groups = {}
    for b, (a, d) in assign.items():
        groups.setdefault(a, []).append((b, d))
    for a, members in groups.items():
        members.sort(key=lambda bd: (bd[1], -old[bd[0]]["mass"]))
        kept, msum = [], 0.0
        for b, d in members:
            if not kept or new[a]["mass"] >= mass_keep*(msum + old[b]["mass"]):
                kept.append(b)
                msum += old[b]["mass"]
            else:
                del assign[b]
        groups[a] = kept

    split_children = {}
    for b in range(no):
        children = np.nonzero(ov[b] > frac*old[b]["mass"])[0]
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


def neighbours(xs, lengths, nn, rcut=None):
    """Bonds from each cluster centre to its nn nearest neighbours (minimum
    image), or to every neighbour closer than rcut if rcut is set. Returns the
    centre index of each bond and the bond vectors x_j - x_i."""
    nc, dim = xs.shape
    d = xs[None, :, :] - xs[:, None, :]
    for a in range(dim):
        d[..., a] = min_image(d[..., a], lengths[a])
    r = np.sqrt((d*d).sum(axis=-1))
    np.fill_diagonal(r, np.inf)
    if rcut:
        i, j = np.nonzero(r < rcut)
    else:
        k = min(nn, nc - 1)
        j = np.argsort(r, axis=1)[:, :k].ravel()
        i = np.repeat(np.arange(nc), k)
    return i, d[i, j]


def ylm(l, ct, ph):
    """Orthonormal spherical harmonics Y_lm(cos theta, phi) for m = 0..l (rows);
    Y_l,-m = (-1)^m conj(Y_lm) gives the rest."""
    dP = np.polynomial.legendre.Legendre.basis(l)
    st = np.sqrt(np.maximum(1.0 - ct*ct, 0.0))
    out = np.empty((l + 1, len(ct)), dtype=complex)
    for m in range(l + 1):
        norm = np.sqrt((2*l + 1)/(4*np.pi)*factorial(l - m)/factorial(l + m))
        out[m] = norm*(-1)**m*st**m*dP(ct)*np.exp(1j*m*ph)
        dP = dP.deriv()
    return out


def steinhardt(clusters, lengths, nn, rcut=None, ls=(4, 6)):
    """Steinhardt bond order of the cluster centres. For each l, returns q_l of
    each cluster, q_l,i = sqrt(4 pi/(2l + 1) sum_m |q_lm,i|^2) with
    q_lm,i = (1/n_i) sum_j Y_lm(bond ij) (NaN without neighbours), and the
    global Q_l from q_lm averaged over all clusters."""
    nc = len(clusters)
    if nc < 2:
        return {l: (np.full(nc, np.nan), np.nan) for l in ls}
    xs = np.array([[c["x"], c["y"], c["z"]] for c in clusters])
    i, d = neighbours(xs, lengths, nn, rcut)
    r = np.sqrt((d*d).sum(axis=1))
    ct = d[:, 2]/r
    ph = np.arctan2(d[:, 1], d[:, 0])
    nb = np.bincount(i, minlength=nc).astype(float)
    has = nb > 0
    nb[~has] = np.nan
    out = {}
    for l in ls:
        y = ylm(l, ct, ph)
        q = np.array([(np.bincount(i, weights=ym.real, minlength=nc)
                       + 1j*np.bincount(i, weights=ym.imag, minlength=nc))/nb for ym in y])
        w = np.full(l + 1, 2.0)
        w[0] = 1.0                                  # m and -m have the same |q_lm|
        norm = 4*np.pi/(2*l + 1)
        ql = np.sqrt(norm*(w[:, None]*np.abs(q)**2).sum(axis=0))
        qbar = q[:, has].mean(axis=1) if has.any() else np.full(l + 1, np.nan)
        out[l] = (ql, np.sqrt(norm*(w*np.abs(qbar)**2).sum()))
    return out


class FieldRDF:
    """Pair correlation g(r) of cell counts, binned over minimum-image displacements."""

    def __init__(self, shape, dx, dr, rmax):
        # displacement indices wrapped into [-n/2, n/2), per direction (z, y, x order)
        iz = np.fft.fftfreq(shape[0])*shape[0]
        iy = np.fft.fftfreq(shape[1])*shape[1]
        ix = np.fft.fftfreq(shape[2])*shape[2]
        r = np.sqrt((iz[:, None, None]*dx[2])**2 + (iy[None, :, None]*dx[1])**2
                    + (ix[None, None, :]*dx[0])**2)
        self.inside = r < rmax
        self.bins = np.floor(r[self.inside]/dr).astype(np.int64)
        nb = int(np.ceil(rmax/dr))
        self.count = np.bincount(self.bins, minlength=nb)
        self.r = np.bincount(self.bins, weights=r[self.inside], minlength=nb)/np.maximum(self.count, 1)
        self.shape = shape
        self.frames = []

    def add(self, n):
        ntot = n.sum()
        fk = np.fft.rfftn(n)
        corr = np.fft.irfftn(fk*np.conj(fk), s=self.shape)
        corr[0, 0, 0] -= ntot   # remove self pairs
        g = corr*n.size/(ntot*(ntot - 1.0))
        gb = np.bincount(self.bins, weights=g[self.inside], minlength=len(self.count))
        self.frames.append(gb/np.maximum(self.count, 1))


class StructureFactor:
    """Shell-averaged static structure factor S(k) = |n_hat(k)|^2 / N_tot of the
    cell counts n, for k != 0 (a Poisson field gives S = 1).

    Wave vectors come from fftfreq, so index i > n/2 is the negative wavenumber
    i - n, and k_d = 2 pi (index)/L_d. Shell j holds the modes with |k| in
    [(j - 1/2) dk, (j + 1/2) dk); only complete shells, up to the smallest
    Nyquist wavenumber pi/dx_d, are kept. Each shell is normalized by its exact
    number of modes."""

    def __init__(self, shape, h, dk):
        # shape and h (cell sizes) in array-axis order
        kmag2 = np.zeros(shape)
        for ax, (n, hh) in enumerate(zip(shape, h)):
            k = 2*np.pi*np.fft.fftfreq(n, d=hh)
            sh = [1]*len(shape)
            sh[ax] = n
            kmag2 = kmag2 + k.reshape(sh)**2
        kmag = np.sqrt(kmag2)
        knyq = min(np.pi/hh for hh in h)
        jmax = int(np.floor(knyq/dk - 0.5))
        bins = np.rint(kmag/dk).astype(np.int64)
        self.valid = (kmag > 0) & (bins >= 1) & (bins <= jmax)
        self.bins = bins[self.valid]
        self.count = np.bincount(self.bins, minlength=jmax + 1)
        self.k = np.bincount(self.bins, weights=kmag[self.valid], minlength=jmax + 1)/np.maximum(self.count, 1)
        self.dk = dk
        self.shape, self.h = tuple(shape), tuple(h)
        self.frames = []

    def add(self, n, keep=True, ntop=6):
        """S(k) of one frame, kept for the shell average if keep. Returns the
        peak-shell statistics of the frame (see peak_shell)."""
        f = np.fft.fftn(n)
        s = (f.real**2 + f.imag**2)/n.sum()
        self.last = s
        sv = s[self.valid]
        sb = np.bincount(self.bins, weights=sv, minlength=len(self.count))/np.maximum(self.count, 1)
        if keep:
            self.frames.append(sb)
        return self.peak_shell(sv, sb, ntop)

    def full_shifted(self):
        """S of the last frame on the full k grid, shifted so that k = 0 is at
        index n//2 (indices above n/2 of the FFT become negative wavenumbers),
        with the k = 0 mode set to 0."""
        s = self.last.copy()
        s.flat[0] = 0.0
        return np.fft.fftshift(s)

    def peak_shell(self, sv, sb, ntop):
        """Bragg spots or a ring: on the shell with the largest mean S, the mean
        |k|, the number of modes n, the participation ratio
        PR = (sum S)^2/(n sum S^2) (about (number of spots)/n for a crystal,
        O(1) for an isotropic ring), and the share of the shell's S in its ntop
        largest modes."""
        j = 1 + int(np.argmax(sb[1:]))
        on = sv[self.bins == j]
        tot = on.sum()
        if tot <= 0.0:                              # a uniform field
            return self.k[j], len(on), np.nan, np.nan
        pr = tot*tot/(len(on)*(on*on).sum())
        top = np.sort(on)[::-1][:ntop].sum()/tot
        return self.k[j], len(on), pr, top


class ClusterRDF:
    """g(r) of cluster centres, accumulated over frames, spherical shells."""

    def __init__(self, lengths, dr, rmax):
        self.lengths = lengths
        self.edges = np.arange(0.0, rmax + 0.5*dr, dr)
        self.edges = self.edges[self.edges <= rmax + 1e-12]
        if self.edges[-1] < rmax - 1e-12:
            self.edges = np.append(self.edges, rmax)
        self.pairs = np.zeros(len(self.edges) - 1)
        self.expected = np.zeros(len(self.edges) - 1)
        self.volume = lengths[0]*lengths[1]*lengths[2]
        self.shell = 4.0*np.pi/3.0*(self.edges[1:]**3 - self.edges[:-1]**3)

    def add(self, clusters):
        nc = len(clusters)
        if nc < 2:
            return
        xs = np.array([[c["x"], c["y"], c["z"]] for c in clusters])
        iu, ju = np.triu_indices(nc, k=1)
        d = xs[iu] - xs[ju]
        for k in range(3):
            d[:, k] = min_image(d[:, k], self.lengths[k])
        r = np.sqrt((d**2).sum(axis=1))
        self.pairs += np.histogram(r, bins=self.edges)[0]
        self.expected += 0.5*nc*(nc - 1)*self.shell/self.volume


def write_sf_full(out, sk, sfac, t, step):
    """Write a shifted full S(k) as a plotfile on the k grid: cell centres at
    k_d = 2 pi (i - n_d//2)/L_d, so cell size 2 pi/L_d. Components S and
    log10 S (0 where S = 0, i.e. at k = 0)."""
    nx = sfac.shape[::-1]
    dk = [2*np.pi/(n*hh) for n, hh in zip(sfac.shape, sfac.h)][::-1]
    plo = [-(n//2 + 0.5)*d for n, d in zip(nx, dk)]
    logs = np.log10(np.where(sk > 0, sk, 1.0))
    write_plotfile(out, np.stack([sk, logs]), ["S", "log10S"], plo, dk, t, step)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfiles", nargs="+", help="3D plotfile directories")
    ap.add_argument("--num-part", type=float, required=True,
                    help="num_part of the run (sets the cell counts for the field g(r))")
    ap.add_argument("--var", default="phi0", help="density variable (default phi0)")
    ap.add_argument("--threshold", type=float, default=5.0,
                    help="cluster cells have phi > threshold * mean(phi) (default 5)")
    ap.add_argument("--mmin", type=float, default=1.e-4,
                    help="minimum cluster mass as a fraction of the total (default 1e-4; "
                         "about 0.05/(number of clusters) separates clusters from noise)")
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
    ap.add_argument("--dk", type=float,
                    help="structure-factor shell width (default: 2 pi / L, the smallest box wavenumber)")
    ap.add_argument("--sf-per-frame", action="store_true",
                    help="also write the shell-averaged S(k) of every frame in the window")
    ap.add_argument("--sf-full", action="store_true",
                    help="write the full S(k) on the k grid, averaged over the frames in the window, "
                         "as the plotfile sf_full_avg (k = 0 at the centre)")
    ap.add_argument("--sf-full-per-frame", action="store_true",
                    help="also write the full S(k) of each frame in the window as sf_full_<step>")
    ap.add_argument("--nn", type=int, default=12,
                    help="bond order: number of nearest neighbours of each cluster (default 12)")
    ap.add_argument("--nn-cut", type=float,
                    help="bond order: use all neighbours closer than this instead of --nn, "
                         "e.g. the first minimum of the cluster g(r)")
    ap.add_argument("--top-modes", type=int, default=12,
                    help="order.txt top_share: number of largest modes on the S(k) peak shell (default 12)")
    ap.add_argument("-o", "--output", default="cluster_diag_3d", help="output directory")
    args = ap.parse_args()
    if args.every < 1 or args.threshold <= 0.0 or args.mmin < 0.0:
        ap.error("--every must be >= 1, --threshold > 0 and --mmin >= 0")

    files = sorted(((read_header(p), p) for p in args.plotfiles), key=lambda hp: hp[0]["step"])
    files = files[::args.every]
    h0 = files[0][0]
    if h0["dim"] != 3:
        sys.exit("only 3D plotfiles are supported; use cluster_diagnostics.py for 2D")
    if args.var not in h0["varnames"]:
        sys.exit("unknown variable %s; available: %s" % (args.var, " ".join(h0["varnames"])))
    ivar = h0["varnames"].index(args.var)
    dlo, dhi = h0["domain"]
    nxyz = [dhi[d] - dlo[d] + 1 for d in range(3)]
    shape = (nxyz[2], nxyz[1], nxyz[0])
    dx = h0["dx"][:3]
    plo = h0["prob_lo"][:3]
    lengths = [h0["prob_hi"][d] - plo[d] for d in range(3)]
    dv = dx[0]*dx[1]*dx[2]
    x = plo[0] + (np.arange(nxyz[0]) + 0.5)*dx[0]
    y = plo[1] + (np.arange(nxyz[1]) + 0.5)*dx[1]
    z = plo[2] + (np.arange(nxyz[2]) + 0.5)*dx[2]
    coords = (x[None, None, :], y[None, :, None], z[:, None, None])

    rcap = 0.5*min(lengths)
    rmax = rcap if args.rmax is None else args.rmax
    if rmax > rcap:
        print("warning: --rmax %g exceeds L/2 = %g; using L/2" % (rmax, rcap), file=sys.stderr)
        rmax = rcap
    if rmax <= 0.0:
        ap.error("--rmax must be positive")
    dr = args.dr if args.dr else min(dx)
    drc = args.dr_clusters if args.dr_clusters else min(lengths)/50.0
    frdf = FieldRDF(shape, dx, dr, rmax)
    sfac = StructureFactor(shape, (dx[2], dx[1], dx[0]), args.dk if args.dk else 2*np.pi/max(lengths))
    crdf = ClusterRDF(lengths, drc, rmax)
    t0, t1 = args.t_range if args.t_range else (-np.inf, np.inf)
    search_radius = args.search_radius if args.search_radius else 0.25*min(lengths)

    os.makedirs(args.output, exist_ok=True)
    fsum = open(os.path.join(args.output, "summary.txt"), "w")
    fcl = open(os.path.join(args.output, "clusters.txt"), "w")
    fev = open(os.path.join(args.output, "events.txt"), "w")
    fsum.write("# t step nclusters mass_fraction largest_mass merges splits forms dissolves\n")
    fcl.write("# t step id mass xc yc zc Rg peak ncells wide q4 q6\n")
    ford = open(os.path.join(args.output, "order.txt"), "w")
    ford.write("# bond order over %s; S(k) peak shell, top_share over the %d largest modes\n"
               % ("neighbours closer than %g" % args.nn_cut if args.nn_cut else
                  "the %d nearest neighbours" % args.nn, args.top_modes))
    ford.write("# t step nclusters q4 q6 Q4_global Q6_global kpeak nmodes PR top_share\n")
    fev.write("# t step type old_ids old_masses new_ids new_masses\n")

    lab_old, old, phi_old, next_id, times = None, [], None, 1, []
    sf_sum, sf_step = None, 0
    warned_mmin = False
    for h, p in files:
        if h["domain"] != h0["domain"] or h["dx"] != h0["dx"]:
            sys.exit("%s: domain or dx differs from %s" % (p, files[0][1]))
        if h["finest_level"] > 0:
            print("warning: %s has finer levels; only level 0 is used" % p, file=sys.stderr)
        phi = read_level0(p, h)[ivar]
        lab, clusters, frac_in, dropped = find_clusters(phi, coords, dv, lengths, args.threshold, args.mmin)
        if dropped > 0.1 and dropped > frac_in and not warned_mmin:
            print("warning: t = %.5g: %.2f of the mass is above the threshold in regions lighter than"
                  " --mmin %g (only %.2f is in clusters); use --mmin of about 0.05/(number of clusters)"
                  % (h["time"], dropped, args.mmin, frac_in), file=sys.stderr)
            warned_mmin = True
        events, next_id, counts = track(lab_old, old, phi_old, lab, clusters, dv, lengths,
                                        search_radius, next_id)
        lab_old, old, phi_old = lab, clusters, phi

        t = h["time"]
        st = steinhardt(clusters, lengths, args.nn, args.nn_cut)
        for ic, c in enumerate(clusters):
            c["q4"], c["q6"] = st[4][0][ic], st[6][0][ic]
        keep = t0 <= t <= t1
        ncount = phi*args.num_part*dv
        kp, nmodes, pr, top = sfac.add(ncount, keep, args.top_modes)
        if keep and (args.sf_full or args.sf_full_per_frame):
            sk = sfac.full_shifted()
            if args.sf_full:
                sf_sum = sk if sf_sum is None else sf_sum + sk
                sf_step = h["step"]
            if args.sf_full_per_frame:
                write_sf_full(os.path.join(args.output, "sf_full_%07d" % h["step"]), sk, sfac, t, h["step"])
        ok = ~np.isnan(st[6][0])
        ford.write("%.10g %d %d %.6f %.6f %.6f %.6f %.8g %d %.6g %.6f\n"
                   % (t, h["step"], len(clusters),
                      st[4][0][ok].mean() if ok.any() else np.nan,
                      st[6][0][ok].mean() if ok.any() else np.nan,
                      st[4][1], st[6][1], kp, nmodes, pr, top))
        largest = max((c["mass"] for c in clusters), default=0.0)
        fsum.write("%.10g %d %d %.6f %.6f %d %d %d %d\n"
                   % (t, h["step"], len(clusters), frac_in, largest, counts["merge"],
                      counts["split"], counts["form"], counts["dissolve"]))
        for c in sorted(clusters, key=lambda c: c["id"]):
            fcl.write("%.10g %d %d %.6g %.8f %.8f %.8f %.6g %.6g %d %d %.6f %.6f\n"
                      % (t, h["step"], c["id"], c["mass"], c["x"], c["y"], c["z"], c["rg"],
                         c["peak"], c["ncells"], int(c["wide"]), c["q4"], c["q6"]))
        for kind, oid, om, nid, nm in events:
            fmt = lambda v: ",".join(str(e) for e in v) if isinstance(v, list) else str(v)
            fmtm = lambda v: ",".join("%.4g" % e for e in v) if isinstance(v, list) else "%.4g" % v
            fev.write("%.10g %d %s %s %s %s %s\n"
                      % (t, h["step"], kind, fmt(oid), fmtm(om), fmt(nid), fmtm(nm)))

        if keep:
            frdf.add(ncount)
            crdf.add(clusters)
            times.append(t)
        print("t = %.5g step %d: %d clusters, %.3f of the mass in clusters, %s"
              % (t, h["step"], len(clusters), frac_in,
                 ", ".join("%d %s" % (v, k) for k, v in counts.items() if v) or "no events"))
    fsum.close()
    fcl.close()
    fev.close()
    ford.close()

    if not times:
        sys.exit("no frames in --t-range")
    if args.sf_full:
        write_sf_full(os.path.join(args.output, "sf_full_avg"), sf_sum/len(times), sfac, times[-1], sf_step)
        print("wrote %s/sf_full_avg (full S(k) averaged over %d frames)" % (args.output, len(times)))
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
    with open(os.path.join(args.output, "sf.txt"), "w") as f:
        f.write("# shell-averaged structure factor S(k) = |n_hat|^2/N over %d frames, t in [%g, %g],"
                " num_part = %g, dk = %g\n" % (len(times), times[0], times[-1], args.num_part, sfac.dk))
        f.write("# k(mean over shell) S nmodes" +
                ("".join(" S(t=%.6g)" % t for t in times) if args.sf_per_frame else "") + "\n")
        sm = np.mean(sfac.frames, axis=0)
        for j in range(1, len(sm)):
            if sfac.count[j] == 0:
                continue
            f.write("%.8g %.8g %d" % (sfac.k[j], sm[j], sfac.count[j]))
            if args.sf_per_frame:
                f.write("".join(" %.8g" % s[j] for s in sfac.frames))
            f.write("\n")
    with open(os.path.join(args.output, "rdf_clusters.txt"), "w") as f:
        f.write("# cluster-centre g(r) over %d frames, t in [%g, %g], dr = %g, rmax = %g\n"
                % (len(times), times[0], times[-1], drc, rmax))
        f.write("# r_lo r_hi g pairs expected\n")
        for b in range(len(crdf.pairs)):
            g = crdf.pairs[b]/crdf.expected[b] if crdf.expected[b] > 0 else 0.0
            f.write("%.8g %.8g %.8g %d %.6g\n" % (crdf.edges[b], crdf.edges[b + 1], g,
                                                 crdf.pairs[b], crdf.expected[b]))
    print("wrote %s/{summary,clusters,events,rdf_field,rdf_clusters,sf,order}.txt" % args.output)


if __name__ == "__main__":
    main()
