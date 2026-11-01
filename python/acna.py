#!/usr/bin/env python3
"""Adaptive common-neighbour analysis and averaged bond order of 3D cluster centres.

Reads clusters.txt from cluster_diagnostics_3d.py (columns t step id mass xc
yc zc ...; unit periodic box), takes one frame (default the last) and
classifies each centre as fcc, hcp, bcc, ico or other by adaptive CNA
(Stukowski, Modelling Simul. Mater. Sci. Eng. 20, 045021, 2012):

  * 12 nearest neighbours, cutoff (1 + sqrt 2)/2 times their mean distance:
    fcc if all 12 bonds have signature (4,2,1), hcp if 6 are (4,2,1) and 6
    are (4,2,2), ico if all 12 are (5,5,5);
  * otherwise 14 nearest neighbours with the bcc cutoff: bcc if 8 bonds are
    (6,6,6) and 6 are (4,4,4).

A signature (ncn, nb, nlcb) is the number of common neighbours of the bond,
the number of bonds among them and the size of the largest connected set of
those bonds. The script also prints the Lechner-Dellago averaged q4bar,
q6bar, w4bar, w6bar per class and the mean distances to the 14 nearest
neighbours, and with --out writes the label of each cluster.

Example:
    python3 acna.py near_3d/diag/clusters.txt --out near_3d/diag/acna.txt
"""

import argparse
from math import factorial, sqrt

import numpy as np

CLASSES = ("bcc", "fcc", "hcp", "ico", "other")


def min_image(d):
    return d - np.round(d)


def frame(path, step=None):
    """Time, step, ids and centres of one frame of clusters.txt (default the last)."""
    rows = [line.split() for line in open(path) if not line.startswith("#")]
    steps = np.array([float(r[1]) for r in rows])
    s = steps.max() if step is None else step
    sel = [r for r, st in zip(rows, steps) if st == s]
    xs = np.array([[float(v) for v in r[4:7]] for r in sel])
    return float(sel[0][0]), s, np.array([int(r[2]) for r in sel]), xs


def longest_chain(bonds):
    """Number of bonds in the largest connected set of bonds (as in OVITO's CNA)."""
    bonds = list(bonds)
    best = 0
    while bonds:
        comp = [bonds.pop()]
        atoms = set(comp[0])
        grew = True
        while grew:
            grew = False
            for b in bonds[:]:
                if b[0] in atoms or b[1] in atoms:
                    comp.append(b)
                    bonds.remove(b)
                    atoms.update(b)
                    grew = True
        best = max(best, len(comp))
    return best


def signatures(vecs, rcut):
    """CNA signatures (ncn, nb, nlcb) of the bonds from the centre to each of vecs."""
    n = len(vecs)
    adj = np.linalg.norm(vecs[:, None] - vecs[None], axis=-1) < rcut
    np.fill_diagonal(adj, False)
    sig = []
    for j in range(n):
        cn = np.flatnonzero(adj[j])
        bonds = [(a, b) for ia, a in enumerate(cn) for b in cn[ia+1:] if adj[a, b]]
        sig.append((len(cn), len(bonds), longest_chain(bonds)))
    return sig


def acna(xs):
    """a-CNA label of each centre, the pair distances and each row's neighbour order."""
    d = min_image(xs[:, None] - xs[None])
    r = np.linalg.norm(d, axis=-1)
    np.fill_diagonal(r, np.inf)
    order = np.argsort(r, axis=1)
    labels = []
    for i in range(len(xs)):
        nb = order[i, :14]
        rr = r[i, nb]
        v = -d[i, nb]   # vectors i -> j
        rc12 = (1 + np.sqrt(2))/2 * rr[:12].mean()
        s12 = signatures(v[:12], rc12)
        n421 = s12.count((4, 2, 1))
        n422 = s12.count((4, 2, 2))
        if n421 == 12:
            labels.append("fcc")
            continue
        if n421 == 6 and n422 == 6:
            labels.append("hcp")
            continue
        if s12.count((5, 5, 5)) == 12:
            labels.append("ico")
            continue
        rc14 = (1 + np.sqrt(2))/2 * (2/np.sqrt(3)*rr[:8].sum() + rr[8:14].sum())/14
        s14 = signatures(v[:14], rc14)
        if s14.count((6, 6, 6)) == 8 and s14.count((4, 4, 4)) == 6:
            labels.append("bcc")
            continue
        labels.append("other")
    return np.array(labels), r, order


def assoc_legendre(l, m, x):
    """Associated Legendre function P_l^m(x), m >= 0, with the Condon-Shortley phase."""
    pmm = np.ones_like(x)
    if m > 0:
        somx2 = np.sqrt((1-x)*(1+x))
        fact = 1.0
        for _ in range(m):
            pmm = -pmm*fact*somx2
            fact += 2
    if l == m:
        return pmm
    pmmp1 = x*(2*m+1)*pmm
    if l == m+1:
        return pmmp1
    for ll in range(m+2, l+1):
        pll = (x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm)/(ll-m)
        pmm, pmmp1 = pmmp1, pll
    return pll


def ylm(l, v):
    """Complex spherical harmonics Y_lm for m = -l..l at unit vectors v (n, 3)."""
    ct = v[:, 2]
    ph = np.arctan2(v[:, 1], v[:, 0])
    out = []
    for m in range(-l, l+1):
        am = abs(m)
        norm = np.sqrt((2*l+1)/(4*np.pi)*factorial(l-am)/factorial(l+am))
        y = norm*assoc_legendre(l, am, ct)*np.exp(1j*am*ph)
        if m < 0:
            y = (-1)**am*np.conj(y)
        out.append(y)
    return np.array(out)


def w3j(j1, j2, j3, m1, m2, m3):
    """Wigner 3j symbol (Racah formula)."""
    if m1 + m2 + m3 != 0:
        return 0.0
    t1 = factorial(j1+j2-j3)*factorial(j1-j2+j3)*factorial(-j1+j2+j3)/factorial(j1+j2+j3+1)
    pre = sqrt(t1*factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2)
               * factorial(j3+m3)*factorial(j3-m3))
    s = 0.0
    for k in range(0, j1+j2+j3+1):
        den = [k, j3-j2+k+m1, j3-j1+k-m2, j1+j2-j3-k, j1-k-m1, j2-k+m2]
        if min(den) < 0:
            continue
        s += (-1)**k/np.prod([factorial(x) for x in den])
    return (-1)**(j1-j2-m3)*pre*s


def wigner3j_sum(l, qlm):
    """sum over m1 + m2 + m3 = 0 of (l l l; m1 m2 m3) q_lm1 q_lm2 q_lm3."""
    w = 0.0
    for m1 in range(-l, l+1):
        for m2 in range(-l, l+1):
            m3 = -m1-m2
            if abs(m3) > l:
                continue
            w = w + w3j(l, l, l, m1, m2, m3)*qlm[m1+l]*qlm[m2+l]*qlm[m3+l]
    return w


def averaged_order(xs, order, nn):
    """Lechner-Dellago averaged (qbar_l, wbar_l) per centre for l = 4, 6, over nn neighbours."""
    d = min_image(xs[:, None] - xs[None])
    q = {}
    for l in (4, 6):
        qlm = np.zeros((len(xs), 2*l+1), complex)
        for i in range(len(xs)):
            v = -d[i, order[i, :nn]]
            v /= np.linalg.norm(v, axis=1)[:, None]
            qlm[i] = ylm(l, v).mean(axis=1)
        qbar = np.array([(qlm[i] + qlm[order[i, :nn]].sum(0))/(nn+1) for i in range(len(xs))])
        mag = np.sqrt(4*np.pi/(2*l+1)*(np.abs(qbar)**2).sum(1))
        wl = np.array([wigner3j_sum(l, qb).real/((np.abs(qb)**2).sum()**1.5) for qb in qbar])
        q[l] = (mag, wl)
    return q


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("clusters", nargs="+", help="clusters.txt files")
    ap.add_argument("--step", type=float, help="frame to analyse (default the last)")
    ap.add_argument("--nn", type=int, default=12,
                    help="neighbours for q4bar/q6bar (12 fcc/hcp, 14 bcc; default 12)")
    ap.add_argument("--out", help="write id, centre, label, q4bar, q6bar of each cluster")
    a = ap.parse_args()
    for path in a.clusters:
        t, step, ids, xs = frame(path, a.step)
        lab, r, order = acna(xs)
        q = averaged_order(xs, order, a.nn)
        rs = np.sort(r, axis=1)
        print(f"== {path}: t={t:.4f} step={step:.0f} N={len(xs)}")
        print("   a-CNA: " + "  ".join(f"{k} {np.mean(lab == k)*100:5.1f}%" for k in CLASSES))
        print("   neighbour distances (mean over clusters): "
              + " ".join(f"{v:.4f}" for v in rs[:, :14].mean(0)))
        for k in CLASSES:
            s = lab == k
            if s.sum() == 0:
                continue
            print(f"   {k:5s} n={s.sum():3d}  q4bar={q[4][0][s].mean():.3f} q6bar={q[6][0][s].mean():.3f} "
                  f"w4bar={q[4][1][s].mean():+.3f} w6bar={q[6][1][s].mean():+.4f}")
        if a.out:
            with open(a.out, "w") as f:
                f.write("# id x y z label q4bar q6bar\n")
                for i in range(len(xs)):
                    f.write(f"{ids[i]} {xs[i, 0]:.6f} {xs[i, 1]:.6f} {xs[i, 2]:.6f} {lab[i]} "
                            f"{q[4][0][i]:.4f} {q[6][0][i]:.4f}\n")


if __name__ == "__main__":
    main()
