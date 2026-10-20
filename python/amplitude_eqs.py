#!/usr/bin/env python3
"""Weakly nonlinear amplitude equations for cluster lattices near D_crit.

For phi = mu0 (1 + u), the only nonlinearity of the deterministic
McKean-Vlasov equation is mu0 div(u grad(U * u)), with Fourier form

  N_k = sum_{p+q=k} B(p,q) u_p u_q,  B(p,q) = -(mu0/2)[(k.q) Uhat(q) + (k.p) Uhat(p)].

Put equal amplitudes A on every vector of a lattice's critical shell
(|k| = k_c). Sum the quadratic term on that shell (the resonant triads), slave
every other generated mode, w_K = -N_K/sigma(K), and feed it back. The result
is

  dA/dt = sigma A + gamma_eff A^2 - G A^3,  sigma = k_c^2 (D_crit - D),

and for gamma_eff != 0 the lattice exists up to the saddle-node

  D_SN = D_crit + gamma_eff^2/(4 G k_c^2),

with amplitude gamma_eff/G at D_crit. For 2D hexagons gamma_eff = gamma and
G = g + 2h, where g (self) and h (modes at 120 degrees) are also printed. For
stripes G = g, and A = sqrt(sigma/g) above threshold.

Lattices: 2D hex, stripes, squares; 3D bcc (12 <110>), fcc (8 <111>),
sc (6 <100>), rods (hexagonal columns), lamellae. A lattice with no triads on
its shell has gamma_eff = 0 (supercritical if G > 0). G < 0 means a nearly
unstable slaved mode; one-shell theory is then not enough (fcc for gema in
3D, whose <200> shell is close to critical).

By default k_c and D_crit = -mu0 Uhat_min are continuum values. Use --k (in
units of 2 pi, e.g. 8 for the (8,0) box mode) and --dcrit to evaluate at the
box mode and the code's discrete D_crit instead.

Example (the 2D and 3D gema runs of exec/dean_kow/interaction):
    python3 amplitude_eqs.py --dim 2 --inputs inputs_fv --k 8 --dcrit 0.5472
    python3 amplitude_eqs.py --dim 3 --inputs inputs_fv_3d
"""

import argparse
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ip_kernels import add_kernel_args, kernel_from_args  # noqa: E402


def shell_vectors(name, kc):
    """Critical-shell wavevectors of a lattice, scaled to |k| = kc."""
    if name in ("hex", "rods"):
        v = [[np.cos(t), np.sin(t)] for t in np.arange(6)*np.pi/3]
    elif name == "stripes" or name == "lamellae":
        v = [[1.0, 0.0], [-1.0, 0.0]]
    elif name == "squares":
        v = [[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]]
    else:
        base = {"bcc": (1, 1, 0), "fcc": (1, 1, 1), "sc": (1, 0, 0)}[name]
        out = set()
        for s in itertools.product((1, -1), repeat=3):
            for p in itertools.permutations(np.array(base)*np.array(s)):
                out.add(tuple(p))
        v = sorted(out)
    v = np.array(v, float)
    if name in ("rods", "lamellae"):
        v = np.hstack([v, np.zeros((len(v), 1))])
    return [kc*w/np.linalg.norm(w) for w in v]


def reduce_lattice(S, Uhat, D, mu0=1.0):
    """(gamma_eff, G) of dA/dt = sigma A + gamma_eff A^2 - G A^3 for equal amplitudes on S."""
    def B(p, q):
        k = p + q
        return -0.5*mu0*(k @ q*Uhat(np.linalg.norm(q)) + k @ p*Uhat(np.linalg.norm(p)))

    def sig(k):
        return -(k @ k)*(D + mu0*Uhat(np.linalg.norm(k)))

    def key(v):
        return tuple(np.round(v, 7))

    shell = {key(v) for v in S}
    n2 = {}
    for p in S:
        for q in S:
            K = p + q
            if np.linalg.norm(K) < 1e-9:
                continue                      # B(p,-p) = 0: the mean is not forced
            n2[key(K)] = n2.get(key(K), 0.0) + B(p, q)
    k1 = S[0]
    gam = n2.get(key(k1), 0.0)
    slaved = {K: -val/sig(np.array(K)) for K, val in n2.items() if K not in shell}
    cubic = 0.0
    for p in S:
        q = key(k1 - p)
        if q in slaved:
            cubic += 2*B(p, np.array(q))*slaved[q]
    return gam, -cubic


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dim", type=int, choices=[2, 3], required=True)
    ap.add_argument("--k", type=float, help="critical |k| in units of 2 pi (box mode); default continuum k_c")
    ap.add_argument("--dcrit", type=float, help="D_crit to use (e.g. the code's discrete value); default -mu0 Uhat(k_c)")
    ap.add_argument("--D", type=float, help="D at which to evaluate the slaved-mode denominators (default D_crit)")
    ap.add_argument("--lattices", nargs="+", help="default: hex stripes squares (2D), bcc rods fcc sc lamellae (3D)")
    add_kernel_args(ap)
    args = ap.parse_args()

    ker = kernel_from_args(args)
    dim = args.dim

    def Uh(k):
        return ker.Uhat(k, dim)

    kc_cont, umin = ker.critical(dim)
    kc = 2*np.pi*args.k if args.k else kc_cont
    dcrit = args.dcrit if args.dcrit else -Uh(kc)
    D = args.D if args.D else dcrit
    print("%dD %s" % (dim, ker.describe()))
    print("Uhat(0) = %.5g   continuum k_c = %.4f (k_c R = %.3f), Uhat_min = %.5g" %
          (Uh(0.0), kc_cont, kc_cont*ker.R, umin))
    print("using k_c = %.4f, D_crit = %.5g, slaved denominators at D = %.5g" % (kc, dcrit, D))
    print("Uhat at k_c, sqrt2 k_c, sqrt3 k_c, 2 k_c, (2/sqrt3) k_c: %s" %
          " ".join("%.5g" % Uh(f*kc) for f in (1, np.sqrt(2), np.sqrt(3), 2, 2/np.sqrt(3))))
    if dcrit <= 0:
        print("D_crit <= 0: the uniform state has no finite-k instability for this kernel")
        return
    lattices = args.lattices or (["hex", "stripes", "squares"] if dim == 2
                                 else ["bcc", "rods", "fcc", "sc", "lamellae"])
    print("\n%-9s %5s %12s %12s %10s %12s %10s" % ("lattice", "modes", "gamma_eff", "G", "A(Dcrit)", "D_SN", "D_SN/Dc"))
    for name in lattices:
        if dim == 2 and name not in ("hex", "stripes", "squares"):
            ap.error(name + " is a 3D lattice")
        if dim == 3 and name in ("hex", "stripes", "squares"):
            ap.error(name + " is a 2D lattice (use rods or lamellae in 3D)")
        S = shell_vectors(name, kc)
        gam, G = reduce_lattice(S, Uh, D)
        if abs(gam) > 1e-9*abs(G) and G > 0:
            dsn = dcrit + gam**2/(4*G*kc**2)
            print("%-9s %5d %12.5g %12.5g %10.4f %12.5g %10.4f" % (name, len(S), gam, G, gam/G, dsn, dsn/dcrit))
        elif abs(gam) > 1e-9*abs(G):
            print("%-9s %5d %12.5g %12.5g %10s %12s %10s   G < 0: one-shell theory fails" % (name, len(S), gam, G, "-", "-", "-"))
        else:
            note = "supercritical, A = sqrt(sigma/G)" if G > 0 else "G < 0: one-shell theory fails"
            print("%-9s %5d %12s %12.5g %10s %12s %10s   %s" % (name, len(S), "0", G, "-", "-", "-", note))
    if dim == 2 and "hex" in lattices:
        S = shell_vectors("hex", kc)
        g = reduce_lattice(shell_vectors("stripes", kc), Uh, D)[1]
        G = reduce_lattice(S, Uh, D)[1]
        print("\n2D hexagons: g = %.5g, h = %.5g (G = g + 2h); g = h = gamma/2 gives D_SN = (7/6) D_crit" % (g, (G - g)/2))


if __name__ == "__main__":
    main()
