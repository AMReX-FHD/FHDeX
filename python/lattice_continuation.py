#!/usr/bin/env python3
"""Follow a cluster-lattice steady state of the deterministic equation in D.

Solves the noise-free McKean-Vlasov equation

  d phi/dt = D lap phi + div(phi grad C),  C = U * phi,  mean(phi) = mu0 = 1,

pseudo-spectrally on the smallest periodic cell of one lattice:

  hex      2D, rectangular cell 4 pi/k x 4 pi/(sqrt3 k), 2 clusters
  stripes  2D, cell 2 pi/k (x) by 2 pi/k (y), one stripe
  bcc      3D cubic cell a = 2 pi sqrt2/k, 2 clusters (|G_110| = k)
  fcc      3D cubic cell a = 2 pi sqrt3/k, 4 clusters (|G_111| = k)

Here k = scale * k_c, with k_c the continuum critical wavenumber, or 2 pi
times --k. U is periodized over the cell, so the cell's Uhat on its own
Fourier modes is the continuum transform (Poisson summation), and the result
is what an infinite lattice of that spacing would do.

Time stepping is semi-implicit: the part linear about phi = 1,
-k^2 (D + Uhat), is implicit and the rest explicit. Starting from Gaussian
clusters (or a cosine for stripes), each D in --D is run for --tmax, starting
from the previous state. So listing D upward follows the branch to its
saddle-node, where it decays to the uniform state; the scan stops there.

For each D this prints phi_max, phi_min, the first-harmonic amplitude A, the
free energy F = int [D phi ln phi + phi C/2] per volume minus the uniform
value Uhat(0)/2, and the max residual |d phi/dt|. A steady branch has a
residual below about 1e-6; a large residual means it is still decaying.

Example (2D hexagons for the gema kernel of inputs_fv, box mode k = 8):
    python3 lattice_continuation.py --lattice hex --inputs inputs_fv --k 8 \\
        --D 0.50 0.55 0.60 0.62 0.64 0.645 0.65
"""

import argparse
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ip_kernels import add_kernel_args, kernel_from_args  # noqa: E402


def cell(lattice, k):
    """(lengths, sites as fractions of the cell, first-harmonic index)."""
    if lattice == "hex":
        return [4*np.pi/k, 4*np.pi/(np.sqrt(3)*k)], np.array([[0, 0], [.5, .5]]), (2, 0)
    if lattice == "stripes":
        L = 2*np.pi/k
        return [L, L], None, (1, 0)
    if lattice == "bcc":
        a = 2*np.pi*np.sqrt(2)/k
        return [a]*3, np.array([[0, 0, 0], [.5, .5, .5]]), (1, 1, 0)
    if lattice == "fcc":
        a = 2*np.pi*np.sqrt(3)/k
        return [a]*3, np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]), (1, 1, 1)
    raise ValueError(lattice)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lattice", required=True, choices=["hex", "stripes", "bcc", "fcc"])
    ap.add_argument("--D", type=float, nargs="+", required=True, help="diffusion coefficients, in order")
    ap.add_argument("--k", type=float, help="lattice |k| in units of 2 pi (default: continuum k_c)")
    ap.add_argument("--scale", type=float, default=1.0, help="multiply the lattice spacing by this (k -> k/scale)")
    ap.add_argument("--n", type=int, help="cells along the first axis (default 96 in 2D, 32 bcc, 40 fcc)")
    ap.add_argument("--dt", type=float, default=2e-6)
    ap.add_argument("--tmax", type=float, default=0.1, help="time run at each D (default 0.1)")
    ap.add_argument("--dcrit", type=float, help="D_crit for the D/D_crit column (default -mu0 Uhat at the lattice k)")
    add_kernel_args(ap)
    args = ap.parse_args()

    ker = kernel_from_args(args)
    dim = 2 if args.lattice in ("hex", "stripes") else 3
    k = 2*np.pi*args.k if args.k else ker.critical(dim)[0]
    k /= args.scale
    lengths, sites, h1 = cell(args.lattice, k)
    n0 = args.n or {"hex": 96, "stripes": 64, "bcc": 32, "fcc": 40}[args.lattice]
    if args.lattice == "hex":
        shape = [n0, max(8, int(round(n0*lengths[1]/lengths[0])))]
    elif args.lattice == "stripes":
        shape = [n0, 8]
    else:
        shape = [n0]*3
    dx = [L/n for L, n in zip(lengths, shape)]

    # grid and periodized potential, axis order (x, y[, z])
    xs = np.meshgrid(*[np.arange(n)*h for n, h in zip(shape, dx)], indexing="ij")
    nim = [int(np.ceil(ker.rmax/L)) + 1 for L in lengths]
    Up = np.zeros(shape)
    for img in itertools.product(*[range(-m, m + 1) for m in nim]):
        r2 = sum((x - np.round(x/L)*L + i*L)**2 for x, L, i in zip(xs, lengths, img))
        Up += ker.U(np.sqrt(r2))
    vol = float(np.prod(dx))
    Uh = np.real(np.fft.fftn(Up))*vol
    ks = np.meshgrid(*[2*np.pi*np.fft.fftfreq(n, h) for n, h in zip(shape, dx)], indexing="ij")
    K2 = sum(kk**2 for kk in ks)
    dcrit = args.dcrit or -Uh[h1]
    funi = 0.5*Uh.flat[0]

    def rhs(ph, D):
        phi = np.real(np.fft.ifftn(ph))
        ch = Uh*ph
        out = -K2*D*ph
        for kk in ks:
            out = out + 1j*kk*np.fft.fftn(phi*np.real(np.fft.ifftn(1j*kk*ch)))
        return out

    if sites is None:
        phi = 1.0 + 0.5*np.cos(k*xs[0])
    else:
        phi = np.full(shape, 0.05)
        width = 0.125*min(lengths)
        for s in sites:
            d2 = 0.0
            for x, L, f in zip(xs, lengths, s):
                d = x - f*L
                d2 = d2 + (d - np.round(d/L)*L)**2
            phi += np.exp(-d2/(2*width**2))
    phi /= phi.mean()
    ph = np.fft.fftn(phi)

    print("%s %dD cell %s, grid %s, |k| = %.4f (scale %.3f), mu0 Uhat(k) = %.5g, Uhat(0) = %.5g"
          % (args.lattice, dim, " x ".join("%.4f" % L for L in lengths), "x".join(map(str, shape)),
             k, args.scale, Uh[h1], Uh.flat[0]), flush=True)
    print("# D  D/Dcrit  phi_max  phi_min  A(k1)  F-F_uni  residual", flush=True)
    npts = float(np.prod(shape))
    for D in args.D:
        L = -K2*(D + Uh)
        for _ in range(int(round(args.tmax/args.dt))):
            ph = (ph + args.dt*(rhs(ph, D) - L*ph))/(1 - args.dt*L)
        phi = np.real(np.fft.ifftn(ph))
        res = np.abs(np.real(np.fft.ifftn(rhs(ph, D)))).max()
        C = np.real(np.fft.ifftn(Uh*ph))
        F = np.mean(D*phi*np.log(np.maximum(phi, 1e-300)) + 0.5*phi*C)
        A = np.abs(ph[h1])/npts
        print("%.5g %.4f %.4g %.4g %.4f %+.5f %.2e" % (D, D/dcrit, phi.max(), phi.min(), A, F - funi, res),
              flush=True)
        if phi.max() < 1.01:
            print("# decayed to the uniform state", flush=True)
            break


if __name__ == "__main__":
    main()
