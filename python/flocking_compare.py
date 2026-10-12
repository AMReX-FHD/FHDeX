#!/usr/bin/env python3
"""Compare a flocking SPDE plotfile (3D: x, y, theta) with a particle plotfile (2D).

The SPDE plotfile holds phi(x, y, theta). The particle plotfile holds one
component per theta bin, phi_t000 ... phi_tNNN, each an estimate of
phi(x, y, theta_k) on the same x, y grid, plus the moments rho, px, py.
Both are normalized to integrate to 1 over (x, y, theta).

Prints relative L2 differences ||a - b|| / ||b|| (b = SPDE) of the fields
below; px and py are both divided by the norm of the SPDE polarization vector
(or of rho when that is zero, as for a theta-symmetric state):
  phi        the full phase-space density
  rho        spatial density   rho(x, y)   = sum_k phi_k dtheta
  px, py     polarization      p(x, y)     = sum_k phi_k (cos, sin)(theta_k) dtheta
  marginal   angular density   psi(theta)  = sum_xy phi dx dy
and the projection of each onto the initial mode
  A = int phi cos(2 pi (kx x/Lx + ky y/Ly)) cos(m (theta - theta0)) dV.

With --variance it also reports, for each plotfile, the mean over cells of
(phi - mean)^2 against the Poisson value mean/(N dV), for a uniform state.

Example:
    python3 flocking_compare.py plt_det_000115 plt_part_000035 --num-part 4e6
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stack_plotfiles_time import read_header, read_level0  # noqa: E402


def load_spde(p):
    h = read_header(p)
    if h["dim"] != 3:
        sys.exit("%s: expected a 3D SPDE plotfile" % p)
    phi = read_level0(p, h)[h["varnames"].index("phi")]   # [ntheta, ny, nx]
    return h, phi


def load_particles(p):
    h = read_header(p)
    if h["dim"] != 2:
        sys.exit("%s: expected a 2D particle plotfile" % p)
    data = read_level0(p, h)
    names = [v for v in h["varnames"] if v.startswith("phi_t")]
    phi = np.stack([data[h["varnames"].index(v)] for v in sorted(names)])  # [ntheta, ny, nx]
    return h, phi


def load_field(p, h, name, dim):
    """A 2D field [ny, nx] from a plotfile, or None if it is not there. In 3D
    the field is the same in every theta plane, so plane 0 is returned."""
    if name not in h["varnames"]:
        return None
    d = read_level0(p, h)[h["varnames"].index(name)]
    return d[0] if dim == 3 else d


def moments(phi, dth):
    nth = phi.shape[0]
    th = (np.arange(nth) + 0.5) * dth
    rho = phi.sum(axis=0) * dth
    px = np.tensordot(np.cos(th), phi, axes=(0, 0)) * dth
    py = np.tensordot(np.sin(th), phi, axes=(0, 0)) * dth
    return th, rho, px, py


def rel_l2(a, b, norm=None):
    """||a - b|| / ||b||, or / norm when given."""
    nb = np.sqrt((b * b).sum()) if norm is None else norm
    return np.sqrt(((a - b) ** 2).sum()) / nb if nb > 0 else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spde", help="SPDE plotfile (3D)")
    ap.add_argument("particles", help="particle plotfile (2D)")
    ap.add_argument("--kx", type=int, default=1, help="mode kx for the projection (default 1)")
    ap.add_argument("--ky", type=int, default=0, help="mode ky (default 0)")
    ap.add_argument("--m", type=int, default=1, help="mode m in theta (default 1)")
    ap.add_argument("--theta0", type=float, default=0.0, help="mode phase theta0 (default 0)")
    ap.add_argument("--num-part", type=float, help="N, for --variance")
    ap.add_argument("--variance", action="store_true",
                    help="report the cell variance of phi against the Poisson value")
    args = ap.parse_args()

    hs, ps = load_spde(args.spde)
    hp, pp = load_particles(args.particles)
    if ps.shape != pp.shape:
        sys.exit("grids differ: SPDE %s, particles %s (ntheta, ny, nx)" % (ps.shape, pp.shape))
    if abs(hs["time"] - hp["time"]) > 1e-9 * max(1.0, abs(hs["time"])):
        print("warning: times differ: SPDE %.10g, particles %.10g" % (hs["time"], hp["time"]))

    dx, dy, dth = hs["dx"]
    lx = hs["prob_hi"][0] - hs["prob_lo"][0]
    ly = hs["prob_hi"][1] - hs["prob_lo"][1]
    nth, ny, nx = ps.shape
    x = hs["prob_lo"][0] + (np.arange(nx) + 0.5) * dx
    y = hs["prob_lo"][1] + (np.arange(ny) + 0.5) * dy

    th, rs, pxs, pys = moments(ps, dth)
    _, rp, pxp, pyp = moments(pp, dth)
    ms = ps.sum(axis=(1, 2)) * dx * dy
    mp = pp.sum(axis=(1, 2)) * dx * dy

    mode = (np.cos(2 * np.pi * (args.kx * (x[None, :] - hs["prob_lo"][0]) / lx
                                + args.ky * (y[:, None] - hs["prob_lo"][1]) / ly))[None]
            * np.cos(args.m * (th - args.theta0))[:, None, None])
    dv = dx * dy * dth

    print("t = %.6g (SPDE)  %.6g (particles)   grid %d x %d x %d" % (hs["time"], hp["time"], nx, ny, nth))
    print("mass                 SPDE %.10f   particles %.10f" % (ps.sum() * dv, pp.sum() * dv))
    print("relative L2 difference (particles vs SPDE):")
    print("  phi(x,y,theta)  %.4e" % rel_l2(pp, ps))
    print("  rho(x,y)        %.4e" % rel_l2(rp, rs))
    pnorm = np.sqrt((pxs * pxs + pys * pys).sum())
    if pnorm < 1e-12 * np.sqrt((rs * rs).sum()):
        pnorm = np.sqrt((rs * rs).sum())
    print("  px(x,y)         %.4e" % rel_l2(pxp, pxs, pnorm))
    print("  py(x,y)         %.4e" % rel_l2(pyp, pys, pnorm))
    print("  psi(theta)      %.4e" % rel_l2(mp, ms))
    rts, rtp = load_field(args.spde, hs, "rhot", 3), load_field(args.particles, hp, "rhot", 2)
    if rts is not None and rtp is not None:
        print("  rhot(x,y)       %.4e   (smoothed density W * rho)" % rel_l2(rtp, rts))
    a_s = (ps * mode).sum() * dv
    a_p = (pp * mode).sum() * dv
    print("mode projection A    SPDE %.6e   particles %.6e   difference %.3e"
          % (a_s, a_p, a_p - a_s))
    if args.num_part:
        # particle sampling standard deviation of A = (1/N) sum_i mode(x_i, theta_i)
        sig = np.sqrt(max((ps * mode * mode).sum() * dv - a_s * a_s, 0.0) / args.num_part)
        print("  particle noise std of A %.3e   difference / std %.2f" % (sig, (a_p - a_s) / sig))

    if args.variance:
        if not args.num_part:
            sys.exit("--variance needs --num-part")
        for name, f in (("SPDE", ps), ("particles", pp)):
            mean = f.mean()
            var = ((f - mean) ** 2).mean()
            print("variance %-9s  %.4e   Poisson mean/(N dV) %.4e   ratio %.4f"
                  % (name, var, mean / (args.num_part * dv), var / (mean / (args.num_part * dv))))


if __name__ == "__main__":
    main()
