#!/usr/bin/env python3
"""Compare inertial Dean-Kawasaki runs (src_deankow/inertial) with each other
and with the exact harmonic-potential solution.

Each plotfile holds phi(x, v) on the phase-space grid (the SPDE field, or the
particle histogram count/(N dx dv)). For every plotfile of each run given
(--runs, one glob per run) the script prints the time, the mass and the
moments <x - x0>, <v>, <(x - x0)^2>, <v^2>, <(x - x0) v> (cell-centre
quadrature, x - x0 by minimum image), next to the exact values.

Exact values: for U = k/2 (x - x0)^2 the means obey the damped oscillator
x' = v, v' = -gamma v - k x, and the covariance obeys the Lyapunov equation
C' = A C + C A^T + B with A = [[0, 1], [-k, -gamma]], B = diag(0, 2 D),
both from the Gaussian initial condition (ic_x0, ic_v0, ic_sx, ic_sv). The
discrete initial state is a sampled Gaussian, so the script starts the exact
solution from the moments of the first plotfile of the first run.

With two or more runs it also compares phi at matching times: the rms of the
difference divided by the Poisson noise sqrt(phi_mean/(N dx dv)), which is
about 1 (or sqrt(2) for two noisy runs) when they agree within sampling noise,
and the L1 distance of each run's last frame from the equilibrium
phi_eq ~ exp(-(v^2/2 + U)/T).

Example:
    python3 inertial_compare.py --runs 'plt_det_*' 'plt_part_*' --num-part 1e6 \\
        --gamma 1 --D 1 --k 100
"""

import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stack_plotfiles_time import read_header, read_level0  # noqa: E402


def load(f):
    h = read_header(f)
    d = read_level0(f, h)
    phi = d[h["varnames"].index("phi")]            # [nv, nx]
    lo, hi = h["prob_lo"], h["prob_hi"]
    nv, nx = phi.shape
    dx = (hi[0] - lo[0])/nx
    dv = (hi[1] - lo[1])/nv
    x = lo[0] + (np.arange(nx) + 0.5)*dx
    v = lo[1] + (np.arange(nv) + 0.5)*dv
    return h["time"], phi, x, v, dx, dv, hi[0] - lo[0]


def moments(phi, x, v, dx, dv, x0, L):
    d = x - x0
    d -= L*np.round(d/L)
    w = phi*dx*dv
    m = w.sum()
    D, V = np.meshgrid(d, v)
    return m, (w*D).sum()/m, (w*V).sum()/m, (w*D*D).sum()/m, (w*V*V).sum()/m, (w*D*V).sum()/m


def exact(t, mean0, cov0, gamma, D, k, nsub=4000):
    """Means and covariance at the times t from (mean0, cov0) at t[0], RK4."""
    A = np.array([[0.0, 1.0], [-k, -gamma]])
    B = np.array([[0.0, 0.0], [0.0, 2.0*D]])
    def f(y):
        m, C = y[:2], y[2:].reshape(2, 2)
        return np.concatenate([A @ m, (A @ C + C @ A.T + B).ravel()])
    y = np.concatenate([mean0, cov0.ravel()])
    out = [y.copy()]
    for t0, t1 in zip(t[:-1], t[1:]):
        n = max(1, int(np.ceil((t1 - t0)*nsub)))
        h = (t1 - t0)/n
        for _ in range(n):
            k1 = f(y); k2 = f(y + 0.5*h*k1); k3 = f(y + 0.5*h*k2); k4 = f(y + h*k3)
            y = y + h*(k1 + 2*k2 + 2*k3 + k4)/6
        out.append(y.copy())
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True, help="one plotfile glob per run")
    ap.add_argument("--num-part", type=float, required=True, help="N, for the Poisson noise level")
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--D", type=float, default=1.0)
    ap.add_argument("--k", type=float, default=100.0)
    ap.add_argument("--x0", type=float, help="potential centre (default: centre of the x domain)")
    args = ap.parse_args()

    runs = []
    for g in args.runs:
        files = sorted(glob.glob(g), key=lambda f: read_header(f)["time"])
        if not files:
            sys.exit("no plotfiles match " + g)
        runs.append((g, [load(f) for f in files]))

    t0, phi0, x, v, dx, dv, L = runs[0][1][0]
    x0 = args.x0 if args.x0 is not None else x[0] - 0.5*dx + 0.5*L
    T = args.D/args.gamma

    # exact solution from the moments of the first frame of the first run
    m0 = moments(phi0, x, v, dx, dv, x0, L)
    mean0 = np.array([m0[1], m0[2]])
    cov0 = np.array([[m0[3] - m0[1]**2, m0[5] - m0[1]*m0[2]], [m0[5] - m0[1]*m0[2], m0[4] - m0[2]**2]])
    times = np.array([fr[0] for fr in runs[0][1]])
    ex = exact(times, mean0, cov0, args.gamma, args.D, args.k)

    print("# exact: damped oscillator / Lyapunov from the first frame of %s; T = D/gamma = %g" % (runs[0][0], T))
    print("# equilibrium: <(x-x0)^2> = T/k = %.6g, <v^2> = T = %.6g" % (T/args.k, T))
    for name, frames in runs:
        print("\n# run %s" % name)
        print("#    t      mass       <x-x0>     <v>       <(x-x0)^2>   <v^2>     <(x-x0)v>  |  exact <x-x0> <v> <(x-x0)^2> <v^2>")
        for n, (t, phi, *_ ) in enumerate(frames):
            m = moments(phi, x, v, dx, dv, x0, L)
            row = "%7.3f %.10f %+.5f %+.5f %.6f %.5f %+.6f" % (t, m[0], m[1], m[2], m[3], m[4], m[5])
            j = np.argmin(abs(times - t))
            if abs(times[j] - t) < 1e-9:
                e = ex[j]
                row += "  |  %+.5f %+.5f %.6f %.5f" % (e[0], e[1], e[2] + e[0]**2, e[5] + e[1]**2)
            print(row)

    # equilibrium distance of each run's last frame
    X, Vv = np.meshgrid(x, v)
    d = X - x0
    d -= L*np.round(d/L)
    eq = np.exp(-(0.5*Vv*Vv + 0.5*args.k*d*d)/T)
    eq /= eq.sum()*dx*dv
    print("\n# L1 distance from phi_eq at the last frame (sum |phi - phi_eq| dx dv)")
    for name, frames in runs:
        t, phi = frames[-1][0], frames[-1][1]
        print("  %-30s t = %.3f  L1 = %.5f" % (name, t, np.abs(phi - eq).sum()*dx*dv))

    # pairwise comparison at matching times, in units of the Poisson noise
    if len(runs) > 1:
        print("\n# rms(phi_a - phi_b) / Poisson noise sqrt(phi_mean/(N dx dv)), cells with phi_mean > 1e-3 max")
        a_name, a_frames = runs[0]
        for b_name, b_frames in runs[1:]:
            bt = {round(fr[0], 6): fr[1] for fr in b_frames}
            for t, pa, *_ in a_frames:
                pb = bt.get(round(t, 6))
                if pb is None:
                    continue
                pm = 0.5*(pa + pb)
                sel = pm > 1e-3*pm.max()
                noise = np.sqrt(pm[sel]/(args.num_part*dx*dv))
                r = np.sqrt(np.mean(((pa - pb)[sel]/noise)**2))
                print("  %s vs %s  t = %.3f  rms/noise = %.3f" % (a_name, b_name, t, r))


if __name__ == "__main__":
    main()
