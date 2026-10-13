#!/usr/bin/env python3
"""Growth of the structure factor on the critical shell, from 2D or 3D plotfiles.

For each plotfile (level 0, periodic) this computes the particle-normalized
structure factor of the density,

  S(k) = |sum_cells n e^{-i k.x}|^2 / N,  n = N phi/sum(phi),

so that an uncorrelated (Poisson) field has S = 1. It reports, per frame, the
mean and maximum of S on the shell |m| in [m0 - w/2, m0 + w/2), with
k = 2 pi m on the unit box, the overall maximum of S, and phi_max/mean.

With --D and --sigma (the "max growth rate sigma" line of the run's startup
output) it compares with linear theory driven by noise, starting from S = 0,

  dS/dt = 2 sigma S + 2 k^2 D,

and reports for each --thresholds value the time the shell mean first exceeds
it, against the linear prediction. The effective growth rate is fitted to
ln S of the shell mean over --fit-range and given as a ratio to sigma. A
ratio of 1 means linear growth. Larger ratios mean the nonlinear (triad)
growth has taken over, as it does near D_crit for cluster crystals.

Example (2D run with D = 0.545, sigma = 5.61 from its startup output):
    python3 sk_growth.py _new_2d_10_vvclose/plt_fv_* --num-part 327680 \\
        --D 0.545 --sigma 5.61 -o vvclose_sk.txt
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stack_plotfiles_time import read_header, read_level0  # noqa: E402


def shell_radius(shape):
    """|m| of every FFT index, for the unit box (k = 2 pi m)."""
    ms = np.meshgrid(*[np.fft.fftfreq(n)*n for n in shape], indexing="ij")
    return np.sqrt(sum(m*m for m in ms))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfiles", nargs="+")
    ap.add_argument("--num-part", type=float, required=True)
    ap.add_argument("--var", default="phi0", help="density variable (default phi0)")
    ap.add_argument("--shell", type=float, help="shell centre |m| (default: |m| of the largest S in the last frame)")
    ap.add_argument("--width", type=float, default=1.0, help="shell width in |m| (default 1)")
    ap.add_argument("--every", type=int, default=1, help="use every n-th plotfile")
    ap.add_argument("--D", type=float, help="diffusion coefficient, for the linear prediction")
    ap.add_argument("--sigma", type=float, help="linear growth rate of the shell (code startup output)")
    ap.add_argument("--thresholds", type=float, nargs="+", default=[100, 1000, 10000])
    ap.add_argument("--fit-range", type=float, nargs=2, default=[30, 3000], metavar=("SLO", "SHI"))
    ap.add_argument("-o", "--output", help="write the per-frame table here")
    args = ap.parse_args()

    frames = []
    for f in args.plotfiles:
        h = read_header(f)
        frames.append((h["step"], f, h))
    frames.sort()
    frames = frames[::args.every]
    N = args.num_part
    M = None
    rows = []
    for step, f, h in frames:
        data = read_level0(f, h)
        phi = data[h["varnames"].index(args.var)]
        if M is None:
            M = shell_radius(phi.shape)
        n = N*phi/phi.sum()
        S = np.abs(np.fft.fftn(n))**2/N
        S.flat[0] = 0.0
        rows.append((h["time"], step, S, phi.max()/phi.mean()))
    m0 = args.shell if args.shell else M.flat[int(np.argmax(rows[-1][2]))]
    sh = (M >= m0 - args.width/2) & (M < m0 + args.width/2)
    table = []
    for t, step, S, pm in rows:
        imax = int(np.argmax(S))
        table.append((t, step, S[sh].mean(), S[sh].max(), S.flat[imax], M.flat[imax], pm))
    table = np.array(table)

    lines = ["# shell |m| in [%.3f, %.3f): %d modes; k = 2 pi m" % (m0 - args.width/2, m0 + args.width/2, sh.sum()),
             "# t step S_shell_mean S_shell_max S_max m_at_max phimax/mean"]
    lines += ["%.6g %d %.6g %.6g %.6g %.4f %.4f" % tuple(r) for r in table]
    if args.output:
        with open(args.output, "w") as fo:
            fo.write("\n".join(lines) + "\n")
    else:
        print("\n".join(lines))

    t, smean = table[:, 0], table[:, 2]
    print("\nsummary: shell |m| = %.3f, %d frames, t = %.4g to %.4g" % (m0, len(t), t[0], t[-1]))
    k2 = (2*np.pi*m0)**2
    for thr in args.thresholds:
        hit = np.nonzero(smean > thr)[0]
        tobs = "%.4g" % t[hit[0]] if hit.size else "never"
        tlin = ""
        if args.D is not None and args.sigma:
            if args.sigma > 0:
                tlin = "  linear %.4g" % (np.log(1 + thr*args.sigma/(k2*args.D))/(2*args.sigma))
            else:
                seq = k2*args.D/(-args.sigma)
                tlin = "  linear: never (stable; S saturates at %.3g)" % seq if seq < thr else ""
        print("  shell mean S > %g: t = %s%s" % (thr, tobs, tlin))
    fit = (smean > args.fit_range[0]) & (smean < args.fit_range[1])
    if fit.sum() >= 2:
        rate = np.polyfit(t[fit], np.log(smean[fit]), 1)[0]/2
        ratio = "  (%.2f x sigma)" % (rate/args.sigma) if args.sigma else ""
        print("  fitted growth rate over S in [%g, %g]: %.4g%s" % (args.fit_range[0], args.fit_range[1], rate, ratio))
    print("  final: shell mean S = %.4g, phi_max/mean = %.3g" % (smean[-1], table[-1, 6]))


if __name__ == "__main__":
    main()
