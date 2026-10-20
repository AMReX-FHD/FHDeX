#!/usr/bin/env python3
"""Bragg spots of the full structure factor of 2D or 3D plotfiles.

Averages |FFT(phi/mean - 1)|^2 over the given plotfiles (level 0, periodic,
unit box so k = 2 pi m) and reports:

  * the shell profile: mean and max S in |m| bins of width --bin;
  * statistics of the peak shell (|m| within --shell-width of the bin with
    the largest mean S, or --shell): number of modes n, max/mean, the
    share of the shell's S in the brightest 12, 24 and 48 modes, and the
    participation ratio PR = (sum S)^2/(n sum S^2);
  * the brightest --top +- pairs on the peak shell, each with its integer
    mode, S relative to the brightest, and the angles (folded to <= 90
    degrees) to the brightest six.

How to read them. An isotropic ring (liquid, amorphous) has exponentially
distributed S on the shell: max/mean about ln n and PR about 0.5. A single
crystal concentrates the shell on its spots: PR about (number of spots)/n
and max/mean about n/(number of spots). The angles identify the lattice:

  2D hexagonal  6 spots, 60 degrees apart
  3D bcc        12 <110> spots, 60 and 90 degrees; next shells at sqrt2, sqrt3, 2
  3D fcc        8 <111> spots, 70.5 degrees; next shell <200> at 1.155 x

Several sets of spots at fixed relative rotations mean grains or twins.

Example:
    python3 sf_spots.py near_3d/plt_fv3d_11*  --every 10
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stack_plotfiles_time import read_header, read_level0  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfiles", nargs="+")
    ap.add_argument("--var", default="phi0", help="density variable (default phi0)")
    ap.add_argument("--every", type=int, default=1, help="use every n-th plotfile")
    ap.add_argument("--last", type=int, help="use only the last n plotfiles (after --every)")
    ap.add_argument("--bin", type=float, default=0.5, help="|m| bin width of the shell profile")
    ap.add_argument("--mmax", type=float, default=20.0, help="largest |m| in the shell profile")
    ap.add_argument("--shell", type=float, help="peak shell centre |m| (default: bin with the largest mean S)")
    ap.add_argument("--shell-width", type=float, default=1.0)
    ap.add_argument("--top", type=int, default=14, help="number of +- spot pairs to list")
    ap.add_argument("--save", help="save the averaged S as a .npy file")
    args = ap.parse_args()

    frames = sorted((read_header(f)["step"], f) for f in args.plotfiles)[::args.every]
    if args.last:
        frames = frames[-args.last:]
    S = None
    for _, f in frames:
        h = read_header(f)
        phi = read_level0(f, h)[h["varnames"].index(args.var)]
        d = phi/phi.mean() - 1.0
        s = np.abs(np.fft.fftn(d))**2
        S = s if S is None else S + s
    S /= len(frames)
    S.flat[0] = 0.0
    if args.save:
        np.save(args.save, S)
    ms = np.meshgrid(*[np.fft.fftfreq(n)*n for n in S.shape], indexing="ij")
    M = np.sqrt(sum(m*m for m in ms))
    t0, t1 = read_header(frames[0][1])["time"], read_header(frames[-1][1])["time"]
    print("%d frames, t = %.5g to %.5g, grid %s" % (len(frames), t0, t1, "x".join(map(str, S.shape))))

    print("\nshell profile: |m| range, modes, mean S, max S")
    best, bestm = -1.0, None
    for lo in np.arange(args.bin, args.mmax, args.bin):
        sh = (M >= lo) & (M < lo + args.bin)
        if not sh.any():
            continue
        mean = S[sh].mean()
        print("  %6.2f-%6.2f %6d %11.4g %11.4g" % (lo, lo + args.bin, sh.sum(), mean, S[sh].max()))
        if mean > best:
            best, bestm = mean, lo + args.bin/2

    m0 = args.shell or bestm
    sh = (M >= m0 - args.shell_width/2) & (M < m0 + args.shell_width/2)
    v = S[sh]
    srt = np.sort(v)[::-1]
    n = v.size
    print("\npeak shell |m| in [%.3f, %.3f): %d modes" % (m0 - args.shell_width/2, m0 + args.shell_width/2, n))
    print("  max/mean %.1f   (ring ~ ln n = %.1f)" % (v.max()/v.mean(), np.log(n)))
    print("  share of shell S in the brightest 12 / 24 / 48 modes: %.3f %.3f %.3f"
          % tuple(srt[:j].sum()/v.sum() for j in (12, 24, 48)))
    print("  participation ratio %.4f   (ring ~ 0.5; 6 spots %.4f, 12 spots %.4f)"
          % (v.sum()**2/(n*(v*v).sum()), 6/n, 12/n))

    idx = np.nonzero(sh.ravel())[0]
    idx = idx[np.argsort(S.ravel()[idx])[::-1]]
    reps = []
    for i in idx:
        m = np.array([mm.flat[i] for mm in reversed(ms)])   # arrays are (z,) y, x
        if any(np.array_equal(m, -r) for r, _ in reps):
            continue
        reps.append((m, S.flat[i]))
        if len(reps) == args.top:
            break
    units = [r/np.linalg.norm(r) for r, _ in reps]
    print("\nbrightest +- pairs on the peak shell (m in axis order x, y[, z]):")
    print("  %-16s %6s %7s   angles (deg, folded to <= 90) to the brightest %d"
          % ("mode", "|m|", "S/max", min(6, len(reps))))
    for (m, s), u in zip(reps, units):
        ang = [np.degrees(np.arccos(min(1.0, abs(u @ w)))) for w in units[:6]]
        print("  %-16s %6.2f %7.3f   %s" % ("(" + ",".join("%d" % c for c in m) + ")", np.linalg.norm(m),
                                         s/reps[0][1], " ".join("%4.0f" % a for a in ang)))


if __name__ == "__main__":
    main()
