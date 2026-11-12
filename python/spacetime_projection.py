#!/usr/bin/env python3
"""Project space-time stacks onto a (space, time) plane and write one PNG.

Input: 3D plotfiles made by stack_plotfiles_time.py, with time along z.
For each stack, a variable is reduced over one space direction (the largest value
over y, by default), so every cluster becomes a track in x-t and a merger is two
tracks joining. The panels sit side by side, with time increasing upward, on one
log10 colour scale, so brightness means the same density in every panel. Only
numpy and the standard library are used.

A text file next to the image (<output> with .txt) gives each panel's pixel
range and its time and space extent, since the PNG has no axes.

Example:
    python3 spacetime_projection.py _ppc1_collapse_n50/st_cn _ppc1_driftmerge_n50/st_dm \\
        -o xt_collapse_vs_driftmerge.png --var measure
    python3 spacetime_projection.py st_cn st_dm -o yt.png --show y --tmin 0.03 --tmax 0.06
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_frames import colourize, write_png  # noqa: E402
from stack_plotfiles_time import read_header, read_level0  # noqa: E402


def slice_times(stack, nz):
    """Time of each z slice from <stack>/times.txt: t_target with --dt, else the plotfile time."""
    with open(os.path.join(stack, "times.txt")) as f:
        head = f.readline().lstrip("#").split()
        rows = np.loadtxt(f, ndmin=2)
    col = head.index("t_target") if "t_target" in head else head.index("time")
    t = rows[:, col]
    if len(t) != nz:
        sys.exit("%s: times.txt has %d rows for %d slices" % (stack, len(t), nz))
    return t


def project(stack, var, show, tmin, tmax):
    h = read_header(stack)
    if h["dim"] != 3:
        sys.exit("%s: spacedim is %d, expected a 3D space-time stack" % (stack, h["dim"]))
    if var is None:
        var = h["varnames"][0]
    if var not in h["varnames"]:
        sys.exit("%s: unknown variable %s; available: %s" % (stack, var, " ".join(h["varnames"])))
    data = read_level0(stack, h)[h["varnames"].index(var)]      # [nt, ny, nx]
    t = slice_times(stack, data.shape[0])
    keep = np.ones(len(t), bool)
    if tmin is not None:
        keep &= t >= tmin
    if tmax is not None:
        keep &= t <= tmax
    if keep.sum() < 2:
        sys.exit("%s: fewer than two slices in the time window" % stack)
    d = 0 if show == "x" else 1
    img = data[keep].max(axis=1 if show == "x" else 2)           # [nt, n_show]
    return img, t[keep], (h["prob_lo"][d], h["prob_hi"][d]), var


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stacks", nargs="+", help="3D space-time plotfiles from stack_plotfiles_time.py")
    ap.add_argument("-o", "--output", default="spacetime.png", help="PNG file")
    ap.add_argument("--var", help="variable (default: the first in each stack)")
    ap.add_argument("--show", choices=("x", "y"), default="x",
                    help="space direction kept; the other is reduced by its maximum (default x)")
    ap.add_argument("--tmin", type=float, help="first time shown")
    ap.add_argument("--tmax", type=float, help="last time shown")
    ap.add_argument("--floor", type=float, default=1.0,
                    help="values below this are drawn darkest (default 1, the mean density)")
    ap.add_argument("--vmax", type=float, help="log10 upper limit (default: largest 99.5th percentile)")
    ap.add_argument("--gap", type=int, default=12, help="pixels between panels")
    ap.add_argument("--scale", type=int, default=1, help="integer pixel magnification")
    args = ap.parse_args()
    if args.floor <= 0.0 or args.gap < 0 or args.scale < 1:
        ap.error("--floor must be > 0, --gap >= 0 and --scale >= 1")

    panels = [project(s, args.var, args.show, args.tmin, args.tmax) for s in args.stacks]
    nt = min(p[0].shape[0] for p in panels)
    if any(p[0].shape[0] != nt for p in panels):
        print("warning: panels have different slice counts; all are cut to the first %d" % nt,
              file=sys.stderr)
    vmin = np.log10(args.floor)
    logs = [np.log10(np.maximum(p[0][:nt], args.floor)) for p in panels]
    vmax = args.vmax if args.vmax is not None else max(np.percentile(v, 99.5) for v in logs)
    if vmax <= vmin:
        sys.exit("vmax %.3g is not above log10(floor) %.3g" % (vmax, vmin))

    parts, cols, x0 = [], [], 0
    for k, v in enumerate(logs):
        if k > 0:
            parts.append(np.ones((nt, args.gap)))
            x0 += args.gap
        parts.append((v - vmin) / (vmax - vmin))
        cols.append((x0, x0 + v.shape[1] - 1))
        x0 += v.shape[1]
    img = np.hstack(parts)[::-1]                                # time increasing upward
    rgb = colourize(img)
    if args.scale > 1:
        rgb = np.repeat(np.repeat(rgb, args.scale, axis=0), args.scale, axis=1)
    write_png(args.output, rgb)

    txt = os.path.splitext(args.output)[0] + ".txt"
    with open(txt, "w") as f:
        f.write("# %s: max over %s of %s, log10 colour scale [%.4g, %.4g], time upward, "
                "%d slices, pixel scale %d\n"
                % (args.output, "y" if args.show == "x" else "x", panels[0][3], vmin, vmax, nt, args.scale))
        f.write("# panel stack columns(first last) %s_lo %s_hi t_first t_last\n" % (args.show, args.show))
        for k, (s, p, c) in enumerate(zip(args.stacks, panels, cols)):
            t = p[1][:nt]
            f.write("%d %s %d %d %.6g %.6g %.6g %.6g\n" % (k, s, c[0], c[1], p[2][0], p[2][1], t[0], t[-1]))
    print("wrote %s (%d x %d) and %s; log10 scale [%.3g, %.3g]"
          % (args.output, rgb.shape[1], rgb.shape[0], txt, vmin, vmax))


if __name__ == "__main__":
    main()
