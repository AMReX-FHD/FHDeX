#!/usr/bin/env python3
"""Render 2D AMReX plotfiles as PNG images of log10(phi), without matplotlib.

Each plotfile becomes <outdir>/frame_<step>.png: level-0 data of one variable,
log10 of max(phi, floor), on a fixed colour scale shared by all frames, with
y up. A colour bar is written to <outdir>/colorbar.png. Only numpy and the
standard library (zlib) are used.

Example:
    python3 plot_frames.py plt_md_* -o plots --every 10
    python3 plot_frames.py plt_md_* -o plots --vmin -1 --vmax 2.5 --scale 2

Make a montage or a movie from the frames with ImageMagick or ffmpeg, e.g.
    montage plots/frame_*.png -tile 6x -geometry +2+2 montage.png
    ffmpeg -framerate 10 -pattern_type glob -i 'plots/frame_*.png' movie.mp4
"""

import argparse
import os
import struct
import sys
import zlib

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stack_plotfiles_time import read_header, read_level0  # noqa: E402

# control points of a perceptually ordered map (dark blue -> teal -> yellow)
_CMAP = np.array([[0.050, 0.030, 0.330], [0.230, 0.320, 0.550], [0.130, 0.570, 0.550],
                  [0.370, 0.780, 0.380], [0.990, 0.910, 0.140]])


def colourize(v):
    """v in [0, 1] -> uint8 RGB by piecewise-linear interpolation of _CMAP."""
    x = np.clip(v, 0.0, 1.0)*(len(_CMAP) - 1)
    i = np.minimum(x.astype(int), len(_CMAP) - 2)
    f = (x - i)[..., None]
    rgb = _CMAP[i]*(1 - f) + _CMAP[i + 1]*f
    return (255*rgb + 0.5).astype(np.uint8)


def write_png(path, rgb):
    """rgb: (ny, nx, 3) uint8, first row at the top."""
    ny, nx, _ = rgb.shape
    raw = b"".join(b"\x00" + rgb[j].tobytes() for j in range(ny))

    def chunk(tag, data):
        c = tag + data
        return struct.pack(">I", len(data)) + c + struct.pack(">I", zlib.crc32(c) & 0xffffffff)

    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n")
        f.write(chunk(b"IHDR", struct.pack(">IIBBBBB", nx, ny, 8, 2, 0, 0, 0)))
        f.write(chunk(b"IDAT", zlib.compress(raw, 6)))
        f.write(chunk(b"IEND", b""))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfiles", nargs="+", help="2D plotfile directories")
    ap.add_argument("-o", "--outdir", default="plots")
    ap.add_argument("--var", default="phi0", help="variable (default phi0)")
    ap.add_argument("--every", type=int, default=1, help="use every n-th plotfile")
    ap.add_argument("--vmin", type=float, help="log10 lower limit (default: log10 of the mean density minus 1)")
    ap.add_argument("--vmax", type=float, help="log10 upper limit (default: max over the frames)")
    ap.add_argument("--scale", type=int, default=1, help="integer pixel magnification")
    args = ap.parse_args()

    frames = sorted((read_header(f)["step"], f) for f in args.plotfiles)[::args.every]
    os.makedirs(args.outdir, exist_ok=True)
    data = []
    for step, f in frames:
        h = read_header(f)
        if h["dim"] != 2:
            sys.exit("%s is not 2D" % f)
        phi = read_level0(f, h)[h["varnames"].index(args.var)]   # shape (ny, nx)
        data.append((step, h["time"], phi))
    mean = np.mean(data[0][2])
    vmin = args.vmin if args.vmin is not None else np.log10(mean) - 1.0
    vmax = args.vmax if args.vmax is not None else max(np.log10(max(d[2].max(), 1e-30)) for d in data)
    floor = 10.0**vmin
    with open(os.path.join(args.outdir, "frames.txt"), "w") as out:
        out.write("# step time max_phi file; colour scale log10(%s) in [%g, %g]\n" % (args.var, vmin, vmax))
        for step, t, phi in data:
            v = (np.log10(np.maximum(phi, floor)) - vmin)/(vmax - vmin)
            img = colourize(v[::-1])                                 # y up
            if args.scale > 1:
                img = img.repeat(args.scale, axis=0).repeat(args.scale, axis=1)
            name = "frame_%07d.png" % step
            write_png(os.path.join(args.outdir, name), img)
            out.write("%d %.8g %.6g %s\n" % (step, t, phi.max(), name))
    bar = colourize(np.linspace(1, 0, 256)[:, None].repeat(24, axis=1))
    write_png(os.path.join(args.outdir, "colorbar.png"), bar)
    print("wrote %d frames to %s; log10(%s) in [%g, %g] (colorbar.png: top = vmax)"
          % (len(data), args.outdir, args.var, vmin, vmax))


if __name__ == "__main__":
    main()
