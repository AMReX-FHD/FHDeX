#!/usr/bin/env python3
"""Stack a sequence of 2D AMReX plotfiles into one 3D plotfile with time as z.

Slice k of the output is the level-0 data of the k-th input plotfile (sorted by
step). The result is a standard single-level 3D plotfile that amrvis, AMReXplorer
and the AMReX Tools/Plotfile utilities can read. Only numpy is required.

A plotfile needs uniform cell spacing, so the z cell size is the mean interval
between outputs, optionally scaled by --tscale. The actual step and time of each
slice are written to <output>/times.txt.

Example:
    python3 stack_plotfiles_time.py plt_hk_* -o spacetime_plt --tscale 10
"""

import argparse
import os
import re
import sys

import numpy as np

BOX_RE = re.compile(r"\(\(([-\d,]+)\)\s*\(([-\d,]+)\)\s*\(([-\d,]+)\)\)")
FAB_RE = re.compile(r"FAB \(\((\d+), \(([\d ]+)\)\),\((\d+), \(([\d ]+)\)\)\)(.*)\s+(\d+)\s*$")


def parse_box(s):
    m = BOX_RE.search(s)
    if m is None:
        raise ValueError("cannot parse box: " + s)
    lo = [int(v) for v in m.group(1).split(",")]
    hi = [int(v) for v in m.group(2).split(",")]
    return lo, hi


def read_header(pltfile):
    """Parse the plotfile Header; returns a dict with level-0 information."""
    with open(os.path.join(pltfile, "Header")) as f:
        lines = [ln.rstrip("\n") for ln in f]
    it = iter(lines)
    h = {"version": next(it)}
    nvars = int(next(it))
    h["varnames"] = [next(it).strip() for _ in range(nvars)]
    h["dim"] = int(next(it))
    h["time"] = float(next(it))
    h["finest_level"] = int(next(it))
    h["prob_lo"] = [float(v) for v in next(it).split()]
    h["prob_hi"] = [float(v) for v in next(it).split()]
    next(it)  # refinement ratios (empty for a single level)
    h["domain"] = parse_box(next(it))
    h["step"] = int(next(it).split()[0])
    h["dx"] = [float(v) for v in next(it).split()]
    for _ in range(h["finest_level"]):
        next(it)  # dx of finer levels
    h["coord_sys"] = int(next(it))
    if h["version"] != "HyperCLaw-V1.1":
        raise ValueError(pltfile + ": unsupported plotfile version " + h["version"])
    return h


def read_cell_h(level_dir):
    """Parse a VisMF Cell_H; returns the list of boxes and (file, offset) per FAB."""
    with open(os.path.join(level_dir, "Cell_H")) as f:
        lines = [ln.strip() for ln in f]
    i = 4  # version, how, ncomp, ngrow
    nbox = int(lines[i].lstrip("(").split()[0])
    boxes = [parse_box(lines[i + 1 + b]) for b in range(nbox)]
    i += nbox + 2  # boxes and ")"
    nfab = int(lines[i])
    fabs = []
    for ln in lines[i + 1:i + 1 + nfab]:
        _, fname, off = ln.split()
        fabs.append((fname, int(off)))
    return boxes, fabs


def read_fab(fname, offset):
    """Read one FAB; returns (lo, hi, data[ncomp, ..., x]) in C order."""
    with open(fname, "rb") as f:
        f.seek(offset)
        header = f.readline().decode("ascii")
        m = FAB_RE.match(header)
        if m is None:
            raise ValueError("cannot parse FAB header in %s: %s" % (fname, header))
        nbytes = int(m.group(1))
        order = m.group(4).split()
        if nbytes not in (4, 8):
            raise ValueError("unsupported real size %d in %s" % (nbytes, fname))
        if order == [str(k) for k in range(nbytes, 0, -1)]:
            endian = "<"
        elif order == [str(k) for k in range(1, nbytes + 1)]:
            endian = ">"
        else:
            raise ValueError("unsupported byte order %s in %s" % (order, fname))
        lo, hi = parse_box(m.group(5))
        ncomp = int(m.group(6))
        shape = [h - l + 1 for l, h in zip(lo, hi)]
        npts = int(np.prod(shape))
        dtype = np.dtype(endian + ("f8" if nbytes == 8 else "f4"))
        data = np.frombuffer(f.read(ncomp * npts * nbytes), dtype=dtype)
    if data.size != ncomp * npts:
        raise ValueError("short read in %s at offset %d" % (fname, offset))
    # Fortran order with x fastest, one component after another
    data = data.reshape([ncomp] + shape[::-1]).astype(np.float64)
    return lo, hi, data


def read_level0(pltfile, h):
    """Read level 0 into a full-domain array [ncomp, ny, nx] (2D) or [ncomp, nz, ny, nx] (3D)."""
    level_dir = os.path.join(pltfile, "Level_0")
    _, fabs = read_cell_h(level_dir)
    dlo, dhi = h["domain"]
    shape = [b - a + 1 for a, b in zip(dlo, dhi)]
    out = np.full([len(h["varnames"])] + shape[::-1], np.nan)
    for fname, off in fabs:
        lo, hi, data = read_fab(os.path.join(level_dir, fname), off)
        sl = tuple(slice(lo[d] - dlo[d], hi[d] - dlo[d] + 1) for d in reversed(range(len(lo))))
        out[(slice(None),) + sl] = data
    if np.isnan(out).any():
        raise ValueError(pltfile + ": level 0 does not cover the domain")
    return out


def box_str(lo, hi):
    return "((%s) (%s) (%s))" % (",".join(map(str, lo)), ",".join(map(str, hi)),
                                 ",".join("0" * len(lo)))


def write_fab(f, lo, hi, data):
    """Write one FAB (little-endian float64) at the current position of f."""
    f.write(("FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (8 7 6 5 4 3 2 1)))%s %d\n"
             % (box_str(lo, hi), data.shape[0])).encode("ascii"))
    f.write(np.ascontiguousarray(data, dtype="<f8").tobytes())


def write_plotfile(out, data, varnames, prob_lo, dx, time=0.0, step=0):
    """Write a single-level, single-box plotfile. data[ncomp, (z,) y, x] is in
    C order (as read_level0 returns it); prob_lo and dx are in x, y, (z) order.
    Existing files in out are overwritten."""
    nc = data.shape[0]
    n = list(data.shape[1:][::-1])
    dim = len(n)
    lo, hi = [0]*dim, [m - 1 for m in n]
    prob_hi = [prob_lo[d] + n[d]*dx[d] for d in range(dim)]
    level_dir = os.path.join(out, "Level_0")
    os.makedirs(level_dir, exist_ok=True)
    with open(os.path.join(level_dir, "Cell_D_00000"), "wb") as f:
        write_fab(f, lo, hi, data)
    flat = data.reshape(nc, -1)
    with open(os.path.join(level_dir, "Cell_H"), "w") as f:
        f.write("1\n1\n%d\n0\n(1 0\n%s\n)\n1\nFabOnDisk: Cell_D_00000 0\n\n" % (nc, box_str(lo, hi)))
        for vals in (flat.min(axis=1), flat.max(axis=1)):
            f.write("1,%d\n%s\n\n" % (nc, "".join("%.17e," % v for v in vals)))
    with open(os.path.join(out, "Header"), "w") as f:
        f.write("HyperCLaw-V1.1\n%d\n" % nc)
        for v in varnames:
            f.write(v + "\n")
        f.write("%d\n%.17g\n0\n" % (dim, time))
        f.write(" ".join("%.17g" % v for v in prob_lo) + "\n")
        f.write(" ".join("%.17g" % v for v in prob_hi) + "\n\n")
        f.write(box_str(lo, hi) + "\n%d\n" % step)
        f.write(" ".join("%.17g" % v for v in dx) + "\n0\n0\n")
        f.write("0 1 %.17g\n%d\n" % (time, step))
        for d in range(dim):
            f.write("%.17g %.17g\n" % (prob_lo[d], prob_hi[d]))
        f.write("Level_0/Cell\n")


def stack(files, out, varnames=None, every=1, tscale=1.0, use_index=False, kchunk=32):
    headers = [(read_header(p), p) for p in files]
    headers.sort(key=lambda hp: hp[0]["step"])
    headers = headers[::every]
    if len(headers) < 2:
        sys.exit("need at least two plotfiles")
    h0 = headers[0][0]

    for h, p in headers:
        if h["dim"] != 2:
            sys.exit("%s: spacedim is %d, expected 2" % (p, h["dim"]))
        for key in ("varnames", "domain", "dx", "prob_lo", "prob_hi", "coord_sys"):
            if h[key] != h0[key]:
                sys.exit("%s: %s differs from %s" % (p, key, headers[0][1]))
        if h["finest_level"] > 0:
            print("warning: %s has %d levels; only level 0 is used"
                  % (p, h["finest_level"] + 1), file=sys.stderr)

    allvars = h0["varnames"]
    varnames = varnames or allvars
    for v in varnames:
        if v not in allvars:
            sys.exit("unknown variable %s; available: %s" % (v, " ".join(allvars)))
    comps = [allvars.index(v) for v in varnames]

    nt = len(headers)
    steps = [h["step"] for h, _ in headers]
    times = np.array([h["time"] for h, _ in headers])
    if use_index:
        dz, zlo = 1.0, 0.0
    else:
        dz = (times[-1] - times[0]) / (nt - 1)
        if dz <= 0.0:
            sys.exit("plotfile times are not increasing; use --index")
        dev = np.max(np.abs(times - (times[0] + dz * np.arange(nt))))
        if dev > 0.01 * dz:
            print("warning: output times deviate from uniform spacing by up to %.3g of the "
                  "mean interval; see times.txt" % (dev / dz), file=sys.stderr)
        dz *= tscale
        zlo = times[0] * tscale - 0.5 * dz
    zhi = zlo + nt * dz

    dlo, dhi = h0["domain"]
    nx, ny = dhi[0] - dlo[0] + 1, dhi[1] - dlo[1] + 1
    nc = len(comps)
    level_dir = os.path.join(out, "Level_0")
    os.makedirs(level_dir, exist_ok=False)

    boxes, fabs, mins, maxs = [], [], [], []
    for c, k0 in enumerate(range(0, nt, kchunk)):
        k1 = min(k0 + kchunk, nt) - 1
        chunk = np.empty((nc, k1 - k0 + 1, ny, nx))
        for k in range(k0, k1 + 1):
            h, p = headers[k]
            chunk[:, k - k0] = read_level0(p, h)[comps]
        lo = [dlo[0], dlo[1], k0]
        hi = [dhi[0], dhi[1], k1]
        fname = "Cell_D_%05d" % c
        with open(os.path.join(level_dir, fname), "wb") as f:
            write_fab(f, lo, hi, chunk)
        boxes.append((lo, hi))
        fabs.append(fname)
        mins.append(chunk.reshape(nc, -1).min(axis=1))
        maxs.append(chunk.reshape(nc, -1).max(axis=1))
        print("wrote slices %d-%d of %d" % (k0, k1, nt))

    def minmax_block(vals):
        return "%d,%d\n" % (len(vals), nc) + "".join(
            "".join("%.17e," % v for v in row) + "\n" for row in vals)

    with open(os.path.join(level_dir, "Cell_H"), "w") as f:
        f.write("1\n1\n%d\n0\n" % nc)
        f.write("(%d 0\n" % len(boxes))
        for lo, hi in boxes:
            f.write(box_str(lo, hi) + "\n")
        f.write(")\n%d\n" % len(fabs))
        for fname in fabs:
            f.write("FabOnDisk: %s 0\n" % fname)
        f.write("\n" + minmax_block(mins) + "\n" + minmax_block(maxs))

    plo = h0["prob_lo"][:2] + [zlo]
    phi = h0["prob_hi"][:2] + [zhi]
    dx = h0["dx"][:2] + [dz]
    tlast = times[-1]
    with open(os.path.join(out, "Header"), "w") as f:
        f.write("HyperCLaw-V1.1\n%d\n" % nc)
        for v in varnames:
            f.write(v + "\n")
        f.write("3\n%.17g\n0\n" % tlast)
        f.write(" ".join("%.17g" % v for v in plo) + "\n")
        f.write(" ".join("%.17g" % v for v in phi) + "\n")
        f.write("\n")
        f.write(box_str([dlo[0], dlo[1], 0], [dhi[0], dhi[1], nt - 1]) + "\n")
        f.write("%d\n" % steps[-1])
        f.write(" ".join("%.17g" % v for v in dx) + "\n")
        f.write("%d\n0\n" % h0["coord_sys"])
        f.write("0 %d %.17g\n%d\n" % (len(boxes), tlast, steps[-1]))
        for lo, hi in boxes:
            for d in range(3):
                f.write("%.17g %.17g\n" % (plo[d] + (lo[d] - (dlo[d] if d < 2 else 0)) * dx[d],
                                           plo[d] + (hi[d] + 1 - (dlo[d] if d < 2 else 0)) * dx[d]))
        f.write("Level_0/Cell\n")

    with open(os.path.join(out, "times.txt"), "w") as f:
        f.write("# k step time z_center\n")
        for k in range(nt):
            f.write("%d %d %.17g %.17g\n" % (k, steps[k], times[k], zlo + (k + 0.5) * dz))

    print("wrote %s: %d x %d x %d cells, z in [%g, %g], dz = %g"
          % (out, nx, ny, nt, zlo, zhi, dz))


def selftest(out, files):
    """Check every z slice of the stacked file against its 2D input, bit for bit."""
    h3 = read_header(out)
    boxes, fabs = read_cell_h(os.path.join(out, "Level_0"))
    slices = {}
    for fname, off in fabs:
        lo, hi, data = read_fab(os.path.join(out, "Level_0", fname), off)
        for k in range(lo[2], hi[2] + 1):
            slices[k] = data[:, k - lo[2]]
    with open(os.path.join(out, "times.txt")) as f:
        rows = [ln.split() for ln in f if not ln.startswith("#")]
    bystep = {}
    for p in files:
        h = read_header(p)
        bystep[h["step"]] = (h, p)
    for k, step, _, _ in rows:
        h, p = bystep[int(step)]
        ref = read_level0(p, h)[[h["varnames"].index(v) for v in h3["varnames"]]]
        if not np.array_equal(ref, slices[int(k)]):
            sys.exit("selftest FAILED at slice %s (%s)" % (k, p))
    print("selftest passed: %d slices match their inputs exactly" % len(rows))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfiles", nargs="+", help="2D plotfile directories")
    ap.add_argument("-o", "--output", default="spacetime_plt", help="output plotfile name")
    ap.add_argument("--vars", nargs="+", help="variables to keep (default: all)")
    ap.add_argument("--every", type=int, default=1, help="use every n-th plotfile")
    ap.add_argument("--tscale", type=float, default=1.0,
                    help="scale factor for the time axis, to set the aspect ratio")
    ap.add_argument("--index", action="store_true",
                    help="use the slice index as z (dz = 1) instead of time")
    ap.add_argument("--kchunk", type=int, default=32, help="z slices per output box")
    ap.add_argument("--selftest", action="store_true",
                    help="after writing, check every slice against its input")
    args = ap.parse_args()
    if args.every < 1 or args.kchunk < 1 or args.tscale <= 0.0:
        ap.error("--every and --kchunk must be >= 1 and --tscale > 0")

    stack(args.plotfiles, args.output, args.vars, args.every, args.tscale,
          args.index, args.kchunk)
    if args.selftest:
        selftest(args.output, args.plotfiles)


if __name__ == "__main__":
    main()
