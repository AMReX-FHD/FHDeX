#!/usr/bin/env python3
"""Free energy of the density field in 2D or 3D interaction plotfiles.

The deterministic part of the src_deankow/interaction model is the gradient flow of

  F[phi] = D int phi ln phi dx + 1/2 int int U(x - y) phi(x) phi(y) dx dy,

with phi normalized to int phi = 1 (MODEL.tex, section 1.2). For each plotfile
this evaluates the discrete form used by the code: the mesh kernel U sampled
at the minimum-image cell offsets (as in mykernel.H), C = U * phi by FFT with
quadrature weight dV, and

  F_ent = D sum phi ln phi dV   (phi ln phi taken as 0 where phi <= 0),
  F_int = 1/2 sum phi C dV,     F = F_ent + F_int.

The uniform reference with the same mass, mu = sum phi dV / V, is
F_uni = D mu ln mu V + 1/2 mu^2 Uhat(0) dV V, the reference of
lattice_continuation.py. F is per unit mass, in units where D is the
temperature; the free energy of the N particles is about N F.

How to read it:
  * Without noise (dorand = 0) F never increases.
  * With noise even a uniform state fluctuates. Equipartition under
    P[phi] ~ exp(-N F/D) gives <F> - F_uni ~ (M - 1) D/(2N) for M cells; this
    noise offset is written in the header. A noisy run starts near it and
    falls below it as clusters or a crystal form.
  * Compare F - F_uni with the crystal branches of lattice_continuation.py,
    which cross 0 at D_F (2D hexagons 0.639, 3D bcc 0.432 for the gema runs).

Output (-o, default free_energy.txt): a header, then one line per plotfile:
  t step F F-F_uni F_ent F_int mean_phi max_phi neg_mass
where neg_mass = sum min(phi, 0) dV, the mass in cells with phi < 0.

Example:
    python3 free_energy.py near_3d/plt_fv3d_* --inputs inputs_fv_3d --D 0.30 \\
        --num-part 2e7 -o near_3d_F.txt
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ip_kernels import add_kernel_args, kernel_from_args, read_inputs  # noqa: E402
from stack_plotfiles_time import read_header, read_level0  # noqa: E402


def mesh_kernel_hat(ker, shape, dx):
    """Uhat dV of the kernel sampled at minimum-image cell offsets; shape and dx
    in array (reversed) axis order."""
    xs = np.meshgrid(*[np.fft.fftfreq(n)*n*h for n, h in zip(shape, dx)], indexing="ij")
    r = np.sqrt(sum(x*x for x in xs))
    return np.real(np.fft.fftn(ker.U(r)))*float(np.prod(dx))


def free_energy(phi, uh_dv, D, dv):
    """F, F_ent, F_int of one frame."""
    C = np.real(np.fft.ifftn(uh_dv*np.fft.fftn(phi)))
    pos = phi > 0
    f_ent = D*np.sum(phi[pos]*np.log(phi[pos]))*dv
    f_int = 0.5*np.sum(phi*C)*dv
    return f_ent + f_int, f_ent, f_int


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfiles", nargs="+")
    ap.add_argument("--var", default="phi0", help="density variable (default phi0)")
    ap.add_argument("--every", type=int, default=1, help="use every n-th plotfile")
    ap.add_argument("--D", type=float, help="diffusion coefficient (default: diff_coeff from --inputs)")
    ap.add_argument("--num-part", type=float, help="num_part, for the noise offset (default: from --inputs)")
    ap.add_argument("-o", "--output", default="free_energy.txt")
    add_kernel_args(ap)
    args = ap.parse_args()

    vals = read_inputs(args.inputs) if args.inputs else {}
    if args.inputs and int(float(vals.get("amr.use_ext_pot", ["0"])[0])):
        sys.exit("amr.use_ext_pot = 1 in %s: the external potential is not included" % args.inputs)
    D = args.D if args.D is not None else float(vals["diff_coeff"][0]) if "diff_coeff" in vals else None
    if D is None:
        sys.exit("give --D, or --inputs with diff_coeff")
    npart = args.num_part if args.num_part is not None else (
        float(vals["num_part"][0]) if "num_part" in vals else None)
    ker = kernel_from_args(args)

    frames = sorted((read_header(f)["step"], f) for f in args.plotfiles)[::args.every]
    h = read_header(frames[0][1])
    dim = h["dim"]
    dlo, dhi = h["domain"]
    shape = [b - a + 1 for a, b in zip(dlo, dhi)][::-1]
    dx = h["dx"][:dim][::-1]
    dv = float(np.prod(dx))
    M = int(np.prod(shape))
    vol = M*dv
    uh_dv = mesh_kernel_hat(ker, shape, dx)
    offset = (M - 1)*D/(2*npart) if npart else float("nan")

    head = ["free energy F = D sum phi ln phi dV + 1/2 sum phi (U*phi) dV, per unit mass",
            "potential: %s; D = %g; num_part = %s; %dD grid %s; Uhat(0) dV = %.6g"
            % (ker.describe(), D, "%g" % npart if npart else "unknown", dim,
               "x".join(map(str, shape[::-1])), uh_dv.flat[0]),
            "noise offset (M-1) D/(2N) = %.6g" % offset,
            "t step F F-F_uni F_ent F_int mean_phi max_phi neg_mass"]
    with open(args.output, "w") as out:
        for line in head:
            out.write("# " + line + "\n")
            print(line)
        for step, f in frames:
            hf = read_header(f)
            phi = read_level0(f, hf)[hf["varnames"].index(args.var)]
            mu = phi.sum()*dv/vol
            funi = (D*mu*np.log(mu) + 0.5*mu*mu*uh_dv.flat[0])*vol
            F, fe, fi = free_energy(phi, uh_dv, D, dv)
            row = "%.10g %d %.10g %.6e %.10g %.10g %.8g %.6g %.4e" % (
                hf["time"], step, F, F - funi, fe, fi, mu, phi.max(), np.minimum(phi, 0).sum()*dv)
            out.write(row + "\n")
            print(row, flush=True)


if __name__ == "__main__":
    main()
