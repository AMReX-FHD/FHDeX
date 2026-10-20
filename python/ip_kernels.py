#!/usr/bin/env python3
"""Interaction potentials of src_deankow/interaction and their Fourier transforms.

Shared by amplitude_eqs.py and lattice_continuation.py. The potentials match
ip_U in src_deankow/interaction/Potentials.H:

  gema   U = eps exp(-(r/R)^alpha)
  hk     U = eps R (1 - r^2/R^2)/2 for r < R, else 0
  morse  U = -eps_att exp(-(r - re)/R_att) + eps_rep exp(-(r - re)/R_rep)
  pp     U = -eps R w(r/R), w(s) = alpha (s^2 - a^2)/2 + beta (a^2 - 1)/2 for s <= a,
                              beta (s^2 - 1)/2 for a < s <= 1, else 0

Uhat(k) is the radial Fourier transform in 2D (2 pi int U J0(kr) r dr) or 3D
(4 pi int U sin(kr)/(kr) r^2 dr), computed by quadrature. With phi normalized
to int phi = 1 on the unit box (mu0 = 1), mu0 Uhat(k) is the quantity in the
code's startup printout, Uhat*cellvol.

Parameters can be read from an inputs file (read_inputs): the last occurrence
of each key wins, as with ParmParse.
"""

import numpy as np


def read_inputs(path):
    """Key -> list of string values for an AMReX inputs file; the last occurrence wins."""
    vals = {}
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if "=" not in line:
                continue
            key, rhs = line.split("=", 1)
            vals[key.strip()] = rhs.split()
    return vals


def bessel_j0(x):
    """J0(x) from the rational and asymptotic approximations of Abramowitz and
    Stegun 9.4.1 and 9.4.3 (absolute error below 1e-7)."""
    x = np.abs(np.asarray(x, float))
    out = np.empty_like(x)
    s = x <= 3.0
    y = (x[s]/3.0)**2
    out[s] = (1.0 + y*(-2.2499997 + y*(1.2656208 + y*(-0.3163866 + y*(0.0444479
              + y*(-0.0039444 + y*0.00021))))))
    xl = x[~s]
    y = 3.0/xl
    f0 = (0.79788456 + y*(-0.00000077 + y*(-0.00552740 + y*(-0.00009512 + y*(0.00137237
          + y*(-0.00072805 + y*0.00014476))))))
    t0 = (xl - 0.78539816 + y*(-0.04166397 + y*(-0.00003954 + y*(0.00262573
          + y*(-0.00054125 + y*(-0.00029333 + y*0.00013558))))))
    out[~s] = f0*np.cos(t0)/np.sqrt(xl)
    return out


class Kernel:
    """One interaction potential, U(r) and its radial transform."""

    def __init__(self, ip_type="gema", eps=666.0, R=0.1, alpha=3.0, eps_att=2.0,
                 eps_rep=1.0, R_att=0.1, R_rep=0.05, re=0.0, pp_a=0.27,
                 pp_alpha=1.0, pp_beta=0.02):
        self.ip_type = ip_type
        self.eps, self.R, self.alpha = eps, R, alpha
        self.eps_att, self.eps_rep, self.R_att, self.R_rep, self.re = eps_att, eps_rep, R_att, R_rep, re
        self.pp_a, self.pp_alpha, self.pp_beta = pp_a, pp_alpha, pp_beta
        if ip_type == "gema":
            self.rmax = 8.0*R
        elif ip_type in ("hk", "pp"):
            self.rmax = R
        elif ip_type == "morse":
            self.rmax = max(re, 0.0) + 40.0*max(R_att, R_rep)
        else:
            raise ValueError("unknown ip_type " + ip_type)
        self._cache = {}

    @classmethod
    def from_inputs(cls, path):
        v = read_inputs(path)

        def get(key, default):
            return float(v[key][0]) if key in v else default

        return cls(ip_type=v.get("amr.ip_type", ["gema"])[0],
                   eps=get("amr.ip_eps", 666.0), R=get("amr.ip_R", 0.1),
                   alpha=get("amr.ip_alpha", 3.0),
                   eps_att=get("amr.ip_eps_att", 2.0), eps_rep=get("amr.ip_eps_rep", 1.0),
                   R_att=get("amr.ip_R_att", 0.1), R_rep=get("amr.ip_R_rep", 0.05),
                   re=get("amr.ip_re", 0.0), pp_a=get("amr.ip_pp_a", 0.27),
                   pp_alpha=get("amr.ip_pp_alpha", 1.0), pp_beta=get("amr.ip_pp_beta", 0.02))

    def describe(self):
        if self.ip_type == "gema":
            return "gema eps=%g R=%g alpha=%g" % (self.eps, self.R, self.alpha)
        if self.ip_type == "hk":
            return "hk eps=%g R=%g" % (self.eps, self.R)
        if self.ip_type == "pp":
            return "pp eps=%g R=%g a=%g alpha=%g beta=%g" % (self.eps, self.R, self.pp_a,
                                                            self.pp_alpha, self.pp_beta)
        return "morse eps_att=%g eps_rep=%g R_att=%g R_rep=%g re=%g" % (
            self.eps_att, self.eps_rep, self.R_att, self.R_rep, self.re)

    def U(self, r):
        r = np.asarray(r, float)
        if self.ip_type == "gema":
            return self.eps*np.exp(-(r/self.R)**self.alpha)
        if self.ip_type == "hk":
            s = r/self.R
            return np.where(r < self.R, 0.5*self.eps*self.R*(1.0 - s*s), 0.0)
        if self.ip_type == "pp":
            s, a = r/self.R, self.pp_a
            w = np.where(s <= a, 0.5*self.pp_alpha*(s*s - a*a) + 0.5*self.pp_beta*(a*a - 1.0),
                         0.5*self.pp_beta*(s*s - 1.0))
            return np.where(r < self.R, -self.eps*self.R*w, 0.0)
        return (-self.eps_att*np.exp(-(r - self.re)/self.R_att)
                + self.eps_rep*np.exp(-(r - self.re)/self.R_rep))

    def Uhat(self, k, dim):
        """Radial Fourier transform at |k| (scalar), cached."""
        key = (dim, round(float(k), 9))
        if key not in self._cache:
            r = np.linspace(0.0, self.rmax, 8001)
            u = self.U(r)
            k = float(k)
            if dim == 3:
                val = 4*np.pi*np.trapz(u*np.sinc(k*r/np.pi)*r*r, r)
            else:
                val = 2*np.pi*np.trapz(u*bessel_j0(k*r)*r, r)
            self._cache[key] = val
        return self._cache[key]

    def critical(self, dim, kmin=1.0, kmax=None, nk=301):
        """(k_c, Uhat_min) from a scan of Uhat(k) over k > 0, refined around the minimum."""
        if kmax is None:
            kmax = 30.0/min(self.R, self.R_rep if self.ip_type == "morse" else self.R)
        ks = np.linspace(kmin, kmax, nk)
        uk = np.array([self.Uhat(k, dim) for k in ks])
        i = int(np.argmin(uk))
        lo, hi = ks[max(i - 1, 0)], ks[min(i + 1, nk - 1)]
        ks = np.linspace(lo, hi, 201)
        uk = np.array([self.Uhat(k, dim) for k in ks])
        i = int(np.argmin(uk))
        return ks[i], uk[i]


def add_kernel_args(ap):
    """Command-line options shared by the scripts that use a Kernel."""
    g = ap.add_argument_group("interaction potential (default: gema eps=666 R=0.1 alpha=3)")
    g.add_argument("--inputs", help="read amr.ip_* from this AMReX inputs file")
    g.add_argument("--ip-type", default="gema", choices=["gema", "hk", "morse", "pp"])
    g.add_argument("--eps", type=float, default=666.0)
    g.add_argument("--R", type=float, default=0.1)
    g.add_argument("--alpha", type=float, default=3.0, help="gema exponent")
    g.add_argument("--eps-att", type=float, default=2.0)
    g.add_argument("--eps-rep", type=float, default=1.0)
    g.add_argument("--R-att", type=float, default=0.1)
    g.add_argument("--R-rep", type=float, default=0.05)
    g.add_argument("--re", type=float, default=0.0)
    g.add_argument("--pp-a", type=float, default=0.27)
    g.add_argument("--pp-alpha", type=float, default=1.0)
    g.add_argument("--pp-beta", type=float, default=0.02)


def kernel_from_args(args):
    if args.inputs:
        return Kernel.from_inputs(args.inputs)
    return Kernel(ip_type=args.ip_type, eps=args.eps, R=args.R, alpha=args.alpha,
                  eps_att=args.eps_att, eps_rep=args.eps_rep, R_att=args.R_att,
                  R_rep=args.R_rep, re=args.re, pp_a=args.pp_a, pp_alpha=args.pp_alpha,
                  pp_beta=args.pp_beta)
