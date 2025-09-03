#ifndef _kernel_functions_K_H_
#define _kernel_functions_K_H_

#include <AMReX.H>

struct Kernel1P
{
    static constexpr int ks = 2;

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real operator() (Real r_in) const noexcept
    {
        Real kernel_1p;

        Real r1 = std::abs(r_in);
        Real r2 = r1 * r1;

        if (r1 <= 0.5)
        {
            kernel_1p = 1;
        }
        else
        {
            kernel_1p = 0.0;
        }

        return kernel_1p;
    }
};

struct Kernel3P
{
    static constexpr int ks = 2;

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real operator() (Real r_in) const noexcept
    {
        Real kernel_3p;

        Real r1 = std::abs(r_in);
        Real r2 = r1 * r1;

        if (r1 <= 0.5)
        {
            kernel_3p = (1 + std::sqrt(1 - (3 * r2))) / 3.0;
        }
        else if (r1 <= 1.5)
        {
            kernel_3p = (5.0 - (3.0 * r1) - std::sqrt(-3 * (1 - r1) * (1 - r1) + 1)) / 6.0;
        }
        else
        {
            kernel_3p = 0.0;
        }

        return kernel_3p;
    }
};

struct Kernel4P
{
    static constexpr int ks = 3;

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real operator() (Real r_in) const noexcept
    {
        Real kernel_4p;

        Real r2 = r_in * r_in;

        if (r_in <= -2)
        {
            kernel_4p = 0.0;
        }
        else if (r_in <= -1)
        {
            kernel_4p = 0.125 * (5.0 + (2.0 * r_in) - std::sqrt(-7.0 - (12.0 * r_in) - (4.0 * r2)));
        }
        else if (r_in <= 0)
        {
            kernel_4p = 0.125 * (3.0 + (2.0 * r_in) + std::sqrt(1.0 - (4.0 * r_in) - (4.0 * r2)));
        }
        else if (r_in <= 1)
        {
            kernel_4p = 0.125 * (3.0 - (2.0 * r_in) + std::sqrt(1.0 + (4.0 * r_in) - (4.0 * r2)));
        }
        else if (r_in <= 2)
        {
            kernel_4p = 0.125 * (5.0 - (2.0 * r_in) - std::sqrt(-7.0 + (12.0 * r_in) - (4.0 * r2)));
        }
        else
        {
            kernel_4p = 0.0;
        }

        return kernel_4p;
    }
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real beta_kernel (Real r)
{
    Real beta;
    Real K = 59.0 / 60 - std::sqrt(29.0) / 20.0;

    // pre-computed ratios
    Real a = 9.0 / 4.0;
    Real b = 3.0 / 2.0;
    Real c = 22.0 / 3.0;
    Real d = 7.0 / 3.0;

    // NOTE: mistake in the paper: b*(K+r**2)*r -> b*(K+r**2)
    // beta = a - b*(K+r**2)*r + (c-7*K)*r - d*r**3
    beta = a - b * (K + (r * r)) + (c - 7.0 * K) * r - d * (r * r * r);

    return beta;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real gamma_kernel (Real r)
{
    Real gamma;
    Real K = 59.0 / 60 - std::sqrt(29.0) / 20.0;

    // pre-computed ratios
    Real a = 11.0 / 32.0;
    Real b = 3.0 / 32.0;
    Real c = 1.0 / 72.0;
    Real d = 1.0 / 18.0;

    gamma = -a * r * r
            + b * (2 * K + (r * r)) * (r * r)
            + c * std::pow(((3 * K - 1) * r + (r * r * r)), 2)
            + d * std::pow(((4 - 3 * K) * r - (r * r * r)), 2);

    return gamma;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real phi1 (Real r)
{
    Real phi1;
    Real K = 59.0 / 60 - std::sqrt(29.0) / 20.0;
    Real sgn = (3.0 / 2.0 - K) / std::abs(3.0 / 2.0 - K);

    // pre-computed ratios
    Real alpha = 28.0;
    Real inv_alpha = 1.0 / (2.0 * alpha);

    phi1 = inv_alpha * ( -beta_kernel(r) + sgn * std::sqrt(beta_kernel(r) * beta_kernel(r) - 4 * alpha * gamma_kernel(r)) );

    return phi1;
}

struct Kernel6P
{
    static constexpr int ks = 4;

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real operator() (Real r_in) const noexcept
    {
        // The 6-point kernel function, based on the paper:
        // >*A Gaussian-like immersed-boundary kernel with three continuous derivatives and improved translational invariance*
        // >Yuanxun Bao, Jason Kaye, Charles S. Peskin, *Journal of Computational Physics* **316**, 139 (2016)
        // >https://dx.doi.org/10.1016/j.jcp.2016.04.024
        // Note also: https://github.com/stochasticHydroTools/IBMethod/blob/master/IBKernels/Kernels.c because the paper above has mistakes (but the repo's OK)

        Real kernel_6p;

        // internal parameters
        Real K = 59.0 / 60 - std::sqrt(29.0) / 20.0;
        Real sgn = (3.0 / 2 - K) / std::abs(3.0 / 2 - K);

        // pre-computed ratios
        Real inv16 = 1.0 / 16.0;
        Real inv8  = 1.0 / 8.0;
        Real inv12 = 1.0 / 12.0;
        Real inv4  = 1.0 / 4.0;
        Real inv6  = 1.0 / 6.0;
        Real rat58 = 5.0 / 8.0;

        // internal variables
        Real r = r_in;

        // compute kernel function
        if (r <= -3)
        {
            kernel_6p = 0.0;
        }
        else if (r <= -2)
        {
            r += 3;
            kernel_6p = phi1(r);
        }
        else if (r <= -1)
        {
            r += 2;
            kernel_6p = -3 * phi1(r) - inv16 + inv8 * (K + (r * r)) + inv12 * ((3 * K - 1) * r + (r * r * r));
        }
        else if (r <= 0)
        {
            r += 1;
            kernel_6p = 2 * phi1(r) + inv4 + inv6 * ((4 - 3 * K) * r - (r * r * r));
        }
        else if (r <= 1)
        {
            kernel_6p = 2 * phi1(r) + rat58 - inv4 * (K + (r * r));
        }
        else if (r <= 2)
        {
            r -= 1;
            kernel_6p = -3 * phi1(r) + inv4 - inv6 * ((4 - 3 * K) * r - (r * r * r));
        }
        else if (r <= 3)
        {
            r -= 2;
            kernel_6p = phi1(r) - inv16 + inv8 * (K + (r * r)) - inv12 * ((3 * K - 1) * r + (r * r * r));
        }
        else
        {
            kernel_6p = 0.0;
        }

        return kernel_6p;
    }
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real integral(Real beta_in, int w_in)
{
    Real beta = beta_in;
    Real w = static_cast<Real>(w_in);
    int nbin = 100;
    Real dxi = w / nbin;
    Real integral = 0.0;

    // mid-point rule for integration
    for (int i = 0; i < nbin; i++)
    {
        integral += dxi * std::exp(beta * (std::sqrt(1 - 4 * (-w * 0.5 + i * dxi + 0.5 * dxi) * (-w * 0.5 + i * dxi + 0.5 * dxi) / w / w) - 1));
    }
    // std::cout << "denominator is " << integral << std::endl;
    return integral;
}

struct KernelES
{
    // KernelES(Real beta_in, int w_in) : beta(beta_in), w(w_in) {}

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real operator() (Real r_in, Real beta, int w) const noexcept
    {
        Real kernel_es;
        Real r = r_in;
        //Real beta = beta_in;
        //Real w = static_cast<Real>(w_in);

        if (r >= -w / 2.0 && r <= w / 2.0)
        {
            kernel_es = std::exp(beta * (std::sqrt(1 - 4 * r * r / w / w) - 1));
        }
        else
        {
            kernel_es = 0.0;
        }

        // std::cout << "kernel_es is " << kernel_es << std::endl;
        return kernel_es;
    }

    // private:
    //     Real beta;
    //     int w;
};

#endif
