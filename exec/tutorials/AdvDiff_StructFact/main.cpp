#include <iostream>
#include <cmath>    // for pow()
#include <random>   // for random number engine and distribution

#include <fftw3.h>

using namespace std;

int main() {

    /*
      SOLVE:

      dphi/dn = -u*dphi/dx + mu*Laplacian(phi) + sqrt(2*mu)*div(W)

    */

///////////////////////////////////////////////////////

    // 0 = centered advection
    // 1 = upwind advection
    // 2 = MOL Godunov
    // 3 = Characteristic Godunov
    int integrator = 0;

    // grid spacing
    double dx = 1.;

    // velocity (constant in space and time)
    double u = 1;

    // viscosity
    double mu = 1.;

    // time step
    double dt = 0.1;

    // # of steps
    int nsteps = 1000000;

///////////////////////////////////////////////////////

    // volume
    double vol = pow(dx,3);

    // store the solution
    // the valid region spans from indices 1 through 64
    // index 0 and 65 are ghost cells (periodic)
    double phiold[66], phinew[66];

    // initialize solution to zero and then fill periodic ghost cells
    for (int i=1; i<=64; ++i) {
        phiold[i] = 0.;
    }
    phiold[0]  = phiold[64];
    phiold[65] = phiold[1];

    // stores the solution update for the predictor and corrector stages
    // we use indices 1 through 64
    double updateold[65], updatenew[65];

    // random number engine and normal distribution (mean 0, standard deviaion 1)
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,1.0);

    // store random numbers on faces
    // we only use indices 1 through 65
    // note the random number on face 65 has the same value as face 1 (periodic)
    double random[66];

    // we use indices 0 through 32 to store the structure factors
    double structFactSum[33];
    double structFactAvg[33];
    int structFactCount = 0;

    // use indices 0 through 63 to store a temporary solution to be passed into fftw
    double fft_in[64];

    for (int istep=1; istep<=nsteps; ++istep) {

        // fill random numbers on faces 1 through 64
        // random number on face 65 has the same value as face 1 (periodic)
        for (int i=1; i<=64; ++i) {
            random[i] = distribution(generator);
        }
        random[65] = random[1];

        if (integrator == 0) {
            // centered advection

            // forward euler predictor
            for (int i=1; i<=64; ++i) {
                updateold[i] =
                    - dt * u * (phiold[i+1] - phiold[i-1]) / (2.*dx)
                    + dt * mu * (phiold[i+1] - 2.*phiold[i] + phiold[i-1]) / (dx*dx)
                    + sqrt(2.*mu*dt/vol) * (random[i+1] - random[i]) / dx;

                phinew[i] = phiold[i] + updateold[i];
            }
            phinew[0]  = phinew[64];
            phinew[65] = phinew[1];

            // trapezoidal corrector
            for (int i=1; i<=64; ++i) {
                updatenew[i] =
                    - dt * u * (phinew[i+1] - phinew[i-1]) / (2.*dx)
                    + dt * mu * (phinew[i+1] - 2.*phinew[i] + phinew[i-1]) / (dx*dx)
                    + sqrt(2.*mu*dt/vol) * (random[i+1] - random[i]) / dx;

                phinew[i] = phiold[i] + 0.5*(updateold[i]+updatenew[i]);
            }
            phinew[0]  = phinew[64];
            phinew[65] = phinew[1];


        }
        else if (integrator == 1) {
            // upwind advection

            // forward euler predictor
            for (int i=1; i<=64; ++i) {
                updateold[i] =
                    - dt * u * (phiold[i] - phiold[i-1]) / dx
                    + dt * mu * (phiold[i+1] - 2.*phiold[i] + phiold[i-1]) / (dx*dx)
                    + sqrt(2.*mu*dt/vol) * (random[i+1] - random[i]) / dx;

                phinew[i] = phiold[i] + updateold[i];
            }
            phinew[0]  = phinew[64];
            phinew[65] = phinew[1];

            // trapezoidal corrector
            for (int i=1; i<=64; ++i) {
                updatenew[i] =
                    - dt * u * (phinew[i] - phinew[i-1]) / dx
                    + dt * mu * (phinew[i+1] - 2.*phinew[i] + phinew[i-1]) / (dx*dx)
                    + sqrt(2.*mu*dt/vol) * (random[i+1] - random[i]) / dx;

                phinew[i] = phiold[i] + 0.5*(updateold[i]+updatenew[i]);
            }
            phinew[0]  = phinew[64];
            phinew[65] = phinew[1];

        }
        else if (integrator == 2) {
            // MOL Godunov

            // compute the slopes in indicies 1 through 64
            // copy slope into ghost cell at index 0
            double slope[66];

            // stores the face values at indicies 1 through 65
            double phiface[66];

            // forward euler predictor

            // 2nd-order centered slopes
            for (int i=1; i<=64; ++i) {
                slope[i] = (phiold[i+1]-phiold[i-1]) / (2.*dx);
            }
            slope[0]  = slope[64];

            // compute face values
            for (int i=1; i<=65; ++i) {
                phiface[i] = phiold[i-1] + 0.5*dx*slope[i-1];
            }

            // update
            for (int i=1; i<=64; ++i) {
                updateold[i] =
                    - dt * u * (phiface[i+1] - phiface[i]) / dx
                    + dt * mu * (phiold[i+1] - 2.*phiold[i] + phiold[i-1]) / (dx*dx)
                    + sqrt(2.*mu*dt/vol) * (random[i+1] - random[i]) / dx;

                phinew[i] = phiold[i] + updateold[i];
            }
            phinew[0]  = phinew[64];
            phinew[65] = phinew[1];

            // trapezoidal corrector

            // 2nd-order centered slopes
            for (int i=1; i<=64; ++i) {
                slope[i] = (phinew[i+1]-phinew[i-1]) / (2.*dx);
            }
            slope[0]  = slope[64];

            // compute face values
            for (int i=1; i<=65; ++i) {
                phiface[i] = phinew[i-1] + 0.5*dx*slope[i-1];
            }

            // update
            for (int i=1; i<=64; ++i) {
                updatenew[i] =
                    - dt * u * (phiface[i+1] - phiface[i]) / dx
                    + dt * mu * (phinew[i+1] - 2.*phinew[i] + phinew[i-1]) / (dx*dx)
                    + sqrt(2.*mu*dt/vol) * (random[i+1] - random[i]) / dx;

                phinew[i] = phiold[i] + 0.5*(updateold[i]+updatenew[i]);
            }
            phinew[0]  = phinew[64];
            phinew[65] = phinew[1];

        }
        else if (integrator == 3) {
            // Characteristic Godunov

            // compute the slopes in indicies 1 through 64
            // copy slope into ghost cell at index 0
            double slope[66];

            // stores the face values at indicies 1 through 65
            double phiface[66];

            // forward euler predictor

            // 2nd-order centered slopes
            for (int i=1; i<=64; ++i) {
                slope[i] = (phiold[i+1]-phiold[i-1]) / (2.*dx);
            }
            slope[0]  = slope[64];

            // compute face values
            for (int i=1; i<=65; ++i) {
                phiface[i] = phiold[i-1] + (0.5*dx - 0.5*dt)*slope[i-1];
            }

            // update
            for (int i=1; i<=64; ++i) {
                updateold[i] =
                    - dt * u * (phiface[i+1] - phiface[i]) / dx
                    + dt * mu * (phiold[i+1] - 2.*phiold[i] + phiold[i-1]) / (dx*dx)
                    + sqrt(2.*mu*dt/vol) * (random[i+1] - random[i]) / dx;

                phinew[i] = phiold[i] + updateold[i];
            }
            phinew[0]  = phinew[64];
            phinew[65] = phinew[1];

            // trapezoidal corrector

            // (resuse same face values as predictor)

            // update
            for (int i=1; i<=64; ++i) {
                updatenew[i] =
                    - dt * u * (phiface[i+1] - phiface[i]) / dx
                    + dt * mu * (phinew[i+1] - 2.*phinew[i] + phinew[i-1]) / (dx*dx)
                    + sqrt(2.*mu*dt/vol) * (random[i+1] - random[i]) / dx;

                phinew[i] = phiold[i] + 0.5*(updateold[i]+updatenew[i]);
            }
            phinew[0]  = phinew[64];
            phinew[65] = phinew[1];

        }

        // copy phinew into phiold (including ghost cells)
        for (int i=0; i<=65; ++i) {
            phiold[i] = phinew[i];
        }

        // FFT - copy into temporary array with 0-based indexing
        for (int i=0; i<64; ++i) {
            fft_in[i] = phinew[i+1];
        }

        // collect stats for structure factor after the first 50% of the run
        if (istep > nsteps/2) {

            // perform the FFT
            int N = 64;
            fftw_complex out[N/2+1];
            fftw_plan p1 = fftw_plan_dft_r2c_1d(N, fft_in, out, FFTW_ESTIMATE);
            fftw_execute(p1);

            // compute structure factor
            structFactCount++;
            for (int i=0; i<33; ++i) {
                structFactSum[i] += (out[i][0]*out[i][0] + out[i][1]*out[i][1]) / 64;
                structFactAvg[i] = structFactSum[i] / structFactCount;
            }

            fftw_destroy_plan(p1);
        }
    }

    for (int i=0; i<33; ++i) {
        std::cout << structFactAvg[i] << std::endl;
    }
}
