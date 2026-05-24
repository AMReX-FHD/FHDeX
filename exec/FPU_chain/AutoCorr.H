#include <vector>
#include <cstddef>
#include <algorithm>
#include <iostream>
#include <random>
#include <cmath>

class AutoCorr {
public:
    // M: number of ensembles, maxLag: K, dt: sampling interval
    AutoCorr(std::size_t M, std::size_t maxLag, double dt)
        : M_(M), K_(maxLag), dt_(dt),
          ensembles_(M),
          sumR_(maxLag + 1, 0.0),
          E_(maxLag + 1, 0)
    {
        for (auto &e : ensembles_) {
            e.N = 0;
            e.S.assign(K_ + 1, 0.0);
            e.last.assign(K_, 0.0);
            e.curR.assign(K_ + 1, 0.0);
        }
    }

    // Add a new sample x for ensemble with index id (0 <= id < M)
    void addSample(std::size_t id, double x) {
        auto &e = ensembles_.at(id);
        std::size_t oldN = e.N;

        // Update cross-terms: S[k] += J[N - k] * x for k=1..min(K, oldN)
        std::size_t maxk = std::min(K_, oldN);
        for (std::size_t k = 1; k <= maxk; ++k) {
            e.S[k] += e.last[k - 1] * x;
        }
        // k = 0 term
        e.S[0] += x * x;

        // Shift ring buffer: last[i] <- previous sample
        double prev = x;
        for (std::size_t i = 0; i < K_; ++i) {
            double tmp = e.last[i];
            e.last[i] = prev;
            prev = tmp;
        }

        e.N = oldN + 1;

        // Recompute R for k = 0..min(K, N-1)
        std::size_t newMaxK = std::min(K_, e.N - 1);

        for (std::size_t k = 0; k <= newMaxK; ++k) {
            double Rnew = dt_ * e.S[k] / static_cast<double>(e.N - k);

            if (oldN > k) {
                // oldR existed; update aggregate
                double delta = Rnew - e.curR[k];
                sumR_[k] += delta;
                e.curR[k] = Rnew;
            } else {
                // k == oldN is new for this ensemble
                sumR_[k] += Rnew;
                e.curR[k] = Rnew;
                if (k <= K_) {
                    E_[k] += 1;
                }
            }
        }
    }

    // Return ensemble-averaged autocorrelation R_bar_avg[k] for k = 0..K
    std::vector<double> getRBarAvg() const {
        std::vector<double> out(K_ + 1, 0.0);
        for (std::size_t k = 0; k <= K_; ++k) {
            if (E_[k] > 0) out[k] = sumR_[k] / static_cast<double>(E_[k]);
        }
        return out;
    }

    std::size_t maxLag() const { return K_; }

private:
    struct EnsembleState {
        std::size_t N;              // samples seen so far for this ensemble
        std::vector<double> S;      // S[k] = sum_j J[j] J[j+k], length K+1
        std::vector<double> last;   // ring buffer of last K samples, length K
        std::vector<double> curR;   // current R_e[k] for k=0..K, length K+1
    };

    const std::size_t M_;
    const std::size_t K_;
    const double dt_;
    std::vector<EnsembleState> ensembles_;

    // Aggregator across ensembles
    std::vector<double> sumR_;        // sum of R_e[k] across ensembles that have data for k
    std::vector<std::size_t> E_;      // number of ensembles with N_e > k
};
