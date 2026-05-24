#ifndef FIX_AVE_CORRELATE_H
#define FIX_AVE_CORRELATE_H

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>

/**
 * FixAveCorrelate — autocorrelation of a single scalar with running average.
 *
 * Computes C(d) = <v(t) * v(t + d*Nevery)> averaged over all time origins.
 *
 * Usage:
 *   1. Construct with (nevery, nrepeat, nfreq).
 *   2. Call end_of_step(timestep, value) every simulation timestep.
 *   3. When it returns true, an output period completed.
 *   4. Use get_correlation(d) to retrieve the running-averaged autocorrelation.
 */
class FixAveCorrelate {
public:
    FixAveCorrelate(int nevery, int nrepeat, int nfreq)
        : nevery_(nevery), nrepeat_(nrepeat), nfreq_(nfreq)
    {
        if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
            throw std::invalid_argument("Nevery, Nrepeat, Nfreq must be > 0");
        if (nfreq % nevery != 0)
            throw std::invalid_argument("Nfreq must be a multiple of Nevery");
        if (static_cast<long long>(nrepeat) * nevery > nfreq)
            throw std::invalid_argument("Nrepeat*Nevery must be <= Nfreq");

        nsample_max_ = nfreq_ / nevery_;
        ring_.resize(nsample_max_, 0.0);
        corr_.assign(nrepeat_, 0.0);
        count_.assign(nrepeat_, 0);
        corr_total_.assign(nrepeat_, 0.0);
        count_total_.assign(nrepeat_, 0);
    }

    bool end_of_step(long long timestep, double value) {
        output_ready_ = false;
        if (timestep % nevery_ != 0) return false;

        int idx = nsample_ % nsample_max_;
        ring_[idx] = value;
        nsample_++;

        int max_delay = std::min(nrepeat_, nsample_);
        for (int d = 0; d < max_delay; d++) {
            int origin = ((nsample_ - 1 - d) % nsample_max_ + nsample_max_) % nsample_max_;
            corr_[d] += ring_[origin] * ring_[idx];
            count_[d]++;
        }

        if (timestep > 0 && timestep % nfreq_ == 0) {
            output_ready_ = true;
            last_output_timestep_ = timestep;

            for (int d = 0; d < nrepeat_; d++) {
                corr_total_[d] += corr_[d];
                count_total_[d] += count_[d];
                corr_[d] = 0.0;
                count_[d] = 0;
            }
            nsample_ = 0;
        }

        return output_ready_;
    }

    double get_correlation(int d) const {
        return (count_total_[d] > 0) ? corr_total_[d] / count_total_[d] : 0.0;
    }

    long long get_count(int d) const { return count_total_[d]; }

    bool output_ready() const { return output_ready_; }

    void print_output(std::ostream& os = std::cout) const {
        os << "# Timestep: " << last_output_timestep_ << "\n";
        os << "# Delay  C(d)  Count\n";
        for (int d = 0; d < nrepeat_; d++) {
            os << std::setw(8) << static_cast<long long>(d) * nevery_
               << "  " << std::scientific << std::setprecision(6)
               << get_correlation(d)
               << "  " << count_total_[d] << "\n";
        }
    }

    void write_output(const std::string& filename) const {
        std::ofstream ofs(filename);
        if (!ofs) throw std::runtime_error("Cannot open file: " + filename);
        print_output(ofs);
    }

private:
    int nevery_, nrepeat_, nfreq_;
    int nsample_max_;
    int nsample_ = 0;

    std::vector<double> ring_;
    std::vector<double> corr_, corr_total_;
    std::vector<long long> count_, count_total_;

    bool output_ready_ = false;
    long long last_output_timestep_ = 0;
};

#endif // FIX_AVE_CORRELATE_H
