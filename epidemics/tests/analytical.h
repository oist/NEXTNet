#pragma once

#include "tests/stdafx.h"
#include "utility.h"

/**
  * @brief The large-population mean-field solution for gamma transmission times
  */
struct meanfield_infpop_gamma {
    /**
     * @brief Create solution for specified mean, variance and R0
     * @param mean mean of the transmission time distribution
     * @param var variance of the transmission time distribution
     * @param r0 average number of subsequent transmissions
     * @param max_terms number of terms used in numeric evaluation
     * @return a solution object
     */
    static meanfield_infpop_gamma mean_variance(double mean, double var, double r0, std::size_t max_terms=1000) {
        const double rho = mean / var;
        const double alpha = mean * rho;
        return meanfield_infpop_gamma(alpha, rho, r0, (unsigned int)max_terms);
    }

    /**
     * @brief Create solution for specified mean, variance and R0
     * @param alpha the shape parameter of the transmission time distribution
     * @param rho the rate parameter of the transmission time distribution
     * @param r0 average number of subsequent transmissions
     * @param max_terms number of terms used in numeric evaluation
     * @return a solution object
     */
    static meanfield_infpop_gamma alpha_rho(double alpha, double rho, double r0, std::size_t max_terms=1000) {
        return meanfield_infpop_gamma(alpha, rho, r0, (unsigned int)max_terms);
    }

    meanfield_infpop_gamma(double alpha_, double rho_, double r0, unsigned int max_terms = 1000)
        :alpha(alpha_), rho(rho_), qp(fraction(1/alpha, max_terms))
        ,w(qp.second, 0), a(qp.second, 0), r(0.0), s(pow(r0, 1.0 / alpha))
    {
        using namespace std::complex_literals;
        // Compute q-th power of p-ths roots w, and poles in the Laplace plane a
        for(unsigned int i=0; i < qp.second; ++i) {
            w.at(i) = std::exp(1.0i * 2.0 * M_PI * (double)i * (double)qp.first / (double)qp.second);
            a.at(i) = rho * (s * w[i] - 1.0);
            r += w.at(i) / a.at(i);
        }
    }

    /**
     * @brief rate of new infections at time t
     * @param t time
     * @return infection rate
     */
    double I(double t) {
        if (t < 0.0) return 0;
        if (t == 0.0) return NAN;

        std::complex<double> y = 0.0;
        for(unsigned int i=0; i < qp.second; ++i)
            y += exp(t * a[i]) * w[i];
        return abs(y) * s * rho / alpha;
    }

    /**
     * @brief total infections at time t
     * @param t
     * @return total number of infections
     */
    double N(double t) {
        if (t < 0.0) return 0;
        if (t == 0.0) return 1.0;

        std::complex<double> y = -r;
        for(unsigned int i=0; i < qp.second; ++i)
            y += exp(t * a[i]) * w[i] / a[i];
        return 1.0 + abs(y) * s * rho / alpha;
    }

private:
    double alpha;
    double rho;
    std::pair<unsigned int, unsigned int> qp;
    std::vector<std::complex<double>> w;
    std::vector<std::complex<double>> a;
    std::complex<double> r;
    double s;
};
