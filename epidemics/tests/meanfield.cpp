#include "tests/stdafx.h"

#include "random.h"
#include "NextReaction.h"
#include "nMGA.h"
#include "utility.h"

#if ENABLE_PLOTTING
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
#endif

namespace {
    struct meanfield_infpop_gamma {
        static meanfield_infpop_gamma mean_variance(double mean, double var, double r0, std::size_t max_terms=1000) {
            const double rho = mean / var;
            const double alpha = mean * rho;
            return meanfield_infpop_gamma(alpha, rho, r0, max_terms);
        }

        static meanfield_infpop_gamma alpha_rho(double alpha, double rho, double r0, std::size_t max_terms=1000) {
            return meanfield_infpop_gamma(alpha, rho, r0, max_terms);
        }

        meanfield_infpop_gamma(double alpha_, double rho_, double r0, unsigned int max_terms = 1000)
            :alpha(alpha_), rho(rho_), qp(fraction(1/alpha, max_terms))
            ,w(qp.second, 0), a(qp.second, 0), r(0.0), s(pow(r0, 1.0 / alpha))
        {
            using namespace std::complex_literals;
            // Compute q-th power of p-ths roots w, and poles in the Laplace plane a
            for(unsigned int i=0; i < qp.second; ++i) {
                w.at(i) = exp(1i * 2.0 * M_PI * (double)i * (double)qp.first / (double)qp.second);
                a.at(i) = rho * (s * w[i] - 1.0);
                r += w.at(i) / a.at(i);
            }
        }

        double I(double t) {
            if (t < 0.0) return 0;
            if (t == 0.0) return NAN;

            std::complex<double> y = 0.0;
            for(unsigned int i=0; i < qp.second; ++i)
                y += exp(t * a[i]) * w[i];
            return abs(y) * s * rho / alpha;
        }

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
        double r0;
        std::pair<unsigned int, unsigned int> qp;
        std::vector<std::complex<double>> w;
        std::vector<std::complex<double>> a;
        std::complex<double> r;
        double s;
    };

    template<typename N, typename S, typename ...Args>
    std::pair<std::vector<absolutetime_t>, std::vector<absolutetime_t>>
    simulate(rng_t& engine, transmission_time& psi, absolutetime_t T, std::size_t Mn, std::size_t Ms, Args&&...args) {
        std::vector<double> t = {};
        for(std::size_t i=0; i < Mn; ++i) {
            N nw(std::forward<Args>(args)..., engine);
            for(std::size_t i=0; i < Ms; ++i) {
                S sim(nw, psi);
                sim.add_infections({ std::make_pair(0, 0.0)});

                while (true) {
                    auto point = sim.step(engine);
                    if (point.second > T)
                        break;
                    t.push_back(point.second);
                }
            }
        }
        std::sort(t.begin(), t.end());

        std::vector<double> y = { (1.0/Mn) * (1.0/Ms) };
        for(std::size_t i=1; i < t.size(); ++i)
            y.push_back(y.back() + (1.0/Mn) * (1.0/Ms) );

        return std::make_pair(t, y);
    }
}

#if ENABLE_PLOTTING
TEST_CASE("Plot large-population mean-field solution for Gamma transmission times on a fully-connected network", "[meanfield]") {
    using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t Mn = 100;
    const std::size_t Ms = 100;
    const std::size_t N1 = 100;
    const std::size_t N2 = 1000;
    const std::size_t N3 = 10000;
    const std::size_t T = 45;
    const double R0 = 2;
    // Every infecteced node has N-1 neighbours, of which in the large-population limit N-2 are susceptible.
    // To trigger subsequent infections amgonst these N-2 susceptible neighbours, we must infect each neighbour
    // with probability R0/(N-2). Note that for the first infected node, this is not strictly speaking correct,
    // since it's infection has no source, and it will thus create R0(N-1)/(N-2) > R0 subsequent infections.
    // This only causes a relative error of |1 - (N-1)/(N-2)| =~= 1/N though.
    const double p1 = R0/(N1-2);
    const double p2 = R0/(N2-2);
    const double p3 = R0/(N3-2);
    const double MEAN = 10;
    const double VARIANCE = 1;
    transmission_time_gamma psi1(MEAN, VARIANCE, 1.0-p1);
    transmission_time_gamma psi2(MEAN, VARIANCE, 1.0-p2);
    transmission_time_gamma psi3(MEAN, VARIANCE, 1.0-p3);

    /* Simulate using next reaction M times */
    auto r1 = simulate<fully_connected, simulate_next_reaction>(engine, psi1, T, Mn, Ms, N1);
    auto r2 = simulate<fully_connected, simulate_next_reaction>(engine, psi2, T, Mn, Ms, N2);
    auto r3 = simulate<fully_connected, simulate_next_reaction>(engine, psi3, T, Mn, Ms, N3);

    /* Evaluate analytical solution */
    std::vector<double> t_analytical;
    std::copy(r1.first.begin(), r1.first.end(), std::back_inserter(t_analytical));
    std::copy(r2.first.begin(), r2.first.end(), std::back_inserter(t_analytical));
    std::copy(r3.first.begin(), r3.first.end(), std::back_inserter(t_analytical));
    std::sort(t_analytical.begin(), t_analytical.end());
    std::vector<double> y_analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for(std::size_t i=0; i < t_analytical.size(); ++i)
        y_analytical.push_back(sol.N(t_analytical[i]));

    plt::figure_size(1600, 1200);
    plt::named_plot("next reaction (N="s + std::to_string(N1) + ")", r1.first, r1.second);
    plt::named_plot("next reaction (N="s + std::to_string(N2) + ")", r2.first, r2.second);
    plt::named_plot("next reaction (N="s + std::to_string(N3) + ")", r3.first, r3.second);
    plt::named_plot("analytical", t_analytical, y_analytical);
    plt::title("Large-population mean-field solution for Gamma transmission times on a fully-connected network");
    plt::legend();
    plt::show();
}
#endif
