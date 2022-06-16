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
            return meanfield_infpop_gamma(alpha, rho, r0, max_terms);
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
        double r0;
        std::pair<unsigned int, unsigned int> qp;
        std::vector<std::complex<double>> w;
        std::vector<std::complex<double>> a;
        std::complex<double> r;
        double s;
    };

    /**
     * @brief simulate and return transmission times
     * @param engine RNG engine
     * @param psi transmission time distribution
     * @param T stopping time
     * @param Mn number of network instances to simulate on
     * @param Ms number of simulations per network instance
     * @param args parameters to pass to the network
     * @return a pair of pairs containing (1) ordered transmission times
     *         and (2) total number of infections
     */
    template<typename N, typename S, typename ...Args>
    std::pair<std::vector<absolutetime_t>, std::vector<absolutetime_t>>
    simulate(rng_t& engine, transmission_time& psi, absolutetime_t T, std::size_t Mn, std::size_t Ms, Args&&...args) {
        std::vector<double> t = {};

        // Create Mn network instances
        for(std::size_t i=0; i < Mn; ++i) {
            N nw(std::forward<Args>(args)..., engine);
            // Simulate Ms times on the current network
            for(std::size_t i=0; i < Ms; ++i) {
                S sim(nw, psi);
                sim.add_infections({ std::make_pair(0, 0.0)});

                // Run simulation, collect transmission times
                while (true) {
                    auto point = sim.step(engine);
                    if (point.second > T)
                        break;
                    t.push_back(point.second);
                }
            }
        }

        // Sort transmission times (only necessary for Mn or Ms > 1)
        std::sort(t.begin(), t.end());

        // Create vector representing the total number of infections
        // at each time point
        std::vector<double> y = { (1.0/Mn) * (1.0/Ms) };
        for(std::size_t i=1; i < t.size(); ++i)
            y.push_back(y.back() + (1.0/Mn) * (1.0/Ms) );

        return std::make_pair(t, y);
    }
}

#if ENABLE_PLOTTING && 0
TEST_CASE("Plot large-population mean-field solution for Gamma transmission times on a fully-connected network", "[meanfield]") {
    using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t Mn = 50;
    const std::size_t Ms = 10;
    const std::size_t N1 = 100;
    const std::size_t N2 = 1000;
    const std::size_t N3 = 10000;
    const std::size_t T = 25;
    const std::size_t X = 400;
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
    std::vector<double> y_analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for(std::size_t i=0; i < X; ++i) {
        t_analytical.push_back((double)T * i / (X-1));
        y_analytical.push_back(sol.N(t_analytical.back()));
    }

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

#if ENABLE_PLOTTING && 0
TEST_CASE("Plot large-population mean-field solution for Gamma transmission times on an Erdös-Reyni network", "[meanfield]") {
    using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t N = 10000;
    const std::size_t M = 100;
    const std::size_t T = 35;
    const std::size_t X = 400;
    const double R0 = 4;
    // Every infecteced node has N-1 neighbours, of which in the large-population limit N-2 are susceptible.
    // To trigger subsequent infections amgonst these N-2 susceptible neighbours, we must infect each neighbour
    // with probability R0/(N-2). Note that for the first infected node, this is not strictly speaking correct,
    // since it's infection has no source, and it will thus create R0(N-1)/(N-2) > R0 subsequent infections.
    // This only causes a relative error of |1 - (N-1)/(N-2)| =~= 1/N though.
    const double MEAN = 10;
    const double VARIANCE = 1;
    transmission_time_gamma psi(MEAN, VARIANCE);

    /* Simulate using next reaction M times */
    auto r1 = simulate<erdos_reyni, simulate_next_reaction>(engine, psi, T, 1, M, N, R0 + 0.9);

    /* Evaluate analytical solution */
    std::vector<double> t_analytical;
    std::vector<double> y_analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for(std::size_t i=0; i < X; ++i) {
        t_analytical.push_back((double)T * i / (X-1));
        y_analytical.push_back(sol.N(t_analytical.back()));
    }

    plt::figure_size(1600, 1200);
    plt::named_plot("next reaction", r1.first, r1.second);
    plt::named_plot("analytical", t_analytical, y_analytical);
    plt::title("Large-population mean-field solution for Gamma transmission times on an Erdös-Reyni network");
    plt::legend();
    plt::show();
}
#endif

#if ENABLE_PLOTTING && 0
TEST_CASE("Plot large-population mean-field solution for Gamma transmission times on an acyclic network", "[meanfield]") {
    using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t M = 1000;
    const std::size_t T = 55;
    const std::size_t X = 400;
    const double R0 = 2;
    // Every infecteced node has N-1 neighbours, of which in the large-population limit N-2 are susceptible.
    // To trigger subsequent infections amgonst these N-2 susceptible neighbours, we must infect each neighbour
    // with probability R0/(N-2). Note that for the first infected node, this is not strictly speaking correct,
    // since it's infection has no source, and it will thus create R0(N-1)/(N-2) > R0 subsequent infections.
    // This only causes a relative error of |1 - (N-1)/(N-2)| =~= 1/N though.
    const double MEAN = 10;
    const double VARIANCE = 1;
    transmission_time_gamma psi(MEAN, VARIANCE);

    /* Simulate using next reaction M times */
    auto r1 = simulate<acyclic, simulate_next_reaction>(engine, psi, T, M, 1, R0+1, true);

    /* Evaluate analytical solution */
    std::vector<double> t_analytical;
    std::vector<double> y_analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for(std::size_t i=0; i < X; ++i) {
        t_analytical.push_back((double)T * i / (X-1));
        y_analytical.push_back(sol.N(t_analytical.back()));
    }

    plt::figure_size(800, 600);
    plt::named_plot("next reaction", r1.first, r1.second);
    plt::named_plot("analytical", t_analytical, y_analytical);
    plt::title("Large-population mean-field solution for Gamma transmission times on an acyclic network");
    plt::legend();
    plt::show();
}
#endif

TEST_CASE("Plot different mean-field large-population limits for Gamma transmission times", "[meanfield]") {
    using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t Mfully_n = 100;
    const std::size_t Mfully_s = 10;
    const std::size_t Macyclic = 10000;
    const std::size_t Merdos = 1000;
    const std::size_t Nfully1 = 100;
    const std::size_t Nfully2 = 1000;
    const std::size_t Nfully3 = 10000;
    const std::size_t Nerdos = 10000;
    const std::size_t T = 45;
    const std::size_t X = 400;
    const double R0 = 2;
    // Every infecteced node has N-1 neighbours, of which in the large-population limit N-2 are susceptible.
    // To trigger subsequent infections amgonst these N-2 susceptible neighbours, we must infect each neighbour
    // with probability R0/(N-2). Note that for the first infected node, this is not strictly speaking correct,
    // since it's infection has no source, and it will thus create R0(N-1)/(N-2) > R0 subsequent infections.
    // This only causes a relative error of |1 - (N-1)/(N-2)| =~= 1/N though.
    const double pfully1 = R0/(Nfully1-2);
    const double pfully2 = R0/(Nfully2-2);
    const double pfully3 = R0/(Nfully3-2);
    const double MEAN = 10;
    const double VARIANCE = 1;
    transmission_time_gamma psifully1(MEAN, VARIANCE, 1.0-pfully1);
    transmission_time_gamma psifully2(MEAN, VARIANCE, 1.0-pfully2);
    transmission_time_gamma psifully3(MEAN, VARIANCE, 1.0-pfully3);
    transmission_time_gamma psi(MEAN, VARIANCE);

    /* Simulate using next reaction M times */
    auto rfully1 = simulate<fully_connected, simulate_next_reaction>(engine, psifully1, T, Mfully_n, Mfully_s, Nfully1);
    auto rfully2 = simulate<fully_connected, simulate_next_reaction>(engine, psifully2, T, Mfully_n, Mfully_s, Nfully2);
    auto rfully3 = simulate<fully_connected, simulate_next_reaction>(engine, psifully3, T, Mfully_n, Mfully_s, Nfully3);
    auto racyclic = simulate<acyclic, simulate_next_reaction>(engine, psi, T, Macyclic, 1, R0+1, true);
    auto rerdos = simulate<erdos_reyni, simulate_next_reaction>(engine, psi, T, Merdos, 1, Nerdos, R0);

    /* Evaluate analytical solution */
    std::vector<double> t_analytical;
    std::vector<double> y_analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for(std::size_t i=0; i < X; ++i) {
        t_analytical.push_back((double)T * i / (X-1));
        y_analytical.push_back(sol.N(t_analytical.back()));
    }

    plt::figure_size(1600, 1200);
    plt::named_plot("next reaction fully-connected (N="s + std::to_string(Nfully1) + ")", rfully1.first, rfully1.second);
    plt::named_plot("next reaction fully-connected (N="s + std::to_string(Nfully2) + ")", rfully2.first, rfully2.second);
    plt::named_plot("next reaction fully-connected (N="s + std::to_string(Nfully3) + ")", rfully3.first, rfully3.second);
    plt::named_plot("next reaction acyclic", racyclic.first, racyclic.second);
    plt::named_plot("next reaction Erdös-Reyni (N="s + std::to_string(Nerdos) + ")", rerdos.first, rerdos.second);
    plt::named_plot("analytical", t_analytical, y_analytical);
    plt::title("Mean-field large-population limits for Gamma transmission times");
    plt::legend();
    plt::show();
}
