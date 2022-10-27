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
    template<typename R>
    std::vector<R> parallel(std::size_t n, rng_t& engine, std::function<R(rng_t&)> body) {
        typedef boost::counting_iterator<std::size_t> count_it;

        std::vector<R> r;
        r.reserve(n);

        sub_rngs rngs(n, engine);
        std::transform(/* std::execution::par, */
                       count_it(0), count_it(n), std::back_inserter(r),
                       [&body, &rngs](std::size_t i) -> R { return body(rngs[i]); });
        return r;
    }

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
    simulate_SIR(rng_t& engine, const transmission_time& psi, absolutetime_t T, std::size_t Mn, std::size_t Ms, Args&&...args) {
        auto ts = parallel<std::vector<double>>(Mn, engine, [&psi, T, Ms, &args...](rng_t& thread_engine){
            std::vector<double> t = {};

            // Simulate Ms times on one network
            N nw(std::forward<Args>(args)..., thread_engine);
            for(std::size_t i=0; i < Ms; ++i) {
                S sim(nw, psi);
                sim.add_infections({ std::make_pair(0, 0.0)});

                // Run simulation, collect transmission times
                while (true) {
                    auto point = sim.step(thread_engine);
                    if (!point || (point->time > T))
                        break;
                    t.push_back(point->time);
                }
            }

            return t;
        });

        // Merge results
        std::vector<double> t;
        for(const auto& tp: ts)
          std::copy(tp.begin(), tp.end(), std::back_inserter(t));

        // Sort transmission times (only necessary for Mn or Ms > 1)
        std::sort(t.begin(), t.end());

        // Create vector representing the total number of infections
        // at each time point
        std::vector<double> y = { (1.0/Mn) * (1.0/Ms) };
        for(std::size_t i=1; i < t.size(); ++i)
            y.push_back(y.back() + (1.0/Mn) * (1.0/Ms) );

        return std::make_pair(t, y);
    }


    /**
     * @brief simulate and return transmission times
     * @param engine RNG engine
     * @param psi transmission time distribution
     * @param T stopping time
     * @param args parameters to pass to the network
     * @return a pair of pairs containing (1) ordered transmission times
     *         and (2) total number of infections
     */
    template<typename N, typename S, typename ...Args>
    void
    simulate_SIS(rng_t& engine, transmission_time& psi, transmission_time& rho,
                 std::vector<absolutetime_t> &times, std::vector<double> &infections, std::vector<double>& infected,
                 absolutetime_t T, Args&&...args)
{
        N nw(std::forward<Args>(args)..., engine);
        S sim(nw, psi, &rho);
        sim.add_infections({ std::make_pair(0, 0.0)});
        int current_infected = 0, total_infected = 0;
        // Run simulation, collect transmission times
        while (true) {
    
            auto point = sim.step(engine);
            if (!point || (point -> time > T))
                break;
            
            switch (point-> kind) {
                case event_kind::infection:
                    ++total_infected;
                    ++current_infected;
                    break;
                case event_kind::reset:
                    --current_infected;
                    break;
                default:
                    throw std::logic_error("unexpected event kind");
            }

            times.push_back(point->time);
            infections.push_back(total_infected);
            infected.push_back(current_infected);
        }
    }
}

#if ENABLE_PLOTTING
TEST_CASE("Plot large-population mean-field SIR solutions for Gamma transmission times on different networks", "[meanfield]") {
    using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t M = 100;
    const std::size_t Nfully1 = 100;
    const std::size_t Nfully2 = 1000;
    const std::size_t Nfully3 = 10000;
    const std::size_t Nerdos = 10000;
    const std::size_t T = 65;
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
    auto rfully1 = simulate_SIR<fully_connected, simulate_next_reaction>(engine, psifully1, T, M, 1, Nfully1);
    auto rfully2 = simulate_SIR<fully_connected, simulate_next_reaction>(engine, psifully2, T, M, 1, Nfully2);
    auto rfully3 = simulate_SIR<fully_connected, simulate_next_reaction>(engine, psifully3, T, M, 1, Nfully3);
    auto racyclic = simulate_SIR<acyclic, simulate_next_reaction>(engine, psi, T, M, 1, R0+1, true);
    auto rerdos = simulate_SIR<erdos_reyni, simulate_next_reaction>(engine, psi, T, M, 1, Nerdos, R0);

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
    plt::title("Large-population SIR mean-field");
    plt::legend();
    plt::save("tests.out/meanfield.sir.pdf");
}
#endif

#if ENABLE_PLOTTING && 0
TEST_CASE("Plot large-population mean-field SIS solution for Gamma transmission and recovery times on an fully connected network", "[meanfield]") {
       using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t N = 5000;
    const std::size_t T = 2000;
    const std::size_t X = 400;
    const double R0 = 5;
    // Every infecteced node has N-1 neighbours, of which in the large-population limit N-2 are susceptible.
    // To trigger subsequent infections amgonst these N-2 susceptible neighbours, we must infect each neighbour
    // with probability R0/(N-2). Note that for the first infected node, this is not strictly speaking correct,
    // since it's infection has no source, and it will thus create R0(N-1)/(N-2) > R0 subsequent infections.
    // This only causes a relative error of |1 - (N-1)/(N-2)| =~= 1/N though.
    const double p = R0/(N-2);

    const double MEAN = 10;
    const double VARIANCE = 50;
    const double MEAN_rho = 70;
    const double VARIANCE_rho = 1;
    transmission_time_gamma psi(MEAN, VARIANCE, 1.0-p);
    transmission_time_gamma rho(MEAN_rho, VARIANCE_rho);

    /* Simulate using next reaction M times */
    std::vector<double> t_sim, y_sim_new, y_sim_total;

    simulate_SIS<fully_connected, simulate_next_reaction>(engine, psi,rho,t_sim, y_sim_new,y_sim_total, T, N);
    /* Evaluate analytical solution */
    std::vector<double> t_analytical;
    std::vector<double> y_analytical;

    for(std::size_t i=0; i < X; ++i) {
        const double t = (double)T * i / (X-1);
        t_analytical.push_back(t);
        y_analytical.push_back(N*(1-1/R0));
    }

    plt::figure_size(800, 600);
    plt::named_plot("next reaction", t_sim, y_sim_total);
    plt::named_plot("analytical", t_analytical, y_analytical);
    plt::title("Large-population SIS mean-field on a fully connected network");
    plt::legend();
    plt::show();
}
#endif

