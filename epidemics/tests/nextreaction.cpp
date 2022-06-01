#include "tests/stdafx.h"

#include "random.h"
#include "NextReaction.h"
#include "nMGA.h"

#if ENABLE_PLOTTING
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
#endif

#if ENABLE_PLOTTING && 0
TEST_CASE("Plot infection times for next reaction scheme on Erdös-Reyni", "[nextreaction]") {
    std::mt19937 engine;

    const double N = 1000;
    const double K = 3;
    const double MEAN = 10;
    const double VARIANCE = 1;

    erdos_reyni network(N, K, engine);
    transmission_time_lognormal psi(MEAN, VARIANCE); 

    /* Simulate using next reaction */
    simulate_next_reaction sim_nextreaction(network, psi);
    sim_nextreaction.add_infections({ std::make_pair(0, 0.0)});

    std::vector<double> times_nextreaction = {};
    while (true) {
        auto point = sim_nextreaction.step(engine);
        if (point.second == INFINITY)
            break;
        REQUIRE((times_nextreaction.empty() || (times_nextreaction.back() <= point.second)));
        times_nextreaction.push_back(point.second);
    }
    std::vector<double> y_nextreaction(times_nextreaction.size(), 0);
    std::iota(y_nextreaction.begin(), y_nextreaction.end(), 1);

    /* Simulate using exact nMGA */
    simulate_nmga sim_ngma_exact(network, psi);
    sim_ngma_exact.approximation_threshold = 10*N;
    sim_ngma_exact.add_infections({ std::make_pair(0, 0.0)});

    std::vector<double> times_nmga_exact = {};
    while (true) {
        auto point = sim_ngma_exact.step(engine);
        if (point.second == INFINITY)
            break;
        REQUIRE((times_nmga_exact.empty() || (times_nmga_exact.back() <= point.second)));
        times_nmga_exact.push_back(point.second);
    }
    std::vector<double> y_nmga_exact(times_nmga_exact.size(), 0);
    std::iota(y_nmga_exact.begin(), y_nmga_exact.end(), 1);

    /* Simulate using approximate nMGA */
    simulate_nmga sim_ngma_approx(network, psi);
    sim_ngma_approx.approximation_threshold = 200;
    sim_ngma_approx.add_infections({ std::make_pair(0, 0.0)});

    std::vector<double> times_nmga_approx = {};
    while (true) {
        auto point = sim_ngma_approx.step(engine);
        if (point.second == INFINITY)
            break;
        REQUIRE((times_nmga_approx.empty() || (times_nmga_approx.back() <= point.second)));
        times_nmga_approx.push_back(point.second);
    }
    std::vector<double> y_nmga_approx(times_nmga_approx.size(), 0);
    std::iota(y_nmga_approx.begin(), y_nmga_approx.end(), 1);

    std::vector<double> y(N, 0);
    std::iota(y.begin(), y.end(), 1);

    plt::figure_size(1600, 1200);
    plt::named_plot("next reaction", times_nextreaction, y_nextreaction);
    plt::named_plot("nMGA exact", times_nmga_exact, y_nmga_exact);
    plt::named_plot("nMGA approx (th. 100)" , times_nmga_approx, y_nmga_approx);
    plt::title("Next reaction vs. nGMA for Erdös-Reyni");
    plt::legend();
    plt::show();
}
#endif
