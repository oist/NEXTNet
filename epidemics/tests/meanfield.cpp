#include "tests/stdafx.h"

#include "random.h"
#include "NextReaction.h"
#include "nMGA.h"

#if ENABLE_PLOTTING
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
#endif

#if ENABLE_PLOTTING
TEST_CASE("Plot mean-field solution for Gamma transmission times on a fully-connected network", "[meanfield]") {
    // TODO: Implement test!
    std::mt19937 engine;

    const std::size_t N = 1000000;
    const std::size_t M = 100;
    const double R0 = 2;
    const double p = R0/N;
    const double MEAN = 10;
    const double VARIANCE = 1;

    fully_connected network(N, engine);
    transmission_time_lognormal psi(MEAN, VARIANCE, 1.0-p);

    /* Simulate using next reaction */
    simulate_next_reaction sim_nextreaction(network, psi);
    sim_nextreaction.add_infections({ std::make_pair(0, 0.0)});

    std::vector<double> times_nextreaction = {};
    while (times_nextreaction.size() < M) {
        auto point = sim_nextreaction.step(engine);
        if (point.second == INFINITY)
            break;
        REQUIRE((times_nextreaction.empty() || (times_nextreaction.back() <= point.second)));
        times_nextreaction.push_back(point.second);
    }
    std::vector<double> y_nextreaction(times_nextreaction.size(), 0);
    std::iota(y_nextreaction.begin(), y_nextreaction.end(), 1);

    plt::figure_size(1600, 1200);
    plt::named_plot("next reaction", times_nextreaction, y_nextreaction);
    plt::title("Next reaction on a fully-connected network");
    plt::legend();
    plt::show();
}
#endif
