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
    std::mt19937 engine;

    const std::size_t M = 100;
    const std::size_t N = 10000;
    const std::size_t T = 45;
    const double R0 = 2;
    const double p = R0/N;
    const double MEAN = 10;
    const double VARIANCE = 1;
    transmission_time_lognormal psi(MEAN, VARIANCE, 1.0-p);

    /* Simulate using next reaction M times */
    std::vector<double> times_nextreaction = {};
    for(std::size_t i=0; i < M; ++i) {
        fully_connected network(N, engine);
        simulate_next_reaction sim_nextreaction(network, psi);
        sim_nextreaction.add_infections({ std::make_pair(0, 0.0)});

        absolutetime_t last = 0.0;
        while (last < T) {
            auto point = sim_nextreaction.step(engine);
            if (point.second == INFINITY)
                break;
            REQUIRE(last <= point.second);
            last = point.second;
            times_nextreaction.push_back(point.second);
        }
    }

    /* Convert list of infection times into a graph */
    std::sort(times_nextreaction.begin(), times_nextreaction.end());
    std::vector<double> y_nextreaction(times_nextreaction.size(), 0);
    std::iota(y_nextreaction.begin(), y_nextreaction.end(), 1);

    plt::figure_size(1600, 1200);
    plt::named_plot("next reaction", times_nextreaction, y_nextreaction);
    plt::title("Next reaction on a fully-connected network");
    plt::legend();
    plt::show();
}
#endif
