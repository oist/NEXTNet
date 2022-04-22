#include "tests/stdafx.h"

#include "random.h"
#include "NextReaction.h"

#if ENABLE_PLOTTING
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
#endif

#if ENABLE_PLOTTING
TEST_CASE("Plot infection times for next reaction scheme on Erdös-Reyni", "[nextreaction]") {
    std::mt19937 engine;

    const double N = 3000;
    const double K = 3;
    const double MEAN = 10;
    const double VARIANCE = 1;

    erdos_reyni network(N, K, engine);
    transmission_time_lognormal psi(MEAN, VARIANCE); 
    simulate_next_reaction simulation(network, psi);
    simulation.add_infections({ std::make_pair(0, 0.0)});

    std::vector<double> times = {};
    while (true) {
        auto point = simulation.step(engine);
        if (point.second == INFINITY)
            break;
        REQUIRE((times.empty() || (times.back() <= point.second)));
        times.push_back(point.second);
    }
    std::vector<double> y(times.size(), 0);
    std::iota(y.begin(), y.end(), 1);

    plt::plot(times, y);
    plt::show();
}
#endif
