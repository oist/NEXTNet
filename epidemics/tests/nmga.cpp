#include "tests/stdafx.h"

#include "random.h"
#include "nMGA.h"

#if ENABLE_PLOTTING
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
#endif

#if ENABLE_PLOTTING && 0
TEST_CASE("Plot infection times for nMGA on Erdös-Reyni", "[nMGA]") {
    std::mt19937 engine;

    const double N = 1000;
    const double K = 3;
    const double MEAN = 10;
    const double VARIANCE = 1;

    erdos_reyni network(N, K, engine);
    transmission_time_lognormal psi(MEAN, VARIANCE);

    simulate_nmga simulation(network, psi);
    simulation.approximation_threshold = N;
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
    plt::title("nMGA for Erdös-Reyni");
    plt::show();
}
#endif
