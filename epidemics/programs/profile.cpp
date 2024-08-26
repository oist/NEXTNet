#include "random.h"
#include "analysis.h"
#include "graph.h"
#include "NextReaction.h"

#include <chrono>

using namespace std::string_literals;

int program_profile(int argc, const char * argv[])
{
    rng_t engine;

    // Parameters, which are network size (N), average number of contacts (K), 
    // mean infection transmission time (Mi) and variance (Vi),
    // mean reset time (Mr) and variance (Vr), number of steps to run (M)
    const int N = 20000;
    const std::size_t M = 100;
    const int K0 = 6;
    const double P_REWIRE = 0.1;
    const double Mi = 5;
    const double Vi = 3;
    const double Mr = 11;
    const double Vr = 5;
    std::cerr << "Running with parameters:" << std::endl;
    std::cerr << "        N: " << N << " (network size)" << std::endl;
    std::cerr << "        K: " << K0 << " (average degree)" << std::endl;
    std::cerr << " P_REWIRE: " << P_REWIRE << " (average degree)" << std::endl;
    std::cerr << "    Mi/Vi: " << Mi << " / " << Vi << " (Infection time mean / variance)" << std::endl;
    std::cerr << "    Mr/Vr: " << Mr << " / " << Vr << " (Reset time mean / variance)" << std::endl;
    std::cerr << "        M: " << M << " (number of steps)" << std::endl;

    // Create contact network graph
    std::cerr << "Creating network" << std::endl;
    auto nw = watts_strogatz(N, K0, P_REWIRE, engine);

    // Create transmission and reset time distributions
    std::cerr << "Creating simulator" << std::endl;
    const auto psi = transmission_time_gamma(Mi, Vi);
    const auto rho = transmission_time_lognormal(Mr, Vr);

    // Run simulation for M steps
    double total_runtime = 0.0;
    for(std::size_t i=0; i < M; ++i) {
        // Create simulation and specifiy initial set of infections
        auto sim = simulate_next_reaction(nw, psi, &rho, false, true, true);
        sim.add_infections({{0, 0.0}});

        const auto start = std::chrono::steady_clock::now();
        while (true) {
            auto r = sim.step(engine);
            if (!r)
                break;
        }
        const auto end = std::chrono::steady_clock::now();
        total_runtime += std::chrono::nanoseconds(end - start).count() / 1e9;
    }

    std::cerr << "Time per simultion run: " << (total_runtime / M) << " seconds " << std::endl;
    return 0;
}
