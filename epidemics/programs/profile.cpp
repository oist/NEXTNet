#include "random.h"
#include "analysis.h"
#include "simulation.h"
#include "graph.h"
#include "NextReaction.h"

using namespace std::string_literals;

int program_profile(int argc, const char * argv[])
{
    rng_t engine;

    // Parameters, which are network size (N), average number of contacts (K), 
    // mean infection transmission time (Mi) and variance (Vi),
    // mean reset time (Mr) and variance (Vr), number of steps to run (M)
    const int N = (argc >= 3) ? atoi(argv[2]) : 1000000;
    const std::size_t M = (argc >= 4) ? atoi(argv[3]) : 1000000;
    const int K = 5;
    const double Mi = 3;
    const double Vi = 30;
    const double Mr = 50;
    const double Vr = 20;
    std::cerr << "Running with parameters:" << std::endl;
    std::cerr << "      N: " << N << " (network size)" << std::endl;
    std::cerr << "      K: " << N << " (average degree)" << std::endl;
    std::cerr << "  Mi/Vi: " << Mi << " / " << Vi << " (Infection time mean / variance)" << std::endl;
    std::cerr << "  Mr/Vr: " << Mr << " / " << Vr << " (Reset time mean / variance)" << std::endl;
    std::cerr << "      M: " << M << " (number of steps)" << std::endl;

    // Create contact network graph
    auto nw = erdos_reyni(N, K, engine);

    // Create transmission and reset time distributions
    const auto psi = transmission_time_lognormal(Mi, Vi);
    const auto rho = transmission_time_lognormal(Mr, Vr);

    // Create simulation and specifiy initial set of infections
    auto sim = simulate_next_reaction(nw, psi, &rho);
    sim.add_infections({{0, 0.0}});

    // Run simulation for M steps
    for(std::size_t i=0; i < M; ++i)
        sim.step(engine);
    
    return 0;
}
