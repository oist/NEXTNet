#include "random.h"
#include "analysis.h"
#include "graph.h"
#include "NextReaction.h"

using namespace std::string_literals;

int program_benchmark_priorityqueue(int argc, const char * argv[])
{
    rng_t engine;

    using namespace std::chrono;
    typedef std::chrono::duration<double, std::micro> usecs;
    typedef std::chrono::duration<double, std::ratio<1>> secs;

    // Parameters, which are network size (N), average number of contacts (K),
    // mean infection transmission time (Mi) and variance (Vi),
    // mean reset time (Mr) and variance (Vr), number of steps to run (M)
    const int BATCH = 10000;
    const int LOG = 100000;
    const int N = (argc >= 3) ? atoi(argv[2]) : 1000000;
    const std::size_t M = (argc >= 4) ? atoi(argv[3]) : 1000000;
    const bool edges_concurrent = (argc >= 5) ? ("true"s == argv[4]) : false;
    const std::optional<std::string> OUTFILE = (argc >= 6) ? std::optional<std::string>(argv[5]) : std::nullopt;
    const int K = 5;
    const double Mi = 3;
    const double Vi = 30;
    const double Mr = 50;
    const double Vr = 20;
    std::cerr << "Running with parameters:" << std::endl;
    std::cerr << "            N: " << N << " (network size)" << std::endl;
    std::cerr << "            K: " << K << " (average degree)" << std::endl;
    std::cerr << "        Mi/Vi: " << Mi << " / " << Vi << " (Infection time mean / variance)" << std::endl;
    std::cerr << "        Mr/Vr: " << Mr << " / " << Vr << " (Reset time mean / variance)" << std::endl;
    std::cerr << "            M: " << M << " (number of steps)" << std::endl;
    std::cerr << "  Conc. edges: " << (edges_concurrent ? "YES"s : "NO"s) << std::endl;
    std::cerr << "Writing results to " << (OUTFILE ? *OUTFILE : "stdout"s) << std::endl;

    // Create contact network graph
    std::cerr << "Creating network" << std::endl;
    auto nw = erdos_reyni(N, K, engine);

    // Create transmission and reset time distributions
    std::cerr << "Creating simulator" << std::endl;
    const auto psi = transmission_time_lognormal(Mi, Vi);
    const auto rho = transmission_time_lognormal(Mr, Vr);

    // Create simulation and specifiy initial set of infections
    auto sim = simulate_next_reaction(nw, psi, &rho, true, edges_concurrent);
    sim.add_infections({{0, 0.0}});

    // Prepare output
    std::ofstream file;
    std::ostream* out;
    if (OUTFILE) {
        file.open(*OUTFILE, std::ios_base::trunc | std::ios_base::out);
        out = &file;
    } else {
        out = &std::cout;
    }
    *out
        << "step\t"
        << "queuestep\t"
        << "runtime\t"
        << "simtime\t"
        << "infected\t"
        << "active_edges\t"
        << "us_per_step\t"
        << "us_per_qstep\t"
        << "\n";
    const auto start = std::chrono::high_resolution_clock::now();
    auto last = start;
    absolutetime_t simtime = NAN;
    std::size_t i = 0;
    std::size_t last_i = 0;
    std::size_t last_qstep = 0;
    auto next_batch = [&i, &last_i, &last_qstep, &start, &last, &sim, &simtime, &out]() {
        const auto cur = std::chrono::high_resolution_clock::now();
        const auto ttime = cur - start;
        const auto btime = cur - last;
        const double btime_us = duration_cast<usecs>(btime).count();
        const double bsteps = i - last_i;
        const std::size_t qsteps = sim.queue_steps_total - last_qstep;
        last = cur;
        last_i = i;
        last_qstep = sim.queue_steps_total;

        if (i == 0)
            return;

        *out
            << i << "\t"
            << sim.queue_steps_total << "\t"
            << duration_cast<secs>(ttime).count() << "\t"
            << simtime << "\t"
            << sim.infected.size() << "\t"
            << sim.active_edges.size() << "\t"
            << btime_us / bsteps << "\t"
            << btime_us / qsteps
            << "\n";
    };

    // Run simulation for M steps
    for(; i < M; ++i) {
        if (i % BATCH == 0)
            next_batch();

        if ((i+1) % LOG == 0)
            std::cerr << "  Step " << (i+1) << " / " << M << std::endl;

        const auto ev = sim.step(engine);
        if (ev)
            simtime = ev->time;
    }
    next_batch();

    out->flush();
    file.close();
    return 0;
}
