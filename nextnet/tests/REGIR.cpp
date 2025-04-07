#include "nextnet/tests/stdafx.h"
#include "nextnet/tests/simulate.h"
#include "nextnet/tests/analytical.h"
#include "nextnet/tests/plot.h"

#include "nextnet/random.h"
#include "nextnet/REGIR.h"
#include "nextnet/NextReaction.h"

TEST_CASE("Stability (REGIR)", "[REGIR]")
{
    using namespace std;

    const int N                     = 530;
    const double MEAN_INFECTION     = 10;
    const double VARIANCE_INFECTION = 1.0;
    const double R0                 = 3;

    transmission_time_gamma psi(MEAN_INFECTION, VARIANCE_INFECTION);

    rng_t engine;
    erdos_reyni network(N, R0, engine);
    simulate_regir::params p;
    p.approximation_threshold = N / 10;
    simulate_regir simulate(network, psi, nullptr, p);
    const int N0 = floor(network.nodes() / 10); // 10% of the network will be infected at t=0;
    for (node_t node = 0; node < N0; node++)
        simulate.add_infections({ std::make_pair(node, 0.0) });

    while (true) {
        auto point = simulate.step(engine);
        if (!point)
            break;
    }
}

#if ENABLE_PLOTTING
TEST_CASE("Plot large-population SIR mean-field (REGIR)", "[REGIR]")
{
    using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t M   = 100;
    const std::size_t T   = 25;
    const std::size_t X   = 400;
    const double R0       = 3;
    const double MEAN     = 10;
    const double VARIANCE = 1;
    transmission_time_gamma psi(MEAN, VARIANCE);

    /* Simulate using REGIR reaction M times */
    simulate_regir::params p;
    auto racyclic = simulate_trajectory<acyclic, simulate_regir>(engine, psi, p, T, M, 1, R0 + 1, true);

    /* Evaluate analytical solution */
    std::pair<std::vector<double>, std::vector<double>> analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for (std::size_t i = 0; i < X; ++i) {
        const double t = (double)T * i / (X - 1);
        analytical.first.push_back(t);
        analytical.second.push_back(sol.N(t));
    }

    plot("regir.sir.mean.pdf", "Large-population SIR mean-field [REGIR]", [&](auto &gp, auto &p) {
        p.add_plot1d(racyclic, "with lines title 'REGIR acyclic'"s);
        p.add_plot1d(analytical, "with lines title 'analytical'");
    });
}
#endif

#if ENABLE_PLOTTING
TEST_CASE("Plot SIS average trajectory (REGIR)", "[REGIR]")
{
    using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t M = 20;
    const int N         = 100;
    const double T      = 15;
    const double R0     = 3;

    const double MEAN         = 2;
    const double VARIANCE     = 2;
    const double MEAN_rho     = 5;
    const double VARIANCE_rho = 3;
    transmission_time_gamma psi(MEAN, VARIANCE);
    transmission_time_gamma rho(MEAN_rho, VARIANCE_rho);

    simulate_regir::params p_regir;
    p_regir.approximation_threshold = 100;

    /* Simulate using REGIR once */
    std::vector<double> t_sim, y_sim_new, y_sim_total;
    average_trajectories(engine, [&](rng_t &engine) {
		struct {
			std::unique_ptr<network> nw;
			std::unique_ptr<simulation_algorithm> simulator;
		} env;
		env.nw.reset(new erdos_reyni(N, R0, engine));
		env.simulator.reset(new simulate_regir(*env.nw, psi, &rho, p_regir));
		return env; }, t_sim, y_sim_new, y_sim_total, T, M);

    /* Simulate using NextReaction once */
    std::vector<double> t_sim_nr, y_sim_new_nr, y_sim_total_nr;
    average_trajectories(engine, [&](rng_t &engine) {
		struct {
			std::unique_ptr<network> nw;
			std::unique_ptr<simulation_algorithm> simulator;
		} env;
		env.nw.reset(new erdos_reyni(N, R0, engine));
		env.simulator.reset(new simulate_next_reaction(*env.nw, psi, &rho));
		return env; }, t_sim_nr, y_sim_new_nr, y_sim_total_nr, T, M);

    plot("regir.sis.average.pdf", "SIS average trajectory [REGIR]", [&](auto &gp, auto &p) {
        p.add_plot1d(std::make_pair(t_sim, y_sim_total), "with lines title 'REGIR'"s);
        p.add_plot1d(std::make_pair(t_sim_nr, y_sim_total_nr), "with lines title 'NextReaction'"s);
    });
}
#endif

TEST_CASE("SIR REGIR on ER Graph", "[REGIR]")
{

    // The SIR can be mapped to a percolation process, which allows to derive some precise results:

    // we call p the probability that an infected nodes transmits the disease to a given neighbour before recovery;
    // At the end of the epidemic (once there are 0 infected in the network, only recovered or susceptibles), the fraction
    // of recovered is equivalent to the size of the Giant component under a percolation process where we remove each edge with probability p.
    // For networks with a poisson distribution, it is possible to have an exact expression of the expected fraction of recovered.

    rng_t engine;

    double MEAN_INFECTION     = 10;
    double VARIANCE_INFECTION = 5;
    double MEAN_RECOVERY      = 14;
    double VARIANCE_RECOVERY  = 7;

    double R0 = 3.0;

    int size = 1000;

    // // for gamma distributions, with the parameters given above, we can calculate p.
    // double p = 0.87852;
    // Pr r that a randomly chosen node belongs to the giant component is:
    double Pg = 0.855518;

    transmission_time_gamma psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_gamma rho(MEAN_RECOVERY, VARIANCE_RECOVERY);

    erdos_reyni network(size, R0, engine);
    simulate_regir::params p;
    p.approximation_threshold = 100;
    p.tau_precision           = 1e-6;
    p.SIR                     = true;
    simulate_regir simulation(network, psi, &rho, p);
    std::uniform_int_distribution<> dis(0, size - 1);

    // initial infected
    for (node_t i = 0; i < 10; i++) {
        node_t rand_node = dis(engine);
        simulation.add_infections({ std::make_pair(rand_node, 0.0) });
    }

    int nb_recovered = 0;
    while (true) {
        auto point = simulation.step(engine);
        if (!point)
            break;

        if (point->kind == epidemic_event_kind::reset)
            nb_recovered++;
    }

    REQUIRE(nb_recovered < size);
    REQUIRE(nb_recovered > 0);
    REQUIRE(std::abs(nb_recovered - Pg * size) / (Pg * size) < 0.1);
}
