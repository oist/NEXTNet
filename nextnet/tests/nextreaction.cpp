#include "nextnet/tests/stdafx.h"
#include "nextnet/tests/simulate.h"
#include "nextnet/tests/analytical.h"
#include "nextnet/tests/statistics.h"
#include "nextnet/tests/plot.h"

#include "nextnet/random.h"
#include "nextnet/NextReaction.h"

#if ENABLE_PLOTTING
TEST_CASE("Plot large-population SIR mean-field (NextReaction)", "[nextreaction]")
{
    using namespace std::string_literals;
    std::mt19937 engine;

    const int M      = 1000;
    const int Nfully = 1000;
    const int Nerdos = 1000;
    const double T   = 25;
    const int X      = 400;
    const double R0  = 3;
    // Every infecteced node has N-1 neighbours, of which in the large-population limit N-2 are susceptible.
    // To trigger subsequent infections amgonst these N-2 susceptible neighbours, we must infect each neighbour
    // with probability R0/(N-2). Note that for the first infected node, this is not strictly speaking correct,
    // since it's infection has no source, and it will thus create R0(N-1)/(N-2) > R0 subsequent infections.
    // This only causes a relative error of |1 - (N-1)/(N-2)| =~= 1/N though.
    const double pfully   = R0 / (Nfully - 2);
    const double MEAN     = 10;
    const double VARIANCE = 1;
    transmission_time_gamma psifully(MEAN, VARIANCE, 1.0 - pfully);
    transmission_time_gamma psi(MEAN, VARIANCE);

    /* Simulate using next reaction M times */
    simulate_next_reaction::params p;
    p.edges_concurrent = false;
    auto rfully        = simulate_trajectory<fully_connected, simulate_next_reaction>(engine, psifully, p, T, M, 1, Nfully);
    auto racyclic      = simulate_trajectory<acyclic, simulate_next_reaction>(engine, psi, p, T, M, 1, R0 + 1, true);
    auto rerdos        = simulate_trajectory<erdos_reyni, simulate_next_reaction>(engine, psi, p, T, M, 1, Nerdos, R0);

    /* Evaluate analytical solution */
    std::pair<std::vector<double>, std::vector<double>> analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for (std::size_t i = 0; i < X; ++i) {
        const double t = (double)T * i / (X - 1);
        analytical.first.push_back(t);
        analytical.second.push_back(sol.N(t));
    }

    plot("nextreaction.sir.mean.pdf", "Large-population SIR mean-field [NextReaction]", [&](auto &gp, auto &p) {
        p.add_plot1d(rfully, "with lines title 'next reaction fully-connected (N="s + std::to_string(Nfully) + ")'"s);
        p.add_plot1d(racyclic, "with lines title 'next reaction acyclic'"s);
        p.add_plot1d(rerdos, "with lines title 'next reaction Erdös-Reyni (N="s + std::to_string(Nerdos) + ")'"s);
        p.add_plot1d(analytical, "with lines title 'analytical'");
    });
}
#endif

#if ENABLE_PLOTTING
TEST_CASE("Plot SIS average trajectory (NextReaction)", "[nextreaction]")
{
    using namespace std::string_literals;
    std::mt19937 engine;

    const int M     = 50;
    const int N     = 500;
    const double T  = 75;
    const double R0 = 3;

    const double MEAN         = 3;
    const double VARIANCE     = 1;
    const double MEAN_rho     = 10;
    const double VARIANCE_rho = 1;
    transmission_time_gamma psi(MEAN, VARIANCE);
    transmission_time_gamma rho(MEAN_rho, VARIANCE_rho);

    /* Simulate using next reaction once, using sequential activation of edges */
    std::vector<double> t_sim_seq, y_sim_new_seq, y_sim_total_seq;
    average_trajectories(engine, [&](rng_t &engine) {
		struct {
			std::unique_ptr<network> nw;
			std::unique_ptr<simulation_algorithm> simulator;
		} env;
		env.nw.reset(new erdos_reyni(N, R0, engine));
		simulate_next_reaction::params p_serial;
		p_serial.edges_concurrent = false;
		env.simulator.reset(new simulate_next_reaction(*env.nw, psi, &rho, p_serial));
		return env; }, t_sim_seq, y_sim_new_seq, y_sim_total_seq, T, M);

    /* Simulate using next reaction once times, using concurrent activation of edges */
    std::vector<double> t_sim_conc, y_sim_new_conc, y_sim_total_conc;
    average_trajectories(engine, [&](rng_t &engine) {
		struct {
			std::unique_ptr<network> nw;
			std::unique_ptr<simulation_algorithm> simulator;
		} env;
		env.nw.reset(new erdos_reyni(N, R0, engine));
		simulate_next_reaction::params p_concurrent;
		p_concurrent.edges_concurrent = true;
		env.simulator.reset(new simulate_next_reaction(*env.nw, psi, &rho, p_concurrent));
		return env; }, t_sim_conc, y_sim_new_conc, y_sim_total_conc, T, M);

    plot("nextreaction.sis.mean.pdf", "SIS average trajectory [NextReaction]", [&](auto &gp, auto &p) {
        p.add_plot1d(std::make_pair(t_sim_seq, y_sim_total_seq), "with lines title 'next reaction (seq. edges)'"s);
        p.add_plot1d(std::make_pair(t_sim_conc, y_sim_total_conc), "with lines title 'next reaction (conc. edges)'"s);
    });
}
#endif

TEST_CASE("Fraction of the recovered nodes for SIR on the Erdos-Renyi graph", "[nextreaction]")
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

    int size = 100000;

    // // for gamma distributions, with the parameters given above, we can calculate p.
    // double p = 0.87852;
    // Pr r that a randomly chosen node belongs to the giant component is:
    double Pg = 0.855518;

    transmission_time_gamma psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_gamma rho(MEAN_RECOVERY, VARIANCE_RECOVERY);

    erdos_reyni network(size, R0, engine);
    simulate_next_reaction::params p;
    p.SIR = true;
    simulate_next_reaction simulation(network, psi, &rho, p);
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

TEST_CASE("Fraction of the recovered nodes for SIR on the weighted Erdos-Renyi graph", "[nextreaction]")
{

    // The SIR can be mapped to a percolation process, which allows to derive some precise results:

    // we call p the probability that an infected nodes transmits the disease to a given neighbour before recovery;
    // At the end of the epidemic (once there are 0 infected in the network, only recovered or susceptibles), the fraction
    // of recovered is equivalent to the size of the Giant component under a percolation process where we remove each edge with probability p.
    // For networks with a poisson distribution, it is possible to have an exact expression of the expected fraction of recovered.

    rng_t engine;

    double MEAN_INFECTION                    = 10;
    double VARIANCE_INFECTION                = 5;
    double MEAN_RECOVERY                     = 14;
    double VARIANCE_RECOVERY                 = 7;
    std::vector<double> WEIGHT_CLASSES       = { 1.0, 5.0 };
    std::vector<double> WEIGHT_PROBABILITIES = { 2.0 / 3.0, 1.0 / 3.0 };

    double R0 = 3.0;

    int size = 100000;

    transmission_time_gamma psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_gamma rho(MEAN_RECOVERY, VARIANCE_RECOVERY);

    weighted_erdos_reyni network(size, R0, WEIGHT_CLASSES, WEIGHT_PROBABILITIES, engine);
    simulate_next_reaction::params p;
    p.SIR = true;
    simulate_next_reaction simulation(network, psi, &rho, p);
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

    /* TODO: Check results */
}

TEST_CASE("Exact reinfection mode", "[nextreaction]")
{
	rng_t engine;
	
	/* Create a two-node network: A <--> B */
	adjacencylist_network nw({{1}, {0}}, true, true);
	transmission_time_periodic psi(1.0/128);
	transmission_time_gamma rho(2, 1);
	const int M = 100;
	const int K = 1;
	
	std::vector<std::size_t> B_infections;
	B_infections.resize(K+2);
	for (int i=0; i < M; ++i) {
		/* Infect node A at time 0 and run simulation until A recovers */
		simulate_next_reaction::params p;
		p.exact_reinfection = true;
		simulate_next_reaction sim(nw, psi, &rho, p);
		sim.add_infections({ std::make_pair(0, 0.0) });
		int B_inf = 0;
		while(true) {
			epidemic_event_t next = *sim.step(engine);
			if ((next.kind == epidemic_event_kind::reset) && (next.node == 0))
				break;
			/* Count number of infections of B */
			if ((next.kind == epidemic_event_kind::infection) && (next.node == 1))
				B_inf++;
			if (B_inf > K)
				break;
		}
		B_infections[std::min(B_inf, K+1)]++;
	}
	
	for(std::size_t i=0; i < K+1; ++i)
		std::cerr << "  " << i << ": " << B_infections[i] << std::endl;
	const double pval1 = ztest(B_infections[0], 1, 0);
	REQUIRE(pval1 >= 0.01);
	const double pval2 = ztest(B_infections[1], sqrt(M)/2.0, M/2.0);
	REQUIRE(pval2 >= 0.01);
}
