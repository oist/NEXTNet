#include "tests/stdafx.h"
#include "tests/simulate.h"
#include "tests/analytical.h"
#include "tests/plot.h"

#include "random.h"
#include "nMGA.h"

TEST_CASE("Stability", "[nMGA]") {
	using namespace std;
	
	const int N = 530;
	const double MEAN_INFECTION = 10;
	const double VARIANCE_INFECTION = 1.0;
	const double MEAN_RECOVERY = 20;
	const double VARIANCE_RECOVERY = 1;
	const double R0 = 3;

	transmission_time_lognormal psi(MEAN_INFECTION, VARIANCE_INFECTION);
	transmission_time_lognormal rho(MEAN_RECOVERY, VARIANCE_RECOVERY);

	rng_t engine;
	erdos_reyni network(N, R0,engine);
	simulate_nmga simulate(network, psi, nullptr, false, N/10);
	const int N0 = floor(network.nodes() / 10); // 10% of the network will be infected at t=0;
	for (node_t node = 0; node < N0; node++)
		simulate.add_infections({ std::make_pair(node, 0.0)});
			
	while (true) {
		auto point = simulate.step(engine);
		if (!point)
			break;
	}
}

#if ENABLE_PLOTTING
TEST_CASE("Plot large-population SIR mean-field (nMGA)", "[nMGA]") {
    using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t M = 1000;
    const std::size_t T = 25;
    const std::size_t X = 400;
    const double R0 = 3;
    const double MEAN = 10;
    const double VARIANCE = 1;
    transmission_time_gamma psi(MEAN, VARIANCE);

    /* Simulate using nMGA reaction M times */
    auto racyclic = simulate_SIR<acyclic, simulate_nmga>(engine, psi, T, M, 1, R0+1, true);

    /* Evaluate analytical solution */
    std::pair<std::vector<double>, std::vector<double>> analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for(std::size_t i=0; i < X; ++i) {
        const double t = (double)T * i / (X-1);
        analytical.first.push_back(t);
        analytical.second.push_back(sol.N(t));
    }

    plot("nmga.sir.mean.pdf", "Large-population SIR mean-field [nMGA]", [&](auto& gp, auto& p) {
        p.add_plot1d(racyclic, "with lines title 'nMGA acyclic'"s);
        p.add_plot1d(analytical, "with lines title 'analytical'");
    });
}
#endif

#if ENABLE_PLOTTING
TEST_CASE("Plot SIS single trajectory (nMGA)", "[nMGA]") {
	using namespace std::string_literals;
	std::mt19937 engine;

	const int N = 10000;
	const double T = 100;
	const double R0 = 3;

	const double MEAN = 3;
	const double VARIANCE = 1;
	const double MEAN_rho = 10;
	const double VARIANCE_rho = 1;
	transmission_time_gamma psi(MEAN, VARIANCE);
	transmission_time_gamma rho(MEAN_rho, VARIANCE_rho);

	/* Simulate using nMGA once */
	std::vector<double> t_sim, y_sim_new, y_sim_total;
	simulate_SIS<erdos_reyni, simulate_nmga>(engine, psi,rho, t_sim, y_sim_new, y_sim_total, T, N, R0);

	plot("nmga.sis.single.pdf", "SIS single trajectory [nMGA]", [&](auto& gp, auto& p) {
		p.add_plot1d(std::make_pair(t_sim, y_sim_total), "with lines title 'nMGA'"s);
	});
}
#endif
