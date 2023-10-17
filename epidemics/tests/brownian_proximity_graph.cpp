#include "tests/stdafx.h"
#include "tests/simulate.h"
#include "tests/plot.h"

#include "brownian_proximity_graph.h"
#include "NextReaction.h"
#include "algorithm.h"

#if 0
TEST_CASE("Epidemic on Brownian proximity graph", "[brownian_proximity_graph]") {
	rng_t engine;

	brownian_proximity_graph g(100, 3.0, 1, 0.1, engine);
	transmission_time_gamma psi(5, 3);
	transmission_time_gamma rho(100, 10);
	simulate_next_reaction nr(g, psi, &rho, false, true);
	simulate_on_dynamic_network sim(nr);
	
	nr.add_infections({ std::make_pair(0, 0.0)} );
	
	while(sim.next(engine) < 1000) {
		const auto ev = *sim.step(engine);
		if (std::holds_alternative<network_event_t>(ev)) {
			const auto nw_ev = std::get<network_event_t>(ev);
			std::cerr << "Network event " << name(nw_ev.kind) << ": "
				<< "time=" << nw_ev.time << ", "
				<< "src=" << nw_ev.source_node <<  ", "
				<< "tgt=" << nw_ev.target_node << std::endl;
		} else if (std::holds_alternative<event_t>(ev)) {
			const auto ep_ev = std::get<event_t>(ev);
			std::cerr << "Epidemic event " << name(ep_ev.kind) << ": "
				<< "time=" << ep_ev.time << ", "
				<< "src=" << ep_ev.source_node << ", "
				<< "tgt=" << ep_ev.node << std::endl;
		}
	}
};
#endif

TEST_CASE("Plot SIS average trajectory on Brownian proximity graph", "[brownian_proximity_graph]")
{
	using namespace std::string_literals;

	const std::size_t M = 1;
	const std::size_t N = 100;
	const double K = 3.0;
	const double RADIUS = 1.0;
	const double D = 0.1;
	const double TAU = 10;
	const double PSI_MEAN = 3;
	const double PSI_VARIANCE = 1;
	const double RHO_MEAN = 10;
	const double RHO_VARIANCE = 1;
	const double TMAX = 500;

	rng_t engine;
	std::vector<double> t_sim, y_sim_new, y_sim_total;
	average_trajectories(engine, [&](rng_t& engine) {
		struct {
			std::unique_ptr<brownian_proximity_graph> g;
			std::unique_ptr<transmission_time_gamma> psi;
			std::unique_ptr<transmission_time_gamma> rho;
			std::unique_ptr<simulate_next_reaction> nr;
			std::unique_ptr<simulate_on_dynamic_network> simulator;
		} env;
		env.g = std::make_unique<brownian_proximity_graph>(N, K, RADIUS, D, engine);
		env.psi = std::make_unique<transmission_time_gamma>(PSI_MEAN, PSI_VARIANCE);
		env.rho = std::make_unique<transmission_time_gamma>(RHO_MEAN, RHO_VARIANCE);
		env.nr = std::make_unique<simulate_next_reaction>(*env.g.get(), *env.psi.get(), env.rho.get(), false, true);
		env.nr->add_infections({ std::make_pair(0, 0.0)});
		env.simulator = std::make_unique<simulate_on_dynamic_network>(*env.nr.get());
		return env;
	}, [](network_or_epidemic_event_t any_ev) {
		/* Translate event into a pair (time, delta) */
		if (std::holds_alternative<event_t>(any_ev)) {
			/* Epidemic event */
			const auto ev = std::get<event_t>(any_ev);
			return std::make_pair(ev.time, delta_infected(ev.kind));
		} else if (std::holds_alternative<network_event_t>(any_ev)) {
			/* Network event */
			const auto ev = std::get<network_event_t>(any_ev);
			return std::make_pair(ev.time, 0);
		} else throw std::logic_error("unknown event type");
	}, t_sim, y_sim_total, y_sim_new, TMAX, M);
	
	std::vector<double> t_sim_erdos, y_sim_erdos_new, y_sim_erdos_total;
	average_trajectories(engine, [&](rng_t& engine) {
		struct {
			std::unique_ptr<dynamic_erdos_reyni> g;
			std::unique_ptr<transmission_time_gamma> psi;
			std::unique_ptr<transmission_time_gamma> rho;
			std::unique_ptr<simulate_next_reaction> nr;
			std::unique_ptr<simulate_on_dynamic_network> simulator;
		} env;
		env.g = std::make_unique<dynamic_erdos_reyni>(N, K, TAU, engine);
		env.psi = std::make_unique<transmission_time_gamma>(PSI_MEAN, PSI_VARIANCE);
		env.rho = std::make_unique<transmission_time_gamma>(RHO_MEAN, RHO_VARIANCE);
		env.nr = std::make_unique<simulate_next_reaction>(*env.g.get(), *env.psi.get(), env.rho.get(), false, true);
		env.nr->add_infections({ std::make_pair(0, 0.0)});
		env.simulator = std::make_unique<simulate_on_dynamic_network>(*env.nr.get());
		return env;
	}, [](network_or_epidemic_event_t any_ev) {
		/* Translate event into a pair (time, delta) */
		if (std::holds_alternative<event_t>(any_ev)) {
			/* Epidemic event */
			const auto ev = std::get<event_t>(any_ev);
			return std::make_pair(ev.time, delta_infected(ev.kind));
		} else if (std::holds_alternative<network_event_t>(any_ev)) {
			/* Network event */
			const auto ev = std::get<network_event_t>(any_ev);
			return std::make_pair(ev.time, 0);
		} else throw std::logic_error("unknown event type");
	}, t_sim_erdos, y_sim_erdos_total, y_sim_erdos_new, TMAX, M);

	
	plot("nextreaction_brownian.sis.mean.pdf", "SIS average trajectory on Brownian proximity graph", [&](auto& gp, auto& p) {
		p.add_plot1d(std::make_pair(t_sim, y_sim_new), "with lines title 'Brownian proximity'"s);
		p.add_plot1d(std::make_pair(t_sim_erdos, y_sim_erdos_new), "with lines title 'Dynamic Erdos'"s);
   });
}