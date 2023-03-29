#include "tests/stdafx.h"
#include "tests/simulate.h"
#include "tests/analytical.h"
#include "tests/plot.h"

#include "random.h"
#include "NextReaction.h"
#include "dynamic_graph.h"
#include "algorithm.h"

namespace {

/**
 * @brief Simple symmetric Z-test (similar to a t-Test but for known variance)
 *
 * @return The symmetric p-value
 */
#if 0
inline double ztest(double mean_obs, double sd_true, double mean_true) {
	using namespace std;
	const double z = (mean_obs - mean_true) / sd_true;
	return 1 - std::erf(abs(z) / sqrt(2));
}
#endif

}

TEST_CASE("Plot SIS average trajectory on dynamic Erdös-Reyni networks", "[nextreaction]")
{
	using namespace std::string_literals;

	const std::size_t M = 100;
	const std::size_t N = 1000;
	const std::size_t K = 3;
	const double TAU = 10;
	const double PSI_MEAN = 3;
	const double PSI_VARIANCE = 1;
	const double RHO_MEAN = 10;
	const double RHO_VARIANCE = 1;
	const double TMAX = 100;

	rng_t engine;
	std::vector<double> t_sim, y_sim_new, y_sim_total;
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
	}, t_sim, y_sim_total, y_sim_new, TMAX, M);g
	
	std::vector<double> t_sim_static, y_sim_new_static, y_sim_total_static;
	average_trajectories(engine, [&](rng_t& engine){
		struct {
			std::unique_ptr<graph> g;
			std::unique_ptr<transmission_time_gamma> psi;
			std::unique_ptr<transmission_time_gamma> rho;
			std::unique_ptr<simulation_algorithm> simulator;
		} env;
		env.g = std::make_unique<erdos_reyni>(N, K, engine);
		env.psi = std::make_unique<transmission_time_gamma>(PSI_MEAN, PSI_VARIANCE);
		env.rho = std::make_unique<transmission_time_gamma>(RHO_MEAN, RHO_VARIANCE);
		env.simulator = std::make_unique<simulate_next_reaction>(*env.g.get(), *env.psi.get(), env.rho.get(), true, false);
		return env;
	}, t_sim_static, y_sim_total_static, y_sim_new_static, TMAX, M);

	
	plot("nextreaction_dynamic.sis.mean.pdf", "SIS average trajectory on dynamic Erdös-Reyni [NextReaction]", [&](auto& gp, auto& p) {
		p.add_plot1d(std::make_pair(t_sim, y_sim_new), "with lines title 'dynamic network'"s);
		p.add_plot1d(std::make_pair(t_sim_static, y_sim_new_static), "with lines title 'static network'"s);
   });
}
