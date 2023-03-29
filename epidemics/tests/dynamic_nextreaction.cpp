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

TEST_CASE("Dynamics on a dynamic Erdös-Reyni network", "[nextreaction]")
{
	const std::size_t N = 1000;
	const std::size_t K = 3;
	const double TAU = 5.0;
	const double PSI_MEAN = 3;
	const double PSI_VARIANCE = 1;
	const double RHO_MEAN = 10;
	const double RHO_VARIANCE = 1;
	const double TMAX = 100;

	rng_t engine;
	dynamic_erdos_reyni g(N, K, TAU, engine);
	transmission_time_gamma psi(PSI_MEAN, PSI_VARIANCE);
	transmission_time_gamma rho(RHO_MEAN, RHO_VARIANCE);
	simulate_next_reaction sim_nr(g, psi, &rho);
	simulate_on_dynamic_network sim(sim_nr);
	
	while (true) {
		/* Determine time of next event */
		const double tnext = sim.next(engine);
		if (tnext > TMAX)
			break;
		
		/* Perform step. Since we're doing SIS, there's always a next event */
		const auto any_ev = *sim.step(engine);
		if (std::holds_alternative<event_t>(any_ev)) {
			/* Epidemic event */
			const auto ev = std::get<event_t>(any_ev);
		} else if (std::holds_alternative<network_event_t>(any_ev)) {
			/* Network event */
			const auto ev = std::get<network_event_t>(any_ev);
		}
	}
}
