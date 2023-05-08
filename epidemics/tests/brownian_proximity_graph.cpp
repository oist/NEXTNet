#include "tests/stdafx.h"

#include "brownian_proximity_graph.h"
#include "NextReaction.h"
#include "algorithm.h"

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
