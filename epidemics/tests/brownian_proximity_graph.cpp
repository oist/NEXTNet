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
	
	while(sim.next(engine) < 1000)
		sim.step(engine);
};
