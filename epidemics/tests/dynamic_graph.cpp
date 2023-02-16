#include "tests/stdafx.h"

#include "dynamic_graph.h"

int edge_index(node_t nodeA, node_t nodeB) {
	/* Translate (src, dst) pair of edge into an index 0 <= i < N*(N - 1)/2
	 * We first translate (src, dist) into (n1, n2) where n1 > n2. Since there
	 * are n1 valid edge descriptors (n1, x) for a node n1, there are
	 * 0 + 1 + ... + (n1 - 1) = n * (n1 - 1) / 2 edge descriptors (n1p, x)
	 * with n1p < n1. The mapping of (n1, n2) to (n1 * (n1 - 1) / 2) + n2
	 * is thus bijective onto [0, E) where E is the number of possible edges.
	 */
	const int n1 = std::max(nodeA, nodeB);
	const int n2 = std::min(nodeA, nodeB);
	REQUIRE(n1 > n2);
	const int e = (n1 * (n1 - 1) / 2) + n2;
	REQUIRE(e >= 0);
	return e;
}

/**
 * @brief Test case to verify `dynamic_erdos_reyni`
 */
TEST_CASE("dynamic Erdös-Reyni", "[dynamic_graph]") {
	const int N = 20;
	const int E = N*(N-1)/2;
	const double P = 0.1;
	const double TAU = 1.0;
	const double TMAX = 10000*TAU;
	
	/* Keep a list of pairs (t, f) for every edge which stores the times
	 * at which the edge appears (t, 1) or disappears (t, 0)
	 */
	std::array<std::vector<std::pair<double, bool>>, E> edge_presence;
	
	std::mt19937 engine;
	dynamic_erdos_reyni g(N, N*P, TAU, engine);
	
	/* Initial edge_presence */
	for(node_t a=0; a < N; ++a) {
		const int k = g.outdegree(a);
		for(int j=0; j < k; ++j) {
			const node_t b = g.neighbour(a, j);
			if (a > b)
				break;
			const int e = edge_index(a, b);
			REQUIRE(edge_presence[e].empty());
			edge_presence[e].emplace_back(0, true);
		}
	}
	for(int e=0; e < E; ++e) {
		if (edge_presence[e].empty())
			edge_presence[e].emplace_back(0, false);
	}
	
	/* Evolve dynamic network */
	while (auto ev = g.step(engine, TMAX)) {
		const network_event_t event = *ev;
		const int e = edge_index(event.source_node, event.target_node);
		
		switch (event.kind) {
			case network_event_kind::neighbour_added:
				edge_presence[e].emplace_back(event.time, true);
				break;
				
			case network_event_kind::neighbour_removed:
				edge_presence[e].emplace_back(event.time, false);
				break;

			default:
				break;
		}
	}
	
	/* Compute average presence and absence times for every edge */
	for(int e=0; e < E; ++e) {
		const auto& p = edge_presence[e];
		if (p.empty())
			continue;
		
		/* Compute present probability and times between events */
		double p_present = 0;
		std::vector<double> time_present;
		std::vector<double> time_absent;

		const auto end = p.end();
		auto curr = p.begin();
		auto prev = curr++;
		for(; curr != end; prev = curr++) {
			REQUIRE(prev->second != curr->second);
			if (curr->second) {
				/* appears */
				time_absent[curr->first - prev->first];
			} else {
				/* disappears */
				time_present[curr->first - prev->first];
				p_present += (curr->first - prev->first) / TMAX;
			}
		}
		if (!prev->second) {
			/* last event was disappear */
			time_absent[TMAX - prev->first];
		} else {
			/* last event was appear */
			time_present[TMAX - prev->first];
			p_present += (TMAX - prev->first) / TMAX;
		}

		
		/* TODO: Check that time_present, time_absent is exponentially distributed */
		std::cout << "p: " << p_present << std::endl;
	}
}
