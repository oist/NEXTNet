#include "tests/stdafx.h"
#include "tests/statistics.h"
#include "dynamic_graph.h"

/**
 * @brief Test case to verify `dynamic_empirical_network`
 */
TEST_CASE("dynamic empirical graph", "[empirical_graph]") {

	double dt = 4;
	dynamic_empirical_network g(TEST_DATA_DIR "/test_empirical_network.txt", dt);
	INFO("adjacencylist size: " << g.adjacencylist.size() << " (Expected: 13)");
	REQUIRE(g.adjacencylist.size()==14);
	INFO("Edges size: " << g.edges.size() << " (Expected: 4)");
	REQUIRE(g.edges.size()==8);

	REQUIRE(g.nodes()==14);
	REQUIRE(g.neighbour(0,0)==13);
	REQUIRE(g.outdegree(0)==2);
	REQUIRE(g.outdegree(10)==0);
}


/**
 * @brief Test case to verify `dynamic_erdos_reyni`
 */
TEST_CASE("dynamic Erdös-Reyni", "[dynamic_graph]") {
	const int N = 20;
	const int E = N*(N-1)/2;
	const double P = 0.1;
	const double TAU = 1;
	const double TMAX = 10000/TAU;
	
	/* Keep a list of pairs (t, f) for every edge which stores the times
	 * at which the edge appears (t, 1) or disappears (t, 0)
	 */
	std::array<std::vector<std::pair<double, bool>>, E> edge_presence;
	
	std::mt19937 engine;
	dynamic_erdos_reyni g(N, (N-1)*P, TAU, engine);
	
	/* Initial edge_presence */
	for(node_t a=0; a < N; ++a) {
		const int k = g.outdegree(a);
		for(int j=0; j < k; ++j) {
			const node_t b = g.neighbour(a, j);
			const int e = edge_index_undirected(a, b);
			if (a > b)
				continue;
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
		const int e = edge_index_undirected(event.source_node, event.target_node);
		
		switch (event.kind) {
			case network_event_kind::neighbour_added:
				if (!edge_presence[e].back().second)
					/* First of fwd/rev pair for this event, add it */
					edge_presence[e].emplace_back(event.time, true);
				else {
					/* Second of fwd/rev pair for this event, should be there already */
					REQUIRE(edge_presence[e].back().first == event.time);
					REQUIRE(edge_presence[e].back().second == true);
				}
				break;
				
			case network_event_kind::neighbour_removed:
				if (edge_presence[e].back().second)
					/* First of fwd/rev pair for this event, add it */
					edge_presence[e].emplace_back(event.time, false);
				else {
					/* Second of fwd/rev pair for this event, should be there already */
					REQUIRE(edge_presence[e].back().first == event.time);
					REQUIRE(edge_presence[e].back().second == false);
				}
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
		std::vector<double> times_present;
		std::vector<double> times_absent;
		double p_present = 0;

		const auto end = p.end();
		auto curr = p.begin();
		auto prev = curr++;
		for(; curr != end; prev = curr++) {
			const double dt = curr->first - prev->first;
			REQUIRE(prev->second != curr->second);
			if (curr->second) {
				/* appeared */
				times_absent.push_back(dt);
			}
			else{
				/* disappeared */
				times_present.push_back(dt);
				p_present += dt / TMAX;
			}
		}
		if (prev->second) {
			/* last event: appeared */
			p_present += (TMAX - prev->first) / TMAX;
		}
		
		const double time_present = std::accumulate(times_present.begin(), times_present.end(), 0.0) / times_present.size();
		const double time_absent = std::accumulate(times_absent.begin(), times_absent.end(), 0.0) / times_absent.size();

		/* Test that we see the exepcted number of events, mean present/absent times and presence probability,
		 *   No of events: Counting process with mean waiting time mu = 1/alpha + 1/beta,
		 *                 variance of waiting time sigma^2 = 1/alpha^2 + 1/beta^2. The mean and
		 *                 variance of the no. of events is then TMAX/mu and TMAX*sigma^2 / mu^3
		 *                 as TMAX goes to infinity (per counting process theory)
		 *   Avg(present): Average of exponentially distributed waiting times with mean 1/(TAU*(1-P))
		 *   Avg(absent) : Average of exponentially distributed waiting times with mean 1/(TAU*P)
		 * We use a simple z-test, i.e. assume normality and known variance of all tested quantities!
		 */
		/* TODO: This test fails for unknown reasons. Figure out why. */
		//const double events = times_present.size() + times_absent.size();
		//const double alpha = P/TAU, beta = 1/TAU - P/TAU;
		//const double mu = 1/alpha + 1/beta;
		//const double sigma2 = 1/pow(alpha, 2) + 1/pow(beta, 2);
		//const double pval_events = ztest(events, sqrt(TMAX*sigma2/pow(mu, 3)), TMAX/mu);
		//REQUIRE(pval_events >= 0.001);

		const double pval_tpresent = ztest(time_present, 1 / (TAU*(1-P)*sqrt(times_present.size())), 1.0 / (TAU*(1-P)));
		REQUIRE(pval_tpresent >= 0.01 / E);
		const double pval_tabsent = ztest(time_absent, 1 / (TAU*P*sqrt(times_absent.size())), 1.0 / (TAU*P));
		REQUIRE(pval_tabsent >= 0.01 / E);
		/* TODO: Figure out the correct mean and variance of p_present */
		//const double pval_ppresent = ztest(p_present, ?, P);
		//REQUIRE(pval_ppresent > 0.01);

		/*
		std::cout << "events: " << events << " (pval " << std::setprecision(3) << pval_events << "), ";
		std::cout << "Avg(present): " << std::setprecision(3) << time_present << " (pval " << std::setprecision(3) << pval_tpresent << "), ";
		std::cout << "Avg(absent): " << std::setprecision(3) << time_absent << " (pval " << std::setprecision(3) << pval_tabsent << "), ";
		std::cout << "P(present): " << std::setprecision(3) << p_present << " (pval " << std::setprecision(3) << pval_ppresent << ")" << std::endl;
		*/
	}
	
	std::cout << "---" << std::endl;
}
