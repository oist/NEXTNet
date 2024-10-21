//
//  dynamic_graph.hpp
//  epidemics
//
//  Created by Florian G. Pflug on 14.02.23.
//

#pragma once

#include "stdafx.h"
#include "types.h"
#include "random.h"
#include "graph.h"

//--------------------------------------
//----------DYNAMIC NETWORKS------------
//--------------------------------------

struct dynamic_network : virtual graph {
	virtual absolutetime_t next(rng_t& engine) = 0;

	virtual std::optional<network_event_t> step(rng_t& engine, absolutetime_t max_time = NAN) = 0;

	virtual void notify_epidemic_event(event_t ev, rng_t& engine);
};

struct graph_mutable : virtual graph
{
	void resize(node_t nodes);

	bool has_edge(node_t src, node_t dst);

	bool add_edge(node_t src, node_t dst);

	bool remove_edge(node_t src, node_t dst);

	virtual node_t nodes();

	virtual node_t neighbour(node_t node, int neighbour_index);

	virtual index_t outdegree(node_t node);

private:
	std::vector<drawable_set<node_t>> adjacencylist;
};

struct dynamic_empirical_network : virtual dynamic_network, virtual graph_mutable
{
	enum edge_duration_kind {
		finite_duration = 1,
		infitesimal_duration = 2
	};

	dynamic_empirical_network(std::string path_to_file, edge_duration_kind contact_type, double dt);

	virtual absolutetime_t next(rng_t& engine);

	virtual std::optional<network_event_t> step(rng_t& engine, absolutetime_t max_time = NAN);

private:
	std::deque<network_event_t> event_queue;

};

// struct activity_driven_network : virtual dynamic_network {

// 	activity_driven_network(std::vector<int> degreelist, double eta, double m, transmission_time& psi, transmission_time* rho, rng_t& engine);

// 	virtual absolutetime_t next(rng_t& engine);

// 	virtual std::optional<network_event_t> step(rng_t& engine, absolutetime_t max_time = NAN);

// 	void add_edge(node_t node, node_t neighbour);
// 	void remove_edge(node_t node, int neighbour_index);

// 	/* Current time (time of last event) */
// 	absolutetime_t current_time = 0.0;
	
// 	/* Time of next edge appearing or vanishing */
// 	absolutetime_t next_time = NAN;
	
// 	/* Unreported event for the reverse edge of the last edge event reported  */
// 	std::optional<network_event_t> reverse_edge_event;

// 	unsigned int edges_absent;
// 	unsigned int edges_present;
	
// 	/* Degree-weighted node distribution */
// 	dyndist::vector_distribution<unsigned> weighted_nodes;
// };

/**
 * @brief Dynamic version of an Erdös-Reyni graph
 *
 * Each edge appears and disappears independently according to a two-state
 * Markov process with rate alpha for an edge appearing and rate beta for
 * the rate disappearing. The dynamic Erdös-Reyni graph is parametrized
 * in terms of the number n of nodes, the average degree k of a node, and
 * the timescale tau on which edges appear and disappear. In terms of these
 * parameter, each edge has probability p_+ = k / (n - 1) to exist at at
 * certain point in time, which is satiesfied for rates of appearance and
 * disappearance of alpha = p_+ / tau and beta = p_- / tau = (1 - p_+) / tau.
 *
 * To see this, consider a two-state Markov process with states + (present)
 * and - (absent) with rates alpha for the transition - -> +, and beta for
 * + -> -. The steady-state probabilites are then p_+ = alpha / (alpha + beta),
 * and p_- = beta / (alpha + beta).
 */
struct dynamic_erdos_reyni : virtual dynamic_network, virtual erdos_reyni {
	dynamic_erdos_reyni(int size, double avg_degree, double timescale, rng_t& engine);

	virtual absolutetime_t next(rng_t& engine);

	virtual std::optional<network_event_t> step(rng_t& engine, absolutetime_t max_time = NAN);

	void add_edge(node_t node, node_t neighbour);
	void remove_edge(node_t node, int neighbour_index);

	/* Probability that a particular edge is present in the steady-state */
	const double edge_probability;
	/* Rate with which absent edges appear */
	const double alpha;
	/* Rate with which present edges vanish */
	const double beta;
	
	/* Current time (time of last event) */
	absolutetime_t current_time = 0.0;
	
	/* Time of next edge appearing or vanishing */
	absolutetime_t next_time = NAN;
	
	/* Unreported event for the reverse edge of the last edge event reported  */
	std::optional<network_event_t> reverse_edge_event;

	unsigned int edges_absent;
	unsigned int edges_present;
	
	/* Degree-weighted node distribution */
	dyndist::vector_distribution<unsigned> weighted_nodes;
};
