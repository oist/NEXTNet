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

	virtual std::optional<network_event_t> step(rng_t& engine, absolutetime_t nexttime = NAN) = 0;

	virtual void notify_epidemic_event(event_t ev);
};

struct dynamic_erdos_reyni : virtual dynamic_network, virtual erdos_reyni {
	dynamic_erdos_reyni(int size, double avg_degree, double timescape, rng_t& engine);

	virtual absolutetime_t next(rng_t& engine);

	virtual std::optional<network_event_t> step(rng_t& engine, absolutetime_t nexttime = NAN);

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
	
	unsigned int edges_absent;
	unsigned int edges_present;
	
	/* Degree-weighted node distribution */
	dyndist::vector_distribution<unsigned> weighted_nodes;
};
