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

struct dynamic_empirical_network : virtual dynamic_network {

	dynamic_empirical_network(std::string path_to_file,double dt);

	virtual absolutetime_t next(rng_t& engine);

	virtual std::optional<network_event_t> step(rng_t& engine, absolutetime_t max_time = NAN);

	/* Current time (time of last event) */
	absolutetime_t current_time = 0.0;
	
	/* Time of next edge appearing or vanishing */
	absolutetime_t next_time = NAN;

	/* timestamps of all appearing or vanishing edges (sorted vector)*/
	std::vector<network_event_t> edges; 

	/* Index of the next network_event*/
	int time_index = 0;

	int max_index = -1;

	node_t nodes();

	node_t neighbour(node_t node, int neighbour_index) ;

	index_t outdegree(node_t node);

    /* Adjacency list of the graph */
    std::vector<std::vector<node_t>>  adjacencylist;
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
