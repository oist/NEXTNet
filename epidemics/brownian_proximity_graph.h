//
//  brownian_proximity_graph.hpp
//  epidemics
//
//  Created by Florian G. Pflug on 24.04.23.
//

#pragma once

#include "stdafx.h"
#include "types.h"
#include "random.h"
#include "dynamic_graph.h"

//--------------------------------------
//----------DYNAMIC NETWORKS------------
//--------------------------------------

struct brownian_proximity_graph : virtual dynamic_network, virtual graph {
	typedef std::unordered_set<node_t> node_set_t;
	typedef std::unordered_map<node_t, index_t> neighbour_map_t;

	typedef std::pair<std::size_t, std::size_t> partition_index_t;

	struct node {
		point position;
		neighbour_map_t neighbour_map;
		std::vector<node_t> neighbours;
	};

	brownian_proximity_graph(node_t N, double R0, double radius, double D, rng_t& engine);

	brownian_proximity_graph(node_t N, double R0, double radius, double D, double dt, rng_t& engine);

	virtual ~brownian_proximity_graph();
	
	virtual node_t nodes();

	virtual node_t neighbour(node_t node, int neighbour_index) = 0;

	virtual int outdegree(node_t node) = 0;

	virtual absolutetime_t next(rng_t& engine);
	
	virtual std::optional<network_event_t> step(rng_t& engine, absolutetime_t nexttime = NAN);

	partition_index_t partition_index(point p) {
		return partition_index_t(std::trunc(p.y / plength), std::trunc(p.x / plength));
	}
	
	node_set_t& partition(partition_index_t i) {
		return partitions.at(i.first * pstride + i.second);
	}
	
	node_set_t& partition(point p) {
		return partition(partition_index(p));
	}
	
	const node_t size;
	const double radius;
	const double length;
	const double diffusivity;
	const double delta_t;
	double current_time = 0.0;
	
	std::vector<node> nodedata;
	
	const double plength;
	const std::size_t pstride;
	std::vector<node_set_t> partitions;
	
	/**
	 * @brief Internal state of the next() member. Stored so that next() can
	 * report events, and continue where it left of the next time it is called.
	 */
	struct state_t {
		bool displacement_done = false;
		node_t node_i = 0;
		bool neighbour_scan_done = false;
		bool neighbour_scan_initialized = false;
		neighbour_map_t::iterator n_it;
		bool range_scan_initialized = false;
		partition_index_t lb;
		partition_index_t ub;
		partition_index_t cur;
		bool partition_scan_initialized = false;
		node_set_t::const_iterator p_it;
	} state;
	
	std::optional<network_event_t> next_event;
};
