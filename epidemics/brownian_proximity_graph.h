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

struct brownian_proximity_graph : virtual dynamic_network, virtual graph, virtual graph_embedding {
	enum node_state_t {
		NONINFECTED = 0,
		INFECTED = 1
	};
	
	struct node_data {
		node_t index;
		point position;
		unsigned int generation;
		node_state_t node_state;
		indexed_set<node_t> neighbours;
	};

	typedef std::vector<node_data> node_vector_t;
	typedef std::pair<unsigned int, unsigned int> partition_index_t;
	typedef std::pair<partition_index_t, unsigned int> partition_node_index_t;

	brownian_proximity_graph(node_t N, double avg_degree, double radius,
							 double D, rng_t& engine);

	brownian_proximity_graph(node_t N, double avg_degree, double radius,
							 double D, double dt, rng_t& engine);
	
	brownian_proximity_graph(node_t N, double avg_degree, double radius,
							 double D0, double D1, double gamma, rng_t& engine);

	brownian_proximity_graph(node_t N, double avg_degree, double radius,
							 double D0, double D1, double gamma, double dt, rng_t& engine);

	virtual ~brownian_proximity_graph();
	
	virtual node_t nodes() override;

	virtual node_t neighbour(node_t node, int neighbour_index) override;

	virtual int outdegree(node_t node) override;

	virtual std::size_t dimensionality() override;

	virtual bool coordinates(const node_t node, std::vector<double>& position) override;

	virtual void bounds(std::vector<double>& a, std::vector<double>& b) override;

	virtual absolutetime_t next(rng_t& engine, absolutetime_t maxtime = INFINITY) override;
	
	virtual std::optional<network_event_t> step(rng_t& engine, absolutetime_t nexttime = NAN) override;

	virtual void notify_epidemic_event(event_t ev, rng_t& engine) override;
	
	double node_diffusivity(const node_data& n);
	
	node_data& node(std::size_t n) {
		return node(node_index.at(n));
	}
	
	node_data& node(partition_node_index_t pni) {
		return partition(pni.first).at(pni.second);
	}
	
	partition_index_t partition_index(std::size_t i) {
		return partition_index_t(i / pstride, i % pstride);
	}
	
	partition_index_t partition_index(point p) {
		return partition_index_t(std::min<float>(std::max<float>(0, std::trunc(p.y / plength)), pstride-1),
								 std::min<float>(std::max<float>(0, std::trunc(p.y / plength)), pstride-1));
	}
	
	node_vector_t& partition(partition_index_t i) {
		assert((i.first >= 0) && (i.first < pstride) && (i.second >= 0) && (i.second < pstride));
		return partitions.at(i.first * pstride + i.second);
	}
	
	node_vector_t& partition(point p) {
		return partition(partition_index(p));
	}
	
	void move_node(node_data& n, partition_index_t pi_old, partition_index_t pi_new);
	
	const node_t size;
	const float radius;
	const float length;
	const double diffusivity_noninfected;
	const double diffusivity_infected;
	const double gamma;
	const double delta_t;
	double current_time = 0.0;
	unsigned int current_generation = 0;
	std::size_t ninfected = 0;
	
	const float plength;
	const std::size_t pstride;

	std::vector<node_vector_t> partitions;
	std::vector<partition_node_index_t> node_index;

	/**
	 * @brief Internal state of the next() member. Stored so that next() can
	 * report events, and continue where it left of the next time it is called.
	 */
	struct state_t {
		bool displacement_done = false;
		std::size_t partition_i = 0;
		bool outer_partition_scan_initialized = false;
		std::size_t outer_partition_node_i;
		bool neighbour_scan_done = false;
		bool neighbour_scan_initialized = false;
		node_t neighbour_idx;
		bool range_scan_initialized = false;
		partition_index_t lb;
		partition_index_t ub;
		partition_index_t cur;
		bool inner_partition_scan_initialized = false;
		std::size_t inner_partition_node_i = 0;
	} state;
	
	std::optional<network_event_t> next_event;
};
