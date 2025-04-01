//
//  weighted_network.hpp
//  epidemics
//
//  Created by Cure Samuel Cyrus on 2022/03/07.
//

#pragma once

#include "stdafx.h"
#include "types.h"
#include "network.h"
#include "random.h"
#include "utility.h"

//--------------------------------------
//--------------WEIGHTED NETWORK--------
//--------------------------------------

/**
 * @brief Extends the abstract interface <network> to support weighted networks
 *
 * The function `neighbour(node_t node, int neighbour_index, double* weight)`
 */
class weighted_network : public virtual network {
public:
	virtual ~weighted_network();

	/**
	 * @brief Returns the target of the i-th outgoing edge of node <node>
	 * @param node node to query
	 * @param neighbour_index index of neighbour to query
	 * @param weight pointer to a double in which the weight of the edge is stored
	 */
	virtual node_t neighbour(node_t node, int neighbour_index, double* weight) = 0;
	
	/**
	 * @brief Returns the target of the i-th outgoing edge of node n.
	 *
	 * Forwards to `neighbour(node_t node, int neighbour_index, double* weight)`
	 */
	virtual node_t neighbour(node_t node, int neighbour_index) override;
};

//--------------------------------------
//------WEIGHTED ADJACENCYLIST GRAPH----
//--------------------------------------

/**
 * @brief Base class for weighted networks defined by an adjacency list.
 *
 * Stores a vector of (neighbour, weight) pairs for every node.
 *
 * Implements functions `neighbour()` and `outdegree()`, the constructor
 * is expected to setup the adjacencylist neighbours.
 */
class weighted_adjacencylist_network : public virtual weighted_network {
public:
	virtual node_t nodes() override;

	virtual node_t neighbour(node_t node, int neighbour_index, double* weight) override;

	virtual index_t outdegree(node_t node) override;

	/* Adjacency list of the graph */
	std::vector<std::vector<std::pair<node_t, double>>>  adjacencylist;
};

//--------------------------------------
//----------WEIGHTED ERDÖS RENYI--------
//--------------------------------------

/**
 * @brief A random Erdös-Reyni network
 */
class weighted_erdos_renyi : public virtual weighted_adjacencylist_network, public virtual network_is_undirected {
public:
	/**
	 * Creates a weighted Erdös-Reyni network in which weights are assigned i.i.d. using
	 * the weightdist function. For each edge, weightdist(engine) is called to randomly generate
	 * a weight.
	 */
	weighted_erdos_renyi(int size, double avg_degree, std::function<double(rng_t&)> weightdist, rng_t& engine);
	
	/**
	 * Creates a weighted Erdös-Reyni network in which weights are assigned i.i.d. by randomly
	 * sampling weights from the vector w with probabilities p.
	 */
	weighted_erdos_renyi(int size, double avg_degree, std::vector<double> w, std::vector<double> p, rng_t& engine);
	
private:
	template<typename ItW, typename ItP>
	static std::function<double(rng_t&)> make_dist(ItW w_begin, ItW w_end, ItP p_begin, ItP p_end) {
		std::vector<double> w(w_begin, w_end);
		std::shared_ptr<std::discrete_distribution<std::size_t>> idx(
			new std::discrete_distribution<std::size_t>(p_begin, p_end));
		if (w.size() - 1 != idx->max())
			throw std::runtime_error("number of weights and probabilities must agree");
		return [idx, w](rng_t& engine) { return w.at((*idx)(engine)); };
	}
};

/* Compatibility with previously miss-spelled name */
typedef weighted_erdos_renyi weighted_erdos_reyni;
