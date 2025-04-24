//
//  weighted_network.hpp
//  epidemics
//
//  Created by Cure Samuel Cyrus on 2022/03/07.
//

#pragma once

#include "nextnet/stdafx.h"
#include "nextnet/types.h"
#include "nextnet/network.h"
#include "nextnet/random.h"
#include "nextnet/utility.h"

//--------------------------------------
//--------------WEIGHTED NETWORK--------
//--------------------------------------

/**
 * @brief Extends the abstract interface <network> to support weighted networks
 *
 * The function `neighbour(node_t node, int neighbour_index, double* weight)`
 */
class weighted_network : public virtual network
{
public:
    virtual ~weighted_network();

    virtual bool is_unweighted();

    /**
     * @brief Returns the target of the i-th outgoing edge of node <node>
     * @param node node to query
     * @param neighbour_index index of neighbour to query
     * @param weight pointer to a double in which the weight of the edge is stored
     */
    using network::neighbour;
    virtual node_t neighbour(node_t node, int neighbour_index, double *weight) = 0;

    /**
     * @brief Returns the target of the i-th outgoing edge of node n.
     *
     * Forwards to `neighbour(node_t node, int neighbour_index, double* weight)`
     */
    virtual node_t neighbour(node_t node, int neighbour_index) override;
};

/**
 * @brief Converts a points to a network into a pointer to a weighted number if
 * the network is actually weighted (i.e. an instance of weighted_network and is_unweighted
 * is false).
 */
weighted_network *as_weighted_network(network *nw);

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
class weighted_adjacencylist_network : public virtual weighted_network
{
public:
	virtual bool is_undirected() override;

	virtual bool is_simple() override;

	virtual node_t nodes() override;

	virtual node_t neighbour(node_t node, int neighbour_index, double *weight) override;

	virtual index_t outdegree(node_t node) override;

protected:
	weighted_adjacencylist_network()
	{}

	weighted_adjacencylist_network(bool undirected_, bool simple_)
		:undirected(undirected_), simple(simple_)
	{}

    /* Adjacency list of the graph */
    std::vector<std::vector<std::pair<node_t, double>>> adjacencylist;

	bool undirected = false;
	bool simple = false;
};

//-------------------------------------------
//--------WEIGHTED EMPIRICAL NETWORK---------
//-------------------------------------------

/**
 * @brief Weighted empirical network read from a file
 *
 * The file must contain lines of the form
 *
 * <nodeA> <csep> <nodeB1> <wsep> <weight1> <csep> <nodeB2> <wsep> <weight2> ...
 *
 * where <csep> is ' ' by default and wsep is ':'. Each such line defines
 * edges A -> B1, A -> B2, ... with weights w1, w2, ... Multiple lines for
 * the same node A are allowed and the edges are merged by summing their weights.
 * This class than thus read adjacencylist and well as edgelist files. Lines
 * starting with '#' are ignored.
 *
 * If undirected is true, the network is assumed to be undirected, i.e.
 * for every edge (a,b) the reverse edge (b,a) is also added. If
 * simplify is true, self-edges are removed. Multi-edges per definition
 * do not exist -- what would be a multi-edge in an unweighted networks
 * is simply an edge with a higher weight here.
 * 
 */
class weighted_empirical_network : public virtual weighted_adjacencylist_network
{
public:
    weighted_empirical_network(
        std::istream& file, bool undirected = true, bool simplify = false,
        node_t idxbase = 1, char csep = ' ', char wsep = ':');
};

//--------------------------------------
//----------WEIGHTED ERDÖS RENYI--------
//--------------------------------------

/**
 * @brief A random Erdös-Reyni network
 */
class weighted_erdos_renyi : public virtual weighted_adjacencylist_network
{
public:
    /**
     * Creates a weighted Erdös-Reyni network in which weights are assigned i.i.d. using
     * the weightdist function. For each edge, weightdist(engine) is called to randomly generate
     * a weight.
     */
    weighted_erdos_renyi(int size, double avg_degree, std::function<double(rng_t &)> weightdist, rng_t &engine);

    /**
     * Creates a weighted Erdös-Reyni network in which weights are assigned i.i.d. by randomly
     * sampling weights from the vector w with probabilities p.
     */
    weighted_erdos_renyi(int size, double avg_degree, std::vector<double> w, std::vector<double> p, rng_t &engine);

private:
    template <typename ItW, typename ItP>
    static std::function<double(rng_t &)> make_dist(ItW w_begin, ItW w_end, ItP p_begin, ItP p_end)
    {
        std::vector<double> w(w_begin, w_end);
        std::shared_ptr<std::discrete_distribution<std::size_t>> idx(
            new std::discrete_distribution<std::size_t>(p_begin, p_end));
        if (w.size() - 1 != idx->max())
            throw std::runtime_error("number of weights and probabilities must agree");
        return [idx, w](rng_t &engine) { return w.at((*idx)(engine)); };
    }
};

/* Compatibility with previously miss-spelled name */
typedef weighted_erdos_renyi weighted_erdos_reyni;
