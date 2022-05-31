//
//  graph.hpp
//  epidemics
//
//  Created by Cure Samuel Cyrus on 2022/03/07.
//

#pragma once

#include "stdafx.h"
#include "types.h"
#include "random.h"

//--------------------------------------
//--------------NETWORKS----------------
//--------------------------------------

/**
 * @brief Abstract interface to a "network", i.e. a (directed) graph.
 * 
 * Provides two functions `neighbour()` and `outdegree()`, which can be used
 * to query the network. `outdegree(n)` must return the number of outgoing edges
 * of node n, and `neighbour(n, i)` must return the target of the i-th outgoing
 * edge.
 */
class graph {
public:    
    /**
     * @brief Returns the target of the i-th outgoing edge of node n.
     */
    virtual node_t neighbour(node_t node, int neighbour_index) = 0;

    /**
     * @brief Returns the number of (outgoing) edges of the given node
     */
    virtual int outdegree(node_t node) = 0;
};


/**
 * @brief A Erdös-Reyni network
 */
class erdos_reyni : public graph {
public:
    erdos_reyni(int size, double avg_degree, rng_t& engine);

    virtual node_t neighbour(node_t node, int neighbour_index);

    virtual index_t outdegree(node_t node);

//private:
    /* Adjacency list of the graph */
    std::vector<std::vector<node_t>> neighbours;
};

class fully_connected : public graph {
public:
    fully_connected(int size, rng_t& engine);

    virtual node_t neighbour(node_t node, int neighbour_index);

    virtual index_t outdegree(node_t node);

    rng_t& engine;
    std::vector<std::vector<node_t>>  neighbours;
};