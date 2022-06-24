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

//--------------------------------------
//----------ADJACENCYLIST GRAPH---------
//--------------------------------------

/**
 * @brief Base class for networks defined by an adjacency list.
 * 
 * Implements functions `neighbour()` and `outdegree()`, the constructor
 * is expected to setup the adjacencylist neighbours.
 */
class graph_adjacencylist : public graph {
public:
    virtual node_t neighbour(node_t node, int neighbour_index);

    virtual index_t outdegree(node_t node);

    /* Adjacency list of the graph */
    std::vector<std::vector<node_t>>  adjacencylist;
};

//--------------------------------------
//--------------ER GRAPH----------------
//--------------------------------------

/**
 * @brief A random Erdös-Reyni network
 */
class erdos_reyni : public graph_adjacencylist {
public:
    erdos_reyni(int size, double avg_degree, rng_t& engine);
};

//--------------------------------------
//-------FULLY CONNECTED----------------
//--------------------------------------

/**
 * @brief A fully-connected network with random edge order
 */
class fully_connected : public graph {
public:
    fully_connected(int size, rng_t& engine);

    virtual node_t neighbour(node_t node, int neighbour_index);

    virtual index_t outdegree(node_t node);

    rng_t& engine;
    std::vector<std::vector<node_t>> neighbours;
};

//--------------------------------------
//---------------TREE-------------------
//--------------------------------------


/**
 * @brief A random acyclic network
 */
class acyclic : public graph {
public:
    static double lambda(double mean, int digits);

    acyclic(double avg_degree, bool reduced_root_degree, rng_t& engine);

    virtual node_t neighbour(node_t node, int neighbour_index);

    virtual index_t outdegree(node_t node);

    rng_t& engine;
    std::poisson_distribution<> degree_distribution;
    bool reduced_root_degree;
    std::deque<std::vector<node_t>> adjacencylist;

private:
    static const node_t incomplete_neighbours = -1;

    void generate_neighbours(node_t node);
};


//--------------------------------------
//------CONFIGURATION MODEL-------------
//--------------------------------------
/**
 * @brief Network from arbitrary degree distribution. 
 */
class config_model : public graph_adjacencylist {
public:
    config_model(std::vector<int> degreelist, rng_t& engine);

    std::size_t selfloops = 0;
    std::size_t multiedges = 0;
};