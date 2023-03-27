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
#include "utility.h"

//--------------------------------------
//--------------NETWORKS----------------
//--------------------------------------

/**
 * @brief Abstract interface to a "network", i.e. a (directed) graph.
 * 
 * Provides three functions `nodes`, `neighbour()` and `outdegree()`, which can be used
 * to query the network. `nodes` returns the number of nodes, `outdegree(n)` must return
 * the number of outgoing edges of node n, and `neighbour(n, i)` must return the target
 * of the i-th outgoing edge.
 */
class graph {
public: 
    virtual ~graph();
    
    /**
     * @brief Return the number of nodes in the graph. If the number is
     * infinite or unknown, -1 is returned.
     */
    virtual node_t nodes();

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
    virtual node_t nodes();

    virtual node_t neighbour(node_t node, int neighbour_index);

    virtual index_t outdegree(node_t node);

    /* Adjacency list of the graph */
    std::vector<std::vector<node_t>>  adjacencylist;
};
//

//--------------------------------------
//--------WATTS STROGATZ GRAPH----------
//--------------------------------------

/**
 * @brief A random Watts-Strogatz network
 */
class watts_strogatz : public graph_adjacencylist {
public:
	watts_strogatz(node_t size, int k, double p, rng_t& engine);

	watts_strogatz(int size, double p, rng_t& engine)
		:watts_strogatz(size, 2, p, engine)
	{}
};

//--------------------------------------
//--------ERDOS REYNI GRAPH-------------
//--------------------------------------

/**
 * @brief A random Erdös-Reyni network
 */
class erdos_reyni : public graph_adjacencylist {
public:
    erdos_reyni(int size, double avg_degree, rng_t& engine);
};

//--------------------------------------
//----------FULLY CONNECTED-------------
//--------------------------------------

/**
 * @brief A fully-connected network with random edge order
 */
class fully_connected : public graph {
public:
    fully_connected(int size, rng_t& engine);

    virtual node_t nodes();

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
//---------CONFIGURATION MODEL----------
//--------------------------------------
/**
 * @brief Network from arbitrary degree distribution. 
 */
class config_model : public graph_adjacencylist {
public:
    config_model(std::vector<int> degreelist, rng_t& engine);

    std::size_t selfloops = 0;
    std::size_t multiedges = 0;
    std::unordered_set<edge_t, pair_hash> edges = {};

};

//------------------------------------------
//--CONFIG MODEL: WITH CORRELATED DEGREES---
//------------------------------------------
/**
 * @brief Network from arbitrary degree distribution. 
 * 
 */
// TODO: Clean this up
#if 0
class config_model_correlated : public graph_adjacencylist {
public:
    config_model_correlated(std::vector<int> degreelist, rng_t& engine, bool assortative);

    std::size_t selfloops = 0;
    std::size_t multiedges = 0;
};
#endif

//--------------------------------------
//--------SCALE FREE NETWORK------------
//--------------------------------------
/**
 * @brief Power Law Network using Barabasi-Albert model.
 *
 * The degree distribution scales with k^-3.
 */
class scale_free : public graph_adjacencylist {
public:
    scale_free(int size, rng_t& engine);

};


//--------------------------------------
//--------IMPORTED NETWORK----------
//--------------------------------------
/**
 * @brief Network generated from an adjacency list
 * 
 * the file must be an adjacency list, i.e., a list of lists;
 */
class imported_network : public graph_adjacencylist {
public:
    imported_network(std::string path_to_file);

private:
    int file_size(std::string path_to_file);

};

//------------------------------------------
//--ADD DEGREE CORRELATION TO THE NETWORK---
//------------------------------------------
/**
 * @brief Add correlation to the network by rewiring its links.
 *
 */
void add_correlation(double r,graph_adjacencylist& nw,rng_t& engine);

// Helper function to verify whether an edge exists or not
bool edge_exists(node_t a, node_t b, const graph_adjacencylist& nw);


//------------------------------------------------
//-----Measure degree correlation in a network----
//------------------------------------------------

/**
 * @brief Average degree of the nearest neighbors of vertices of degree k.
 * 
 * A measure of degree correlation: when the network is uncorrelated,
 * knn(k) should be independent of k.
 *  
 */
std::vector<double> knn(graph_adjacencylist& nw);


/**
 * @brief Pearson correlation to measure the assortativity of a network
 *
 * assortativity a = num/den
 * num = sum_kk' [ w(k,k') * k * k' ] - ( sum_k [ w(k) * k ] ) ^ 2
 * den = sum_k [ w(k) * k ^ 2] - ( sum_k [ w(k) * k ] ) ^ 2
 *
 */
double assortativity(graph_adjacencylist& nw);


/**
 * @brief fraction of non-self links in the network that connect a node of degree k
 * to a node of degree k prime
 *
 */
std::vector<std::vector<double>> Wkk(graph_adjacencylist& nw);
