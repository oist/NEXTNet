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
class graph_adjacencylist : public virtual graph {
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
class watts_strogatz : public virtual graph_adjacencylist {
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
class erdos_reyni : public virtual graph_adjacencylist {
public:
    erdos_reyni(int size, double avg_degree, rng_t& engine);
};

//--------------------------------------
//----------FULLY CONNECTED-------------
//--------------------------------------

/**
 * @brief A fully-connected network with random edge order
 */
class fully_connected : public virtual graph {
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
class acyclic : public virtual graph {
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
class config_model : public virtual graph_adjacencylist {
public:
    config_model(std::vector<int> degreelist, rng_t& engine);

    std::size_t selfloops = 0;
    std::size_t multiedges = 0;
    std::unordered_set<edge_t, pair_hash> edges = {};

};

// some specific networks using the configuration model:
// TODO: Make these class methods

std::vector<int> lognormal_degree_list(double mean, double variance, int size, rng_t& engine);
std::vector<int> powerlaw_degree_list(double exponent, int size, rng_t& engine);

//--------------------------------------
//----- CLUSTERED CONFIGURATION MODEL---
//--------------------------------------

/**
 * @brief Network from arbitrary degree distribution.
 *
 * Based on the algorithm described by Serrano & Boguna, 2005.
 */
class config_model_clustered_serrano : public virtual graph_adjacencylist {
public:
	/**
	 * @brief Converts c(k) into a number of triangles per degree class.
	 *
	 * Given c(k) and given the degrees of all nodes, we generate a vector containing a triangle count per degree.
	 * The vector satisfies c(k) *on average*, see footnote in Serrano & Boguna, 2005 p. 3.
	 *
	 * @param ck the function c(k) which returns the probability that two neighbours of a node of degree k are themselves neighbours
	 * @param degrees list of node degrees
	 * @param engine random number generator
	 * @return a list of triangles per nodes of degree 0, 1, 2, ...
	 */
	static std::vector<int> triangles_binomial(std::function<double(int)> ck, std::vector<int> degrees, rng_t& engine);
	
	/**
	 * @brief Generates a random network with nodes of the given degrees given numbers of triangles per degree class
	 *
	 * This is the most general form where any expected number of triangles in any degree class k can be prescribed.
	 * However, some combinations of degrees and triangles may not allow a network to be constructed.
	 *
	 * @param degrees list of node degrees
	 * @param triangles number of triangles overlapping nodes of degree 0, 1, 2, ...,
	 * @param beta parameter that defines degree class probabilities P(k)
	 * @param engine random number generator
	 */
	config_model_clustered_serrano(std::vector<int> degrees, std::vector<int> triangles, double beta, rng_t& engine);

    /**
     * @brief Generates a random network with nodes of the given degrees and clustering given by c(k)
     *
     * c(k) defines the probability for a node of degree k that two randomly chosen neighbours
     * are themselves neighbours, i.e. form a triangle.
     *
     * @param degrees list of node degrees
     * @param ck the function c(k) which returns the probability that two neighbours of a node are themselves neighbours
     * @param beta parameter that defines degree class probabilities P(k)
     * @param engine random number generator
     */
    config_model_clustered_serrano(std::vector<int> degrees, std::function<double(int)> ck, double beta, rng_t& engine);

	/**
	 * @brief Generates a random network with nodes of the given degrees and triangles specified by alpha
	 *
	 * The number of triangles per degree class is chosen such that the probability c(k) that
	 * two randomly chosen neighbours of a node are themselves neighbours is c(k)=0.5*(k-1)^alpha
	 * where k is the degree of the node.
	 *
	 * @param degrees list of node degrees
	 * @param alpha parameter of c(k)
	 * @param beta parameter that defines degree class probabilities P(k)
	 * @param engine random number generator
	 */
	config_model_clustered_serrano(std::vector<int> degree, double alpha, double beta, rng_t& engine);
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
class config_model_correlated : public virtual graph_adjacencylist {
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
class scale_free : public virtual graph_adjacencylist {
public:
    scale_free(int size, rng_t& engine,int m = 1);

};

//--------------------------------------
//--------LATTICE-----------------------
//--------------------------------------

template<unsigned int D>
class cubic_lattice : public virtual graph {
public:
    const static unsigned int dimension = D;
    const static std::size_t max_coordinate_digits = std::numeric_limits<node_t>::digits / D;
    const static std::size_t max_edge_length = 1 << max_coordinate_digits;

    typedef long coordinate_t;
    const static std::size_t coordinate_type_digits = std::numeric_limits<coordinate_t>::digits;
    typedef std::array<coordinate_t, dimension> coordinates_t;

    const std::size_t edge_length;
    const std::size_t coordinate_digits;
    const std::size_t coordinate_mask;
    const coordinate_t coordinate_max;
    const coordinate_t coordinate_min;

    cubic_lattice()
        :cubic_lattice(max_edge_length)
    {}

    cubic_lattice(const std::size_t _edge_length)
        :edge_length(_edge_length)
        ,coordinate_digits(std::ceil(std::log2(edge_length)))
        ,coordinate_mask(((std::size_t)1 << coordinate_digits) - 1)
        ,coordinate_min(-(edge_length/2))
        ,coordinate_max((edge_length-1)/2)
    {
        if (edge_length > max_edge_length)
            throw std::runtime_error("edge length exceeds maximal supported length");
        assert(coordinate_type_digits >= coordinate_digits);
        assert(coordinate_max - coordinate_min + 1 == edge_length);
    }

    node_t nodes() {
        return std::pow(edge_length, dimension);
    }

    node_t neighbour(node_t n, int neighbour_index) {
        // Decode neighbour index into coordinate components
        coordinates_t c = coordinates(n);
        // Count the number of extremal (i.e min or max) and non-extermal components
        std::size_t e = 0;
        std::array<std::size_t, dimension> ec = { -1 };
        for(std::size_t i = 0; i < dimension; ++i) {
            if ((c[i] != coordinate_max) && (c[i] != coordinate_min))
                continue;
            ec[e++] = i;
        }
        // For ever non-extremal component there are two neighbours
        // for every extermal component one neighbour.
        const std::size_t ne = dimension - e;
        if (neighbour_index < 0) {
            return -1;
        } else if (neighbour_index < ne) {
            // Increment non-extermal component
            c[neighbour_index] += 1;
        } else if (neighbour_index < 2*ne) {
            // Decrement non-extermal component
            c[neighbour_index - ne] -= 1 ;
        } else if (neighbour_index < 2*ne + e) {
            // Increment/decrement extremal component
            const std::size_t ei = neighbour_index - 2*ne;
            const std::size_t ci = ec[ei];
            c[ci] += (c[ci] == coordinate_max) ? -1 : 1;
        } else {
           return -1;
        }
        return node(c);
    }

    int outdegree(node_t node) {
        // Decode neighbour index into coordinate components
        coordinates_t c = coordinates(node);
        // Count the number of extremal (i.e min or max) and non-extermal components
        std::size_t e = 0;
        for(std::size_t i = 0; i < dimension; ++i) {
            if ((c[i] != coordinate_max) && (c[i] != coordinate_min))
                continue;
            e++;
        }
        // For ever non-exterminal component there is one neighbour, otherwise two
        return 2*dimension - e;
    }

    coordinate_t extract_coordinate(node_t node, std::size_t i) const {
        assert(node >= 0);
        // Extract relevat bits
        const std::size_t bits = (((std::size_t)node) >> (i*coordinate_digits)) & coordinate_mask;
        // Convert into a signed coordinate. Have to take care to sign-extend properly
        const std::size_t shift = coordinate_type_digits - coordinate_digits + 1;
        const coordinate_t c = (coordinate_t)(bits << shift) >> shift;
        assert(c >= coordinate_min);
        assert(c <= coordinate_max);
        return c;
    }

    node_t embedd_coordinate(std::size_t i, coordinate_t c) const {
        // Chheck that all bits beyond coordinate_digits are either all zero or all one.
        assert(((c >> coordinate_digits) == -1) || ((c >> coordinate_digits) == 0));
        return (((std::size_t)c) & coordinate_mask) << (i*coordinate_digits);
    }

    template<std::size_t... i>
    auto extract_coordinates(node_t node, std::index_sequence<i...>) const {
        return std::array<coordinate_t, sizeof...(i)>{{extract_coordinate(node, i)...}};
    }

    /**
     * @brief Convert node index into d-dimensional coordinate vector
     * @param node node index
     * @return d-dimensional coordinate vector as a std::array<coordinate_t, D>
     */
    coordinates_t coordinates(node_t node) {
        return extract_coordinates(node, std::make_index_sequence<dimension>{});
    }

    /**
     * @brief Convert d-dimensional coordindate vector into node index
     * @param p d-dimensional coordinate vector as a std::array<coordinate_t, D>
     * @return node index
     */
    node_t node(const coordinates_t& p) const {
        std::size_t bits = 0;
        for(std::size_t i=0; i < D; ++i)
            bits |= embedd_coordinate(i, p[i]);
        return (node_t)bits;
    }
};

template<unsigned int D>
const unsigned int cubic_lattice<D>::dimension;

template<unsigned int D>
const std::size_t cubic_lattice<D>::max_coordinate_digits;

template<unsigned int D>
const std::size_t cubic_lattice<D>::max_edge_length;

template<unsigned int D>
const std::size_t cubic_lattice<D>::coordinate_type_digits;

typedef cubic_lattice<2> cubic_lattice_2d;
typedef cubic_lattice<3> cubic_lattice_3d;
typedef cubic_lattice<4> cubic_lattice_4d;
typedef cubic_lattice<5> cubic_lattice_5d;
typedef cubic_lattice<6> cubic_lattice_6d;
typedef cubic_lattice<7> cubic_lattice_7d;
typedef cubic_lattice<8> cubic_lattice_8d;

//--------------------------------------
//--------IMPORTED NETWORK----------
//--------------------------------------
/**
 * @brief Network generated from an adjacency list
 * 
 * the file must be an adjacency list, i.e., a list of lists;
 */
class imported_network : public virtual graph_adjacencylist {
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
