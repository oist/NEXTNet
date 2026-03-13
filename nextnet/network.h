//
//  graph.hpp
//  epidemics
//
//  Created by Cure Samuel Cyrus on 2022/03/07.
//

#pragma once

#include "nextnet/stdafx.h"
#include "nextnet/types.h"
#include "nextnet/random.h"
#include "nextnet/utility.h"

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
class network
{
public:
    virtual ~network();

    /**
     * @brief Whether the graph is undirected, i.e. whether for every
     * edge (i,j) there is also an edge (j,i)
     */
    virtual bool is_undirected() = 0;
    
    /**
     * @brief Whether the graph is weighted, i.e. each edge has an
     * associated weight. On unweighted networks, each edge has
     * implicit weight 1.0.
     */
    virtual bool is_unweighted() = 0;

    /**
     * @brief Whether the network is layered, i.e. whether edges have
     * an associated layer. On non-layered networks, all edges exist on
     * layer 0.
     */
    virtual bool is_layered() = 0;


    /**
     * @brief Whether the graph is simple, i.e. does not contain
     * self edges (i,i) or multi-edges (i.e. multipel edges (i,j) for the same
     * nodes i and j).
     */
    virtual bool is_simple() = 0;
    
    /**
     * @brief Whether the individual layers are simple, i.e. there are no
     * self edges (i,i) and each layer contains each edge (i,j) at most once.
     */
    virtual bool has_simple_layers();
    
    /**
     * @brief Return the number of nodes in the graph. If the number is
     * infinite or unknown, -1 is returned.
     */
    virtual node_t nodes();

    /**
     * @brief Returns the target of the i-th outgoing edge of node n, the edge's layer in `layer` and the edge's weight in `weight`
     */
    virtual node_t neighbour(node_t node, int neighbour_index, edgelayer_t* layer, double* weight) = 0;

    /**
     * @brief Returns the target of the i-th outgoing edge of node n
     * Forwards to the 4-parameter version by default
     */
    virtual node_t neighbour(node_t node, int neighbour_index);

    /**
     * @brief Returns the number of (outgoing) edges of the given node
     */
    virtual int outdegree(node_t node) = 0;
};

/**
 * @brief Convenience class to mark a graph as undirected
 *
 * This avoid having to override is_undirected() in all undirected graphs to
 * return false. Instead, it suffices to additionally inherit from network_is_undirected
 */
class network_is_undirected : public virtual network
{
    virtual bool is_undirected();
};

/**
 * @brief Convenience class to mark a graph as simple
 *
 * This avoid having to override is_simple() in all undirected graphs to
 * return false. Instead, it suffices to additionally inherit from network_is_simple
 */
class network_is_simple : public virtual network
{
    virtual bool is_simple();
};

/**
 * @brief Convenience class to mark a graph as not having layers
 *
 * This avoid having to override is_layered() in all non-layered graphs to
 * return false. Instead, it suffices to additionally inherit from network_is_not_layered.
 *
 * In addition to returning true from is_layered(), `neighbour(node, index, layer)` is
 * overriden to alway set layer to zero and to call `neighbour(node, index)`.
 */
class network_is_not_layered_and_not_weighted : public virtual network
{
    virtual bool is_layered() override;

    virtual bool is_unweighted() override;

    /**
     * @brief Returns the target of the i-th outgoing edge of node n and the edge's layer in `layer`
     * Forwards to `neighbour(node, index)` by default and sets layer to zero
     */
    virtual node_t neighbour(node_t node, int neighbour_index, edgelayer_t* layer, double* weight) override;

    /**
     * @brief Returns the target of the i-th outgoing edge of node n and the edge's layer in `layer`
     */
    virtual node_t neighbour(node_t node, int neighbour_index) override = 0;
};

/**
 * @brief Interface for embedded graphs, i.e. graphs whose nodes are located in d-dimensional space
 *
 * Provides function to query the number of dimensions of the embedding space, and the
 * coordinates of individual nodes.
 */
class network_embedding
{
public:
    virtual ~network_embedding();

    /**
     * @brief Returns the dimensionality of the embedding space
     * @return number of dimensions
     */
    virtual std::size_t dimensionality() = 0;

    /**
     * @brief Returns the embedding coordinates of the specified node
     * @param node node index
     * @param position reference to a vector, first d entries are filled with the coordinates
     */
    virtual bool coordinates(const node_t node, std::vector<double> &position) = 0;

    /**
     * @brief Returns the bounding box of the embedding
     * @param a lower-left corner of vertex bounding box
     * @param b upper-right corner of vertex bounding box
     */
    virtual void bounds(std::vector<double> &a, std::vector<double> &b) = 0;
};

//--------------------------------------
//----------ADJACENCYLIST NETWEORK------
//--------------------------------------

namespace adjacencylist_types {
    template<bool Layered>
    struct member_layer {
        member_layer(edgelayer_t l) :layer(l) {}
        edgelayer_t layer;
    };
    
    template<>
    struct member_layer<false> {
        member_layer(edgelayer_t l) {}
        constexpr static edgelayer_t layer = 0;
    };
    
    template<bool Weighted>
    struct member_weight {
        member_weight(double w) :weight(w) {}
        double weight;
    };
    
    template<>
    struct member_weight<false> {
        member_weight(double w) {}
        constexpr static double weight = 1.0;
    };

    template<bool Layered, bool Weighted>
    struct initializer{};

    template<> struct initializer<false, false> {
        typedef node_t type;
        static constexpr node_t node(type v) { return v; }
        static constexpr edgelayer_t layer(type v) { return 0; }
        static constexpr double weight(type v) { return 1.0; }
    };
    template<> struct initializer<false, true> {
        typedef std::tuple<node_t, double> type;
        static constexpr node_t node(type v) { return std::get<0>(v); }
        static constexpr edgelayer_t layer(type v) { return 0; }
        static constexpr double weight(type v) { return std::get<1>(v); }
    };
    template<> struct initializer<true, false> {
        typedef std::tuple<node_t, edgelayer_t> type;
        static constexpr node_t node(type v) { return std::get<0>(v); }
        static constexpr edgelayer_t layer(type v) { return std::get<1>(v); }
        static constexpr double weight(type v) { return 1.0; }
    };
    template<> struct initializer<true, true> {
        typedef std::tuple<node_t, edgelayer_t, double> type;
        static constexpr node_t node(type v) { return std::get<0>(v); }
        static constexpr edgelayer_t layer(type v) { return std::get<1>(v); }
        static constexpr double weight(type v) { return std::get<2>(v); }
    };
}

/**
 * @brief Base class for networks defined by an adjacency list.
 *
 * Implements functions `neighbour()` and `outdegree()`, the constructor
 * is expected to setup the adjacencylist neighbours.
 */
template<bool Layered, bool Weighted>
struct adjacencylist_network : public virtual network
{
    typedef adjacencylist_types::member_layer<Layered> layer_member_t;
    typedef adjacencylist_types::member_weight<Weighted> weight_member_t;
    
    struct adjacencylist_entry_t: public layer_member_t, weight_member_t
    {
        node_t node;
        
        typedef adjacencylist_types::initializer<Layered, Weighted> initializer_t;
        typedef typename initializer_t::type initializer_value_t;

        adjacencylist_entry_t(initializer_value_t v)
            :node(initializer_t::node(v))
            ,layer_member_t(initializer_t::layer(v))
            ,weight_member_t(initializer_t::weight(v))
        {}
    };
    
    typedef std::vector<std::vector<adjacencylist_entry_t>> adjacencylist_t;

    adjacencylist_network(bool undirected_, bool simple_)
        : undirected(undirected_)
        , simple(simple_)
        , layers_simple(simple_)
    {
    }

    adjacencylist_network(bool undirected_, bool simple_, bool layers_simple_)
        : undirected(undirected_)
        , simple(simple_)
        , layers_simple(layers_simple_)
    {
    }

    
    adjacencylist_network(adjacencylist_t &&al, bool undirected_, bool simple_)
        : adjacencylist(std::move(al))
        , undirected(undirected_)
        , simple(simple_)
        , layers_simple(simple_)
    {
    }

    adjacencylist_network(adjacencylist_t &&al, bool undirected_, bool simple_, bool layers_simple_)
        : adjacencylist(std::move(al))
        , undirected(undirected_)
        , simple(simple_)
        , layers_simple(layers_simple_)
    {
    }

    virtual bool is_undirected() { return undirected; }

    virtual bool is_unweighted() { return !Weighted; }

    virtual bool is_layered() { return Layered; }

    virtual bool is_simple() { return simple; }

    virtual bool has_simple_layers() { return layers_simple; }

    virtual node_t nodes() { return (node_t)adjacencylist.size(); }

    virtual int outdegree(node_t node) { return (index_t)adjacencylist.at(node).size(); }

    using network::neighbour;
    virtual node_t neighbour(node_t node, int neighbour_index, edgelayer_t* layer, double* weight) {
        const auto &n = adjacencylist.at(node);
        if ((neighbour_index < 0) || (n.size() <= (unsigned int)neighbour_index)) return -1;
        const auto& e = n[neighbour_index];
        if (layer) *layer = e.layer;
        if (weight) *weight = e.weight;
        return e.node;
    }

protected:
    adjacencylist_t adjacencylist;
    bool undirected    = false;
    bool simple        = false;
    bool layers_simple = false;
};

typedef adjacencylist_network<false, false> unlayered_unweighted_adjacencylist_network;
typedef adjacencylist_network<false, true> unlayered_weighted_adjacencylist_network;
typedef adjacencylist_network<true, false> layered_unweighted_adjacencylist_network;
typedef adjacencylist_network<true, true> layered_weighted_adjacencylist_network;

//--------------------------------------
//--------- PROXY NETWORK --------------
//--------------------------------------

struct proxy_network : public virtual network
{
    proxy_network(std::unique_ptr<network>&& nw)
        :impl(std::move(nw))
    {}
    
    virtual bool is_undirected() override { return impl->is_undirected(); }
    virtual bool is_simple() override { return impl->is_simple(); }
    virtual bool has_simple_layers() override { return impl->has_simple_layers(); }
    virtual bool is_layered() override { return impl->is_layered(); }
    virtual node_t nodes() override{ return impl->nodes(); }
    virtual node_t neighbour(node_t n, int i, edgelayer_t* l, double* w) override { return impl->neighbour(n, i, l, w); }
    virtual int outdegree(node_t n) override { return impl->outdegree(n); }
    
protected:
    std::unique_ptr<network> impl;
};

//--------------------------------------
//--------WATTS STROGATZ GRAPH----------
//--------------------------------------

/**
 * @brief A random Watts-Strogatz network
 */
class watts_strogatz : public virtual unlayered_unweighted_adjacencylist_network
{
public:
    watts_strogatz(node_t size, int k, double p, rng_t &engine);

    watts_strogatz(int size, double p, rng_t &engine)
        : watts_strogatz(size, 2, p, engine)
    {
    }
};

//--------------------------------------
//--------ERDOS RENYI GRAPH-------------
//--------------------------------------

/**
 * @brief A random Erdös-Reyni network
 */
class erdos_renyi : public virtual unlayered_unweighted_adjacencylist_network
{
public:
    erdos_renyi() = delete;    
    erdos_renyi(int size, double avg_degree, rng_t &engine);
};

/* Compatibility with miss-spelled previous name */
typedef erdos_renyi erdos_reyni;

//--------------------------------------
//----------FULLY CONNECTED-------------
//--------------------------------------

/**
 * @brief A fully-connected network with random edge order
 */
class fully_connected : public virtual network
    , public virtual network_is_undirected
    , public virtual network_is_simple
    , public virtual network_is_not_layered_and_not_weighted
{
public:
    fully_connected(int size, rng_t &engine);

    virtual node_t nodes() override;

    virtual node_t neighbour(node_t node, int neighbour_index) override;

    virtual index_t outdegree(node_t node) override;

    rng_t &engine;
    std::vector<std::vector<node_t>> neighbours;
};

//--------------------------------------
//---------------TREE-------------------
//--------------------------------------

/**
 * @brief A random acyclic network
 */
class acyclic : public virtual network
    , public virtual network_is_undirected
    , public virtual network_is_simple
    , public virtual network_is_not_layered_and_not_weighted
{
public:
    static double lambda(double mean, int digits);

    acyclic(double avg_degree, bool reduced_root_degree, rng_t &engine);

    virtual node_t neighbour(node_t node, int neighbour_index) override;

    virtual index_t outdegree(node_t node) override;

    rng_t &engine;
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
class config_model : public virtual unlayered_unweighted_adjacencylist_network
{
public:
    config_model(std::vector<int> degreelist, rng_t &engine);

private:
    std::size_t selfloops  = 0;
    std::size_t multiedges = 0;
};

// some specific networks using the configuration model:
// TODO: Make these class methods

std::vector<int> lognormal_degree_list(double mean, double variance, int size, rng_t &engine);
std::vector<int> powerlaw_degree_list(double exponent, int size, rng_t &engine);

//--------------------------------------
//----- CLUSTERED CONFIGURATION MODEL---
//--------------------------------------

/**
 * @brief Network from arbitrary degree distribution.
 *
 * Based on the algorithm described by Serrano & Boguna, 2005.
 */
class config_model_clustered_serrano : public virtual unlayered_unweighted_adjacencylist_network
{
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
    static std::vector<int> triangles_binomial(std::function<double(int)> ck, std::vector<int> degrees, rng_t &engine);

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
    config_model_clustered_serrano(std::vector<int> degrees, std::vector<int> triangles, double beta, rng_t &engine);

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
    config_model_clustered_serrano(std::vector<int> degrees, std::function<double(int)> ck, double beta, rng_t &engine);

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
    config_model_clustered_serrano(std::vector<int> degree, double alpha, double beta, rng_t &engine);

    /**
     * @brief True if the algorithm was unable to satisfied the requested number of triangles
     */
    bool triangles_unsatisfied = false;
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
class config_model_correlated : public virtual unlayered_unweighted_adjacencylist_network {
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
class barabasi_albert : public virtual unlayered_unweighted_adjacencylist_network
{
public:
    barabasi_albert(int size, rng_t &engine, int m = 1);
};

//--------------------------------------
//--------LATTICE-----------------------
//--------------------------------------

template <unsigned int D>
class cubic_lattice : public virtual network
    , public virtual network_embedding
    , public virtual network_is_undirected
    , public virtual network_is_simple
    , public virtual network_is_not_layered_and_not_weighted
{
public:
    const static unsigned int dimension = D;

    typedef std::ptrdiff_t coordinate_t;
    typedef std::array<coordinate_t, dimension> coordinates_t;

    const std::size_t edge_length;
    const std::size_t total_nodes;
    const coordinate_t coordinate_min;
    const coordinate_t coordinate_max;

    cubic_lattice()
        : cubic_lattice(std::floor(std::pow(2.0, std::log2((double)std::numeric_limits<node_t>::max() + 1.0) / D)))
    {
    }

    cubic_lattice(const std::size_t _edge_length)
        : edge_length(_edge_length)
        , total_nodes(std::pow<double>(edge_length, dimension))
        , coordinate_min(-((edge_length - 1) / 2))
        , coordinate_max((edge_length) / 2)
    {
        if (std::pow<double>(edge_length, dimension) > (double)std::numeric_limits<node_t>::max() + 1.0)
            throw std::runtime_error("edge length exceeds maximal supported length");
    }

    virtual ~cubic_lattice() {}

    virtual node_t nodes() override
    {
        return (node_t)total_nodes;
    }

    virtual node_t neighbour(node_t nidx, int neighbour_index) override
    {
        if ((nidx < 0) || ((std::size_t)nidx >= total_nodes))
            return -1;
        // Decode neighbour index into coordinate components
        coordinates_t coords = coordinates(nidx);

        for (std::size_t i = 0; i < dimension; ++i) {
            if (coords[i] == coordinate_min) {
                if (neighbour_index == 0)
                    coords[i] += 1;
                neighbour_index -= 1;
            } else if (coords[i] == coordinate_max) {
                if (neighbour_index == 0)
                    coords[i] -= 1;
                neighbour_index -= 1;
            } else {
                if (neighbour_index == 0)
                    coords[i] += 1;
                else if (neighbour_index == 1)
                    coords[i] -= 1;
                neighbour_index -= 2;
            }
        }
        return node(coords);
    }

    virtual int outdegree(node_t node) override
    {
        if ((node < 0) || ((std::size_t)node >= total_nodes))
            return -1;
        // Decode neighbour index into coordinate components
        coordinates_t c = coordinates(node);
        // Count the number of extremal (i.e min or max) and non-extermal components
        std::size_t e = 0;
        for (std::size_t i = 0; i < dimension; ++i) {
            if ((c[i] != coordinate_max) && (c[i] != coordinate_min))
                continue;
            e++;
        }
        // For ever non-exterminal component there is one neighbour, otherwise two
        return (int)(2 * dimension - e);
    }

    virtual std::size_t dimensionality() override
    {
        return dimension;
    }

    virtual bool coordinates(const node_t node, std::vector<double> &position) override
    {
        if ((node < 0) || ((std::size_t)node >= total_nodes))
            return false;
        const auto c = coordinates(node);
        position.resize(dimension);
        for (std::size_t i = 0; i < dimension; ++i)
            position[i] = c[i];
        return true;
    }

    virtual void bounds(std::vector<double> &a, std::vector<double> &b) override
    {
        a.resize(dimension);
        b.resize(dimension);
        for (std::size_t i = 0; i < dimension; ++i) {
            a[i] = coordinate_min;
            b[i] = coordinate_max;
        }
    }

    /**
     * @brief Convert node index into d-dimensional coordinate vector
     * @param node node index
     * @return d-dimensional coordinate vector as a std::array<coordinate_t, D>
     */
    coordinates_t coordinates(node_t node)
    {
        if ((node < 0) || ((std::size_t)node >= total_nodes))
            throw std::range_error("invalid node index");
        coordinates_t c;
        for (std::size_t i = 0; i < dimension; ++i) {
            c[i] = coordinate_min + (coordinate_t)(node % edge_length);
            node /= edge_length;
        }
        return c;
    }

    /**
     * @brief Convert d-dimensional coordindate vector into node index
     * @param p d-dimensional coordinate vector as a std::array<coordinate_t, D>
     * @return node index
     */
    node_t node(const coordinates_t &c) const
    {
        node_t node = 0;
        for (std::ptrdiff_t i = dimension - 1; i >= 0; --i) {
            if ((c[i] < coordinate_min) || (c[i] > coordinate_max))
                throw std::range_error("node coordinates are out of range");
            node *= edge_length;
            node += (c[i] - coordinate_min);
        }
        return node;
    }
};

template <unsigned int D>
const unsigned int cubic_lattice<D>::dimension;

typedef cubic_lattice<2> cubic_lattice_2d;
typedef cubic_lattice<3> cubic_lattice_3d;
typedef cubic_lattice<4> cubic_lattice_4d;
typedef cubic_lattice<5> cubic_lattice_5d;
typedef cubic_lattice<6> cubic_lattice_6d;
typedef cubic_lattice<7> cubic_lattice_7d;
typedef cubic_lattice<8> cubic_lattice_8d;

//--------------------------------------
//-------- EMPIRICAL NETWORK -----------
//--------------------------------------

/**
 * @brief Weighted empirical network read from a file
 *
 * The file must contain lines of the form
 *
 * <nodeA> <csep> <nodeB1> <csep> <nodeB2> ...
 *
 * where <csep> is ' ' by default. Each such line defines edges A -> B1,
 * A -> B2, ... with weights w1, w2, ... Multiple lines for the same node A
 * are allowed and the edges are merged; this class than thus read
 * adjacencylist and well as edgelist files. Lines starting with '#' are
 * ignored.
 *
 * If undirected is true, the network is assumed to be undirected, i.e.
 * for every edge (a,b) the reverse edge (b,a) is also added. If
 * simplify is true, self-edges and multi-edges are removed.
 *
 */
class empirical_network : public virtual unlayered_unweighted_adjacencylist_network
{
public:
    empirical_network(
        std::istream &file, bool undirected = true, bool simplify = false,
        node_t idxbase = 1, char sep = ' ');
};


//---------------------------------------------------
//-----Compute the reproduction_matrix matrix --------
//---------------------------------------------------

std::vector<std::vector<double>> reproduction_matrix(
    network &nw, double *out_r, double *out_c, double *out_k1, double *out_k2,
    double *out_k3, double *out_m1, double *out_m2,
    double *out_R0, double *out_R_r, double *R_pert);

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
std::vector<double> knn(network &nw);

/**
 * @brief Pearson correlation to measure the assortativity of a network
 *
 * assortativity a = num/den
 * num = sum_kk' [ w(k,k') * k * k' ] - ( sum_k [ w(k) * k ] ) ^ 2
 * den = sum_k [ w(k) * k ^ 2] - ( sum_k [ w(k) * k ] ) ^ 2
 *
 */
double assortativity(network &nw);

/**
 * @brief fraction of non-self links in the network that connect a node of degree k
 * to a node of degree k prime
 *
 */
std::vector<std::vector<double>> Wkk(network &nw);

#if 0
//------------------------------------------
//--ADD DEGREE CORRELATION TO THE NETWORK---
//------------------------------------------
/**
 * @brief Add correlation to the network by rewiring its links.
 *
 */
void add_correlation(double r, adjacencylist_network &nw, rng_t &engine);
#endif
