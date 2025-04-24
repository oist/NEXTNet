//
//  graph.cpp
//  epidemics
//
//  Created by Cure Samuel Cyrus on 2022/03/07.
//

#include "nextnet/stdafx.h"
#include "nextnet/types.h"
#include "nextnet/utility.h"
#include "nextnet/weighted_network.h"

using boost::math::erfc;

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*------ WEIGHTED NETWORK: BASE CLASS ----------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

weighted_network *as_weighted_network(network *nw)
{
    if (weighted_network *wnw = dynamic_cast<weighted_network *>(nw))
        return wnw->is_unweighted() ? nullptr : wnw;
    else
        return nullptr;
}

bool weighted_network::is_unweighted()
{
    return false;
}

node_t weighted_network::neighbour(node_t node, int neighbour_index)
{
    double dummy;
    return this->neighbour(node, neighbour_index, &dummy);
}

weighted_network::~weighted_network()
{
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*------ NETWORK: WEIGHTED ADJACENCY LIST ------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

bool weighted_adjacencylist_network::is_undirected()
{
	return undirected;
}

bool weighted_adjacencylist_network::is_simple()
{
	return simple;
}

node_t weighted_adjacencylist_network::nodes()
{
    return (node_t)adjacencylist.size();
}

node_t weighted_adjacencylist_network::neighbour(node_t node, int neighbour_index, double *weight)
{
    const auto &n = adjacencylist.at(node);
    if ((neighbour_index < 0) || (n.size() <= (unsigned int)neighbour_index))
        return -1;
    auto neighbour_weight = n[neighbour_index];
    if (weight != nullptr)
        *weight = neighbour_weight.second;
    return neighbour_weight.first;
}

index_t weighted_adjacencylist_network::outdegree(node_t node)
{
    return (index_t)adjacencylist.at(node).size();
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*------ NETWORK: WEIGHTED ERDÖS-RENYI ---------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

weighted_erdos_renyi::weighted_erdos_renyi(int size, double avg_degree,
                                           std::function<double(rng_t &)> weightdist, rng_t &engine)
    :weighted_adjacencylist_network(true, true)
{
    /*--------------Initialisation--------------

     Construct erdos-Reyni graph and for each link we add an infection time:
     => A[i][j] is the time node i (or j) takes to infect node j (or i) once it has itself been infected by someone else.

     */
    const double p = avg_degree / (size - 1); // probability of an edge: if size ->infty and degree-> fixed then we get Poisson Graph.

    /* Using geometric distribution */
    std::geometric_distribution<> skip_edge(p); // comment: equals 0 with prob. p

    adjacencylist.resize(size);
    for (int i = 0; i < size; i++) {
        int s = skip_edge(engine);
        // for very low probabilities p, skip_edge can return very large numbers
        // and so we have to be carefull not to overflow. Hence the slightly weird
        // coding of this loop; we test that we don't skip past the end of the nodes
        // before we attempt to update j, that way j cannot overflow.
        for (int j = -1; (s >= 0) && (s < (i - j - 1)); s = skip_edge(engine)) {
            j += 1 + s;
            assert(j < i);
            const double w = weightdist(engine);
            adjacencylist[i].emplace_back(j, w);
            adjacencylist[j].emplace_back(i, w);
        }
    }
}

weighted_erdos_renyi::weighted_erdos_renyi(int size, double avg_degree,
                                           std::vector<double> w, std::vector<double> p,
                                           rng_t &engine)
    : weighted_erdos_renyi(size, avg_degree, make_dist(w.begin(), w.end(), p.begin(), p.end()), engine)
{
}
