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

//--------------------------------------
//--------IMPORTED NETWORK----------
//--------------------------------------

weighted_empirical_network::weighted_empirical_network(
    std::istream& file, bool undirected, bool simplify,
    node_t idxbase, char csep, char wsep)
	: weighted_adjacencylist_network(undirected, true)
{
	// Read adjacencylist / edgelist file and create edge multi-set
	const bool csep_is_space = std::isspace(csep);
	const bool wsep_is_space = std::isspace(wsep);
	std::unordered_map<edge_t, double, pair_hash> edges = {};
    std::string line;
    node_t size = 0;
    std::size_t l = 0;
    while (std::getline(file, line)) {
        // Skip empty and commented lined
        if ((line.size() == 0) || (line[0] == '#'))
            continue;
        
        // Read line from string buffer
		std::stringstream is(line);
        ++l;
		
		// Read source node
    	node_t n1; 
    	if (!(is >> n1)) {
            // Skip first non-commented line if it looks like a header
            if (l == 1) continue;
    		throw std::runtime_error("invalid line: " + line);
        }
        if (n1 <= 0)
			throw std::runtime_error("invalid node " + std::to_string(n1));
    	
        // Translate to zero-based node indices and keep track of size
        if (n1 < idxbase)
			throw std::runtime_error("invalid node " + std::to_string(n1));
        n1 -= idxbase;
        size = std::max(size, n1 + 1);
        
    	// Read list of outgoing edges separated by csep
    	while (true) {
    		// Consume sep unless it's whitespace, break if end of line
    		char c = csep;
			if (!csep_is_space && !(is >> c)) break;
			if (csep != c)
				throw std::runtime_error("invalid line: " + line);
    		
    		// Read next target node, break if end of line
    		node_t n2;
    		if (!(is >> n2)) {
                if (csep_is_space) break;
                // There was a non-ws csep but no following node, complain
                throw std::runtime_error("invalid line: " + line);
    		}
            
            // Read wsep and next weight, can't hit end of line here
    		c = wsep;
			if (!wsep_is_space && !(is >> c)) break;
			if (wsep != c)
				throw std::runtime_error("invalid line: " + line);
            double w;
            if (!(is >> w))
				throw std::runtime_error("invalid line: " + line);
                
            // Translate to zero-based node indices and keep track of size
            if (n2 < idxbase)
				throw std::runtime_error("invalid node " + std::to_string(n2));
            n2 -= idxbase;
            size = std::max(size, n2 + 1);
            
            // Validate weight
            if ((w < 0) || !std::isfinite(w))
				throw std::runtime_error("invalid weight " + std::to_string(w));
                        
			// Collect edge
			const edge_t e = (undirected
			                  ? edge_t { std::min(n1, n2), std::max(n1, n2) }
			                  : edge_t { n1, n2 } );
			edges[e] += w;
    	}
    	
    	// Must read whole line
    	if (!is.eof())
    		throw std::runtime_error("invalid line: " + line);
	}
	
	// Convert edges into adjacency list
    adjacencylist.resize(size);
	for(const auto& v: edges) {
		edge_t e = v.first;
		
		// Skip self-edges if a simple graph was requested
		const bool selfedge = (e.first == e.second);
		if (selfedge && simplify)
			continue;
		
		// Keep track of whether the network is simple;
		simple = simple && !selfedge;
		
		// For undirected graphs, add both forward- and reverse-edges
		for(int r=0; r < ((undirected && !selfedge) ? 2 : 1 ); ++r) {
			if (r) std::swap(e.first, e.second);
			
			// Get adjacencylist for source node
			auto& al = adjacencylist[e.first];
			
			// Retain mulitplicity of edge unless simplify is true
			al.emplace_back(e.second, v.second);
		}
	}
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
