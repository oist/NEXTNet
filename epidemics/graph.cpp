//
//  graph.cpp
//  epidemics
//
//  Created by Cure Samuel Cyrus on 2022/03/07.
//

#include "stdafx.h"
#include "types.h"
#include "graph.h"

erdos_reyni::erdos_reyni(int size, double avg_degree, rng_t& engine ){
    
    /*--------------Initialisation--------------

     Construct erdos-Reyni graph and for each link we add an infection time:
     => A[i][j] is the time node i (or j) takes to infect node j (or i) once it has itself been infected by someone else.
     
     */
    
    const double p = avg_degree/size; // probability of an edge: if size ->infty and degree-> fixed then we get Poisson Graph.
    
    /* Using geometric distribution */
    std::geometric_distribution<> skip_edge(p);// comment: equals 0 with prob. p

    neighbours.resize(size);
    for (int i=0; i<size; i++) {
        for (int j=skip_edge(engine); j<i; j += 1 + skip_edge(engine)) {
            neighbours[i].push_back(j);
            neighbours[j].push_back(i);
        }
    }
}

node_t erdos_reyni::neighbour(node_t node, int neighbour_index) {
    const auto& n = neighbours.at(node);
    if ((neighbour_index < 0) || (n.size() <= (unsigned int)neighbour_index))
        return -1;
    //std::cout << node<< "neighbour: "<< n[neighbour_index]<<"\n";
    return n[neighbour_index];
}

index_t erdos_reyni::outdegree(node_t node) {
    return (index_t) neighbours.at(node).size();
}

fully_connected::fully_connected(int size, rng_t& engine_)
    :engine(engine_)
    ,neighbours(size, std::vector<node_t>({}))
{}

node_t fully_connected::neighbour(node_t node, int neighbour_index) {
    // index has to be in the range [0, size - 1]
    if ((neighbour_index < 0) || (neighbour_index >= neighbours.size() - 1))
        return -1;
    // get neighbours of node <node>
    const auto& n = neighbours.at(node);
    if (!n.empty())
        return n.at(neighbour_index);
    // neighbours not yet generated, generate it now
    std::vector<node_t> np;
    np.reserve(neighbours.size() - 1);
    // add nodes [0,...,node-1,node+1,...,size] as neighbours
    for(int i=0; i < node; ++i)
        np.push_back(i);
    for(int i=node+1, m=neighbours.size(); i < m; ++i)
        np.push_back(i);
    // shuffle
    std::shuffle(np.begin(), np.end(), engine);
    // store and return
    neighbours[node] = np;        
    return np.at(neighbour_index);
}

index_t fully_connected::outdegree(node_t node) {
    // network is fully connected, every node is every node's neighbour
    return neighbours.size() - 1;
}
