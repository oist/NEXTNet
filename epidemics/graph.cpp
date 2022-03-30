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
    
//    std::vector<std::vector<node_t>> neighbours;
    
    const double p = avg_degree/size; // probability of an edge: if size ->infty and degree-> fixed then we get Poisson Graph.
    
    
     //Using naive Bernoulli implementation

//    std::bernoulli_distribution has_edge(p);
//
//    neighbours.resize(size);
//    for (int i=0; i<size; i++) {
//        for (int j=0; j<i; j++) {
//            if (!has_edge(engine)) {
//                continue;
//            }
//            const double tau = infection_distribution.sample(engine);
//            neighbours[i].push_back(std::make_pair(j, tau));
//            neighbours[j].push_back(std::make_pair(i, tau));
//        }
//    }
//
    /* Using geometric distribution */
    std::geometric_distribution<> skip_edge(p);// comment: equals 0 with prob. p

    neighbours.resize(size);
    for (int i=0; i<size; i++) {
        for (int j=skip_edge(engine); j<i; j += 1 + skip_edge(engine)) {
            neighbours[i].push_back(i);
            neighbours[j].push_back(j);
        }
    }
    
    for (int i =0; i< size; i++) {
        sort(neighbours[i].begin(), neighbours[i].end(),
             [](const auto& a, const auto& b) { return a.second < b.second; });
    }
}

node_t erdos_reyni::neighbour(node_t node, int neighbour_index) {
    const auto& n = neighbours.at(node);
    if ((neighbour_index < 0) || (n.size() <= (unsigned int)neighbour_index))
        return -1;
    return n[neighbour_index];
}

int erdos_reyni::outdegree(node_t node) {
    return neighbours.at(node).size();
}
