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

class graph {
public:    
    /*
     * returns:
     * - the ith neighbour of node (in order of infection times)
     * - the infection time from node to the ith neigbhour of node.
     */
    virtual node_t neighbour(node_t node, int neighbour_index) = 0;

    /**
     * Returns the number of (outgoing) edges of the given node
     */
    virtual int outdegree(node_t node) = 0;
};


class erdos_reyni : public graph {
public:
    erdos_reyni(int size, double avg_degree, rng_t& engine);

    virtual node_t neighbour(node_t node, int neighbour_index);

    virtual int outdegree(node_t node);

//private:
    /* Adjacency list of the graph */
    std::vector<std::vector<node_t>> neighbours;
};
