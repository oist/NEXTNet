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
    virtual std::pair<node_t, interval_t> neighbour(node_t node, int neighbour_index) = 0;
};


class erdos_reyni : public graph {
public:
    erdos_reyni(int size, double avg_degree, const beta& infection_distribution, rng_t& engine);

    virtual std::pair<node_t, interval_t> neighbour(node_t node, int neighbour_index);

//private:
    /* Adjacency list of the graph */
    std::vector<std::vector<std::pair<node_t,interval_t>>> neighbours;
};
