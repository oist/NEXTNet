#pragma once

#include "stdafx.h"
#include "types.h"
#include "graph.h"

//--------------------------------------
//-------------SIMULATION---------------
//--------------------------------------

class simulate_next_reaction {
public:
    graph& network;
    std::unordered_set<node_t> infected;
    
    simulate_next_reaction(class graph& nw)
        :network(nw)
    {}
    
    void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v);
    
    std::pair<node_t, absolutetime_t> step();
    
private:
    struct infectiontimes_entry {
        /*
         * Absolute time of infection
         */
        double time;
        
        /*
         * Node that is (putatively) infected
         */
        node_t node;
        
        /*
         * Source node that causes the node's infection, it's
         * original infection time, and the node's neighbour
         * index within the souce node
         */
        double source_time = INFINITY;
        node_t source_node = -1;
        int neighbour_index = -1;
        
        bool operator< (const infectiontimes_entry& o) const { return time < o.time; }
        bool operator<= (const infectiontimes_entry& o) const { return time <= o.time; }
        bool operator== (const infectiontimes_entry& o) const { return time == o.time; }
        bool operator!= (const infectiontimes_entry& o) const { return time != o.time; }
        bool operator>= (const infectiontimes_entry& o) const { return time >= o.time; }
        bool operator> (const infectiontimes_entry& o) const { return time > o.time; }
    };
    
    std::priority_queue<infectiontimes_entry, std::deque<infectiontimes_entry>,
                        std::greater<infectiontimes_entry>>
      infectiontimes;
};




/* Simulates path */

std::vector<double> simulatePath(std::vector<double>& infection_times, int n_max, const lognormal_beta& infection_distribution, rng_t engine);





void print_matrix(std::vector<std::vector<double>>& A);

void simulateManyPaths(int nb_paths, rng_t& engine);
