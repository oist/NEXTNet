#pragma once

#include "stdafx.h"
#include "types.h"
#include "graph.h"

class simulate_next_reaction {
public:
    graph& network;
    transmission_time& psi;
    std::unordered_set<node_t> infected;
    
    simulate_next_reaction(class graph& nw, class transmission_time& psi_)
        :network(nw), psi(psi_)
    {}
    
    void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v);
    
    std::pair<node_t, absolutetime_t> step(rng_t& engine);
    
private:
    struct active_edges_entry {
        /*
         * Absolute time of infection
         */
        absolutetime_t time;
        
        /*
         * Node that is (putatively) infected
         */
        node_t node;
        
        /*
         * Source node that causes the node's infection, it's
         * original infection time, and the node's neighbour
         * index within the souce node
         */
        absolutetime_t source_time = INFINITY;
        node_t source_node = -1;
        index_t neighbour_index = -1;
        index_t neighbours_remaining = 0;
        
        bool operator< (const active_edges_entry& o) const { return time < o.time; }
        bool operator<= (const active_edges_entry& o) const { return time <= o.time; }
        bool operator== (const active_edges_entry& o) const { return time == o.time; }
        bool operator!= (const active_edges_entry& o) const { return time != o.time; }
        bool operator>= (const active_edges_entry& o) const { return time >= o.time; }
        bool operator> (const active_edges_entry& o) const { return time > o.time; }
    };
    
    std::priority_queue<active_edges_entry, std::deque<active_edges_entry>,
                        std::greater<active_edges_entry>>
      active_edges;
};
