#pragma once

#include "stdafx.h"
#include "types.h"
#include "permutation.h"
#include "algorithm.h"
#include "random.h"
#include "graph.h"

class simulate_next_reaction : public simulation_algorithm {
public:
    simulate_next_reaction(graph& nw, const class transmission_time& psi_,
                           const class transmission_time* rho_ = nullptr,
                           bool shuffle_neighbours_ = true)
        :network(nw), psi(psi_), rho(rho_), shuffle_neighbours(shuffle_neighbours_)
    {}

    virtual graph& get_network() const;

    virtual const class transmission_time& transmission_time() const;

    virtual const class transmission_time* reset_time() const;

    virtual void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v);
    
    virtual std::optional<event_t> step(rng_t& engine);

    virtual bool is_infected(node_t) const;
    
private:
    graph& network;
    const class transmission_time& psi;
    const class transmission_time* rho = nullptr;
    bool shuffle_neighbours = true;
    std::unordered_set<node_t> infected;
    
    struct active_edges_entry {
        /*
         * Event kind represented by this edge (infection or reset)
         * Reset events are self-loops, a consequently they obey
         * source_node=-1, neighbour_index=-1, neighbours_remaining=0.
         */
        event_kind kind;

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
         * original infection time, it's reset time, and the
         * node's neighbour index within the souce node
         */
        absolutetime_t source_time = INFINITY;
        node_t source_node = -1;
        absolutetime_t source_reset = INFINITY;
        permutation<node_t> source_permutation;
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
	
	std::optional<event_t> step_infection(const active_edges_entry& next, rng_t& engine);
	
	std::optional<event_t> step_reset(const active_edges_entry& next, rng_t& engine);
};
