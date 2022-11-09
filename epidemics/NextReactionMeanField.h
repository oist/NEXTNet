#pragma once

#include "stdafx.h"
#include "types.h"
#include "algorithm.h"
#include "random.h"

class simulate_next_reaction_mean_field : public simulation_algorithm {
public:
    simulate_next_reaction_mean_field(int N_, double R0_, const class transmission_time& psi_,
                                      const class transmission_time* rho_ = nullptr)
        : N(N_), R0(R0_), psi(psi_), rho(rho_)
    {}
    
    virtual const class transmission_time& transmission_time() const;

    virtual const class transmission_time* reset_time() const;

    virtual void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v);
    
    virtual std::optional<event_t> step(rng_t& engine);

    virtual bool is_infected(node_t) const;

    const int N;
    
    const double R0;
    
    virtual graph& get_network() const;
    
private:
    const double p = R0 / (N - 1);
    const class transmission_time& psi;
    const class transmission_time* rho = nullptr;
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
