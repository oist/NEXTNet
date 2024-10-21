#pragma once

#include "stdafx.h"
#include "types.h"
#include "permutation.h"
#include "algorithm.h"
#include "random.h"
#include "graph.h"

class simulate_next_reaction : public simulation_algorithm {
public:
	struct params {
		params() noexcept {} ;
		
		bool shuffle_neighbours = true;
		bool edges_concurrent = true;
		bool SIR = false;
	};
	
	simulate_next_reaction(graph& nw, const class transmission_time& psi_,
                           const class transmission_time* rho_ = nullptr,
						   params p_ = params())
		:network(nw), psi(psi_), rho(rho_),
		 p(p_), shuffle_neighbours(p.shuffle_neighbours && !p.edges_concurrent)
    {}

    virtual graph& get_network() const;

    virtual const class transmission_time& transmission_time() const;

    virtual const class transmission_time* reset_time() const;

    virtual void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v);
    
    // thermal_infections manually marks the nodes as infected and manually adds the active edges only if the firing time is later than thermal.
    virtual void add_thermal_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v,double Lambda,rng_t& engine);
    
	virtual absolutetime_t next(rng_t& engine);

	virtual std::optional<event_t> step(rng_t& engine, absolutetime_t maxtime = INFINITY,
										event_filter_t event_filter = std::nullopt);
	
	virtual void notify_infected_node_neighbour_added(network_event_t event, rng_t& engine);

	virtual void notify_infected_contact(network_event_t event, rng_t& engine);

    virtual bool is_infected(node_t) const;

    struct infected_state_t {
        infected_state_t(absolutetime_t inf, absolutetime_t res)
            :infection_time(inf), reset_time(res)
        {}

        absolutetime_t infection_time;
        absolutetime_t reset_time;
    };
    typedef std::unordered_map<node_t, infected_state_t> infected_nodes_t;
    infected_nodes_t infected;
    
    graph& network;
    const class transmission_time& psi;
    const class transmission_time* rho = nullptr;
	const params p;
	const bool shuffle_neighbours;

	int removed = 0; // number of nodes that have recovered and cannot be re infected. (only active in the SIR case).

    int current_nb_of_infected(){
        return (int) infected.size() - removed;
    }

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
    
#if NEXT_REACTION_QUEUE == STD_PRIORITY_QUEUE_DEQUE
    std::priority_queue<active_edges_entry, std::deque<active_edges_entry>,
                        std::greater<active_edges_entry>>
      active_edges;
	
	const active_edges_entry& top_edge() { return active_edges.top(); };

	void pop_edge() { active_edges.pop(); };
	
	void push_edge(active_edges_entry e) { active_edges.push(e); };

#elif NEXT_REACTION_QUEUE == EXT_PRIO_QUEUE
	rollbear::prio_queue<32, absolutetime_t, active_edges_entry>
	  active_edges;
	
	const active_edges_entry& top_edge() { return active_edges.top().second; };

	void pop_edge() { active_edges.pop(); };
	
	void push_edge(active_edges_entry e) { active_edges.push(e.time, e); };
#endif

    std::size_t queue_steps_total = 0;
	
	std::optional<event_t> step_infection(const active_edges_entry& next, event_filter_t evf, rng_t& engine);
	
	std::optional<event_t> step_reset(const active_edges_entry& next, event_filter_t evf, rng_t& engine);
};

