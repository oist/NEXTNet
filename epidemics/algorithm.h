//
//  algorithm.h
//  epidemics_update
//
//  Created by Cure Samuel Cyrus on 2022/02/19.
//

#pragma once

#include "stdafx.h"
#include "types.h"

#include "random.h"
#include "graph.h"
#include "dynamic_graph.h"

/**
 * TODO: Rename to epidemic_simulation_algorithm
 */
struct simulation_algorithm {
    virtual ~simulation_algorithm() {};
    
    virtual graph& get_network() const = 0;

    virtual const class transmission_time& transmission_time() const = 0;

    virtual const class transmission_time* reset_time() const = 0;
	
	virtual absolutetime_t next(rng_t& engine) = 0;

	virtual std::optional<event_t> step(rng_t& engine, absolutetime_t maxtime = INFINITY,
										event_filter_t event_filter = std::nullopt) = 0;
	
	virtual void notify_infected_node_neighbour_added(network_event_t event, rng_t& engine) = 0;

	virtual void notify_infected_contact(network_event_t event, rng_t& engine);

    virtual void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v) = 0;

    virtual bool is_infected(node_t) const = 0;
};

struct simulate_on_dynamic_network {
	simulate_on_dynamic_network(simulation_algorithm& sim);
	
	absolutetime_t next(rng_t& engine);

	std::optional<network_or_epidemic_event_t> step(rng_t& engine, absolutetime_t maxtime = INFINITY) ;
	
	bool simulation_event_filter(event_t ev);
	
	dynamic_network* network;
	simulation_algorithm& simulation;
	
	enum class neighbour_state_t : unsigned char {
		admissible = 0,
		masked = 1,
		transmitted = 2
	};
	
	typedef std::unordered_map<node_t, neighbour_state_t> neighbours_states_t;
	std::unordered_map<node_t, neighbours_states_t> infected_neighbour_state;
};
