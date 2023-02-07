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

/**
 * TODO: Rename to epidemic_simulation_algorithm
 */
struct simulation_algorithm {
    virtual ~simulation_algorithm() {};
    
    virtual graph& get_network() const = 0;

    virtual const class transmission_time& transmission_time() const = 0;

    virtual const class transmission_time* reset_time() const = 0;
	
	virtual absolutetime_t next() = 0;

    virtual std::optional<event_t> step(rng_t& engine, absolutetime_t nexttime = NAN,
										event_filter_t event_filter = std::nullopt) = 0;
	
	virtual void notify_infected_node_neighbour_added(network_event_t event);
	
    virtual void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v) = 0;

    virtual bool is_infected(node_t) const = 0;
};

struct epidemic_on_dynamic_network_simulation {
	absolutetime_t next();

	std::optional<network_or_epidemic_event_t> step(rng_t& engine, absolutetime_t nexttime = NAN) ;
	
	bool simulation_event_filter(event_t ev);
	
	dynamic_network* network;
	simulation_algorithm* simulation;
	
	absolutetime_t network_next = NAN;
	absolutetime_t simulation_next = NAN;
	
	typedef std::unordered_map<node_t, bool> neighbour_state_t;
	std::unordered_map<node_t, neighbour_state_t> infected_neighbour_state;
};
