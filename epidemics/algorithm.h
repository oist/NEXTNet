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

struct simulation_algorithm {
    virtual ~simulation_algorithm() {};
    
    virtual graph& get_network() = 0;

    virtual class transmission_time& transmission_time() = 0;

    virtual class transmission_time* reset_time() = 0;

    virtual std::optional<event_t> step(rng_t& engine) = 0;

    virtual void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v) = 0;

    virtual bool is_infected(node_t) = 0;
};
