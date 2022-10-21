//
//  types.hpp
//  epidemics
//
//  Created by Florian G. Pflug on 2022/03/08.
//

#pragma once

#include "stdafx.h"

typedef std::mt19937 rng_t;
typedef int node_t;
typedef int index_t;
typedef double interval_t;
typedef double absolutetime_t;
typedef std::pair<node_t,node_t> edge_t;

enum class event_kind {
    none = 0, infection = 1, reset = 2
};

struct event_t {
    event_kind kind;
    node_t node;
    absolutetime_t time;
};
