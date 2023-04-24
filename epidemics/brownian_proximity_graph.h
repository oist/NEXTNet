//
//  brownian_proximity_graph.hpp
//  epidemics
//
//  Created by Florian G. Pflug on 24.04.23.
//

#pragma once

#include "stdafx.h"
#include "types.h"
#include "random.h"
#include "dynamic_graph.h"

//--------------------------------------
//----------DYNAMIC NETWORKS------------
//--------------------------------------

struct brownian_proximity_graph : virtual dynamic_network, virtual erdos_reyni {
	using bg = boost::geometry;

	brownian_proximity_graph(int N, double R0, double r, double D, rng_t& engine);

	virtual absolutetime_t next(rng_t& engine);

	virtual std::optional<network_event_t> step(rng_t& engine, absolutetime_t nexttime = NAN);

	const double length;
	const double subdiv_length;
	const int subdivs;
	const double diffusivity;
	const double contact_radius;
	const int size;
	const double dt;

	struct point {
		double x, y;
	};

	typedef std::vector<point> points_t;
	std::vector<points_t> subdivs_t;
};
