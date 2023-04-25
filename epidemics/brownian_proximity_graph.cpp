//
//  brownian_proximity_graph.cpp
//  epidemics
//
//  Created by Florian G. Pflug on 14.02.23.
//

#include "brownian_proximity_graph.h"

brownian_proximity_graph::brownian_proximity_graph(node_t N, double R0, double r, double D, double dt, rng_t& engine)
	:size(N)
	,radius(r)
	,length(std::sqrt(size * M_PI * std::pow(radius, 2) / R0))
	,diffusivity(D)
	,delta_t(dt)
	,plength(2*radius)
	,pstride(std::ceil(length / plength))
	,partitions(pstride * pstride)
{
	nodedata.resize(size);
	
	/* Uniformly distribute nodes */
	std::uniform_real_distribution<float> coord(0, length);
	for(node_t i=0; i < size; ++i) {
		const point p = { coord(engine), coord(engine) };
		nodedata[i].position = p;
		partition(p).insert(i);
	}
	
	/* Determine neighbour relationships */
	for(node_t i=0; i < size; ++i) {
		node& n = nodedata[i];
		/* Iterate over partitions that contain possible neighbours */
		const partition_index_t lb = partition_index(point(n.position.x - radius, n.position.y - radius));
		const partition_index_t ub = partition_index(point(n.position.x + radius, n.position.y + radius));
		for(partition_index_t cur = lb; cur.first  < ub.first; ++cur.first) {
			for(; cur.second  < ub.second; ++cur.second) {
				const auto& p = partition(state.cur);
				/* Scan nodes in partition */
				for(auto it = p.begin(), end = p.end(); it != end; ++it) {
					const node_t j = *state.p_it;
					/* Check if node is a neighbour */
					if (distance(n.position, nodedata.at(j).position) <= radius) {
						/* Add to neighbour map and vector */
						auto r = n.neighbour_map.insert(std::make_pair(j, (index_t)n.neighbours.size()));
						n.neighbours.push_back(j);
						assert(r.second);
					}
				}
				/* Re-start scan for next partition */
				state.partition_scan_initialized = false;
			}
		}
	}
}
			
			
node_t brownian_proximity_graph::nodes() {
	return size;
}

node_t brownian_proximity_graph::neighbour(node_t n, int neighbour_index) {
	const node& nn = nodedata.at(n);
	if ((neighbour_index < 0) || (neighbour_index >= nn.neighbours.size()))
		return -1;
	return nn.neighbours[neighbour_index];
}

index_t brownian_proximity_graph::outdegree(node_t node) {
	return (index_t)nodedata.at(node).neighbours.size();
}

absolutetime_t brownian_proximity_graph::next(rng_t& engine)
{
	using std::swap;

	/* If the next event was already determined, report it */
	if (next_event)
		return next_event->time;
	
	/* Otherwise, step time until something happens */
    while (true) {
		/* Move nodes unless already done for the current time step */
		if (!state.displacement_done) {
			std::normal_distribution<float> delta(0, sqrt(2*diffusivity*delta_t));
			for(std::size_t i=0; i < size; ++i) {
				/* Get node */
				node& n = nodedata.at(i);
				/* Displace node */
				point& p = n.position;
				const partition_index_t ppi_old = partition_index(p);
				p.x += std::remainder(delta(engine), 2*length);
				p.y += std::remainder(delta(engine), 2*length);
				/* Enforce reflective boundaries */
				if (p.x < 0)
					p.x = -p.x;
				if (p.x >= length)
					p.x = 2*length - p.x;
				if (p.y < 0)
					p.y = -p.y;
				if (p.y >= length)
					p.y = 2*length - p.y;
				/* Update partitions */
				const partition_index_t ppi_new = partition_index(p);
				if (ppi_old != ppi_new) {
					partition(ppi_old).erase(i);
					partition(ppi_new).insert(i);
				}
			}
			/* Update time, mark time step as completed */
			current_time += delta_t;
			state.displacement_done = true;
		}
		
		/* Update neighbour relationships */
		for(; state.node_i < size; ++state.node_i) {
			node& n = nodedata.at(state.node_i);
			
			/* Remove neighbours that are out of range */
			if (!state.neighbour_scan_done) {
				/* It might seem that neighbour_scan_done is unnecessary, because the for loop would
				 * immediatly exit upon re-entering anway. However, that assumes that the value of the
				 * end interator never changes, which is not guaranteed and probably not true. We thus
				 * do need neighbour_scan_done to ensure we don't re-enter even if the neighbours
				 * set is modified below.
				 */
				if (!state.neighbour_scan_initialized)  {
					state.n_it = n.neighbour_map.begin();
					state.neighbour_scan_initialized = true;
				}
				for(const auto end = n.neighbour_map.end(); state.n_it != end; ++state.n_it) {
					const node_t j = state.n_it->first;
					if (distance(n.position, nodedata.at(j).position) > radius) {
						/* Remove neighbour
						 * First, find its index in the ordered neighbours vector */
						const std::size_t jidx = state.n_it->second;
						/* Swap it with the last element in the neighbours vector,
						 * take care to update that element's index in the map */
						node_t& e = n.neighbours.back();
						n.neighbour_map[e] = jidx;
						swap(n.neighbours[jidx], e);
						/* Finally, remove from neighbours vector */
						n.neighbours.pop_back();
						/* And from the map, moving the iterator forward since we'll skip the loop's ++ */
						state.n_it = n.neighbour_map.erase(state.n_it);
						/* Create event */
						next_event = network_event_t {
							.kind = network_event_kind::neighbour_removed,
							.source_node = state.node_i, .target_node = j, .time = current_time
						};
						/* And report it */
						return next_event->time;
					}
				}
				state.neighbour_scan_done = true;
			}
			
			/* Add neighbours that moved into range
			 * Interate over all partition that contain possible neighbours
			 */
			if (!state.range_scan_initialized) {
				state.lb = partition_index(point(n.position.x - radius, n.position.y - radius));
				state.ub = partition_index(point(n.position.x + radius, n.position.y + radius));
				state.cur = state.lb;
				state.range_scan_initialized = true;
			}
			for(; state.cur.first  < state.ub.first; ++state.cur.first) {
				for(; state.cur.second  < state.ub.second; ++state.cur.second) {
					/* Scan nodes in partition
					 * Keeping the iterator around is safe here because partitions are not
					 * modified while the scan is in progress. Partitions are only modified
					 * during the displacement step
					 */
					const auto& p = partition(state.cur);
					if (!state.partition_scan_initialized) {
						state.p_it = p.begin();
						state.partition_scan_initialized = true;
					}
					for(const auto end = p.end(); state.p_it != end; ++state.p_it) {
						/* Get node */
						const node_t j = *state.p_it;
						/* Check if node is in range */
						if (distance(n.position, nodedata.at(j).position) <= radius) {
							/* Try to add new neighbour to the map, will fail if already a neighbour   */
							auto r = n.neighbour_map.insert(std::make_pair(j, (index_t)n.neighbours.size()));
							if (r.second) {
								/* Found a new neighbour, add to neighbours vector, update index */
								n.neighbours.push_back(j);
								assert(n.neighbours.at(r.first->second) == j);
								/* Continue with the next neighbour upon the next call */
								++state.p_it;
								/* Create event and return */
								next_event = network_event_t {
									.kind = network_event_kind::neighbour_added,
									.source_node = state.node_i, .target_node = j, .time = current_time
								};
								return next_event->time;
							}
						}
					}
					/* Re-start scan for next partition */
					state.partition_scan_initialized = false;
				}
			}
		}
		
		/* Time step completed, reset state */
		state = state_t();
    }
}

std::optional<network_event_t> brownian_proximity_graph::step(rng_t& engine, absolutetime_t max_time) {
	/* Generate next event */
	if (!next_event)
		next(engine);
	
	/* Return if event lies past max_time */
	if (max_time < next_event->time)
		return std::nullopt;
	
	/* Process event */
	const auto ev = next_event;
	next_event = std::nullopt;
	return ev;
}
