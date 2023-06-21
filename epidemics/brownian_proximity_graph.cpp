//
//  brownian_proximity_graph.cpp
//  epidemics
//
//  Created by Florian G. Pflug on 14.02.23.
//

#include "brownian_proximity_graph.h"

using namespace std::literals;

brownian_proximity_graph::brownian_proximity_graph(node_t N, double avg_degree, double r, double D, rng_t& engine)
	:brownian_proximity_graph(N, avg_degree, r, D, std::pow(r, 2) / (100 * 2 * D),engine)
{}

brownian_proximity_graph::brownian_proximity_graph(node_t N, double avg_degree, double r, double D, double dt, rng_t& engine)
	:size(N)
	,radius(r)
	,length(std::sqrt((double)size * M_PI * std::pow(radius, 2) / avg_degree))
	,diffusivity(D)
	,delta_t(dt)
	,current_time(0.0)
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
		for(partition_index_t cur = lb; cur.first  <= ub.first; ++cur.first) {
			for(cur.second = lb.second; cur.second  <= ub.second; ++cur.second) {
				const auto& p = partition(cur);
				/* Scan nodes in partition */
				for(auto it = p.begin(), end = p.end(); it != end; ++it) {
					const node_t j = *it;
					/* Check if node is a neighbour */
					if (distance(n.position, nodedata.at(j).position) <= radius) {
						/* Add to neighbour map and vector */
						auto r = n.neighbour_map.insert(std::make_pair(j, (index_t)n.neighbours.size()));
						n.neighbours.push_back(j);
						assert(r.second);
						_unused(r);
					}
				}
			}
		}
	}
}
			
brownian_proximity_graph::~brownian_proximity_graph()
{
}
			
node_t brownian_proximity_graph::nodes() {
	return size;
}

node_t brownian_proximity_graph::neighbour(node_t n, int neighbour_index) {
	const node& nn = nodedata.at(n);
	if ((neighbour_index < 0) || (neighbour_index >= (int) nn.neighbours.size()))
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
			for(std::size_t i=0; i < (size_t) size; ++i) {
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
					state.neighbour_idx = 0;
					state.neighbour_scan_initialized = true;
				}
				for(; state.neighbour_idx < (int) n.neighbours.size(); ++state.neighbour_idx) {
					const node_t j = n.neighbours[state.neighbour_idx];
					if (distance(n.position, nodedata.at(j).position) > radius) {
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
				state.neighbour_scan_initialized = false;
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
			for(; state.cur.first <= state.ub.first; ++state.cur.first) {
				for(; state.cur.second  <= state.ub.second; ++state.cur.second) {
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
						/* Get putative neighbour */
						const node_t j = *state.p_it;
						/* Check if node is in range */
						if ((state.node_i != j) && (distance(n.position, nodedata.at(j).position) <= radius)) {
							/* Check if node is already a neigbhour, if not report new neighbour */
							if (n.neighbour_map.find(j) == n.neighbour_map.end()) {
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
				/* Re-start inner loop */
				state.cur.second = state.lb.second;
			}
			state.range_scan_initialized = false;

			/* Re-de neighbour scan for next node */
			state.neighbour_scan_done = false;
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
	const auto ev = *next_event;
	switch (ev.kind) {
		case network_event_kind::neighbour_added: {
			/* Get node */
			node& n = nodedata.at(state.node_i);
			assert(state.node_i == ev.source_node);
			/* Get putative neighbour */
			const node_t j = *state.p_it;
			assert(state.range_scan_initialized);
			assert(state.partition_scan_initialized);
			/* Add neighbour */
			auto r = n.neighbour_map.insert(std::make_pair(j, (index_t)n.neighbours.size()));
			assert(r.second);
			n.neighbours.push_back(j);
			assert(n.neighbours.at(r.first->second) == j);
			/* Continue with the next putative neighbour upon the next call */
			++state.p_it;
			_unused(r);
			break;
		}
		case network_event_kind::neighbour_removed: {
			/* Get node */
			node& n = nodedata.at(state.node_i);
			assert(state.node_i == ev.source_node);
			assert(state.neighbour_scan_initialized);
			assert(!state.neighbour_scan_done);
			/* Find neighbour in neighbour map */
			node_t& j = n.neighbours[state.neighbour_idx];
			assert(j == ev.target_node);
			auto nmap_it = n.neighbour_map.find(j);
			/* Prepare to swap with last element in the neighbour list */
			node_t& e = n.neighbours.back();
			/* Update neighbour map with new indices */
			n.neighbour_map[e] = state.neighbour_idx;
			n.neighbour_map.erase(nmap_it);
			/* Swap elements in neighbour list and truncate */
			using std::swap;
			swap(j, e);
			n.neighbours.pop_back();
			break;
		}
		default:
			throw std::logic_error("unknown network even type"s + name(ev.kind));
	}
	
	
	next_event = std::nullopt;
	return ev;
}
