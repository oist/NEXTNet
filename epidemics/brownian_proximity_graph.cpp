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
	,node_index(N)
{
	/* Uniformly distribute nodes */
	std::uniform_real_distribution<float> coord(0, length);
	for(node_t i=0; i < size; ++i) {
		/* Create node */
		const point p = { coord(engine), coord(engine) };
		node_data node_tmp;
		node_tmp.index = i;
		node_tmp.position = p;
		const partition_index_t pi = partition_index(p);
		
		/* Insert node into partition and node index */
		node_vector_t& pn = partition(pi);
		node_index[i] = std::make_pair(pi, pn.size());
		pn.push_back(node_tmp);
	}
	
	/* Determine neighbour relationships */
	/* Interare over all nodes, partition-wise */
	for(node_vector_t& pn: partitions) {
		for(node_data& n1: pn) {
			/* Iterate over partitions that contain possible neighbours */
			const partition_index_t lb = partition_index(point(n1.position.x - radius, n1.position.y - radius));
			const partition_index_t ub = partition_index(point(n1.position.x + radius, n1.position.y + radius));
			for(partition_index_t cur = lb; cur.first  <= ub.first; ++cur.first) {
				for(cur.second = lb.second; cur.second  <= ub.second; ++cur.second) {
					const auto& p = partition(cur);
					/* Scan nodes in partition */
					for(const node_data& n2: p) {
						/* Check if node is a neighbour */
						if (distance(n1.position, n2.position) <= radius) {
							/* Add to neighbour map and vector */
							auto r = n1.neighbour_map.insert(std::make_pair(n2.index, (index_t)n1.neighbours.size()));
							n1.neighbours.push_back(n2.index);
							assert(r.second);
							_unused(r);
						}
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
	if ((n < 0) || (n >= node_index.size()))
		return -1;
	const node_data& nn = node(n);
	if ((neighbour_index < 0) || (neighbour_index >= (int) nn.neighbours.size()))
		return -1;
	return nn.neighbours[neighbour_index];
}

index_t brownian_proximity_graph::outdegree(node_t n) {
	if ((n < 0) || (n >= node_index.size()))
		return -1;
	const node_data& nn = node(n);
	return (index_t)nn.neighbours.size();
}

std::size_t brownian_proximity_graph::dimensionality() {
	return 2;
}

void brownian_proximity_graph::bounds(std::vector<double>& a, std::vector<double>& b) {
	a.resize(2);
	b.resize(2);
	a[0] = b[0] = 0.0;
	a[1] = b[1] = length;
}

bool brownian_proximity_graph::coordinates(const node_t n, std::vector<double>& position) {
	if ((n < 0) || (n >= node_index.size()))
		return -1;
	const node_data& nn = node(n);
    position.resize(2);
    position[0] = nn.position.x;
    position[1] = nn.position.y;
    return true;
}

void brownian_proximity_graph:: move_node(node_data& n, partition_index_t pi_old, partition_index_t pi_new)
{
	/* Add to new partition and update node index */
	node_vector_t& p_new = partition(pi_new);
	const partition_node_index_t pin_new(pi_new, p_new.size());
	node_index[n.index] = pin_new;
	p_new.emplace_back(std::move(n));
	/* Find index in old partition
	 * Note: n has been moved from already!
	 */
	node_vector_t& p_old = partition(pi_old);
	const std::ptrdiff_t ni = &n - &p_old.front();
	assert((0 <= ni) && (ni < p_old.size()));
	/* Move last node of current partition to positing of removed node */
	node_data& last = p_old.back();
	if (&n != &last) {
		using std::swap;
		const partition_node_index_t pin_old(pi_old, ni);
		assert(node_index[last.index].first == pi_old);
		node_index[last.index] = pin_old;
		n = std::move(last);
	}
	/* Finally, remove from old partition */
	p_old.pop_back();
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
			for(node_vector_t& p: partitions) {
				for(std::size_t i=0; i < p.size();) {
					/* Get node */
					node_data& n = p[i];
					/* Skip nodes which were already displaced & moved into this partition */
					if (n.generation == current_generation + 1) {
						++i;
						continue;
					}
					assert(n.generation == current_generation);
					/* Displace node, increment generation counter */
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
					/* Update generation to prevent further moves during this round */
					n.generation++;
					/* Update partitions */
					const partition_index_t ppi_new = partition_index(p);
					if (ppi_old != ppi_new) {
						move_node(n, ppi_old, ppi_new);
						/* Continue with same node index since current node was moved */
					} else {
						/* Continue with next node index */
						++i;
					}
				}
			}
			/* Update time, mark time step as completed */
			current_time += delta_t;
			current_generation++;
			state.displacement_done = true;
		}
		
		/* Update neighbour relationships */
		/* Scan nodes, partition-wise */
		for(; state.partition_i < partitions.size(); ++state.partition_i) {
			node_vector_t p = partitions[state.partition_i];
			if (!state.outer_partition_scan_initialized) {
				state.outer_partition_node_i = 0;
				state.outer_partition_scan_initialized = true;
			}
			for(; state.outer_partition_node_i < p.size(); ++state.outer_partition_node_i) {
				/* Current node */
				node_data& n = p[state.outer_partition_node_i];
				
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
						const node_data& n2 = node(n.neighbours[state.neighbour_idx]);
						assert(n2.index == n.neighbours[state.neighbour_idx]);
						if (distance(n.position, n2.position) > radius) {
							/* Create event */
							next_event = network_event_t {
								.kind = network_event_kind::neighbour_removed,
								.source_node = n.index, .target_node = n2.index, .time = current_time
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
						if (!state.inner_partition_scan_initialized) {
							state.inner_partition_node_i = 0;
							state.inner_partition_scan_initialized = true;
						}
						for(; state.inner_partition_node_i < p.size(); ++state.inner_partition_node_i) {
							/* Get putative neighbour */
							const node_data& n2 = p[state.inner_partition_node_i];
							/* Check if node is in range */
							if ((&n != &n2) && (distance(n.position, n2.position) <= radius)) {
								/* Check if node is already a neigbhour, if not report new neighbour */
								if (n.neighbour_map.find(n2.index) == n.neighbour_map.end()) {
									/* Create event and return */
									next_event = network_event_t {
										.kind = network_event_kind::neighbour_added,
										.source_node = n.index, .target_node = n2.index, .time = current_time
									};
									return next_event->time;
								}
							}
						}
						/* Re-start scan for next partition */
						state.inner_partition_scan_initialized = false;
					}
					/* Re-start inner loop */
					state.cur.second = state.lb.second;
				}
				state.range_scan_initialized = false;
				
				/* Re-de neighbour scan for next node */
				state.neighbour_scan_done = false;
			}
			/* Re-start scan for next partition */
			state.outer_partition_scan_initialized = false;
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
			const partition_node_index_t pin(partition_index(state.partition_i), state.outer_partition_node_i);
			node_data& n = node(pin);
			assert(n.index == ev.source_node);
			/* Get putative neighbour */
			const node_vector_t p = partition(state.cur);
			const node_data& n2 = p[state.inner_partition_node_i];
			assert(state.range_scan_initialized);
			assert(state.inner_partition_scan_initialized);
			/* Add neighbour */
			auto r = n.neighbour_map.insert(std::make_pair(n2.index, (index_t)n.neighbours.size()));
			assert(r.second);
			n.neighbours.push_back(n2.index);
			assert(n.neighbours.at(r.first->second) == n2.index);
			/* Continue with the next putative neighbour upon the next call */
			++state.inner_partition_node_i;
			_unused(r);
			break;
		}
		case network_event_kind::neighbour_removed: {
			/* Get node */
			const partition_node_index_t pin(partition_index(state.partition_i), state.outer_partition_node_i);
			node_data& n = node(pin);
			assert(n.index == ev.source_node);
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
