//
//  dynamic_graph.cpp
//  epidemics
//
//  Created by Florian G. Pflug on 14.02.23.
//

#include "dynamic_graph.h"

void dynamic_network::notify_epidemic_event(event_t ev, rng_t& engine) {
	/* Do nothing by default */
}

dynamic_erdos_reyni::dynamic_erdos_reyni(int size, double avg_degree, double timescale, rng_t& engine)
	:erdos_reyni(size, avg_degree, engine)
	,edge_probability(avg_degree / (size - 1))
	,alpha(edge_probability * timescale)
	,beta((1.0 - edge_probability) * timescale)
	,edges_present(0)
{
	/* Initial degree-weights node distribution and present/absent edge counters */
	for(node_t i = 0; (std::size_t)i < this->adjacencylist.size(); ++i) {
		const unsigned k = this->adjacencylist[i].size();
		weighted_nodes.push_back(k);
		edges_present += k;
	}
	assert(edges_absent % 2 == 0);
	edges_present /= 2;
	edges_absent = size * (size - 1) / 2 - edges_present;
}

absolutetime_t dynamic_erdos_reyni::next(rng_t& engine) {
	if (!std::isnan(next_time))
		return next_time;
	
	/* Gillespie algorith: Draw time of next event */
	const double rate = edges_absent * alpha + edges_present * beta;
	next_time = current_time + std::exponential_distribution<>(rate)(engine);
	return next_time;
}

std::optional<network_event_t> dynamic_erdos_reyni::step(rng_t& engine, absolutetime_t max_time) {
	/* Determine time of next event if necessary, return if after max_time */
	if (std::isnan(next_time))
		next(engine);
	if (max_time < next_time)
		return std::nullopt;
	
	/* Get time of next event */
	const double event_time = next_time;
	next_time = NAN;
	
	/* Update current time */
	current_time = event_time;
	
	/* Draw type of event, edge appearing or disappearing */
	const double p = alpha * edges_absent / (alpha*edges_absent + beta*edges_present);
	const bool edge_appears = std::bernoulli_distribution(p)(engine);
	if (edge_appears) {
		/* We have to draw uniformly from the set of absent edges.
		 * We assume that most edges are absent, and thus just draw
		 * uniformly from all possible edges until we pick one that
		 * is actually absent.
		 * TODO: This is inefficient if the graph is almost complete.
		 */
		node_t src, dst;
		while (true) {
			std::uniform_int_distribution<node_t> uniform_node(0, this->adjacencylist.size()-1);
			/* Draw source node */
			src = uniform_node(engine);
			/* Draw destination node that isn't the same as the source node */
			do {
				dst = uniform_node(engine);
			} while (src == dst);
			/* Stop if the edge is indeed currently absent */
			const  std::vector<node_t>& al = this->adjacencylist.at(src);
			if (std::find(al.begin(), al.end(), dst) == al.end())
				break;
		};
		/* Add edge and return event */
		add_edge(src, dst);
		return network_event_t {
			.kind = network_event_kind::neighbour_added,
			.source_node = src, .target_node = dst, .time = event_time
		};
	} else {
		/* We have to draw uniformly from the set of present edges.
		 * We do so by drawing a random node, with nodes weighted by
		 * their degree, and then drawing a random neighbour of that node
		 */
		const node_t src = weighted_nodes(engine);
		const std::vector<node_t>& al = this->adjacencylist.at(src);
		std::uniform_int_distribution<node_t> uniform_neighbour(0, al.size()-1);
		const int neighbour_index = uniform_neighbour(engine);
		/* Remove edge */
		const node_t dst = al[neighbour_index];
		remove_edge(src, neighbour_index);
		return network_event_t {
			.kind = network_event_kind::neighbour_removed,
			.source_node = src, .target_node = dst, .time = event_time
		};
	}
}

void dynamic_erdos_reyni::add_edge(node_t node, node_t neighbour) {
	/* add forward edge */
	std::vector<node_t>& al_node = this->adjacencylist.at(node);
	al_node.push_back(neighbour);
	weighted_nodes[node] += 1;

	/* add reverse edge */
	std::vector<node_t>& al_neighbour = this->adjacencylist.at(neighbour);
	al_neighbour.push_back(node);
	weighted_nodes[neighbour] += 1;
	
	/* update counters */
	edges_absent -= 1;
	edges_present += 1;
}

void dynamic_erdos_reyni::remove_edge(node_t node, int neighbour_index) {
	using std::swap;
	
	/* Remove forward edge */
	assert(weighted_nodes[node] > 0);
	std::vector<node_t>& al_node = this->adjacencylist.at(node);
	assert(!al_node.empty());
	const node_t neighbour = al_node[neighbour_index];
	swap(al_node[neighbour_index], al_node.back());
	al_node.pop_back();
	weighted_nodes[node] -= 1;
	
	/* Remove reverse edge */
	std::vector<node_t>& al_neighbour = this->adjacencylist.at(neighbour);
	auto i = std::find(al_neighbour.begin(), al_neighbour.end(), node);
	swap(*i, al_neighbour.back());
	al_neighbour.pop_back();
	weighted_nodes[neighbour] -= 1;
	
	/* update counters */
	edges_absent += 1;
	edges_present -= 1;
}

