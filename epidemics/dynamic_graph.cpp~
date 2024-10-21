//
//  dynamic_graph.cpp
//  epidemics
//
//  Created by Florian G. Pflug on 14.02.23.
//

#include "dynamic_graph.h"

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- DYNAMIC_NETWORK ------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

void dynamic_network::notify_epidemic_event(event_t ev, rng_t& engine) {
	/* Do nothing by default */
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- MUTABLE_GRAPH --------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

void graph_mutable::resize(node_t nodes)
{
	adjacencylist.resize((std::size_t)nodes);
}

bool graph_mutable::has_edge(node_t src, node_t dst)
{
	auto& al = adjacencylist.at(src);
	return (al.find(dst) != al.end());
}

bool graph_mutable::add_edge(node_t src, node_t dst)
{
	return adjacencylist.at(src).insert(dst).second;
}

bool graph_mutable::remove_edge(node_t src, node_t dst) {
	return (adjacencylist.at(src).erase(dst) > 0);
}

node_t graph_mutable::nodes() {
	return adjacencylist.size();
}

node_t graph_mutable::neighbour(node_t node, int neighbour_index) {
	const auto& al = adjacencylist.at(node);
	if ((neighbour_index < 0) || ((std::size_t)neighbour_index >= al.size()))
		return -1;
	return al[neighbour_index];
}

index_t graph_mutable::outdegree(node_t node) {
	return (index_t) adjacencylist.at(node).size();
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- DYNAMIC NETWORK: EMPIRICAL -==----------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

dynamic_empirical_network::dynamic_empirical_network(std::string path_to_file, dynamic_empirical_network::edge_duration_kind contact_type, interval_t dt)
{
	node_t max_node = 0;

	/* read whitespace-separated file */
	std::ifstream file(path_to_file);
	if (!file.is_open())
		throw std::runtime_error("unable to open file: " + path_to_file);
	std::string line;
	while (std::getline(file, line)) {
		/* read line of the form: src <space> dst <space> time */
		std::istringstream iss(line);
		int src,dst;
		double time;
		if (!(iss >> src >> dst >> time) || (src < 0) || (dst < 0))
			throw std::runtime_error("unable to parse line: " + line);
		max_node = std::max(max_node, std::max(src, dst));

		/* queue event */
		switch (contact_type) {
			case finite_duration: {
				network_event_t ev1 {
					.kind = network_event_kind::neighbour_added,
					.source_node = src,
					.target_node = dst,
					.time = time - dt/2
				};
				event_queue.push_back(ev1);

				network_event_t ev2 {
					.kind = network_event_kind::neighbour_removed,
					.source_node = src,
					.target_node = dst,
					.time = time + dt/2
				};
				event_queue.push_back(ev2);
				break;
			}
			case infitesimal_duration: {
				network_event_t ev {
					.kind = network_event_kind::instantenous_contact,
					.source_node = src,
					.target_node = dst,
					.time = time,
					.infitesimal_duration = dt
				};
				event_queue.push_back(ev);
				break;
			}
		}
	}

	/* make sure events are sorted */
	auto cmp = [](const network_event_t& a, const network_event_t& b) {
		return a.time < b.time;
	};
	if (!std::is_sorted(event_queue.begin(), event_queue.end(), cmp))
		std::sort(event_queue.begin(), event_queue.end(), cmp);

	/* Allocate adjacency list */
	resize(max_node + 1);
}

absolutetime_t dynamic_empirical_network::next(rng_t& engine)
{
	/* find next actual event, skipping over events that do nothing */
	for(network_event_t& event: event_queue) {
		switch (event.kind) {
			case network_event_kind::neighbour_added:
				if (!has_edge(event.source_node, event.target_node))
					return event.time;
				break;
			case network_event_kind::neighbour_removed:
				if (has_edge(event.source_node, event.target_node))
					return event.time;
				break;
			default:
				return event.time;
		}
	}
	return INFINITY;
}

std::optional<network_event_t> dynamic_empirical_network::step(rng_t& engine, absolutetime_t max_time)
{
	/* loop until we find an event, necessary because add/remove may be skipped if edge already exists */
	while (!event_queue.empty()) {
		/* process event if not later than max_time */
		network_event_t event = event_queue.front();
		if (!std::isnan(max_time) && (event.time > max_time))
			break;
		event_queue.pop_front();

		/* update graph
		 * NOTE: The event-skipping logic below must be kept in sync with next()
		 */
		switch (event.kind) {
			case network_event_kind::neighbour_added:
				/* add edge, report event unless edge already existed */
				if (add_edge(event.source_node, event.target_node))
					return event;
				break;
			case network_event_kind::neighbour_removed:
				/* remove edge, report event if edge didn't exist */
				if (remove_edge(event.source_node, event.target_node))
					return event;
				break;
			default:
				/* never skipped */
				return event;
		}
	}
	return std::nullopt;
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- DYNAMIC NETWORK: ERDÖS REYNI -----------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

dynamic_erdos_reyni::dynamic_erdos_reyni(int size, double avg_degree, double timescale, rng_t& engine)
	:erdos_reyni(size, avg_degree, engine)
	,edge_probability(avg_degree / (size - 1))
	,alpha(edge_probability / timescale)
	,beta((1.0 - edge_probability) / timescale)
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
	
	/* No unreported reverse-edge event should exit */
	assert(!reverse_edge_event);

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
	
	/* Report last reported event again but for the reverse edge */
	if (reverse_edge_event) {
		assert(current_time == next_time);
		assert(reverse_edge_event->time == next_time);
		const network_event_t r = *reverse_edge_event;
		reverse_edge_event = std::nullopt;
		next_time = NAN;
		return r;
	}

	/* Update current time
	 * We don't reset next_time here because we'll only report the forward event below.
	 * The reverse event will be reported upon the next invocation of step(), and which
	 * time next_time will be reset to NAN. Only then will future calls to next() draw
	 * a new random next_time
	 */
	current_time = next_time;
	
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
		reverse_edge_event = network_event_t {
			.kind = network_event_kind::neighbour_added,
			.source_node = dst, .target_node = src, .time = next_time
		};
		return network_event_t {
			.kind = network_event_kind::neighbour_added,
			.source_node = src, .target_node = dst, .time = next_time
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
		reverse_edge_event = network_event_t {
			.kind = network_event_kind::neighbour_removed,
			.source_node = dst, .target_node = src, .time = next_time
		};
		return network_event_t {
			.kind = network_event_kind::neighbour_removed,
			.source_node = src, .target_node = dst, .time = next_time
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

