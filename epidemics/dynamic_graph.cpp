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

dynamic_empirical_network::dynamic_empirical_network(std::string path_to_file,double dt){

    std::ifstream file(path_to_file); // Open the CSV file
    
	if (!file.is_open()) {
        throw std::runtime_error("Error: Unable to open file: " + path_to_file);
    }
	std::set<int> set_active_edges;

	int nb_nodes = 0;

    std::string line;
    while (std::getline(file, line)) { // Read each line of the CSV file
        std::vector<int> row;
        std::istringstream iss(line);
		int src,dst;
		double next_time;
		if (iss >> src >> dst >> next_time){
			bool duplicate = false;


			for (auto it = edges.rbegin(); it != edges.rend(); ++it) {

				/* this edge already existed at some point*/
				if ((it->source_node == src) && (it->target_node ==dst)){
					
					/* we need to update and extend the edge's life duration*/
					if (it -> time + dt >= next_time) {
						it -> time = next_time + dt;
						assert(it-> kind == network_event_kind::neighbour_removed);
						duplicate=true;
					}
					/*break from the vector loop*/
					break;
				}
			}
			if (!duplicate){
				edges.push_back( network_event_t{
					.kind = network_event_kind::neighbour_added,
					.source_node = src, .target_node = dst, .time = next_time
				});

				edges.push_back( network_event_t{
					.kind = network_event_kind::neighbour_removed,
					.source_node = src, .target_node = dst, .time = next_time + dt
				});
			}
			nb_nodes = std::max({nb_nodes, src, dst});
			
		}
	}



	adjacencylist.resize(nb_nodes+1);
	// Sort the vector in increasing order of 'time'
    std::sort(edges.begin(), edges.end(), [](const network_event_t& a, const network_event_t& b) {
        return a.time < b.time;
    });

	for (node_t i = 0; i < (int) edges.size(); i++)
	{
		const network_event_t edge = edges[i];
		if (edge.time > 0)
			break;
		adjacencylist[edge.source_node].push_back(edge.target_node);
		// adjacencylist[edge.target_node].push_back(edge.source_node);
	}
	

	max_index = (int) edges.size();

}

node_t dynamic_empirical_network::nodes() {
    return (node_t)adjacencylist.size();
}    

node_t dynamic_empirical_network::neighbour(node_t node, int neighbour_index) {
    const auto& n = adjacencylist.at(node);
    if ((neighbour_index < 0) || (n.size() <= (unsigned int)neighbour_index))
        return -1;
    return n[neighbour_index];
}

index_t dynamic_empirical_network::outdegree(node_t node) {
    return (index_t) adjacencylist.at(node).size();
}


absolutetime_t dynamic_empirical_network::next(rng_t& engine) {
	// if (!std::isnan(next_time))
	// 	return next_time;
	
	/* No unreported reverse-edge event should exit */
	// assert(!reverse_edge_event);

	if (time_index >= (int) edges.size()){
		return INFINITY;
	}
	next_time = edges[time_index].time;
	return next_time;
}

int dynamic_empirical_network::present_edges(double t){

	
	int nb_edges = 0;
	for (auto edge : edges) {
		// if (edge.time < t)
		// 	continue;
		
		if (edge.time > t)
			break;

		if (edge.kind==network_event_kind::neighbour_added){
			nb_edges++;
		} else if (edge.kind==network_event_kind::neighbour_removed){
			nb_edges--;
		}

	
	}
	return nb_edges;
}



double dynamic_empirical_network::average_degree(double t){

	std::unordered_map<int, int> degree;
	int nb_edges = 0;


	for (auto edge : edges) {
		// if (edge.time < t)
		// 	continue;
		
		if (edge.time > t)
			break;

		const int src = edge.source_node;
		const int dst = edge.target_node;

		if (edge.kind==network_event_kind::neighbour_added){
			nb_edges++;
			degree[src]++;
			degree[dst]++;
		} else if (edge.kind==network_event_kind::neighbour_removed){
			degree[src]--;
			degree[dst]--;
			nb_edges--;
		}
	}

	int active_nodes = 0;
	// Iterate over the map and count nodes having at least degree 1
    for (const auto& pair : degree) {
        if (pair.second >= 1) {
            active_nodes++;
        }
    }
	return (double) 2*nb_edges/active_nodes;
}

std::optional<network_event_t> dynamic_empirical_network::step(rng_t& engine, absolutetime_t max_time) {
	
	/* Determine time of next event if necessary, return if after max_time */
	if (std::isnan(next_time))
		next(engine);
	if (max_time < next_time)
		return std::nullopt;

	if (time_index == max_index)
		return std::nullopt;

	// /* Report last reported event again but for the reverse edge */
	// if (reverse_edge_event) {
	// 	assert(current_time == next_time);
	// 	assert(reverse_edge_event->time == next_time);
	// 	const network_event_t r = *reverse_edge_event;
	// 	reverse_edge_event = std::nullopt;
	// 	next_time = NAN;
	// 	return r;
	// }


	/* Update current time
	 * We don't reset next_time here because we'll only report the forward event below.
	 * The reverse event will be reported upon the next invocation of step(), and which
	 * time next_time will be reset to NAN. Only then will future calls to next() draw
	 * a new random next_time
	 */
	current_time = next_time;

	network_event_t next_event = edges[time_index];
	time_index++;

	if (next_event.kind==network_event_kind::neighbour_removed){

		auto it = std::find(adjacencylist[next_event.source_node].begin(), adjacencylist[next_event.source_node].end(), next_event.target_node);
    	if (it != adjacencylist[next_event.source_node].end()) adjacencylist[next_event.source_node].erase(it);

		// it = std::find(adjacencylist[next_event.target_node].begin(), adjacencylist[next_event.target_node].end(), next_event.source_node);
    	// if (it != adjacencylist[next_event.target_node].end()) adjacencylist[next_event.target_node].erase(it);

		
	} else if (next_event.kind==network_event_kind::neighbour_added){

		auto it = std::find(adjacencylist[next_event.source_node].begin(), adjacencylist[next_event.source_node].end(), next_event.target_node);
    	if (it == adjacencylist[next_event.source_node].end())
			adjacencylist[next_event.source_node].push_back(next_event.target_node);
		// adjacencylist[next_event.target_node].push_back(next_event.source_node);
	}

	// reverse_edge_event = network_event_t {
	// 	.kind = next_event.kind,
	// 	.source_node = next_event.target_node, .target_node =  next_event.source_node, .time = next_time
	// };
	return next_event;


}


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

