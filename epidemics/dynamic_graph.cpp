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

std::vector<std::vector<double>> dynamic_empirical_network::compute_number_of_edges(rng_t& engine){

	std::vector<std::vector<double>> average_degree;
	int number_of_edges=0;

	while (auto ev = step(engine)) {
		const network_event_t event = *ev;

		switch (event.kind){
			case network_event_kind::neighbour_added: 
				number_of_edges++;
				break;
			case network_event_kind::neighbour_removed:
				number_of_edges--;
				break;
			case network_event_kind::instantenous_contact:
				number_of_edges++;
				break;
			default: throw std::logic_error("invalid event kind");
		}
		const double time = event.time;
		average_degree.push_back({time, (double) number_of_edges});
		     
	}


	return average_degree;
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- DYNAMIC NETWORK: SIRX NETWORK ----------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

dynamic_sirx_network::dynamic_sirx_network(graph& network_, double kappa0_, double kappa_)
	:network(network_), network_is_undirected(network.is_undirected()), network_size(network.nodes())
	,kappa0(kappa0_), kappa(kappa_), queue_next_flipped(is_undirected())
{
	for(node_t n = 0; n < network_size; ++n)
		nonremoved.insert(n);
}

dynamic_sirx_network::node_state_t dynamic_sirx_network::state(node_t node)
{
	if (is_infected(node)) {
		if (is_removed(node))
			return X;
		else
			return I;
	} else {
		if (is_removed(node))
			return R;
		else
			return S;
	}
}


bool dynamic_sirx_network::is_undirected()
{
	return network_is_undirected;
}

node_t dynamic_sirx_network::nodes()
{
	return network_size;
}

node_t dynamic_sirx_network::neighbour(node_t node, int neighbour_index)
{
	if (is_removed(node))
		return -1;
	return network.neighbour(node, neighbour_index);
}

index_t dynamic_sirx_network::outdegree(node_t node)
{
	if (is_removed(node))
		return 0;
	return network.outdegree(node);
}

void dynamic_sirx_network::notify_epidemic_event(event_t ev, rng_t& engine)
{
	assert(current_time <= ev.time);
	assert(std::isnan(next_time) || ev.time <= next_time);
	switch(ev.kind) {
		case event_kind::infection:
		case event_kind::outside_infection:
			/* mark node as infected */
			infected.insert(ev.node);
			/* keep track of infected, non-removed nodes. if the nw is undirected,
			 * removed nodes cannot be infected, so the node must be non-removed
			 */
			assert(!is_undirected() || !is_removed(ev.node));
			if (is_undirected() || !is_removed(ev.node))
				infected_nonremoved.insert(ev.node);
			/* clear previously generated next_time since the rates have changed */
			if (next_time > ev.time) {
				next_time = NAN;
				queue.clear();
			}
			break;
		case event_kind::reset:
			/* mark node as no longer infected */
			infected.erase(ev.node);
			infected_nonremoved.erase(ev.node);
			/* clear previously generated next_time since the rates have changed */
			if (next_time > ev.time) {
				next_time = NAN;
				queue.clear();
			}
			break;
		default:
			break;
	}
}

absolutetime_t dynamic_sirx_network::next(rng_t& engine)
{
	/* generate next time unless already done previously */
	double base_time = current_time;
	while (std::isnan(next_time)) {
		/* should not have queued events */
		assert(queue.empty());

		/* I. find the time of the next node removal */
		/* removal rate is kappa0 for all non-removed nodes,
		 * and additionally kappa for all non-removed infected nodes */
		const double r0 = nonremoved.size() * kappa0;
		const double r = infected_nonremoved.size() * kappa;
		if (r + r0 == 0.0) {
			/* total rates are both zero, next event at infinity */
			next_time = INFINITY;
			return next_time;
		}
		std::exponential_distribution dist(r0 + r);
		next_time = base_time + dist(engine);

		/* II. find the removed node
		 * Determine whether to remove non-specifically or specifically an infected node.
		 * Total rate of non-specific removal is r0 = size * kappa0,
		 * total rate of specific removal of infected nodes is r = infected.size * kappa,
		 * so the probability of a non-specific removal is p0 = r0 / (r + r0).
		 */
		node_t n = -1;
		const double p0 = r0 / (r0 + r);
		if (std::bernoulli_distribution(p0)(engine)) {
			/* non-specific removal */
			n = *nonremoved(engine);
		} else {
			/* specific removal of infected node */
			n = *infected_nonremoved(engine);
		}

		/* III. queue events for the removal of all (outgoing) edges of the node
		 * For undirected networks, we remove incoming and outgoing nodes, for directed
		 * networks only all outgoing nodes (so the node can still be infected!)
		 * NOTE: for directed networks, it would be great to remove also incoming edges,
		 * but there's currently no way to find them all. So instead we keep incoming edges,
		 * meaning that removed nodes can still be infected for directed networks,
		 * but can't infect others.
		 */
		for(int i=0, nn=0; (nn = network.neighbour(n, i)) >= 0; ++i) {
			queue.push_back(network_event_t {
				.kind = network_event_kind::neighbour_removed,
				.source_node = n,
				.target_node = nn,
				.time = next_time
			});
		}

		/* if the node has no neighbours, mark as removed and move to next event */
		if (queue.empty()) {
			nonremoved.erase(n);
			infected_nonremoved.erase(n);
			base_time = next_time;
			next_time = NAN;
		}
	}

	return next_time;
}

std::optional<network_event_t> dynamic_sirx_network::step(rng_t& engine, absolutetime_t max_time)
{
	/* make sure the next event was generated, do nothing if past max_time */
	next(engine);
	if (next_time > max_time)
		return std::nullopt;

	/* get next event off queue, it's time must be next_time */
	assert(!queue.empty());
	network_event_t ev = queue.front();
	const node_t node = ev.source_node;
	assert(ev.time == next_time);
	assert(current_time <= next_time);
	current_time = next_time;
	/* for undirected networks, report each edge twice, flipped and unflipped */
	if (queue_next_flipped) {
		/* flip edge this time, next time report unflipped edge */
		assert(is_undirected());
		std::swap(ev.source_node, ev.target_node);
		queue_next_flipped = false;
	} else {
		/* report unflipped edge and dequeue */
		queue.pop_front();
		queue_next_flipped = is_undirected();
		/* and clear next_time only once we've emptied the queue */
		if (queue.empty())
			next_time = NAN;
	}

	/* when first reporting an edge, update network */
	if (!queue_next_flipped) {
		nonremoved.erase(node);
		infected_nonremoved.erase(node);
	}

	/* return event */
	return ev;
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


/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- DYNAMIC NETWORK: ACTIVITY DRIVEN -------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

activity_driven_network::activity_driven_network(std::vector<double> activity_rates, double eta, double m, double recovery_rate, rng_t& engine)
    : activity_rates(activity_rates), eta(eta), m(m), recovery_rate(recovery_rate), engine(engine)
{
    const int N = (int) activity_rates.size(); 
	resize(N);

	active_nodes.resize(N, false);

    std::uniform_int_distribution<> distr(0, N - 1);

    for (node_t node = 0; node < N; node++)
    {
	    const double ai = activity_rates[node];

        const double activation_time = std::exponential_distribution<>(eta * ai)(engine);

		active_edges_entry e;
		e.kind = network_event_kind::neighbour_added;
		e.source_node = node;
		e.target_node = -1;
		e.time = activation_time;
		e.activity_event = activity_event_kind::activate;
		push_edge(e);
        
    }
}

absolutetime_t activity_driven_network::next(rng_t& engine){

	if (active_edges.empty())
		return INFINITY;
	
	auto next = top_edge();
	
	while( next.activity_event != activity_event_kind::none){

		switch (next.activity_event)
		{
		case activity_event_kind::activate:
			/* Dequeue event */
			pop_edge();
			activate_node(next.source_node,next.time);
			next = top_edge();
			break;
		case activity_event_kind::deactivate:
			/* Dequeue event */
			pop_edge();
			deactivate_node(next.source_node,next.time);
			next = top_edge();
			break;
		default:
			throw std::runtime_error("invalid event kind");
			break;
		}

	}

	return next.time;


}

void activity_driven_network::activate_node(node_t node,double time){

	std::uniform_int_distribution<> distr(0, (int) activity_rates.size() - 1);

	if (active_nodes[node])
		throw std::logic_error("node shouldnt be already active");

	active_nodes[node] = true;

	for (int i = 0; i < m; i++)
	{

		node_t neigh = node;
		while(neigh==node)
			neigh = distr(engine);

		if (has_edge(node, neigh) && has_edge(neigh,node))
			continue;
		if (node==neigh)
			continue;
		 
		active_edges_entry e;
		e.kind = network_event_kind::neighbour_added;
		e.source_node = node;
		e.target_node = neigh;
		e.time = time;
		e.activity_event = activity_event_kind::none;
		push_edge(e);

		active_edges_entry e2;
		e2.kind = network_event_kind::neighbour_added;
		e2.source_node = neigh;
		e2.target_node = node;
		e2.time = time;
		e2.activity_event = activity_event_kind::none;
		push_edge(e2);
	}

	active_edges_entry e;
	e.kind = network_event_kind::neighbour_removed;
	e.source_node = node;
	e.target_node = -1;
	e.time =time + std::exponential_distribution<>(recovery_rate)(engine);
	e.activity_event = activity_event_kind::deactivate;
	push_edge(e);

}

void activity_driven_network::deactivate_node(node_t node,double time){

	std::uniform_int_distribution<> distr(0, (int) activity_rates.size() - 1);

	if (!active_nodes[node])
		throw std::logic_error("node shouldnt be inactive");
	active_nodes[node]=false;


	const int k = outdegree(node);
	for (int ki = 0; ki < k ; ki++){
		const node_t nei = neighbour(node,ki);

		if ((!has_edge(node, nei)) && (!has_edge(nei, node)))
			throw std::logic_error("there should be a link between nei and node");

		if (node ==nei)
			throw std::logic_error("no self links were allowed in the first place");

		active_edges_entry e;
		e.kind = network_event_kind::neighbour_removed;
		e.source_node = node;
		e.target_node = nei;
		e.time = time;
		e.activity_event = activity_event_kind::none;
		push_edge(e);

		active_edges_entry e2;
		e2.kind = network_event_kind::neighbour_removed;
		e2.source_node = nei;
		e2.target_node = node;
		e2.time = time;
		e2.activity_event = activity_event_kind::none;
		push_edge(e2);		
	}

	active_edges_entry e;
	e.kind = network_event_kind::neighbour_added;
	e.source_node = node;
	e.target_node = -1;
	e.time =time + std::exponential_distribution<>(eta * activity_rates[node])(engine);
	e.activity_event = activity_event_kind::activate;
	push_edge(e);

}



std::optional<network_event_t> activity_driven_network::step(rng_t& engine, absolutetime_t max_time){

	absolutetime_t next_time = next(engine);

	if (next_time > max_time)
		return std::nullopt;

	/* If there are no more infection times, stop */
	if (active_edges.empty())
		return std::nullopt;


	auto next = top_edge();
	if (next.time > max_time)
		return std::nullopt;
	
	pop_edge();

	std::uniform_int_distribution<> distr(0, (int) activity_rates.size() - 1);

	const node_t src = next.source_node;
	const node_t dst = next.target_node;
	// Handle event before returning the edge
	switch (next.kind) {
		case network_event_kind::neighbour_added:{
				
			add_edge(src,dst);

			break;
		}
		case network_event_kind::neighbour_removed:{

			const bool r1 =remove_edge(src,dst);
			// const bool r2 = remove_edge(dst,src);
			if (!r1)
				throw std::logic_error("removing edeges fialed");
			break;	
		}
		default: 
			throw std::runtime_error("invalid event kind");
	}


	return network_event_t{
					.kind = next.kind,
					.source_node = next.source_node, .target_node = next.target_node, .time = next.time
				};
}
