//
//  dynamic_graph.cpp
//  epidemics
//
//  Created by Florian G. Pflug on 14.02.23.
//

#include "nextnet/temporal_network.h"

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- TEMPORAL_NETWORK -----------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

void temporal_network::notify_epidemic_event(epidemic_event_t ev, rng_t &engine)
{
    /* Do nothing by default */
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- MUTABLE_NETWORK ------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

void mutable_network::resize(node_t nodes)
{
    adjacencylist.resize((std::size_t)nodes);
}

bool mutable_network::has_edge(node_t src, node_t dst)
{
    auto &al = adjacencylist.at(src);
    return (al.find(dst) != al.end());
}

bool mutable_network::add_edge(node_t src, node_t dst)
{
    return adjacencylist.at(src).insert(dst).second;
}

bool mutable_network::remove_edge(node_t src, node_t dst)
{
    return (adjacencylist.at(src).erase(dst) > 0);
}

node_t mutable_network::nodes()
{
    return (node_t)adjacencylist.size();
}

node_t mutable_network::neighbour(node_t node, int neighbour_index)
{
    const auto &al = adjacencylist.at(node);
    if ((neighbour_index < 0) || ((std::size_t)neighbour_index >= al.size()))
        return -1;
    return al[neighbour_index];
}

index_t mutable_network::outdegree(node_t node)
{
    return (index_t)adjacencylist.at(node).size();
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- MUTABLE_WEIGHTED_NETWORK ---------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

void mutable_weighted_network::resize(node_t nodes)
{
    adjacencylist.resize((std::size_t)nodes);
}

bool mutable_weighted_network::has_edge(node_t src, node_t dst)
{
    auto &al = adjacencylist.at(src);
    return (al.find(dst) != al.end());
}

bool mutable_weighted_network::has_edge(node_t src, node_t dst, double *weight)
{
    auto &al     = adjacencylist.at(src);
    const auto i = al.find(dst);
    if (i == al.end()) {
        *weight = NAN;
        return false;
    }
    *weight = i->second;
    return true;
}

double mutable_weighted_network::edge_weight(node_t src, node_t dst)
{
    auto &al = adjacencylist.at(src);
    auto i   = al.find(dst);
    if (i == al.end())
        return NAN;
    return i->second;
}

void mutable_weighted_network::add_edge(node_t src, node_t dst, double weight)
{
    if ((weight < 0) || !std::isfinite(weight))
        throw std::range_error("weights must be positive and finite");
    if (weight == 0.0)
        remove_edge(src, dst);
    else
        adjacencylist.at(src)[dst] = weight;
}

bool mutable_weighted_network::remove_edge(node_t src, node_t dst)
{
    return (adjacencylist.at(src).erase(dst) > 0);
}

node_t mutable_weighted_network::nodes()
{
    return (node_t)adjacencylist.size();
}

node_t mutable_weighted_network::neighbour(node_t node, int neighbour_index)
{
    const auto &al = adjacencylist.at(node);
    if ((neighbour_index < 0) || ((std::size_t)neighbour_index >= al.size()))
        return -1;
    return al[neighbour_index].first;
}

node_t mutable_weighted_network::neighbour(node_t node, int neighbour_index, double *weight)
{
    const auto &al = adjacencylist.at(node);
    if ((neighbour_index < 0) || ((std::size_t)neighbour_index >= al.size()))
        return -1;
    const auto &n = al[neighbour_index];
    if (weight != nullptr)
        *weight = n.second;
    return n.first;
}

index_t mutable_weighted_network::outdegree(node_t node)
{
    return (index_t)adjacencylist.at(node).size();
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- NEXT_REACTION_NETWORK ------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

next_reaction_network::next_reaction_network(network_kind kind_)
    : kind(kind_)
{
}

bool next_reaction_network::is_undirected()
{
    return !(kind & directed_kind);
}

bool next_reaction_network::is_unweighted()
{
    return !(kind & weighted_kind);
}

void next_reaction_network::queue_add_edge(node_t src, node_t dst, double w, absolutetime_t time, tag_t tag)
{
    queue_network_event(network_event_t{
                            .kind        = network_event_kind::neighbour_added,
                            .source_node = src,
                            .target_node = dst,
                            .weight      = w,
                            .time        = time },
                        tag);
}

void next_reaction_network::queue_remove_edge(node_t src, node_t dst, double w, absolutetime_t time, tag_t tag)
{
    queue_network_event(network_event_t{
                            .kind        = network_event_kind::neighbour_removed,
                            .source_node = src,
                            .target_node = dst,
                            .weight      = w,
                            .time        = time },
                        tag);
}

void next_reaction_network::queue_instantenous_contact(node_t src, node_t dst, double w, absolutetime_t time, tag_t tag)
{
    queue_network_event(network_event_t{
                            .kind        = network_event_kind::instantenous_contact,
                            .source_node = src,
                            .target_node = dst,
                            .weight      = w,
                            .time        = time },
                        tag);
}

void next_reaction_network::queue_network_event(network_event_t ev, tag_t tag)
{
    if (ev.time < current_time)
        throw std::range_error("cannot queue past events");

    push_event(queue_entry{ .time = ev.time, .tag = tag, .event = ev });

    /* Forget cached next_time */
    if (ev.time < next_time)
        next_time = NAN;
}

void next_reaction_network::queue_callback(absolutetime_t time, event_callback_t cb, tag_t tag)
{
    if (time < current_time)
        throw std::range_error("cannot queue past events");

    push_event(queue_entry{ .time = time, .tag = tag, .event = cb });

    /* Forget cached next_time */
    if (time < next_time)
        next_time = NAN;
}

void next_reaction_network::clear_tag(tag_t tag)
{
    auto &queue_by_tag = event_queue.get<by_tag>();
    auto r             = queue_by_tag.equal_range(tag);
    queue_by_tag.erase(r.first, r.second);
	
	/* Forget cached next_time */
	next_time = NAN;
}

absolutetime_t next_reaction_network::next(rng_t &engine, absolutetime_t maxtime)
{
    /* If there's a cached next_time, return it */
    if (!std::isnan(next_time))
        return next_time;

    /* Find time of next event that step() will report */
    while (true) {
		if (event_queue.empty()) {
			next_time = INFINITY;
			break;
		}

        /* Fetch next event, return nothing if past max_time */
        const queue_entry &top = top_event();
		if (top.time > maxtime) {
			next_time = INFINITY;
			break;
		}

        /* Check whether step() would skip the event or not */
        if (const network_event_t *nwev = std::get_if<network_event_t>(&top.event)) {
            assert(nwev->time == top.time);
            network_event_t ev = *nwev;
            switch (nwev->kind) {
                case network_event_kind::neighbour_added:
                    double w;
                    if (has_edge(ev.source_node, ev.target_node, &w) &&
                        (!(kind & weighted_kind) || (w == ev.weight))) {
                        pop_event();
                        continue;
                    }
                    break;
                case network_event_kind::neighbour_removed:
                    if (!has_edge(ev.source_node, ev.target_node)) {
                        pop_event();
                        continue;
                    }
                    break;
                default:
                    break;
            };
        } else if (const event_callback_t *cbev = std::get_if<event_callback_t>(&top.event)) {
            /* Execute callback. Since this is not a reported event, we continue */
            cbev->operator()();
            pop_event();
            continue;
        }

        next_time = top.time;
		break;
    }

	return next_time;
}

std::optional<network_event_t> next_reaction_network::step(rng_t &engine, absolutetime_t max_time)
{
    while (true) {
        if (event_queue.empty())
            return std::nullopt;

        /* Fetch next event, return nothing if past max_time */
        const queue_entry &top = top_event();
        if (top.time > max_time)
            return std::nullopt;

        assert(current_time <= top.time);
        current_time = top.time;

        /* Decode and handle event */
        assert(std::isnan(next_time) || (next_time >= top.time));
        if (const network_event_t *nwev = std::get_if<network_event_t>(&top.event)) {
            /* Network event */
            assert(nwev->time == top.time);

            /* Get event and flip edge for second pass on undirected networks */
            network_event_t ev = *nwev;
            if (flipped_edge_next)
                std::swap(ev.source_node, ev.target_node);

            /* Time of the event should agree with cached next_time */
            assert(std::isnan(next_time) || (next_time == ev.time));

            /* Execute event
             * NOTE: The logic of which events we skip below must match next(),
             * otherwise we break the rule that if next() reports a certain
             * time t for the next event, step(t) must not return nullopt.
             */
            const double weight = (kind & weighted_kind) ? ev.weight : 1.0;
            switch (ev.kind) {
                case network_event_kind::neighbour_added:
                    double w;
                    if (has_edge(ev.source_node, ev.target_node, &w)) {
                        /* Skip if edge exists with correct weight */
                        if (!(kind & weighted_kind) || (w == weight)) {
                            next_time = NAN;
                            pop_event();
                            continue;
                        }
                        /* Remove existing edge with incorrect weight */
                        if ((kind & weighted_kind) && (w != weight)) {
                            remove_edge(ev.source_node, ev.target_node);
                            return network_event_t{
                                .kind        = network_event_kind::neighbour_removed,
                                .source_node = ev.source_node,
                                .target_node = ev.target_node,
                                .weight      = w,
                                .time        = ev.time
                            };
                        }
                    }
                    add_edge(ev.source_node, ev.target_node, weight);
                    break;
                case network_event_kind::neighbour_removed:
                    /* Skip event if the edge does not exist */
                    if (!has_edge(ev.source_node, ev.target_node)) {
                        next_time = NAN;
                        pop_event();
                        continue;
                    }
                    remove_edge(ev.source_node, ev.target_node);
                    break;
                default:
                    /* do nothing */
                    break;
            };

            /* On undirected networks, if this was the first pass for an edge event,
             * do a second pass for the flipped version of the edge
             */
            flipped_edge_next = (!flipped_edge_next &&
                                 !(kind & directed_kind) &&
                                 ((ev.kind == network_event_kind::neighbour_added) ||
                                  (ev.kind == network_event_kind::neighbour_removed)));

            /* Dequeue event (unless we're doing the flipped version next) and return it */
            if (!flipped_edge_next) {
                next_time = NAN;
                pop_event();
            }

            return ev;
        } else if (const event_callback_t *cbev = std::get_if<event_callback_t>(&top.event)) {
            /* Execute callback. Since this is not a reported event, we continue */
            cbev->operator()();
            pop_event();
        } else {
            throw std::logic_error("unknown event in queue");
        }
    }
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- EMPIRICAL_CONTACT_NETWORK --------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

empirical_contact_network::empirical_contact_network(
    std::string path_to_file, network_kind kind,
    edge_duration_kind contact_type, interval_t dt)
    : next_reaction_network(kind)
{
    current_time    = INFINITY;
    node_t max_node = 0;

    /* read whitespace-separated file */
    std::ifstream file(path_to_file);
    if (!file.is_open())
        throw std::runtime_error("unable to open file: " + path_to_file);
    std::string line;
    while (std::getline(file, line)) {
        /* read line of the form: src <space> dst <space> time */
        std::istringstream iss(line);
        int src, dst;
        double time;
        if (!(iss >> src >> dst >> time) || (src < 0) || (dst < 0))
            throw std::runtime_error("unable to parse line: " + line);
        max_node = std::max(max_node, std::max(src, dst));

        /* queue event */
        switch (contact_type) {
            case finite_duration:
                current_time = std::min(current_time, time - dt / 2.0);
                queue_add_edge(src, dst, 1.0, time - dt / 2.0);
                queue_remove_edge(src, dst, 1.0, time + dt / 2.0);
                break;

            case infitesimal_duration:
                current_time = std::min(current_time, time);
                queue_instantenous_contact(src, dst, dt, time);
                break;
        }
    }

    /* Allocate adjacency list */
    resize(max_node + 1);
}

std::vector<std::vector<double>> empirical_contact_network::compute_number_of_edges(rng_t &engine)
{
    std::vector<std::vector<double>> average_degree;
    int number_of_edges = 0;

    while (auto ev = step(engine)) {
        const network_event_t event = *ev;

        switch (event.kind) {
            case network_event_kind::neighbour_added:
                number_of_edges++;
                break;
            case network_event_kind::neighbour_removed:
                number_of_edges--;
                break;
            case network_event_kind::instantenous_contact:
                number_of_edges++;
                break;
            default:
                throw std::logic_error("invalid event kind");
        }
        const double time = event.time;
        average_degree.push_back({ time, (double)number_of_edges });
    }

    return average_degree;
}


/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- DYNAMIC NETWORK: ERDÖS RENYI -----------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

temporal_erdos_renyi::temporal_erdos_renyi(int size, double avg_degree, double timescale, rng_t &engine)
    : erdos_renyi(size, avg_degree, engine)
    , edge_probability(avg_degree / (size - 1))
    , alpha(edge_probability / timescale)
    , beta((1.0 - edge_probability) / timescale)
    , edges_present(0)
{
    /* Initial degree-weights node distribution and present/absent edge counters */
    for (node_t i = 0; (std::size_t)i < this->adjacencylist.size(); ++i) {
        const unsigned k = (unsigned)this->adjacencylist[i].size();
        weighted_nodes.push_back(k);
        edges_present += k;
    }
    assert(edges_present % 2 == 0);
    edges_present /= 2;
    edges_absent = size * (size - 1) / 2 - edges_present;
}

absolutetime_t temporal_erdos_renyi::next(rng_t &engine, absolutetime_t)
{
    if (!std::isnan(next_time))
        return next_time;

    /* No unreported reverse-edge event should exit */
    assert(!reverse_edge_event);

    /* Gillespie algorith: Draw time of next event */
    const double rate = edges_absent * alpha + edges_present * beta;
    next_time         = current_time + std::exponential_distribution<>(rate)(engine);
    return next_time;
}

std::optional<network_event_t> temporal_erdos_renyi::step(rng_t &engine, absolutetime_t max_time)
{
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
        reverse_edge_event      = std::nullopt;
        next_time               = NAN;
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
    const double p          = alpha * edges_absent / (alpha * edges_absent + beta * edges_present);
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
            std::uniform_int_distribution<node_t> uniform_node(0, (node_t)this->adjacencylist.size() - 1);
            /* Draw source node */
            src = uniform_node(engine);
            /* Draw destination node that isn't the same as the source node */
            do {
                dst = uniform_node(engine);
            } while (src == dst);
            /* Stop if the edge is indeed currently absent */
            const std::vector<node_t> &al = this->adjacencylist.at(src);
            if (std::find(al.begin(), al.end(), dst) == al.end())
                break;
        };
        /* Add edge and return event */
        add_edge(src, dst);
        reverse_edge_event = network_event_t{
            .kind        = network_event_kind::neighbour_added,
            .source_node = dst,
            .target_node = src,
            .time        = next_time
        };
        return network_event_t{
            .kind        = network_event_kind::neighbour_added,
            .source_node = src,
            .target_node = dst,
            .time        = next_time
        };
    } else {
        /* We have to draw uniformly from the set of present edges.
         * We do so by drawing a random node, with nodes weighted by
         * their degree, and then drawing a random neighbour of that node
         */
        const node_t src              = (node_t)weighted_nodes(engine);
        const std::vector<node_t> &al = this->adjacencylist.at(src);
        std::uniform_int_distribution<node_t> uniform_neighbour(0, (node_t)al.size() - 1);
        const int neighbour_index = uniform_neighbour(engine);
        /* Remove edge */
        const node_t dst = al[neighbour_index];
        remove_edge(src, neighbour_index);
        reverse_edge_event = network_event_t{
            .kind        = network_event_kind::neighbour_removed,
            .source_node = dst,
            .target_node = src,
            .weight      = 1.0,
            .time        = next_time
        };
        return network_event_t{
            .kind        = network_event_kind::neighbour_removed,
            .source_node = src,
            .target_node = dst,
            .weight      = 1.0,
            .time        = next_time
        };
    }
}

void temporal_erdos_renyi::add_edge(node_t node, node_t neighbour)
{
    /* add forward edge */
    std::vector<node_t> &al_node = this->adjacencylist.at(node);
    al_node.push_back(neighbour);
    weighted_nodes[node] += 1;

    /* add reverse edge */
    std::vector<node_t> &al_neighbour = this->adjacencylist.at(neighbour);
    al_neighbour.push_back(node);
    weighted_nodes[neighbour] += 1;

    /* update counters */
    edges_absent -= 1;
    edges_present += 1;
}

void temporal_erdos_renyi::remove_edge(node_t node, int neighbour_index)
{
    using std::swap;

    /* Remove forward edge */
    assert(weighted_nodes[node] > 0);
    std::vector<node_t> &al_node = this->adjacencylist.at(node);
    assert(!al_node.empty());
    const node_t neighbour = al_node[neighbour_index];
    swap(al_node[neighbour_index], al_node.back());
    al_node.pop_back();
    weighted_nodes[node] -= 1;

    /* Remove reverse edge */
    std::vector<node_t> &al_neighbour = this->adjacencylist.at(neighbour);
    auto i                            = std::find(al_neighbour.begin(), al_neighbour.end(), node);
    swap(*i, al_neighbour.back());
    al_neighbour.pop_back();
    weighted_nodes[neighbour] -= 1;

    /* update counters */
    edges_absent += 1;
    edges_present -= 1;
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- TEMPORAL_SIRX_NETWORK ------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

temporal_sirx_network::network_kind temporal_sirx_network::kind(network* nw)
{
	network_kind r = (network_kind)0;
	r = (network_kind) (r | (nw->is_undirected() ? (network_kind)0 : directed_kind));
	r = (network_kind) (r | (dynamic_cast<weighted_network*>(nw) ? weighted_kind : (network_kind)0));
	return r;
}

temporal_sirx_network::temporal_sirx_network(network& nw, double kappa0_, double kappa_, rng_t& engine)
	: next_reaction_network(kind(&nw))
	, kappa0(kappa0_)
	, kappa(kappa_)
	, state(nw.nodes(), node_state{ .removed = false, .infected = false })
{
	weighted_network* wnw = dynamic_cast<weighted_network*>(&nw);
	
	/* Copy network structure */
	resize((node_t)nw.nodes());
	for(node_t n = 0, N = nw.nodes(); n < N; ++n) {
		node_t nn;
		double w = 1.0;
		for(int i = 0; (nn = (wnw ? wnw->neighbour(n, i, &w) : nw.neighbour(n, i))) >= 0; i++)
			add_edge(n, nn, w);
	}
	
	/* Queue removals */
	for(node_t n = 0, N = nw.nodes(); n < N; ++n)
		queue_removal(n, state[n], 0, engine);
}

void temporal_sirx_network::queue_removal(node_t node, node_state s, absolutetime_t time, rng_t& engine)
{
	/* Compute deactivation rate */
	const double k = s.infected ? (kappa0 + kappa) : kappa0;
	/* Deactivate immediately if rate is infinite, otherwise generate waiting time and queue */
	if (k == INFINITY)
		remove_node(node, time);
	else if (k > 0) {
		const double t = time + std::exponential_distribution<>(k)(engine);
		queue_callback(t, [this, node, t]() { this->remove_node(node, t); }, node);
	}
}

void temporal_sirx_network::remove_node(node_t node, absolutetime_t time)
{
	/* Fetch and update node state */
	node_state &s = state[node];
	if (s.removed)
		throw std::logic_error("node is already removed");
	s.removed = true;

	node_t nn;
	for (int i = 0; (nn = neighbour(node, i)) >= 0; ++i)
		queue_remove_edge(node, nn, 1.0, time);
}

void temporal_sirx_network::notify_epidemic_event(epidemic_event_t ev, rng_t &engine)
{
	/* Fetch node state */
	node_state &s = state[ev.node];

	/* Update infected state of node, re-genereate waiting times if necessary */
	switch (ev.kind) {
		case epidemic_event_kind::outside_infection:
		case epidemic_event_kind::infection:
			if (s.infected)
				throw std::logic_error("node is already infected");
			s.infected = true;

			/* if the remove rate of the node changed, re-generate waiting time */
			if (!s.removed && (kappa != 0)) {
				clear_tag(ev.node);
				queue_removal(ev.node, s, ev.time, engine);
			}
			break;

		case epidemic_event_kind::reset:
			if (!s.infected)
				throw std::logic_error("node is not infected");
			s.infected = false;

			/* if the remove rate of the node changed, re-generate waiting time */
			if (!s.removed && (kappa != 0)) {
				clear_tag(ev.node);
				queue_removal(ev.node, s, ev.time, engine);
			}
			break;

		default:
			throw std::logic_error("unknown event type");
			break;
	}
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- ACTIVITY_DRIVEN_NETWORK ----------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

activity_driven_network::activity_driven_network(std::vector<double> activity_, std::size_t m_,
                                                 double eta_sus_, double eta_inf_, double b_sus_, double b_inf_,
                                                 rng_t &engine_)
    : next_reaction_network()
    , activity(std::move(activity_))
    , eta_sus(eta_sus_)
    , eta_inf(eta_inf_)
    , b_sus(b_sus_)
    , b_inf(b_inf_)
    , m(m_)
    , engine(engine_)
    , state(activity.size(), node_state{ .active = false, .infected = false })
{
    /* Allocate space in adjacency lists for N nodes */
    resize((node_t)activity.size());

    /* Initialze nodes */
    for (node_t n = 0; n < (node_t)activity.size(); ++n) {
        /* Compute activation rate ai, and steady-state probability p_active */
        const double ai       = eta_sus * activity[n];
        const double p_active = ai / (ai + b_sus);
        /* Set state according to steady-state probability */
        const bool active = std::bernoulli_distribution(p_active)(engine);
        if (active)
            activate_node(n, 0);
        else if (ai > 0)
            queue_activation(n, state[n], 0);
    }

    /* Execute the network events queued for time zero to finish initialization */
    while (next(engine, 0) == 0)
        step(engine, 0);
}

void activity_driven_network::queue_activation(node_t node, node_state s, absolutetime_t time)
{
    /* Compute activation rate */
    const double eta = s.infected ? eta_inf : eta_sus;
    const double ai  = eta * activity[node];
    /* Activate immediately if rate is infinite, otherwise generate waiting time and queue */
    if (ai == INFINITY)
        activate_node(node, time);
    else if (ai > 0) {
        const double t = time + std::exponential_distribution<>(ai)(engine);
        queue_callback(t, [this, node, t]() { this->activate_node(node, t); }, node);
    }
}

void activity_driven_network::queue_deactivation(node_t node, node_state s, absolutetime_t time)
{
    /* Compute deactivation rate */
    const double b = s.infected ? b_inf : b_sus;
    /* Deactivate immediately if rate is infinite, otherwise generate waiting time and queue */
    if (b == INFINITY)
        deactivate_node(node, time);
    else if (b > 0) {
        const double t = time + std::exponential_distribution<>(b)(engine);
        queue_callback(t, [this, node, t]() { this->deactivate_node(node, t); }, node);
    }
}

void activity_driven_network::activate_node(node_t node, absolutetime_t time)
{
    /* Fetch and update node state */
    node_state &s = state[node];
    if (s.active)
        throw std::logic_error("node is already active");
    s.active = true;

    /* Sample m nodes from possible neighbours {1,...,n-1,n+1,...,N} */
    for (std::size_t nn : sample_without_replacement(nodes() - 1, m, engine)) {
        /* Skip node itself */
        if ((node_t)nn >= node)
            ++nn;
        queue_add_edge(node, (node_t)nn, 1.0, time);
    }

    /* Queue next deactivaton */
    queue_deactivation(node, s, time);
}

void activity_driven_network::deactivate_node(node_t node, absolutetime_t time)
{
    /* Fetch and update node state */
    node_state &s = state[node];
    if (!s.active)
        throw std::logic_error("node is already inactive");
    s.active = false;

    node_t nn;
    for (int i = 0; (nn = neighbour(node, i)) >= 0; ++i)
        queue_remove_edge(node, nn, 1.0, time);

    /* Queue next activation */
    queue_activation(node, s, time);
}

void activity_driven_network::notify_epidemic_event(epidemic_event_t ev, rng_t &engine)
{
    /* Fetch node state */
    node_state &s = state[ev.node];

    /* Update infected state of node, re-genereate waiting times if necessary */
    switch (ev.kind) {
        case epidemic_event_kind::outside_infection:
        case epidemic_event_kind::infection:
            if (s.infected)
                throw std::logic_error("node is already infected");
            s.infected = true;

            /* if the deactivation rate of the node changed, re-generate waiting time */
            if (s.active && (b_sus != b_inf)) {
                clear_tag(ev.node);
                queue_deactivation(ev.node, s, ev.time);
            }
            break;

        case epidemic_event_kind::reset:
            if (!s.infected)
                throw std::logic_error("node is not infected");
            s.infected = false;

            /* if the deactivation rate of the node changed, re-generate waiting time */
            if (s.active && (eta_sus != eta_inf)) {
                clear_tag(ev.node);
                queue_activation(ev.node, s, ev.time);
            }
            break;

        default:
            throw std::logic_error("unknown event type");
            break;
    }
}
