//
//  dynamic_graph.cpp
//  epidemics
//
//  Created by Florian G. Pflug on 14.02.23.
//

#include "temporal_network.h"

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
}

absolutetime_t next_reaction_network::next(rng_t &engine, absolutetime_t maxtime)
{
    /* If there's a cached next_time, return it */
    if (!std::isnan(next_time))
        return next_time;

    /* Find time of next event that step() will report */
    while (true) {
        if (event_queue.empty())
            return INFINITY;

        /* Fetch next event, return nothing if past max_time */
        const queue_entry &top = top_event();
        if (top.time > maxtime)
            return INFINITY;

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
        return next_time;
    }
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
/*----------- TEMPORAL_SIRX_NETWORK ------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

temporal_sirx_network::temporal_sirx_network(class network &network_, double kappa0_, double kappa_)
    : network(network_)
    , network_is_undirected(network.is_undirected())
    , network_size(network.nodes())
    , kappa0(kappa0_)
    , kappa(kappa_)
    , queue_next_flipped(is_undirected())
{
    for (node_t n = 0; n < network_size; ++n)
        nonremoved.insert(n);
}

temporal_sirx_network::node_state_t temporal_sirx_network::state(node_t node)
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

bool temporal_sirx_network::is_undirected()
{
    return network_is_undirected;
}

node_t temporal_sirx_network::nodes()
{
    return network_size;
}

node_t temporal_sirx_network::neighbour(node_t node, int neighbour_index)
{
    if (is_removed(node))
        return -1;
    return network.neighbour(node, neighbour_index);
}

index_t temporal_sirx_network::outdegree(node_t node)
{
    if (is_removed(node))
        return 0;
    return network.outdegree(node);
}

void temporal_sirx_network::notify_epidemic_event(epidemic_event_t ev, rng_t &engine)
{
    assert(current_time <= ev.time);
    assert(std::isnan(next_time) || ev.time <= next_time);
    switch (ev.kind) {
        case epidemic_event_kind::infection:
        case epidemic_event_kind::outside_infection:
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
        case epidemic_event_kind::reset:
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

absolutetime_t temporal_sirx_network::next(rng_t &engine, absolutetime_t)
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
        const double r  = infected_nonremoved.size() * kappa;
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
        node_t n        = -1;
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
        for (int i = 0, nn = 0; (nn = network.neighbour(n, i)) >= 0; ++i) {
            queue.push_back(network_event_t{
                .kind        = network_event_kind::neighbour_removed,
                .source_node = n,
                .target_node = nn,
                .weight      = 1.0,
                .time        = next_time });
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

std::optional<network_event_t> temporal_sirx_network::step(rng_t &engine, absolutetime_t max_time)
{
    /* make sure the next event was generated, do nothing if past max_time */
    next(engine);
    if (next_time > max_time)
        return std::nullopt;

    /* get next event off queue, it's time must be next_time */
    assert(!queue.empty());
    network_event_t ev = queue.front();
    const node_t node  = ev.source_node;
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
    assert(edges_absent % 2 == 0);
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
/*----------- ACTIVITY_DRIVEN_NETWORK ----------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

activity_driven_network::activity_driven_network(std::vector<double> activity, double eta, double m, double b, rng_t &engine)
    : next_reaction_network()
    , activity(activity)
    , eta(eta)
    , m(m)
    , b(b)
    , engine(engine)
    , active(activity.size(), false)
{
    resize(activity.size());

    for (node_t n = 0; n < (node_t)activity.size(); ++n) {
        const double ai       = eta * activity[n];
        const double p_active = ai / (ai + b);
        const bool active     = std::bernoulli_distribution(p_active)(engine);
        if (active && (b > 0)) {
            activate_node(n, 0);
        } else if (ai > 0) {
            const double t = std::exponential_distribution<>(ai)(engine);
            queue_callback(t, [this, n, t]() { this->activate_node(n, t); });
        }
    }

    /* Execute the events queued for time zero */
    while (next(engine, 0) == 0)
        step(engine, 0);
}

void activity_driven_network::activate_node(node_t node, absolutetime_t time)
{
    if (active[node])
        throw std::logic_error("node is already active");

    active[node] = true;

    /* Sample m nodes from possible neighbours {1,...,n-1,n+1,...,N} */
    for (std::size_t nn : sample_without_replacement(nodes() - 1, m, engine)) {
        /* Skip node itself */
        if ((node_t)nn >= node)
            ++nn;
        queue_add_edge(node, (node_t)nn, 1.0, time);
    }

    /* Queue next deactivaton */
    const double t = time + std::exponential_distribution<>(b)(engine);
    queue_callback(t, [this, node, t]() { this->deactivate_node(node, t); });
}

void activity_driven_network::deactivate_node(node_t node, absolutetime_t time)
{
    if (!active[node])
        throw std::logic_error("node is already inactive");

    active[node] = false;

    node_t nn;
    for(int i = 0; (nn = neighbour(node, i)) >= 0; ++i)
        queue_remove_edge(node, nn, 1.0, time);

    /* Queue next activation */
    const double ai = eta * activity[node];
    const double t  = time + std::exponential_distribution<>(ai)(engine);
    queue_callback(t, [this, node, t]() { this->activate_node(node, t); });
}
