//
//  dynamic_graph.hpp
//  epidemics
//
//  Created by Florian G. Pflug on 14.02.23.
//

#pragma once

#include "nextnet/stdafx.h"
#include "nextnet/types.h"
#include "nextnet/utility.h"
#include "nextnet/random.h"
#include "nextnet/network.h"
#include "nextnet/weighted_network.h"

/**
 * @brief Temporal network interface, i.e. networks which evolve over time
 *
 * Implementations must provide next() which determines the time of the next event,
 * and step() which performs the corresponding change. Temporal networks may react
 * to the epidemic state, see notify_epidemic_event().
 */
struct temporal_network : public virtual network
{
    virtual absolutetime_t next(rng_t &engine, absolutetime_t maxtime = INFINITY) = 0;

    virtual std::optional<network_event_t> step(rng_t &engine, absolutetime_t max_time = NAN) = 0;

    virtual void notify_epidemic_event(epidemic_event_t ev, rng_t &engine);
};

/**
 * @brief Helper class for implementing mutable networks
 *
 * Provides basic mutators like add_edge() and remove_edge(). This avoids having
 * to implement these in every temporal networks, since the basic temporal network
 * interface is agnostic to network represenations and does not provide mutators.
 */
struct mutable_network : public virtual network
{
    void resize(node_t nodes);

    bool has_edge(node_t src, node_t dst);

    bool add_edge(node_t src, node_t dst);

    bool remove_edge(node_t src, node_t dst);

    virtual node_t nodes();

    virtual node_t neighbour(node_t node, int neighbour_index);

    virtual index_t outdegree(node_t node);

private:
    std::vector<indexed_set<node_t>> adjacencylist;
};

/**
 * @brief Helper class for implementing mutable weighted networks
 *
 * Provides basic mutators like add_edge() and remove_edge(). This avoids having
 * to implement these in every temporal networks, since the basic temporal network
 * interface is agnostic to network represenations and does not provide mutators.
 */
struct mutable_weighted_network : public virtual network
    , public virtual weighted_network
{
    void resize(node_t nodes);

    bool has_edge(node_t src, node_t dst);

    bool has_edge(node_t src, node_t dst, double *weight);

    double edge_weight(node_t src, node_t dst);

    void add_edge(node_t src, node_t dst, double weight);

    bool remove_edge(node_t src, node_t dst);

    virtual node_t nodes();

    virtual node_t neighbour(node_t node, int neighbour_index);

    virtual node_t neighbour(node_t node, int neighbour_index, double *weight);

    virtual index_t outdegree(node_t node);

private:
    std::vector<indexed_map<node_t, double>> adjacencylist;
};

struct by_time
{
};
struct by_tag
{
};

/**
 * @brief Base class for nextreaction-based temporal networks.
 *
 * Network events can be queued and are reported in order of ascending time.
 * In addition to network events, arbitrary callbacks can be queued and are executed
 * at the appropriate time.
 */
struct next_reaction_network : public virtual network
    , public virtual weighted_network
    , public virtual temporal_network
    , virtual mutable_weighted_network
{
    /**
     * @brief Bitfield describing the network properties (all constants must be powers of 2!)
     */
    enum network_kind {
        directed_kind = 1,
        weighted_kind = 2
    };

    typedef std::function<void()> event_callback_t;

    typedef std::uint64_t tag_t;

protected:
    struct queue_entry
    {
        absolutetime_t time;
        tag_t tag;
        std::variant<
            network_event_t,
            event_callback_t>
            event;

        bool operator<(const queue_entry &o) const { return time < o.time; }
        bool operator<=(const queue_entry &o) const { return time <= o.time; }
        bool operator==(const queue_entry &o) const { return time == o.time; }
        bool operator!=(const queue_entry &o) const { return time != o.time; }
        bool operator>=(const queue_entry &o) const { return time >= o.time; }
        bool operator>(const queue_entry &o) const { return time > o.time; }
    };

public:
    next_reaction_network()
        : next_reaction_network((network_kind)0)
    {
    }

    next_reaction_network(network_kind kind_);

    /**
     * @brief Adds an edge from src to dist with weight w at time t
     *
     * If the network is undirected, this will also add the reverse edge
     * from dst to src with the same weight at the same time.
     */
    void queue_add_edge(node_t src, node_t dst, double w, absolutetime_t time, tag_t tag = 0);

    /**
     * @brief Removed the edge from src to dist at time t
     *
     * If the network is undirected, this will also remove the reverse edge
     * from dst to src at the same time.
     */
    void queue_remove_edge(node_t src, node_t dst, double w, absolutetime_t time, tag_t tag = 0);

    /**
     * @brief Instantenous contact from src to dst with weight w at time t
     *
     * If the network is undirected, this will also create a contact in the
     * reverse direction from dst to src at the same time
     */
    void queue_instantenous_contact(node_t src, node_t dst, double w, absolutetime_t time, tag_t tag = 0);

    /**
     * @brief Queue network event
     *
     * This queues arbitrary network events
     */
    void queue_network_event(network_event_t ev, tag_t tag = 0);

    /**
     * @brief Queue callback event
     *
     * This queues arbitrary callbacks which are executed when next() or step()
     * encounters them in the queue.
     */
    void queue_callback(absolutetime_t time, event_callback_t cb, tag_t tag = 0);

    void clear_tag(tag_t tag);

    virtual bool is_undirected() override;

    virtual bool is_unweighted() override;

    virtual absolutetime_t next(rng_t &engine, absolutetime_t maxtime = INFINITY) override;

    virtual std::optional<network_event_t> step(rng_t &engine, absolutetime_t max_time = NAN) override;

    const network_kind kind;

    double next_time = NAN;

    double current_time = 0.0;

    bool flipped_edge_next = false;

private:
    typedef bmi::multi_index_container<
        queue_entry,
        bmi::indexed_by<
            bmi::ordered_non_unique<
                bmi::tag<by_time>,
                bmi::member<queue_entry, absolutetime_t, &queue_entry::time>>,
            bmi::hashed_non_unique<
                bmi::tag<by_tag>,
                bmi::member<queue_entry, uint64_t, &queue_entry::tag>>>>
        queue_type;

    queue_type event_queue;

    const queue_entry &top_event() { return *event_queue.get<by_time>().begin(); }

    void pop_event()
    {
        auto &i = event_queue.get<by_time>();
        i.erase(i.begin());
    }

    void push_event(queue_entry e) { event_queue.insert(e); };

#if 0
    std::priority_queue<queue_entry, std::deque<queue_entry>,
                        std::greater<queue_entry>>
        event_queue;

    const queue_entry &top_event() { return event_queue.top(); };

    void pop_event() { event_queue.pop(); };

    void push_event(queue_entry e) { event_queue.push(e); };

#elif 0
    rollbear::prio_queue<32, absolutetime_t, queue_entry>
        event_queue;

    const queue_entry &top_event() { return event_queue.top().second; };

    void pop_event() { event_queue.pop(); };

    void push_event(queue_entry e) { event_queue.push(e.time, e); };
#endif
};

/**
 * @brief Temporal network defined by short contacts at pre-determined times
 *
 * Contact times are read from a file, and either represent instantenous contacts,
 * or contacts of some finite (usually short) duration dt.
 */
struct empirical_contact_network : public virtual next_reaction_network
{
    enum edge_duration_kind {
        finite_duration      = 1,
        infitesimal_duration = 2
    };

    empirical_contact_network(std::string file,
                              edge_duration_kind contact_type, double dt)
        : empirical_contact_network(file, (network_kind)0, contact_type, dt)
    {
    }

    empirical_contact_network(std::string file, network_kind,
                              edge_duration_kind contact_type, double dt);

    std::vector<std::vector<double>> compute_number_of_edges(rng_t &engine);
};

/**
 * @brief Temporal network modeled after the SIRX network of Maier & Brockmann, 2020
 *
 * Nodes in a arbitrary underlying network are removed with a rate based on their
 * infection state. Non-infected nodes have removal rate kappa0, infected nodes have
 * elevated rate kappa + kappa0.
 */
struct temporal_sirx_network : public virtual temporal_network
{
    enum node_state_t { S = 1,
                        I = 2,
                        R = 3,
                        X = 4 };

    temporal_sirx_network(network &network, double kappa0, double kappa);

    virtual bool is_undirected() override;

    virtual node_t nodes() override;

    virtual node_t neighbour(node_t node, int neighbour_index) override;

    virtual index_t outdegree(node_t node) override;

    virtual void notify_epidemic_event(epidemic_event_t ev, rng_t &engine) override;

    virtual absolutetime_t next(rng_t &engine, absolutetime_t maxtime = INFINITY) override;

    virtual std::optional<network_event_t> step(rng_t &engine, absolutetime_t max_time = NAN) override;

    bool is_removed(node_t node) { return (nonremoved.find(node) == nonremoved.end()); }

    bool is_infected(node_t node) { return (infected.find(node) != infected.end()); }

    node_state_t state(node_t node);

    network &network;

    const bool network_is_undirected;

    const node_t network_size;

    const double kappa0;

    const double kappa;

    indexed_set<node_t> nonremoved;

    indexed_set<node_t> infected_nonremoved;

    std::unordered_set<node_t> infected;

    std::deque<network_event_t> queue;

    bool queue_next_flipped = false;

    double current_time = 0.0;

    double next_time = NAN;
};

/**
 * @brief Dynamic version of an Erdös-Reyni graph
 *
 * Each edge appears and disappears independently according to a two-state
 * Markov process with rate alpha for an edge appearing and rate beta for
 * the rate disappearing. The dynamic Erdös-Reyni graph is parametrized
 * in terms of the number n of nodes, the average degree k of a node, and
 * the timescale tau on which edges appear and disappear. In terms of these
 * parameter, each edge has probability p_+ = k / (n - 1) to exist at at
 * certain point in time, which is satiesfied for rates of appearance and
 * disappearance of alpha = p_+ / tau and beta = p_- / tau = (1 - p_+) / tau.
 *
 * To see this, consider a two-state Markov process with states + (present)
 * and - (absent) with rates alpha for the transition - -> +, and beta for
 * + -> -. The steady-state probabilites are then p_+ = alpha / (alpha + beta),
 * and p_- = beta / (alpha + beta).
 */
struct temporal_erdos_renyi : public virtual temporal_network
    , public virtual erdos_renyi
{
    temporal_erdos_renyi(int size, double avg_degree, double timescale, rng_t &engine);

    virtual absolutetime_t next(rng_t &engine, absolutetime_t maxtime = INFINITY) override;

    virtual std::optional<network_event_t> step(rng_t &engine, absolutetime_t max_time = NAN) override;

    void add_edge(node_t node, node_t neighbour);
    void remove_edge(node_t node, int neighbour_index);

    /* Probability that a particular edge is present in the steady-state */
    const double edge_probability;
    /* Rate with which absent edges appear */
    const double alpha;
    /* Rate with which present edges vanish */
    const double beta;

    /* Current time (time of last event) */
    absolutetime_t current_time = 0.0;

    /* Time of next edge appearing or vanishing */
    absolutetime_t next_time = NAN;

    /* Unreported event for the reverse edge of the last edge event reported  */
    std::optional<network_event_t> reverse_edge_event;

    unsigned int edges_absent;
    unsigned int edges_present;

    /* Degree-weighted node distribution */
    dyndist::vector_distribution<unsigned> weighted_nodes;
};

/* Compatibility with previously miss-spelled name */
typedef temporal_erdos_renyi temporal_erdos_reyni;

/**
 * @brief Activity-driven network model of Cai, Nie & Holme, 2024.
 *
 * Here, nodes are initially inactive and have degree zero. Node i activates with
 * rate a[i] * eta and upon activation connects to m other uniformly chosen nodes
 * (which are not necessarily active). Active nodes inactivate with constant rate
 * recover_rate.
 */
struct activity_driven_network : virtual public next_reaction_network
{
    enum activity_event_kind {
        none       = 0,
        activate   = 1,
        deactivate = 2
    };

    /* bitfield used to store node states */
    struct node_state
    {
        unsigned char active : 1;
        unsigned char infected : 1;
    };

    activity_driven_network(std::vector<double> activity_rates, double eta, double m, double b, rng_t &engine)
        : activity_driven_network(std::move(activity_rates), m, eta, eta, b, b, engine)
    {
    }

    activity_driven_network(std::vector<double> activity_rates, std::size_t m,
                            double eta_sus, double eta_inf, double b_sus, double b_inf,
                            rng_t &engine);

    virtual void notify_epidemic_event(epidemic_event_t ev, rng_t &engine) override;

    void activate_node(node_t node, absolutetime_t time);
    void deactivate_node(node_t node, absolutetime_t time);

    void queue_activation(node_t node, node_state s, absolutetime_t time);
    void queue_deactivation(node_t node, node_state s, absolutetime_t time);

    const std::vector<double> activity;
    const double eta_sus;
    const double eta_inf;
    const double b_sus;
    const double b_inf;
    const std::size_t m;

    rng_t &engine;

    std::vector<node_state> state;
};
