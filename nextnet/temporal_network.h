//
//  dynamic_graph.hpp
//  epidemics
//
//  Created by Florian G. Pflug on 14.02.23.
//

#pragma once

#include "stdafx.h"
#include "types.h"
#include "random.h"
#include "network.h"

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
 * @brief Temporal network defined by short contacts at pre-determined times
 *
 * Contact times are read from a file, and either represent instantenous contacts,
 * or contacts of some finite (usually short) duration dt.
 */
struct empirical_contact_network : public virtual temporal_network
    , public virtual mutable_network
{
    enum edge_duration_kind {
        finite_duration      = 1,
        infitesimal_duration = 2
    };

    empirical_contact_network(std::string path_to_file, edge_duration_kind contact_type, double dt);

    virtual absolutetime_t next(rng_t &engine, absolutetime_t maxtime = INFINITY) override;

    virtual std::optional<network_event_t> step(rng_t &engine, absolutetime_t max_time = NAN) override;

    std::vector<std::vector<double>> compute_number_of_edges(rng_t &engine);

private:
    std::deque<network_event_t> event_queue;
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
struct activity_driven_network : virtual temporal_network
    , virtual mutable_network
{

    // activity_driven_network(std::vector<int> degreelist, double eta, double m, transmission_time& psi, transmission_time* rho,rng_t& engine);

    enum activity_event_kind {
        none       = 0,
        activate   = 1,
        deactivate = 2
    };

    std::vector<double> activity_rates;
    double eta;
    double m;
    double recovery_rate;
    rng_t engine;

    // Constructor to initialize the variables
    activity_driven_network(std::vector<double> activity_rates, double eta, double m, double recovery_rate, rng_t &engine);

    virtual absolutetime_t next(rng_t &engine, absolutetime_t maxtime = INFINITY);

    virtual std::optional<network_event_t> step(rng_t &engine, absolutetime_t max_time = NAN);

    void activate_node(node_t node, double time);
    void deactivate_node(node_t node, double time);

    // void to_equilibrium(rng_t& engine,absolutetime_t max_time = NAN);
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> advance_time(rng_t &engine, absolutetime_t max_time = NAN);

    std::vector<bool> active_nodes;

    int nb_edges = 0;

    /* Time of next edge appearing or vanishing */
    absolutetime_t next_time = NAN;

    struct active_edges_entry
    {
        /*
         * Event kind r (activation or desactivation)
         */
        network_event_kind kind;

        /*
         * Absolute time of event
         */
        absolutetime_t time;

        /*
         * Node that the edge is pointing to
         */
        node_t target_node;

        /*
         * Node that the edge is leaving from
         */
        node_t source_node;

        activity_event_kind activity_event;

        bool operator<(const active_edges_entry &o) const { return time < o.time; }
        bool operator<=(const active_edges_entry &o) const { return time <= o.time; }
        bool operator==(const active_edges_entry &o) const { return time == o.time; }
        bool operator!=(const active_edges_entry &o) const { return time != o.time; }
        bool operator>=(const active_edges_entry &o) const { return time >= o.time; }
        bool operator>(const active_edges_entry &o) const { return time > o.time; }
    };

#if NEXT_REACTION_QUEUE == STD_PRIORITY_QUEUE_DEQUE
    std::priority_queue<active_edges_entry, std::deque<active_edges_entry>,
                        std::greater<active_edges_entry>>
        active_edges;

    const active_edges_entry &top_edge() { return active_edges.top(); };

    void pop_edge() { active_edges.pop(); };

    void push_edge(active_edges_entry e) { active_edges.push(e); };

#elif NEXT_REACTION_QUEUE == EXT_PRIO_QUEUE
    rollbear::prio_queue<32, absolutetime_t, active_edges_entry>
        active_edges;

    const active_edges_entry &top_edge() { return active_edges.top().second; };

    void pop_edge() { active_edges.pop(); };

    void push_edge(active_edges_entry e) { active_edges.push(e.time, e); };
#endif
};
