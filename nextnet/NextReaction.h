#pragma once

#include "nextnet/stdafx.h"
#include "nextnet/types.h"
#include "nextnet/permutation.h"
#include "nextnet/algorithm.h"
#include "nextnet/random.h"
#include "nextnet/network.h"
#include "nextnet/weighted_network.h"

class simulate_next_reaction : public simulation_algorithm
{
public:
    struct params
    {
        params() noexcept {};

        bool shuffle_neighbours = true;
        bool edges_concurrent   = true;
        bool SIR                = false;
    };

    simulate_next_reaction(network &nw_, const class transmission_time &psi_,
                           const class transmission_time *rho_ = nullptr,
                           params p_                           = params())
        : nw(nw_)
        , nw_weighted(as_weighted_network(&nw))
        , psi(psi_)
        , rho(rho_)
        , p(p_)
        , shuffle_neighbours(p.shuffle_neighbours && !p.edges_concurrent)
    {
        if (!p.edges_concurrent && (nw_weighted != nullptr))
            throw std::runtime_error("sequential edges mode is not supported for concurrent networks");
    }

    virtual network &get_network() const override;

    virtual const class transmission_time &transmission_time() const override;

    virtual const class transmission_time *reset_time() const override;

    virtual void add_infections(const std::vector<std::pair<node_t, absolutetime_t>> &v) override;

    virtual absolutetime_t next(rng_t &engine) override;

    virtual std::optional<epidemic_event_t> step(rng_t &engine, absolutetime_t maxtime = INFINITY,
                                                 event_filter_t event_filter = std::nullopt) override;

    virtual void notify_infected_node_neighbour_added(network_event_t event, rng_t &engine) override;

    virtual void notify_infected_contact(network_event_t event, rng_t &engine) override;

    virtual bool is_infected(node_t) const override;

    struct infected_state_t
    {
        infected_state_t(absolutetime_t inf, absolutetime_t res)
            : infection_time(inf)
            , reset_time(res)
        {
        }

        absolutetime_t infection_time;
        absolutetime_t reset_time;
    };
    typedef std::unordered_map<node_t, infected_state_t> infected_nodes_t;
    infected_nodes_t infected;

    network &nw;
    weighted_network *nw_weighted;
    const class transmission_time &psi;
    const class transmission_time *rho = nullptr;
    const params p;
    const bool shuffle_neighbours;

    int removed = 0; // number of nodes that have recovered and cannot be re infected. (only active in the SIR case).

    int current_nb_of_infected()
    {
        return (int)infected.size() - removed;
    }

    struct active_edges_entry
    {
        /*
         * Event kind represented by this edge (infection or reset)
         * Reset events are self-loops, a consequently they obey
         * source_node=-1, neighbour_index=-1, neighbours_remaining=0.
         */
        epidemic_event_kind kind;

        /*
         * Absolute time of infection
         */
        absolutetime_t time;

        /*
         * Node that is (putatively) infected
         */
        node_t node;

        /*
         * Source node that causes the node's infection, it's
         * original infection time, it's reset time, and the
         * node's neighbour index within the souce node
         */
        absolutetime_t source_time  = INFINITY;
        node_t source_node          = -1;
        absolutetime_t source_reset = INFINITY;
        permutation<node_t> source_permutation;
        index_t neighbour_index      = -1;
        index_t neighbours_remaining = 0;
        
        /*
         * Whether the edge exited only instantaneously
         */
        bool instantaneous_edge = false;

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

    std::size_t queue_steps_total = 0;

    std::optional<epidemic_event_t> step_infection(const active_edges_entry &next, event_filter_t evf, rng_t &engine);

    std::optional<epidemic_event_t> step_reset(const active_edges_entry &next, event_filter_t evf, rng_t &engine);
};
