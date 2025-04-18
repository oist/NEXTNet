#pragma once

#include "nextnet/stdafx.h"
#include "nextnet/types.h"
#include "nextnet/algorithm.h"
#include "nextnet/random.h"

class simulate_next_reaction_mean_field : public simulation_algorithm
{
public:
    simulate_next_reaction_mean_field(int N_, double R0_, const class transmission_time &psi_,
                                      const class transmission_time *rho_ = nullptr, bool SIR_ = false)
        : N(N_)
        , R0(R0_)
        , psi(psi_)
        , rho(rho_)
        , SIR(SIR_)
    {
    }

    virtual const class transmission_time &transmission_time() const override;

    virtual const class transmission_time *reset_time() const override;

    virtual void add_infections(const std::vector<std::pair<node_t, absolutetime_t>> &v) override;

    virtual absolutetime_t next(rng_t &engine) override;

    virtual std::optional<epidemic_event_t> step(rng_t &engine, absolutetime_t maxtime = INFINITY,
                                                 event_filter_t event_filter = std::nullopt) override;

    virtual void notify_infected_node_neighbour_added(network_event_t event, rng_t &engine) override;

    virtual bool is_infected(node_t) const override;

    const int N;

    const double R0;

    virtual network &get_network() const override;

    int current_nb_of_infected()
    {
        return (int)infected.size() - removed;
    }

private:
    const double p = R0 / (N - 1);
    const class transmission_time &psi;
    const class transmission_time *rho = nullptr;
    std::unordered_set<node_t> infected;
    bool SIR    = false;
    int removed = 0; // number of individuals that have recovered and cannot be re-infected. (only active in the SIR case).

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

        bool operator<(const active_edges_entry &o) const { return time < o.time; }
        bool operator<=(const active_edges_entry &o) const { return time <= o.time; }
        bool operator==(const active_edges_entry &o) const { return time == o.time; }
        bool operator!=(const active_edges_entry &o) const { return time != o.time; }
        bool operator>=(const active_edges_entry &o) const { return time >= o.time; }
        bool operator>(const active_edges_entry &o) const { return time > o.time; }
    };

    std::priority_queue<active_edges_entry, std::deque<active_edges_entry>,
                        std::greater<active_edges_entry>>
        active_edges;

    std::optional<epidemic_event_t> step_infection(const active_edges_entry &next, rng_t &engine);

    std::optional<epidemic_event_t> step_reset(const active_edges_entry &next, rng_t &engine);
};
