#pragma once

#include "nextnet/stdafx.h"
#include "nextnet/types.h"
#include "nextnet/algorithm.h"
#include "nextnet/network.h"

struct all_rates_zero : std::underflow_error
{
    all_rates_zero()
        : underflow_error("all active edges report a hazardrate of zero"){};
};

class simulate_nmga : public simulation_algorithm
{
private:
    struct active_edges_entry
    {
        /*
         * Event kind represented by this edge (infection or reset)
         * Reset events are self-loops, a consequently they obey
         * source_node=-1, neighbour_index=-1, neighbours_remaining=0.
         */
        epidemic_event_kind kind = epidemic_event_kind::none;

        absolutetime_t source_time = INFINITY;
        node_t source              = -1;
        node_t target              = -1;
        double lambda              = NAN;
    };

    struct active_edges_hasher
    {
        std::size_t operator()(const active_edges_entry &e) const
        {
            return hash_combine(0, e.kind, e.source, e.target);
        }
    };

    struct active_edges_cmp
    {
        bool operator()(const active_edges_entry &a, const active_edges_entry &b) const
        {
            return (a.kind == b.kind) & (a.source == b.source) && (a.target == b.target);
        }
    };

    typedef std::unordered_set<active_edges_entry, active_edges_hasher, active_edges_cmp>
        active_edges_t;

    struct outside_infections_entry
    {
        absolutetime_t time;
        node_t node;

        bool operator<(const outside_infections_entry &o) const { return time < o.time; }
        bool operator<=(const outside_infections_entry &o) const { return time <= o.time; }
        bool operator==(const outside_infections_entry &o) const { return time == o.time; }
        bool operator!=(const outside_infections_entry &o) const { return time != o.time; }
        bool operator>=(const outside_infections_entry &o) const { return time >= o.time; }
        bool operator>(const outside_infections_entry &o) const { return time > o.time; }
    };

    typedef std::priority_queue<outside_infections_entry, std::deque<outside_infections_entry>,
                                std::greater<outside_infections_entry>>
        outside_infections_t;

    double current_time = NAN;
    double lambda_total = NAN;
    outside_infections_t outside_infections;
    active_edges_t active_edges;
    std::vector<double> active_edges_cumulative_finite_lambdas;
    std::vector<std::size_t> active_edges_infinite_lambdas;
    std::uniform_real_distribution<double> unif01_dist;
    std::optional<epidemic_event_t> next_event;

public:
    struct params
    {
        params(){};

        int approximation_threshold = 100;
        double maximal_dt           = NAN;
        double tau_precision        = 1e-6;
        bool SIR                    = false;
    };

    network &nw;
    const class transmission_time &psi;
    const class transmission_time *rho = nullptr;
    const params p;
    const double max_dt;
    std::unordered_set<node_t> infected;

    simulate_nmga(class network &nw_, const class transmission_time &psi_,
                  const class transmission_time *rho_ = nullptr,
                  params p_                           = params())
        : nw(nw_)
        , psi(psi_)
        , rho(rho_)
        , p(p_)
        , max_dt((std::isfinite(p.maximal_dt) && (p.maximal_dt > 0)) ? p.maximal_dt : find_maximal_dt(psi))

    {
    }

    virtual class network &get_network() const override;

    virtual const class transmission_time &transmission_time() const override;

    virtual const class transmission_time *reset_time() const override;

    virtual bool is_infected(node_t) const override;

    virtual void add_infections(const std::vector<std::pair<node_t, absolutetime_t>> &v) override;

    virtual absolutetime_t next(rng_t &engine) override;

    virtual std::optional<epidemic_event_t> step(rng_t &engine, absolutetime_t maxtime = INFINITY,
                                                 event_filter_t event_filter = std::nullopt) override;

    virtual void notify_infected_node_neighbour_added(network_event_t event, rng_t &engine) override;

private:
    int removed = 0;

    static double find_maximal_dt(const class transmission_time &psi);

    void add_active_edge(const active_edges_entry &e)
    {
        active_edges.insert(e);
    }

    void remove_active_edge(active_edges_t::iterator &it)
    {
        it = active_edges.erase(it);
    }

    double phi(absolutetime_t t, interval_t tau);

    interval_t invphi(absolutetime_t t, double u);

    void update_active_edge_lambdas(double time);

    active_edges_t::iterator draw_active_edge(rng_t &engine);

    interval_t next_time_exact(rng_t &engine);

    interval_t next_time_approximation(rng_t &engine);
};
