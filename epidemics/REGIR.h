#pragma once

#include "stdafx.h"
#include "types.h"
#include "algorithm.h"
#include "graph.h"

class simulate_regir : public simulation_algorithm  {
private:
    struct active_edges_entry {
		/*
		 * Event kind represented by this edge (infection or reset)
		 * Reset events are self-loops, a consequently they obey
		 * source_node=-1, neighbour_index=-1, neighbours_remaining=0.
		 */
		event_kind kind = event_kind::none;

		absolutetime_t source_time = INFINITY;
        node_t source = -1;
        node_t target = -1;
        double lambda = NAN;
    };

	/*
    struct active_edges_hasher {
        std::size_t operator()(const active_edges_entry& e) const {
            return hash_combine(0, e.kind, e.source, e.target);
        }
    };

    struct active_edges_cmp {
        bool operator()(const active_edges_entry& a, const active_edges_entry& b) const {
            return (a.kind == b.kind) & (a.source == b.source) && (a.target == b.target);
        }
    };

    typedef std::unordered_set<active_edges_entry, active_edges_hasher, active_edges_cmp>
        active_edges_t;
	 */
	
	typedef std::deque<active_edges_entry> active_edges_t;

    struct outside_infections_entry {
        absolutetime_t time;
        node_t node;

        bool operator< (const outside_infections_entry& o) const { return time < o.time; }
        bool operator<= (const outside_infections_entry& o) const { return time <= o.time; }
        bool operator== (const outside_infections_entry& o) const { return time == o.time; }
        bool operator!= (const outside_infections_entry& o) const { return time != o.time; }
        bool operator>= (const outside_infections_entry& o) const { return time >= o.time; }
        bool operator> (const outside_infections_entry& o) const { return time > o.time; }
    };

    typedef std::priority_queue<outside_infections_entry, std::deque<outside_infections_entry>,
                                std::greater<outside_infections_entry>>
        outside_infections_t;

	const double lambda_max;
    double current_time = NAN;
	double lambda_total = NAN;
    outside_infections_t outside_infections;
    active_edges_t active_edges;
    std::vector<double> active_edges_cumulative_finite_lambdas;
    std::vector<std::size_t> active_edges_infinite_lambdas;
    std::uniform_real_distribution<double> unif01_dist;
    std::optional<event_t> next_event;
	std::optional<std::size_t> next_event_edge_pos;
    
public:
	struct params {
		params() {};
		
		int approximation_threshold = 100;
		double tau_precision = 1e-6;
		bool SIR = false;
	};
	
    graph& network;
    const class transmission_time& psi;
	const class transmission_time* rho = nullptr;
	const params p;
	std::unordered_set<node_t> infected;
	
	simulate_regir(graph& nw, const class transmission_time& psi_,
				   const class transmission_time* rho_ = nullptr,
				   params p_ = params())
        :lambda_max(std::max(psi_.hazardbound(INFINITY),
							 (rho_ ? rho_->hazardbound(INFINITY) : 0)))
	    ,network(nw), psi(psi_), rho(rho_)
        ,p(p_)
    {}

	virtual graph& get_network() const;

	virtual const class transmission_time& transmission_time() const;

	virtual const class transmission_time* reset_time() const;

	virtual bool is_infected(node_t) const;
	
	virtual void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v);
	
	virtual absolutetime_t next(rng_t& engine);

	virtual std::optional<event_t> step(rng_t& engine, absolutetime_t maxtime = INFINITY,
										event_filter_t event_filter = std::nullopt);

	virtual void notify_infected_node_neighbour_added(network_event_t event, rng_t& engine);

private:
	int removed = 0;

	void add_active_edge(const active_edges_entry& e);
	
	void remove_active_edge(active_edges_t::iterator& it);
	
	double phi(absolutetime_t t, interval_t tau);

	interval_t invphi(absolutetime_t t, double u);

	double edge_hazard_rate(const active_edges_entry& e, double time);

	void update_active_edge_lambdas(double time);
	
	active_edges_t::iterator draw_active_edge_hazardrates(rng_t& engine);

	active_edges_t::iterator draw_active_edge_uniform(rng_t& engine);

	absolutetime_t next_time_outside_infection(absolutetime_t max_time);
	
	interval_t next_time_exact(rng_t& engine);

	interval_t next_time_exponential(rng_t& engine);
};
