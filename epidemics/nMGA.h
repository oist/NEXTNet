#pragma once

#include "stdafx.h"
#include "types.h"
#include "algorithm.h"
#include "graph.h"

class simulate_nmga : public simulation_algorithm  {
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

    double current_time = NAN;

    double lambda_total = NAN;

    std::priority_queue<outside_infections_entry, std::deque<outside_infections_entry>,
                        std::greater<outside_infections_entry>>
        outside_infections;

    std::vector<active_edges_entry> active_edges;

    std::vector<double> active_edges_cumulative_finite_lambdas;
    std::vector<std::size_t> active_edges_infinite_lambdas;
    std::uniform_real_distribution<double> unif01_dist;
    
public:
    graph& network;
    const class transmission_time& psi;
	const class transmission_time* rho = nullptr;
	bool shuffle_neighbours = true;
    int approximation_threshold = 100;
    double maximal_dt = NAN;
    double tau_precision = 1e-6;
    std::unordered_set<node_t> infected;
    
	simulate_nmga(graph& nw, const class transmission_time& psi_,
				  const class transmission_time* rho_ = nullptr,
				  bool shuffle_neighbours_ = true,
				  int threshold = 100, double max_dt = NAN,
				  double tauprec = 1e-6)
        :network(nw), psi(psi_), rho(rho_)
        ,shuffle_neighbours(shuffle_neighbours_)
        ,approximation_threshold(threshold)
        ,maximal_dt((std::isfinite(max_dt) && (max_dt > 0)) ?
                    max_dt : find_maximal_dt(psi_))
        ,tau_precision(tauprec)
    {}

	virtual graph& get_network() const;

	virtual const class transmission_time& transmission_time() const;

	virtual const class transmission_time* reset_time() const;

	virtual bool is_infected(node_t) const;
	
	virtual void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v);
	
	virtual std::optional<event_t> step(rng_t& engine);

private:
	static double find_maximal_dt(const class transmission_time& psi);
	
	void add_active_edge(const active_edges_entry& e) {
		active_edges.push_back(e);
	}
	
	void remove_active_edge(std::vector<active_edges_entry>::iterator& it) {
		/* We remove the element pointed to by <it> by swapping it with
		 * the last element and then removing the last element. This is
		 * faster then moving all element to the right of <it>. */
		if (it != active_edges.end() - 1) {
			/* Removing any element that is not the last */
			*it = active_edges.back();
			active_edges.pop_back();
		} else {
			/* Removing the last element */
			active_edges.pop_back();
			it = active_edges.end();
		}
	}
	
	double phi(absolutetime_t t, interval_t tau);

	interval_t invphi(absolutetime_t t, double u);

	void update_active_edge_lambdas();

	std::vector<active_edges_entry>::iterator draw_active_edge(rng_t& engine);

	interval_t next_time_exact(rng_t& engine);

	interval_t next_time_approximation(rng_t& engine);
};

