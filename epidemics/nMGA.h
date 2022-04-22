#pragma once

#include "stdafx.h"
#include "types.h"
#include "graph.h"

class simulate_nmga {
private:
    struct active_edges_entry {
        absolutetime_t source_time;
        node_t source;
        node_t target;
        double lambda;
    };

    double current_time;

    double lambda_total;

    std::vector<active_edges_entry> active_edges;

    std::vector<double> active_edges_cumulative_finite_lambdas;
    std::vector<std::size_t> active_edges_infinite_lambdas;
    std::uniform_real_distribution<double> unif01_dist;
    
    void add_active_edge(const active_edges_entry& e) {
        active_edges.push_back(e);
    }
    
    void remove_active_edge(std::vector<active_edges_entry>::iterator it) {
        /* We remove the element pointed to by <it> by swapping it with
         * the last element and then removing the last element. This is
         * faster then moving all element to the right of <it>. */
        *it = active_edges.back();
        active_edges.pop_back();
    }
    
    double phi(absolutetime_t t, interval_t tau);

    interval_t invphi(absolutetime_t t, double u);

    void update_active_edge_lambdas();

    std::vector<active_edges_entry>::iterator draw_active_edge(rng_t& engine);
    
public:
    graph& network;
    transmission_time& psi;
    int approximation_threshold = 100;
    double tau_precision = 1e-6;
    std::unordered_set<node_t> infected;
    
    simulate_nmga(class graph& nw, class transmission_time& psi_, int threshold = 100, double tauprec = 1e-6)
        :network(nw), psi(psi_)
        ,approximation_threshold(threshold), tau_precision(tauprec)
    {}
    
    void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v);
    
    std::pair<node_t, absolutetime_t> step(rng_t& engine);

    interval_t next_time_exact(rng_t& engine);

    interval_t next_time_approximation(rng_t& engine);
};

