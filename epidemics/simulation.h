#pragma once

#include "stdafx.h"
#include "types.h"
#include "graph.h"

//--------------------------------------
//-------------SIMULATION---------------
//--------------------------------------

class simulate_next_reaction {
public:
    graph& network;
    std::unordered_set<node_t> infected;
    
    simulate_next_reaction(class graph& nw)
        :network(nw)
    {}
    
    void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v);
    
    std::pair<node_t, absolutetime_t> step(rng_t& engine);
    
private:
    struct infectiontimes_entry {
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
         * original infection time, and the node's neighbour
         * index within the souce node
         */
        absolutetime_t source_time = INFINITY;
        node_t source_node = -1;
        int neighbour_index = -1;
        
        bool operator< (const infectiontimes_entry& o) const { return time < o.time; }
        bool operator<= (const infectiontimes_entry& o) const { return time <= o.time; }
        bool operator== (const infectiontimes_entry& o) const { return time == o.time; }
        bool operator!= (const infectiontimes_entry& o) const { return time != o.time; }
        bool operator>= (const infectiontimes_entry& o) const { return time >= o.time; }
        bool operator> (const infectiontimes_entry& o) const { return time > o.time; }
    };
    
    std::priority_queue<infectiontimes_entry, std::deque<infectiontimes_entry>,
                        std::greater<infectiontimes_entry>>
      infectiontimes;
};

class simulate_nmga {
private:
    struct active_edges_entry {
        absolutetime_t source_time;
        node_t source;
        node_t target;
    };
    double current_time;
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
    double lambda_total(std::vector<double>* lambda_finite_cumulative = nullptr,
                        std::vector<std::size_t>* lambda_infinite = nullptr);
    
public:
    graph& network;
    beta& infection_time_distribution;
    int approximation_threshold;
    double tau_precision;
    std::unordered_set<node_t> infected;
    
    simulate_nmga(class graph& nw, class beta& dist, int threshold = 100, double tauprec = 1e-6)
        :network(nw), infection_time_distribution(dist)
        ,approximation_threshold(threshold), tau_precision(tauprec)
    {}
    
    void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v);
    
    std::pair<node_t, absolutetime_t> step(rng_t& engine);
    interval_t next_time_exact(rng_t& engine);
    interval_t next_time_approximation(rng_t& engine);
};

/* Simulates path */

std::vector<double> simulatePath(std::vector<double>& infection_times, int n_max, const lognormal_beta& infection_distribution, rng_t engine);





void print_matrix(std::vector<std::vector<double>>& A);

void simulateManyPaths(int nb_paths, rng_t& engine);
