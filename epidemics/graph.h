//
//  graph.hpp
//  epidemics
//
//  Created by Cure Samuel Cyrus on 2022/03/07.
//

#pragma once

#include "stdafx.h"
#include "random.h"

typedef int node_t;
typedef double interval_t;
typedef double absolutetime_t;

class beta {
public:
    /*
     * sample() returns a single (iid) sample from the
     * (normalized versio of) the distribution.
     */
    virtual interval_t sample(rng_t& engine) const = 0;

    /*
     * sample_next() returns the next single sample from the
     * the distribution condtioned on the last infection time of the individual.
     * /!\ returning infinity is possible and represnets no further infection.
     */
    virtual interval_t sample_next(interval_t last, rng_t& engine) const = 0;
    
    /*
     * sample_next_conditional() returns same as above but takes into consideration the
     * current number of healthy neighbours.
     * i.e if there are no more healthy neighbours
     * it raises an error as no more infections are possible.
     */

    virtual interval_t sample_next_conditional(interval_t last, int healthy, rng_t& engine) const = 0;
};

class graph {
public:
    
    /*
     * returns:
     * - the ith neighbour of node (in order of infection times)
     * - the infection time from node to the ith neigbhour of node.
     */
    virtual std::pair<node_t, interval_t> neighbour(node_t node, int neighbour_index) = 0;
};

class simulator {
public:
    graph& network;
    std::unordered_set<node_t> infected;
    
    simulator(class graph& nw)
        :network(nw)
    {}
    
    void add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v);
    
    std::pair<node_t, absolutetime_t> step();
    
private:
    struct infectiontimes_entry {
        /*
         * Absolute time of infection
         */
        double time;
        
        /*
         * Node that is (putatively) infected
         */
        node_t node;
        
        /*
         * Source node that causes the node's infection, it's
         * original infection time, and the node's neighbour
         * index within the souce node
         */
        double source_time = INFINITY;
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

class lognormal_beta : public beta {
public:
    const double mean;
    const double variance;
    const double r0;
    const double mu = 2 * log(mean) - 0.5 * log( pow(mean,2.0)+ variance );
    const double sigma = sqrt( log( 1 + variance/pow(mean,2.0)));

public:
    lognormal_beta(double _mean, double _variance, double _r0)
        :mean(_mean), variance(_variance), r0(_r0)
    {}
    
    virtual interval_t sample(rng_t& engine) const;

    virtual interval_t sample_next(interval_t last, rng_t& engine) const;
    
    virtual interval_t sample_next_conditional(interval_t last, int healthy, rng_t& engine) const;

private:
    mutable std::lognormal_distribution<interval_t> log_distribution = std::lognormal_distribution<interval_t>(mu,sigma);
};

class erdos_reyni : public graph {
public:
    erdos_reyni(int size, double avg_degree, const beta& infection_distribution, rng_t& engine);

    virtual std::pair<node_t, interval_t> neighbour(node_t node, int neighbour_index);

private:
    /* Adjacency list of the graph */
    std::vector<std::vector<std::pair<node_t,interval_t>>> neighbours;
};

