//
//  graph.cpp
//  epidemics
//
//  Created by Cure Samuel Cyrus on 2022/03/07.
//

#include "graph.h"
#include "random.h"

std::pair<node_t, absolutetime_t> simulator::step() {
    while (true) {
        /* If there are no more infection times, stop */
        if (infectiontimes.empty())
            return std::make_pair(-1, INFINITY);

        /* Fetch the next putatively infected node */
        const auto next = infectiontimes.top();
        infectiontimes.pop();
        
        /* Lazily enqueue next sibling of the putatively infected node if necessary */
        if (next.neighbour_index >= 0) {
            const auto sibling = network.neighbour(next.source_node, next.neighbour_index+1);
            if ((sibling.first >= 0) && (std::isfinite(sibling.second))) {
                /* Create sibling's infection times entry and add to queue */
                infectiontimes_entry e;
                e.time = next.source_time + sibling.second;
                e.node = sibling.first;
                e.source_time = next.source_time;
                e.source_node = next.source_node;
                e.neighbour_index = next.neighbour_index + 1;
                infectiontimes.push(e);
            }
        }
        
        /* Check if the putatively infected node is already infected, if so we're done */
        if (infected.find(next.node) != infected.end())
            continue;
        
        /* Mark the node as infected */
        infected.insert(next.node);
        
        /* Add the infecte node's first neigbhour to the infection times queue */
        const auto neighbour = network.neighbour(next.node, 0);
        if ((neighbour.first >= 0) && (std::isfinite(neighbour.second))) {
            infectiontimes_entry e;
            e.time = next.time + neighbour.second;
            e.node = neighbour.first;
            e.source_time = next.time;
            e.source_node = next.node;
            e.neighbour_index = 0;
            infectiontimes.push(e);
        }
        
        return std::make_pair(next.node, next.time);
    }
}

void simulator::add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v) {
    for(const auto& ve: v) {
        infectiontimes_entry e;
        e.time = ve.second;
        e.node = ve.first;
        infectiontimes.push(e);
    }
}

interval_t lognormal_beta::sample(rng_t& engine) const {
    return log_distribution(engine);
}

interval_t lognormal_beta::sample_next(interval_t last, rng_t& engine) const {
    throw std::logic_error("not implemented yet");
}

interval_t lognormal_beta::sample_next_conditional(interval_t last, int healthy, rng_t& engine) const {
    throw std::logic_error("not implemented yet");
}

erdos_reyni::erdos_reyni(int size, double avg_degree, const beta& infection_distribution, rng_t& engine ){
    
    /*--------------Initialisation--------------

     Construct erdos-Reyni graph and for each link we add an infection time:
     => A[i][j] is the time node i (or j) takes to infect node j (or i) once it has itself been infected by someone else.
     
     */
    
//    std::vector<std::vector<node_t>> neighbours;
    
    const double p = avg_degree/size; // probability of an edge: if size ->infty and degree-> fixed then we get Poisson Graph.
    
    std::bernoulli_distribution has_edge(p);

    neighbours.resize(size);
    for (int i=0; i<size; i++) {
        for (int j=0; j<i; j++) {
            if (!has_edge(engine)) {
                continue;
            }
            const double tau = infection_distribution.sample(engine);
            neighbours[i].push_back(std::make_pair(j, tau));
            neighbours[j].push_back(std::make_pair(i, tau));
        }
    }
    
    for (int i =0; i< size; i++) {
        sort(neighbours[i].begin(), neighbours[i].end(),
             [](const auto& a, const auto& b) { return a.second < b.second; });
    }
}

std::pair<node_t, interval_t>
erdos_reyni::neighbour(node_t node, int neighbour_index) {
    const auto& n = neighbours.at(node);
    if ((neighbour_index < 0) || (n.size() <= neighbour_index))
        return std::make_pair(-1, INFINITY);
    return n[neighbour_index];
}

