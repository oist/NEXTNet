#include "stdafx.h"
#include "NextReaction.h"
#include "random.h"
#include "graph.h"
#include "types.h"
#include "utility.h"


std::pair<node_t, absolutetime_t> simulate_next_reaction::step(rng_t& engine) {
    while (true) {
        /* If there are no more infection times, stop */
        if (active_edges.empty())
            return std::make_pair(-1, INFINITY);

        /* Fetch the next putatively infected node */
        const auto next = active_edges.top();
        active_edges.pop();
        
        /* Lazily enqueue next sibling of the putatively infected node if necessary */
        if (next.neighbours_remaining > 0) {
            const node_t sibling = network.neighbour(next.source_node, next.neighbour_index + 1);
            if (sibling < 0)
                throw std::logic_error(std::string("neighbour ") + std::to_string(next.neighbour_index + 1) +
                                       " of node " + std::to_string(next.source_node) + " is invalid");
            /* Create sibling's infection times entry and add to queue */
            const double tau = psi.sample(engine, next.time - next.source_time,
                                          next.neighbours_remaining);
            active_edges_entry e;
            e.time = next.source_time + tau;
            e.node = sibling;
            e.source_time = next.source_time;
            e.source_node = next.source_node;
            e.neighbour_index = next.neighbour_index + 1;
            e.neighbours_remaining = next.neighbours_remaining - 1;
            active_edges.push(e);
        }
        
        /* Check if the putatively infected node is already infected, if so we're done */
        if (infected.find(next.node) != infected.end())
            continue;
        
        /* Mark the node as infected */
        infected.insert(next.node);
        
        /* Add the infecte node's first neigbhour to the infection times queue */
        const node_t neighbour = network.neighbour(next.node, 0);
        if (neighbour >= 0) {
            /* Create node's infection times entry and add to queue */
            const int neighbours_remaining = network.outdegree(neighbour);
            const double tau = psi.sample(engine, 0, neighbours_remaining);
            active_edges_entry e;
            e.time = next.time + tau;
            e.node = neighbour;
            e.source_time = next.time;
            e.source_node = next.node;
            e.neighbour_index = 0;
            e.neighbours_remaining = neighbours_remaining;
            active_edges.push(e);
        }
        
        return std::make_pair(next.node, next.time);
    }
}

void simulate_next_reaction::add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v) {
    for(const auto& ve: v) {
        active_edges_entry e;
        e.time = ve.second;
        e.node = ve.first;
        active_edges.push(e);
    }
}
