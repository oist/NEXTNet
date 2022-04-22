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

        /* Fetch the next infection time, i.e. the time where the next edge fires
         * From here on, the variables have the following meaing:
         *   next.time: The current time, i.e. the time at which the edge fired
         *   next.node: The node that just got infected (putatively, since it might already be)
         *   next.source_node: The node that caused the infection
         */
        const auto next = active_edges.top();
        active_edges.pop();
        
        /* Activate the next sibling edges of the edge that just fired.
         * This is necessary because we don't activate all outgoing edges of a node when it
         * becomes infected, but rather just the first one. When that first edge fires, we
         * activate the next outgoing edge with the same originating node, called a *sibling*
         * edge.
         */
        if (next.neighbours_remaining > 0) {
            const node_t sibling = network.neighbour(next.source_node, next.neighbour_index + 1);
            if (sibling < 0) {
                /* This should never happen unless the graph reported the wrong number of outgoing edges */
                throw std::logic_error(std::string("neighbour ") + std::to_string(next.neighbour_index + 1) +
                                       " of node " + std::to_string(next.source_node) + " is invalid");
            }
            /* Create sibling's infection times entry and add to queue
             * When sampling transmission times, we always work in a frame of reference where
             * the infection occured at time 0, therefore we translate the current time accordingly
             * when sampling tau here.
             */
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
        
        /* Activate the first outgoing edge of the newly infected node
         * Further edges are activated iteratively after the previously one fired,
         * see the code block above labelled "Activate the next sibling edges of the
         * edge that just fired".
         */
        const node_t neighbour = network.neighbour(next.node, 0);
        if (neighbour >= 0) {
            /* Create target node's infection times entry and add to queue */
            const int neighbours_total = network.outdegree(next.node);
            const double tau = psi.sample(engine, 0, neighbours_total);
            active_edges_entry e;
            e.time = next.time + tau;
            e.node = neighbour;
            e.source_time = next.time;
            e.source_node = next.node;
            e.neighbour_index = 0;
            e.neighbours_remaining = neighbours_total - 1;
            active_edges.push(e);
        }
        
        /* Return the infected node and it's infection time */
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
