#include "stdafx.h"
#include "NextReaction.h"
#include "random.h"
#include "graph.h"
#include "types.h"
#include "utility.h"

std::optional<event_t> simulate_next_reaction::step(rng_t& engine) {
    while (true) {
        /* If there are no more infection times, stop */
        if (active_edges.empty())
			return std::nullopt;

        /* Fetch the next infection/reset time, i.e. the time where the next edge fires */
        const auto next = active_edges.top();
        active_edges.pop();
        
        /* Perform event */
        std::optional<event_t> result;
        switch (next.kind) {
			case event_kind::infection:
				result = step_infection(next, engine);
				break;
			case event_kind::reset:
				result = step_reset(next, engine);
				break;
            default: throw std::logic_error("invalid event kind");
        }

        /* Return event (unless skipped, in which case we continue) */
        if (result)
            return result;        
    }
}

std::optional<event_t> simulate_next_reaction::step_infection(const active_edges_entry& next, rng_t& engine) {
    /*
     * Here, the variables have the following meaning:
     *   next.time: The current time, i.e. the time at which the edge fired
     *   next.node: The node that just got infected (putatively, since it might already be)
     *   next.source_node: The node that caused the infection
     */

    /* The infection's source node cannot yet have resetted */
    assert(next.source_reset > next.time);

    /*
     * Activate the next sibling edges of the edge that just fired.
     * This is necessary because we don't activate all outgoing edges of a node when it
     * becomes infected, but rather just the first one. When that first edge fires, we
     * activate the next outgoing edge with the same originating node, called a *sibling*
     * edge.
     */
    if (next.neighbours_remaining > 0) {
        /* This only occurs for infection edges, reset self-loops have no neighbours */
        const node_t sibling = network.neighbour(next.source_node, next.neighbour_index + 1);
        if (sibling < 0) {
            /* This should never happen unless the graph reported the wrong number of outgoing edges */
            throw std::logic_error(std::string("neighbour ") + std::to_string(next.neighbour_index + 1) +
                                    " of node " + std::to_string(next.source_node) + " is invalid");
        }
        /* Create sibling's infection times entry and add to queue
         * When sampling transmission times, we pass the "current time" in a frame of reference where
         * the infection occured at time 0. The returned transmission time, however, is relative to the
         * current time.
         */
        const double tau = psi.sample(engine, next.time - next.source_time,
                                      next.neighbours_remaining);
        if (std::isnan(tau) || (tau < 0))
            throw std::logic_error("transmission times must be non-negative");
        const double t = next.time + tau;
        /* Only queue if the infection occurs before the infeceting node's reset */
        if (t < next.source_reset) {
            assert(std::isfinite(t));
            active_edges_entry e;
            e.time = t;
            e.node = sibling;
            e.source_time = next.source_time;
            e.source_node = next.source_node;
            e.source_reset = next.source_reset;
            e.neighbour_index = next.neighbour_index + 1;
            e.neighbours_remaining = next.neighbours_remaining - 1;
            active_edges.push(e);
        }
    }
    
    /* Check if the putatively infected node is already infected, if so we're done */
    if (is_infected(next.node))
        return std::nullopt;
    
    /* Node becomes infected.
     * Mark the node as infected, generate reset time, insert recovery self-loop
     */
    infected.insert(next.node);
    const interval_t tau_reset = (rho ? rho->sample(engine, 0, 1) : INFINITY);
    const absolutetime_t node_reset_time = next.time + tau_reset;
    if (tau_reset < INFINITY) {
        active_edges_entry e;
        e.kind = event_kind::reset;
        e.time = node_reset_time;
        e.node = next.node;
        active_edges.push(e);
    }
    
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
        if (std::isnan(tau) || (tau < 0))
            throw std::logic_error("transmission times must be non-negative");
        const double t = next.time + tau;
        /* Only queue if the infection occurs before the infeceting node's reset */
        if (t < node_reset_time) {
            assert(std::isfinite(t));
            active_edges_entry e;
            e.time = t;
            e.node = neighbour;
            e.source_time = next.time;
            e.source_node = next.node;
            e.source_reset = node_reset_time;
            e.neighbour_index = 0;
            e.neighbours_remaining = neighbours_total - 1;
            active_edges.push(e);
        }
    }
    
    /* Return the infection event */
    return event_t { .kind = event_kind::infection, .node = next.node, .time = next.time };
}

std::optional<event_t> simulate_next_reaction::step_reset(const active_edges_entry& next, rng_t& engine) {
    /* The node cannot yet have resetted */
    assert(is_infected(next.node));
    
    /* Mark node as no longer infected */
    infected.erase(next.node);

    /* Return the reset event */
    return event_t { .kind = event_kind::reset, .node = next.node, .time = next.time };
}

void simulate_next_reaction::add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v) {
    for(const auto& ve: v) {
        active_edges_entry e;
        e.time = ve.second;
        e.node = ve.first;
        active_edges.push(e);
    }
}

bool simulate_next_reaction::is_infected(node_t node) {
    return (infected.find(node) != infected.end());
}

graph& simulate_next_reaction::get_network() {
    return network;
}

transmission_time& simulate_next_reaction::transmission_time() {
    return psi;
}

transmission_time* simulate_next_reaction::reset_time() {
    return rho;
}
