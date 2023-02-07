#include "stdafx.h"
#include "NextReaction.h"
#include "random.h"
#include "graph.h"
#include "types.h"
#include "utility.h"

absolutetime_t simulate_next_reaction::next()
{
	if (active_edges.empty())
		return INFINITY;
	
	/* Fetch the next infection/reset time, i.e. the time where the next edge fires */
	const auto next = top_edge();
	return next.time;
}

std::optional<event_t>
simulate_next_reaction::step(rng_t& engine, absolutetime_t nexttime, event_filter_t evf)
{
    while (true) {
        /* If there are no more infection times, stop */
        if (active_edges.empty())
			return std::nullopt;

        /* Fetch the next infection/reset time, i.e. the time where the next edge fires.
		 * Return before handling the event if it's time is larger than nexttime
		 */
        const auto next = top_edge();
		if (next.time > nexttime)
			return std::nullopt;
		
		/* Dequeue event */
        pop_edge();
        ++queue_steps_total;
        
        /* Perform event */
        std::optional<event_t> result;
        switch (next.kind) {
            case event_kind::infection:
            case event_kind::outside_infection:
				result = step_infection(next, evf, engine);
				break;
			case event_kind::reset:
				result = step_reset(next, evf, engine);
				break;
            default: throw std::logic_error("invalid event kind");
        }

        /* Return event (unless skipped, in which case we continue) */
        if (result)
            return result;        
    }
}

void simulate_next_reaction::notify_infected_node_neighbour_added(network_event_t event)
{
	throw std::logic_error("unimplemented");
}

std::optional<event_t> simulate_next_reaction::step_infection(const active_edges_entry& next, event_filter_t evf, rng_t& engine)
{
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
		/* Should never happen if we're making sibling edges active concurrently */
		assert(!edges_concurrent);
		
        /* This only occurs for infection edges, reset self-loops have no neighbours */
        const node_t neighbour_id = next.source_permutation[next.neighbour_index + 1];
        const node_t sibling = network.neighbour(next.source_node, neighbour_id);
        if (sibling < 0) {
            /* This should never happen unless the graph reported the wrong number of outgoing edges */
            throw std::logic_error(std::string("neighbour ") + std::to_string(neighbour_id) +
                                   " (index " + std::to_string(next.neighbour_index + 1) + ") " +
                                   std::to_string(next.source_node) + " is invalid");
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
        /* Only queue if the infection occurs before the infecting node's reset */
        if (t < next.source_reset) {
            assert(std::isfinite(t));
            active_edges_entry e;
            e.kind = event_kind::infection;
            e.time = t;
            e.node = sibling;
            e.source_time = next.source_time;
            e.source_node = next.source_node;
            e.source_reset = next.source_reset;
            e.source_permutation = std::move(next.source_permutation);
            e.neighbour_index = next.neighbour_index + 1;
            e.neighbours_remaining = next.neighbours_remaining - 1;
			push_edge(e);
        }
    }

	/* Create event */
	const event_t ev { .kind = next.kind, .node = next.node, .source_node = next.source_node, .time = next.time };

    /* Check if event is blocked or putatively infected node is already infected, if so we're done */
    if (is_event_blocked(ev, evf) || is_infected(next.node))
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
		push_edge(e);
    }
    
    /* Activate the first outgoing edge of the newly infected node
     * Further edges are activated iteratively after the previously one fired,
     * see the code block above labelled "Activate the next sibling edges of the
     * edge that just fired".
     */

    /* First, create a permutation to shuffle the neighbours if necessary */
    const int neighbours_total = network.outdegree(next.node);
    permutation<node_t> p;
    if (rho && shuffle_neighbours) {
		/* We should never shuffle neighbours if sibling edges are made active concurrently */
		assert(!edges_concurrent);

		if (neighbours_total < 0)
            throw std::runtime_error("cannot shuffle neighbours if nodes have undefined or infinite degree");
        p = permutation<node_t>(neighbours_total, engine);
    }

	/* We now add either only the first (in our permutation) edge to the list
	 * active edges (if !edges_concurrent), or we add all outgoing edges at the same
	 * time (if edges_concurrent). For simplicity, we traverse the neighbours in
	 * the order defined by the permutation in both cases, but note that if we add
	 * edges concurrently, shuffle_neighbours is always false and so the permutation
	 * is actually the identity map.
	 */
	int r;
	if (edges_concurrent)
		/* Interate over all neighbours */
		r = neighbours_total;
	else if (neighbours_total >= 1)
		/* Only handle the first neighbour here */
		r = 1;
	else
		/* No neighbours */
		r = 0;
	for(int neighbour_i = 0; neighbour_i < r; ++neighbour_i) {
		/* Get i-th neighbour according to the permutation */
		const node_t neighbour = network.neighbour(next.node, p[neighbour_i]);
		
		/* This should never happen unless the graph reported the wrong number
		 * of outgoing edges */
		if (neighbour < 0) {
			throw std::logic_error(std::string("neighbour ") + std::to_string(neighbour) +
								   " (index " + std::to_string(neighbour_i + 1) +
								   ") " + std::to_string(next.node) + " is invalid");
		}

		/* Create target node's infection times entry and add to queue.
		 * If we're adding edges sequentially, the distribution we sample
		 * from is the *minimum* time taken over the remaining neighbours.
		 * If we're adding edges concurrently, however, we simply sample
		 * from the original psi distribution. The minimum in this case is
		 * the implicit result of putting all the times into the reaction queue,
		 * and processing it in order.
		 */
		const int m = edges_concurrent ? 1 : neighbours_total;
		const double tau = psi.sample(engine, 0, m);
        if (std::isnan(tau) || (tau < 0))
            throw std::logic_error("transmission times must be non-negative");
        const double t = next.time + tau;
        /* Only queue if the infection occurs before the infecting node's reset */
        if (t < node_reset_time) {
            assert(std::isfinite(t));
            active_edges_entry e;
            e.kind = event_kind::infection;
            e.time = t;
            e.node = neighbour;
            e.source_time = next.time;
            e.source_node = next.node;
            e.source_reset = node_reset_time;
            e.source_permutation = std::move(p);
            e.neighbour_index = 0;
			e.neighbours_remaining = m - 1;
			push_edge(e);
        }
    }

    /* Return the infection event */
	return ev;
}

std::optional<event_t> simulate_next_reaction::step_reset(const active_edges_entry& next, event_filter_t evf, rng_t& engine) {
    /* The node cannot yet have resetted */
    assert(is_infected(next.node));

	/* Create event and query filter */
	const event_t ev { .kind = event_kind::reset, .node = next.node, .source_node = next.source_node, .time = next.time };
	if (is_event_blocked(ev, evf))
		return std::nullopt;
	
    /* if SIR, increase counter of removed/recovered nodes, and leave the node as infected so that 
    * the node cannot be reinfected. 
    * ( this is just a trick to avoid having to create a new state, the node is not actually infected anymore and does not generate new infections.)
    * else, Mark node as no longer infected. In that case the node simply returns to the susceptible state and can get reinfected again later. */
    if (SIR){
        removed += 1;
    } else {
        infected.erase(next.node);
    }

    /* Return the reset event */
	return ev;
}

void simulate_next_reaction::add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v) {
    for(const auto& ve: v) {
        active_edges_entry e;
        e.kind = event_kind::outside_infection;
        e.time = ve.second;
        e.node = ve.first;
		push_edge(e);
    }
}

bool simulate_next_reaction::is_infected(node_t node) const {
    return (infected.find(node) != infected.end());
}

graph& simulate_next_reaction::get_network() const {
    return network;
}

const transmission_time& simulate_next_reaction::transmission_time() const {
    return psi;
}

const transmission_time* simulate_next_reaction::reset_time() const {
    return rho;
}


