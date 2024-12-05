#include "stdafx.h"
#include "NextReaction.h"
#include "random.h"
#include "network.h"
#include "types.h"
#include "utility.h"

absolutetime_t simulate_next_reaction::next(rng_t& engine)
{
	if (active_edges.empty())
		return INFINITY;
	
	/* Fetch the next infection/reset time, i.e. the time where the next edge fires */
	const auto next = top_edge();
	return next.time;
}

std::optional<epidemic_event_t>
simulate_next_reaction::step(rng_t& engine, absolutetime_t maxtime, event_filter_t evf)
{
    if (std::isnan(maxtime))
        throw std::range_error("maxtime must be finite or +INFINITY");

    while (true) {
        /* If there are no more infection times, stop */
        if (active_edges.empty())
			return std::nullopt;

        /* Fetch the next infection/reset time, i.e. the time where the next edge fires.
		 * Return before handling the event if it's time is larger than nexttime
		 */
        const auto next = top_edge();
        if (next.time > maxtime)
			return std::nullopt;
		
		/* Dequeue event */
        pop_edge();
        ++queue_steps_total;
        
        /* Perform event */
        std::optional<epidemic_event_t> result;
        switch (next.kind) {
            case epidemic_event_kind::infection:
            case epidemic_event_kind::outside_infection:
				result = step_infection(next, evf, engine);
				break;
			case epidemic_event_kind::reset:
				result = step_reset(next, evf, engine);
				break;
            default: throw std::logic_error("invalid event kind");
        }

        /* Return event (unless skipped, in which case we continue) */
        if (result)
            return result;        
    }
}

void simulate_next_reaction::notify_infected_node_neighbour_added(network_event_t event, rng_t& engine)
{
	/* A neighbour was added to an already infected node. We have to add an active edge
	 * that corresponds to the new neighbour. The transmission time is generated to be
	 * larger than the time since infection of the node.
	 */

	/* Query state of source node */
	const auto source_state = infected.find(event.source_node);
	if (source_state == infected.end())
		throw std::logic_error("failed to query state of infected node");

	/* Generate firing time of the newly added edge */
	assert(event.time >= source_state->second.infection_time);
	const double tau = psi.sample(engine, event.time - source_state->second.infection_time, 1);
	if (std::isnan(tau) || (tau < 0))
		throw std::logic_error("transmission time must be positive");
	const double t = event.time + tau;

	/* If the edge fires after the node has resetted (or never) it has no effect */
	if (t >= source_state->second.reset_time)
		return;

	/* Add edge */
	assert(std::isfinite(t));
	active_edges_entry e;
	e.kind = epidemic_event_kind::infection;
	e.time = t;
	e.node = event.target_node;
	e.source_time = source_state->second.infection_time;
	e.source_node = event.source_node;
	e.source_reset = source_state->second.reset_time;
	e.neighbour_index = 0;
	e.neighbours_remaining = 0;
	push_edge(e);
}

void simulate_next_reaction::notify_infected_contact(network_event_t event, rng_t& engine)
{
	/* An instantenous contact occured, check if it leads to a transmission.
	 * Event must occur now, not in the future
	 */
	assert(active_edges.empty() || (top_edge().time >= event.time));

	/* Query state of source node */
	const auto source_state = infected.find(event.source_node);
	if (source_state == infected.end())
		throw std::logic_error("notify_infected_contact called for non-infected source node");

	/* Computing tranmission probability */
	assert(event.time >= source_state->second.infection_time);
	const double p = psi.hazardrate(event.time - source_state->second.infection_time) * event.infitesimal_duration;
	if (!std::bernoulli_distribution(p)(engine))
		return;

	/* Queue infection, this will be the next event that occurs (see assert above) */
	active_edges_entry e;
	e.kind = epidemic_event_kind::infection;
	e.time = event.time;
	e.node = event.target_node;
	e.source_time = source_state->second.infection_time;
	e.source_node = event.source_node;
	e.source_reset = source_state->second.reset_time;
	e.neighbour_index = -1;
	e.neighbours_remaining = 0;
	push_edge(e);
}

std::optional<epidemic_event_t> simulate_next_reaction::step_infection(const active_edges_entry& next, event_filter_t evf, rng_t& engine)
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
		assert(!p.edges_concurrent);
		
        /* This only occurs for infection edges, reset self-loops have no neighbours */
        const node_t neighbour_id = next.source_permutation[next.neighbour_index + 1];
        const node_t sibling = nw.neighbour(next.source_node, neighbour_id);
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
            e.kind = epidemic_event_kind::infection;
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
	const epidemic_event_t ev { .kind = next.kind, .source_node = next.source_node, .node = next.node, .time = next.time };

    /* Check if event is blocked or putatively infected node is already infected, if so we're done */
    if (is_event_blocked(ev, evf) || is_infected(next.node))
        return std::nullopt;
    
    /* Node becomes infected.
     * Mark the node as infected, generate reset time, insert recovery self-loop
     */
    const interval_t tau_reset = (rho ? rho->sample(engine, 0, 1) : INFINITY);
    const absolutetime_t node_reset_time = next.time + tau_reset;
    if (tau_reset < INFINITY) {
        active_edges_entry e;
        e.kind = epidemic_event_kind::reset;
        e.time = node_reset_time;
        e.node = next.node;
		push_edge(e);
    }
    infected.emplace(next.node, infected_state_t(next.time, node_reset_time));

    /* Activate the first outgoing edge of the newly infected node
     * Further edges are activated iteratively after the previously one fired,
     * see the code block above labelled "Activate the next sibling edges of the
     * edge that just fired".
     */

    /* First, create a permutation to shuffle the neighbours if necessary */
    const int neighbours_total = nw.outdegree(next.node);
    permutation<node_t> pi;
    if (shuffle_neighbours) {
		/* We should never shuffle neighbours if sibling edges are made active concurrently */
		assert(!p.edges_concurrent);

		if (neighbours_total < 0)
            throw std::runtime_error("cannot shuffle neighbours if nodes have undefined or infinite degree");
        pi = permutation<node_t>(neighbours_total, engine);
    }

	/* We now add either only the first (in our permutation) edge to the list
	 * active edges (if !edges_concurrent), or we add all outgoing edges at the same
	 * time (if edges_concurrent). For simplicity, we traverse the neighbours in
	 * the order defined by the permutation in both cases, but note that if we add
	 * edges concurrently, shuffle_neighbours is always false and so the permutation
	 * is actually the identity map.
	 */
	int r;
	if (p.edges_concurrent)
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
		const node_t neighbour = nw.neighbour(next.node, pi[neighbour_i]);
		
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
		const int m = p.edges_concurrent ? 1 : neighbours_total;
		const double tau = psi.sample(engine, 0, m);
        if (std::isnan(tau) || (tau < 0))
            throw std::logic_error("transmission times must be non-negative");
        const double t = next.time + tau;
        /* Only queue if the infection occurs before the infecting node's reset */
        if (t < node_reset_time) {
            assert(std::isfinite(t));
            active_edges_entry e;
            e.kind = epidemic_event_kind::infection;
            e.time = t;
            e.node = neighbour;
            e.source_time = next.time;
            e.source_node = next.node;
            e.source_reset = node_reset_time;
            if (!p.edges_concurrent)
                e.source_permutation = std::move(pi);
            e.neighbour_index = 0;
			e.neighbours_remaining = m - 1;
			push_edge(e);
        }
    }

    /* Return the infection event */
	return ev;
}

std::optional<epidemic_event_t> simulate_next_reaction::step_reset(const active_edges_entry& next, event_filter_t evf, rng_t& engine) {
    /* The node cannot yet have resetted */
    assert(is_infected(next.node));

	/* Create event and query filter */
	const epidemic_event_t ev { .kind = epidemic_event_kind::reset, .source_node = next.source_node, .node = next.node, .time = next.time };
	if (is_event_blocked(ev, evf))
		return std::nullopt;
	
	/* In SIR mode, reset events do not make nodes susceptible again, but only terminate the infections
	 * phase early. We implement that via the following hack that leaves the node marked as "infected"
	 * (of which a more appropriate name is this mode would be recovered).
	 * NOTE: This hack should eventually be removed, and be replaced by a separate class that uses
	 * event filters to implement SIR mode. Through an appropriate filter, nodes would be added to a
	 * "removed" set upon receiving a reset event, and that set would be queried before allowing
	 * infections to proceed.
	 */
	if (p.SIR) {
		/* SIR mode, just could the number of removed nodes */
		removed += 1;
	} else {
		/* SIS mode, mark node as not infected */
		infected.erase(ev.node);
	}

    /* Return the reset event */
	return ev;
}

void simulate_next_reaction::add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v) {
    for(const auto& ve: v) {
        active_edges_entry e;
        e.kind = epidemic_event_kind::outside_infection;
        e.time = ve.second;
        e.node = ve.first;
		push_edge(e);
    }
}


void simulate_next_reaction::add_thermal_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v,double Lambda, rng_t& engine) {
    std::uniform_real_distribution<> ran(0.0,1.0);
    for(const auto& ve: v) {

        node_t node = ve.first;
        time_t time = ve.second;

        // mark the node as infected
        infected.emplace(node, infected_state_t(time, INFINITY));
        assert(is_infected(node));
        //manually add the infections if infection happpens after their age
        for (int i =0; i < nw.outdegree(node) ; i++ ){
            node_t neighbour = nw.neighbour(node,i);
            const double tau = psi.sample(engine, 0, 1);
            const double thermal =-log(1-ran(engine))/Lambda;

            if (tau<thermal){
                continue;
            }
            active_edges_entry e;
            e.kind = epidemic_event_kind::infection;
            e.time = time+tau;
            e.node = neighbour;
            e.source_time = time;
            e.source_node = node;
            push_edge(e);
        }
    }
}

bool simulate_next_reaction::is_infected(node_t node) const {
    return (infected.find(node) != infected.end());
}

network& simulate_next_reaction::get_network() const {
    return nw;
}

const transmission_time& simulate_next_reaction::transmission_time() const {
    return psi;
}

const transmission_time* simulate_next_reaction::reset_time() const {
    return rho;
}
