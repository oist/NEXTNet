#include "stdafx.h"
#include "NextReactionMeanField.h"
#include "random.h"
#include "types.h"
#include "utility.h"


std::optional<event_t> simulate_next_reaction_mean_field::step(rng_t& engine) {
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
            case event_kind::outside_infection:
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

std::optional<event_t> simulate_next_reaction_mean_field::step_infection(const active_edges_entry& next, rng_t& engine) {
    /*
     * Here, the variables have the following meaning:
     *   next.time: The current time, i.e. the time at which the edge fired
     *   next.node: The node that just got infected (putatively, since it might already be)
     */
    
    /* Check if the putatively infected node is already infected, if so we're done */
    if (is_infected(next.node))
        return std::nullopt;

    /* Node becomes infected.
     * Mark the node as infected, generate reset time, insert recovery
     */
    infected.insert(next.node);
    const interval_t tau_reset = (rho ? rho->sample(engine, 0, 1) : INFINITY);
    if (tau_reset < INFINITY) {
        active_edges_entry e;
        e.kind = event_kind::reset;
        e.time = next.time + tau_reset;
        e.node = next.node;
        active_edges.push(e);
    }

    /* Instead of generating k=Poi(R0) and then sampling without replacement k nodes,
    we can sample them using the geomertic distribution.*/
    std::geometric_distribution<> skip_edge(p); // comment: equals 0 with prob. p

    for (node_t node=0; node < N; node++) {
        int s = skip_edge(engine);
        if (s + node >= N)
            break;
        node += (node == next.node ? 1 + s : s);
        const interval_t tau_inf = psi.sample(engine, 0, 1); // sample simple r.v from psi distribution with age 0
        if (tau_inf > tau_reset) // transmission will not happen. (Note: in our simulations we often choose rho s.t Prob[tau_inf > tau_reset] ~ 1).
            continue;
        active_edges_entry e;
        e.kind = event_kind::infection;
        e.time = next.time + tau_inf;
        e.node = node;
        active_edges.push(e);
    }

    /* Return the infection event */
    return event_t { .kind = next.kind, .node = next.node, .time = next.time };
}

std::optional<event_t> simulate_next_reaction_mean_field::step_reset(const active_edges_entry& next, rng_t& engine) {
    /* The node cannot yet have resetted */
    assert(is_infected(next.node));
    
    /* if SIR, increase counter of removed/recovered nodes, and leave the node as infected so that 
    * the node cannot be reinfected. (just a trick to avoid creating a new state, the node is not actually infected anymore.)
    * else, Mark node as no longer infected. In that case the node simply returns to the susceptible state and can get reinfected again later. */
    if (SIR){
        removed += 1;
    } else {
        infected.erase(next.node);
    }

    /* Return the reset event */
    return event_t { .kind = event_kind::reset, .node = next.node, .time = next.time };
}

void simulate_next_reaction_mean_field::add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v) {
    for(const auto& ve: v) {
        active_edges_entry e;
        e.kind = event_kind::outside_infection;
        e.time = ve.second;
        e.node = ve.first;
        active_edges.push(e);
    }
}

bool simulate_next_reaction_mean_field::is_infected(node_t node) const {
    return (infected.find(node) != infected.end());
}

 graph& simulate_next_reaction_mean_field::get_network() const {
     throw std::logic_error("get_network() is unsupported for mean field simulations");
 }

const transmission_time& simulate_next_reaction_mean_field::transmission_time() const {
    return psi;
}

const transmission_time* simulate_next_reaction_mean_field::reset_time() const {
    return rho;
}


