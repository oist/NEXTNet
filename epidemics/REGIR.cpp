#include "stdafx.h"
#include "REGIR.h"
#include "random.h"
#include "graph.h"
#include "types.h"
#include "utility.h"

graph& simulate_regir::get_network() const {
	return network;
}

const class transmission_time& simulate_regir::transmission_time() const {
	return psi;
}

const class transmission_time* simulate_regir::reset_time() const {
	return rho;
}

bool simulate_regir::is_infected(node_t node) const {
	return (infected.find(node) != infected.end());
}

void simulate_regir::add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v)
{
	for(const auto& ve: v) {
		/* Add infection to sorted list of outside infections */
		outside_infections_entry inf;
		inf.node = ve.first;
		inf.time = ve.second;
		outside_infections.push(inf);
	}
}

absolutetime_t simulate_regir::next(rng_t& engine)
{
    /* If the next event was already determined, just return it */
    if (next_event)
        return next_event->time;

    /* If the current time wasn't yet set, start at the earliest time at which
     * and edge becomes active. */
    absolutetime_t base_time = current_time;
    if (std::isnan(base_time)) {
        base_time = outside_infections.empty() ? INFINITY : outside_infections.top().time;
        for(const active_edges_entry& e: active_edges)
            base_time = std::min(base_time, e.source_time);
    }

    /* If there are no active edges and no outside infections, there is no next event */
    if (active_edges.empty() && outside_infections.empty())
        base_time = INFINITY;
    if (std::isinf(base_time))
        return INFINITY;

	/* Use exact algorithm if number of active edgesd is small */
    const bool use_exact_algorithm = ((p.approximation_threshold < 0) ||
                                      (active_edges.size() <= (unsigned int)p.approximation_threshold));

    /* Generate next event */
    double tau = NAN;
	absolutetime_t next_time;
	active_edges_t::iterator edge_i;
	
    if (active_edges.empty()) {
        /* No active edges, but there are outside infections (see below) */

        tau = INFINITY;
		
		/* Check for outside infections. next_time_outside_infection fills next_event */
		return next_time_outside_infection(INFINITY);
    } else if (use_exact_algorithm) {
        /* Exact version */

        /* Note: The exact version does not use the harard rates lambda, we
         * thus do not have to update them before drawing tau
         */
        tau = next_time_exact(engine);
		next_time = base_time + tau;
		
		/* Check for outside infections. next_time_outside_infection fills next_event */
		const absolutetime_t tau_outside = next_time_outside_infection(base_time + tau);
		if (!std::isnan(tau_outside))
			return tau_outside;
		
		/* Draw random edge according to current hazard rates */
		update_active_edge_lambdas(next_time);
		edge_i = draw_active_edge_hazardrates(engine);
	} else {
		/* REGIR version */
		
		do {
			/* Generate candidate time */
			tau = next_time_exponential(engine);
			
			/* Draw uniformly from all active edges */
			edge_i = draw_active_edge_uniform(engine);
			
			/* Rejection step
			 * Reject with probability lambda / lambda_max, where
			 * lambda is the hazard rate of the selected edge
			 */
			const double q = unif01_dist(engine);
			const double lambda = edge_hazard_rate(*edge_i, base_time + tau);
			if (lambda > lambda_max)
				throw std::logic_error("hazard rate exceeds presumed maximum");
			if (q > lambda/lambda_max)
				continue;
		} while (false);
		next_time = base_time + tau;
			
		/* Check for outside infections. next_time_outside_infection fills next_event */
		const absolutetime_t tau_outside = next_time_outside_infection(base_time + tau);
		if (!std::isnan(tau_outside))
			return tau_outside;
	}

    /* Found next event, store & return its time */
	next_event_edge_pos = edge_i - active_edges.begin();
    assert(edge_i->kind != event_kind::outside_infection);
    next_event = event_t {
            .kind = edge_i->kind,
            .source_node = edge_i->source, .node = edge_i->target,
            .time = next_time
    };
    return next_event->time;
}

std::optional<event_t> simulate_regir::step(rng_t& engine, absolutetime_t maxtime, event_filter_t evf)
{
    if (std::isnan(maxtime))
        throw std::range_error("maxtime must be finite or +INFINITY");

    while (true) {
        /* Determine next event */
        if (!next_event)
            next(engine);

        /* If there is no event or we'd move past max time */
        if (!next_event || (next_event->time > maxtime))
            return std::nullopt;

		/* Event will be handled or skipped, so update current_time and reset next_event */
		const event_t ev = *next_event;
		current_time = next_event->time;
		next_event = std::nullopt;

		/* Remove corresponding active edge or outside infection */
		switch (ev.kind) {
			case event_kind::outside_infection: {
				if (outside_infections.empty() || (outside_infections.top().time != ev.time) || (outside_infections.top().node != ev.node))
					throw std::logic_error("next outside infection changed between next() and step()");
				outside_infections.pop();
				break;
			}

			case event_kind::infection:
			case event_kind::reset: {
				/* Re-find active edge.
				 * NOTE: It would be more efficient to pass the iterator from next() to here. However,
				 * that would make it impossible to check whether the edge is actually still part of the
				 * active edges set, which is a reasonable safe-guard against coding mistakes
				 */
				active_edges_t::iterator edge_i = active_edges.begin() + *next_event_edge_pos;
				remove_active_edge(edge_i);
				next_event_edge_pos = std::nullopt;
				break;
			}

			default:
				throw std::logic_error("unknown event kind");
		}

		/* Query filter and skip event if indicated */
		if (is_event_blocked(ev, evf))
			continue;

		/* Handle event */
		switch (ev.kind) {
			case event_kind::infection:
			case event_kind::outside_infection: {
				/* Infection event */

				/* If the node is already infected, ignore the event */
				if (infected.find(ev.node) != infected.end())
					continue;

				/* Mark node as infected */
				infected.insert(ev.node);

				/* Make recovery self-loop if there's a reset time distribution ... */
				if (rho) {
					active_edges_entry e;
					e.kind = event_kind::reset;
					e.source = ev.node;
					e.source_time = current_time;
					e.target = ev.node;
					add_active_edge(e);
				}

				/* ... and outgoing edges active */
				const int neighbours = network.outdegree(ev.node);
				for(int j=0; j < neighbours; ++j) {
					const node_t neighbour = network.neighbour(ev.node, j);
					if (neighbour < 0) {
						/* This should never happen unless the graph reported the wrong number of outgoing edges */
						throw std::logic_error(std::string("neighbour ") + std::to_string(j + 1) +
												" of node " + std::to_string(ev.node) + " is invalid");
					}
					active_edges_entry e;
					e.kind = event_kind::infection;
					e.source = ev.node;
					e.source_time = current_time;
					e.target = neighbour;
					add_active_edge(e);
				}
				break;
			}

			case event_kind::reset: {
				/* Reset event */

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

				/* Remove active edges originating from the resetted node */
				for(auto i = active_edges.begin(); i != active_edges.end();) {
					/* Skip edges originating anywhere else */
					if (i->source != ev.node) {
						++i ;
						continue;
					}

					/* Remove edge, interator points to next element afterwards */
					remove_active_edge(i);
				}
				break;
			}

			default:
				throw std::logic_error("unknown event kind");
		}

        /* Return handled event */
        return ev;
    }
}

absolutetime_t simulate_regir::next_time_outside_infection(absolutetime_t max_time) {
	/* Check if we would jump over any outside infections. If there is an
	 * outside infection, we discard the tau we just generated but that's OK.
	 */
	while (!outside_infections.empty()) {
		/* Check if the outside infection occurs before the generated tau */
		const outside_infections_entry inf = outside_infections.top();
		if (inf.time > max_time)
			return NAN;

		/* Skip and remove if the node is already infected */
		if (infected.find(inf.node) != infected.end()) {
			outside_infections.pop();
			continue;
		}

		/* Found next event, store & return its time */
		next_event = event_t {
				.kind = event_kind::outside_infection,
				.source_node = -1, .node = inf.node,
				.time = inf.time
		};
		return next_event->time;
	}
	
	return NAN;
}

void simulate_regir::notify_infected_node_neighbour_added(network_event_t event, rng_t& engine)
{
    throw std::logic_error("REGIR currently does not support dynamic networks");
}

interval_t simulate_regir::next_time_exact(rng_t& engine) {
    /* Determine time of next event by inverting the global survival function phi */
    const double u = unif01_dist(engine);
    const double tau = invphi(current_time, u);
    return tau;
}

interval_t simulate_regir::next_time_exponential(rng_t& engine)
{
    /* Determine time of next event */
    return std::exponential_distribution<double>(lambda_max * active_edges.size())(engine);
}

double simulate_regir::phi(absolutetime_t t, interval_t tau) {
	/* This implements our version of equation 4 from the Boguna paper */
	if (std::isinf(tau))
		return 0;
	double r = 1;
	for(const active_edges_entry& e: active_edges) {
		/* Translate t into the edge's frame of reference, i.e. into the
		 * time since the edge's process started. Note that te will be
		 * negative for edges which aren't active yet.
		 */
		const double te = t - e.source_time;
		/* Skip edges which are inactive *and* which are still inactive
		 * at time t + tau. These edges play no role for times [t, t + tau].
		 */
		if (te + tau < 0)
			continue;
		/* Compute the probability that the edge does not fire within
		 * time [te, te + tau] given that it hasn't fired in [0, te],
		 * where te represents t in the edge's frame of reference. Note
		 * that it is possible for an edge to be inactive at time t
		 * (meaning te < 0), but active at time t + tau (i.e. te + tau >= 0).
		 * For such edges, we avoid passing a negative time for the condition
		 * by setting the time to zero instead, which produces the desired
		 * result because it makes the denominator in the conditonal probability
		 * equal to one.
		 */
		const double tp = std::max(te, 0.0);
		/* Pick correct distribution to compute the survival probability with */
		double p;
		switch (e.kind) {
			case event_kind::outside_infection:
			case event_kind::infection:
				p = psi.survivalprobability(te + tau - tp, tp, 1);
				break;
			case event_kind::reset:
				assert(rho);
				p = rho->survivalprobability(te + tau - tp, tp, 1);
				break;
			default:
				throw std::logic_error("unknown event kind");
		}
		/* Update result */
		r *= p;
	}
	return r;
}

interval_t simulate_regir::invphi(absolutetime_t t, double u) {
	return inverse_survival_function(u, p.tau_precision, [&,t] (double tau) { return phi(t, tau); });
}

double simulate_regir::edge_hazard_rate(const active_edges_entry& e, double time) {
	/* Translate t into the edge's frame of reference */
	const double te = time - e.source_time;
	/* For edges not yet active, lambda is zero */
	if (te >= 0) {
		/* Compute lambda, check that it's valid and update */
		double lambda = NAN;
		switch (e.kind) {
			case event_kind::outside_infection:
			case event_kind::infection:
				lambda = psi.hazardrate(te);
				break;
			case event_kind::reset:
				assert(rho);
				lambda = rho->hazardrate(te);
				break;
			default:
				throw std::logic_error("unknown event kind");
		}
		if ((!std::isfinite(lambda) || (lambda < 0)))
			throw std::domain_error("hazardrates must be non-negative and finite");
		return lambda;
	} else {
		return 0;
	}
}
				
void simulate_regir::update_active_edge_lambdas(double time)
{
	double total = 0;
	/* Recompute lambda for every active edge */
	for(const active_edges_entry& ec: active_edges) {
		/* Sets only provide const iterators since elements in a set are immutable.
		 * However, it's actually safe to modify elements provided that the modifications
		 * affect neither the hash value (as computed by the hasher specified in the set's
		 * type) nor equality with other set elements (as defined by the comparator specied
		 * in the set's type). Therefore, the following const cast is safe, provided that we
		 *
		 * DO NOT MODIFIY ANY FIELD USED BY active_edges_hash OR active_edges_cmp!
		 *
		 * Since we only update lambda below, we're OK.
		 */
		active_edges_entry& e = const_cast<active_edges_entry&>(ec);
		e.lambda = edge_hazard_rate(e, time);
		total += e.lambda;
	}
	/* Update sum over all lambdas */
	lambda_total = total;
	if (!std::isfinite(lambda_total))
		throw std::overflow_error("sum over all hazardrates is infinite or NAN");
	if (lambda_total <= 0.0)
		throw std::underflow_error("sum over all hazardrates is zero");
}

auto simulate_regir::draw_active_edge_hazardrates(rng_t& engine) -> active_edges_t::iterator
{
	const double q = unif01_dist(engine) * lambda_total;
	double l = 0.0;
	for(auto i = active_edges.begin(); i != active_edges.end(); ++i) {
		l += i->lambda;
		if (l >= q)
			return i;
	}
	throw std::logic_error("inconsistent state, failed to draw an edge");
}

auto simulate_regir::draw_active_edge_uniform(rng_t& engine) -> active_edges_t::iterator
{
	if (active_edges.empty())
		throw std::logic_error("no active edges state, cannot draw");

	std::uniform_int_distribution<std::size_t> idx_dist(0, active_edges.size() - 1);
	return active_edges.begin() + idx_dist(engine);
}

void simulate_regir::add_active_edge(const active_edges_entry& e) {
	active_edges.push_back(e);
}

void simulate_regir::remove_active_edge(active_edges_t::iterator& it) {
	/* Shouldn't try to remove the next event edge */
	const std::size_t pos = it - active_edges.begin();
	assert(next_event_edge_pos && (*next_event_edge_pos != pos));
	
	/* Swap with last element if not already the last element */
	if (&*it != &active_edges.back()) {
		/* Update next event edge position if necessary */
		if (next_event_edge_pos && (*next_event_edge_pos == active_edges.size()))
			next_event_edge_pos =  it - active_edges.begin();
		/* Exchange */
		using std::swap;
		swap(*it, active_edges.back());
	}
	
	/* Edge to be removed is now the last one */
	active_edges.pop_back();
	
	/* Let iterator point to the "next" element */
	it = active_edges.begin() + pos;
}
