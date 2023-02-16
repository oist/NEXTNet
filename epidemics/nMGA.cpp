#include "stdafx.h"
#include "nMGA.h"
#include "random.h"
#include "graph.h"
#include "types.h"
#include "utility.h"

graph& simulate_nmga::get_network() const {
	return network;
}

const class transmission_time& simulate_nmga::transmission_time() const {
	return psi;
}

const class transmission_time* simulate_nmga::reset_time() const {
	return rho;
}

bool simulate_nmga::is_infected(node_t node) const {
	return (infected.find(node) != infected.end());
}

void simulate_nmga::add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v)
{
	for(const auto& ve: v) {
		/* Add infection to sorted list of outside infections */
		outside_infections_entry inf;
		inf.node = ve.first;
		inf.time = ve.second;
		outside_infections.push(inf);

		/* Now add outgoing nodes to active edges */
		const int neighbours = network.outdegree(ve.first);
		for(int j=0; j < neighbours; ++j) {
			const node_t neighbour = network.neighbour(ve.first, j);
			if (neighbour < 0) {
				/* This should never happen unless the graph reported the wrong number of outgoing edges */
				throw std::logic_error(std::string("neighbour ") + std::to_string(j + 1) +
										" of node " + std::to_string(ve.first) + " is invalid");
			}
			/* Add future active edge */
			active_edges_entry e;
			e.kind = event_kind::outside_infection;
			e.source_time = ve.second;
			e.source = ve.first;
			e.target = neighbour;
			active_edges.push_back(e);
		}
	}
}

absolutetime_t simulate_nmga::next(rng_t& engine)
{
	/* First, make sure current_time is setup */
	if (std::isnan(current_time)) {
		/* If not set, set to the earliest infection time of any node */
		double t = INFINITY;
		for(const active_edges_entry& e: active_edges)
			t = std::min(t, e.source_time);
		current_time = t;
	}
	
	/* Second, if the current time is infinite, we're done */
	if (std::isinf(current_time))
		return INFINITY;
	
	/* Third, if we have already determined the time of the next event, return it */
	if (!std::isnan(next_time))
		return next_time;

	/* Now, find the time until the next event */
	double tau = NAN;
	const bool use_exact_algorithm = ((approximation_threshold < 0) ||
									  (active_edges.size() <= (unsigned int)approximation_threshold));

	/* First, draw the time until the next event */
	if (use_exact_algorithm) {
		/* Exact version */

		/* Note: The exact version does not use the harard rates lambda, we
		 * thus do not have to update them before drawing tau
		 */
		tau = next_time_exact(engine);
	} else {
		/* Approximate version */

		/*
		 * The following loop fixes an issue in the original NMGA algorithm
		 * If we fail to draw a next event time, either because its unreasonably
		 * large or because all the hazard rates are zero, we skip ahead a bit
		 * and try again.
		 */
		for(;; current_time += maximal_dt) {
			try {
				/* First, update hazard rates lambda and lambda_total */
				update_active_edge_lambdas();
				
				/* Then, draw the time of the next event
				 * Note: Here, this draw *does* depend on the hazard rates
				 * Only accept time increments that dont exceed the maximum
				 * allowed time step!
				 */
				tau = next_time_approximation(engine);
				if (tau <= maximal_dt)
					break;
			} catch (const all_rates_zero& e) {
				/* All rates were zero. This typically happens if
				 * the age distribution hasn't convereged when we switch
				 * to the approximate algorithm. Since the rates are zero,
				 * we assume it's going to be a while since the next event
				 * occurs, and skip ahead maximal_dt time units.
				 */
			}
		}
	}

	/* Check if we would jump over any outside infections,
	 * if return the time of the earliest one instead.
	 */
	while (!outside_infections.empty()) {
		/* Check if the outside infection occurs before the generated tau */
		const outside_infections_entry inf = outside_infections.top();
		if (inf.time > current_time + tau)
			break;

		/* Check if the node is not already infected */
		outside_infections.pop();
		if (infected.find(inf.node) != infected.end())
			break;

		/* Forget previously generated tau, return time of outside infection instead */
		next_time = inf.time;
		return inf.time;
	}
	
	/* Return time of next infection */
	next_time = current_time + tau;
	return next_time;
}

std::optional<event_t> simulate_nmga::step(rng_t& engine, absolutetime_t max_time, event_filter_t evf)
{
    while (true) {
		/* Find time of next event unless already done */
		if (std::isnan(next_time))
			next(engine);
		
		/* If we'd run past max_time, don't to anything */
		if (next_time > max_time)
			return std::nullopt;
		
        /* Check if the next event is an outside infection */
        while (!outside_infections.empty()) {
            /* Check if the outside infection occurs before the generated tau */
            const outside_infections_entry inf = outside_infections.top();
            if (inf.time > next_time)
                break;

            /* Check if the node is not already infected */
            outside_infections.pop();
            if (infected.find(inf.node) != infected.end())
                break;

            /* Advance time to time of outside infection, forget next_time
			 * so that we generate a new one next time. Note that the neighbours
			 * of the node that gets infected here were already made active by add_infections()
             */
			assert(current_time <= inf.time);
			assert(next_time == inf.time);
            current_time = inf.time;
			next_time = NAN;
            infected.insert(inf.node);

            /* Report event */
            return event_t { .kind = event_kind::outside_infection,
                             .node = inf.node, .source_node = -1,
							 .time = current_time };
        }

        /* Now determine which event takes place.*/
		const bool use_exact_algorithm = ((approximation_threshold < 0) ||
										  (active_edges.size() <= (unsigned int)approximation_threshold));
		std::vector<active_edges_entry>::iterator edge_i;
        if (use_exact_algorithm) {
            /* Exact version */

            /* Then, update the current time, reset next_time */
            current_time = next_time;
			next_time = NAN;

            /* And update lambdas and lambda_total since we didn't do so before. */
            update_active_edge_lambdas();

            /* Now determine which event takes place */
			edge_i = draw_active_edge(engine);
        } else {
            /* Approximate version */
            
           /* Note that we do not recompute the hazard rates that we
            * updated before choosing tau above, which amounts
            * to doing the draw *before* updating the current time
            */
			edge_i = draw_active_edge(engine);

            /* Finally, update current time, reset next_time */
            current_time = next_time;
			next_time = NAN;
        }
        
        /* Copy selected edge and remove from active edge list */
		const active_edges_entry edge = *edge_i;
        remove_active_edge(edge_i);
		
		/* Create event and query filter */
		event_t ev { .kind = edge.kind, .node = edge.target, .source_node = edge.source, .time = current_time };
		if (is_event_blocked(ev, evf))
			continue;
		
		/* Handle event */
		switch (edge.kind) {
			case event_kind::infection:
			case event_kind::outside_infection: {
				/* If the node is already infected, ignore the event */
				if (infected.find(edge.target) != infected.end())
					continue;
				
				/* Mark node as infected */
				infected.insert(edge.target);
				
				/* Make recovery self-loop if there's a reset time distribution */
				if (rho) {
					active_edges_entry e;
					e.kind = event_kind::reset;
					e.source = edge.target;
					e.source_time = current_time;
					e.target = edge.target;
					add_active_edge(e);
				}
				
				/* and outgoing edges active */
				const int neighbours = network.outdegree(edge.target);
				for(int j=0; j < neighbours; ++j) {
					const node_t neighbour = network.neighbour(edge.target, j);
					if (neighbour < 0) {
						/* This should never happen unless the graph reported the wrong number of outgoing edges */
						throw std::logic_error(std::string("neighbour ") + std::to_string(j + 1) +
												" of node " + std::to_string(edge.target) + " is invalid");
					}
					active_edges_entry e;
					e.kind = event_kind::infection;
					e.source = edge.target;
					e.source_time = current_time;
					e.target = neighbour;
					add_active_edge(e);
				}
				break;
			}
				
			case event_kind::reset: {
				/* Reset event */
				
				/* if SIR, increase counter of removed/recovered nodes, and leave the node as infected so that 
				* the node cannot be reinfected. 
				* ( this is just a trick to avoid having to create a new state, the node is not actually infected anymore and does not generate new infections.)
				* else, Mark node as no longer infected. In that case the node simply returns to the susceptible state and can get reinfected again later. */
				if (SIR){
					removed += 1;
				} else { // SIS
					/* Mark node as not infected */
					infected.erase(edge.target);
				}


				
				/* Remove active edges originating from the resetted node */
				for(auto i = active_edges.begin(); i != active_edges.end();) {
					/* Skip edges originating anywhere else */
					if (i->source != edge.target) {
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
		
        /* Return event */
		return ev;
    }
}

void simulate_nmga::notify_infected_node_neighbour_added(network_event_t event)
{
	throw std::logic_error("unimplemented");
}

interval_t simulate_nmga::next_time_exact(rng_t& engine) {
    /* Determine time of next event by inverting the global survival function phi */
    const double u = unif01_dist(engine);
    const double tau = invphi(current_time, u);
    return tau;
}

interval_t simulate_nmga::next_time_approximation(rng_t& engine)
{
    /* Determine time of next event using equation (8) of Boguna. */
    return std::exponential_distribution<double>(lambda_total)(engine);
}

double simulate_nmga::find_maximal_dt(const class transmission_time& psi) {
	// We find dt such that F(t + dt) - F(t) < dp,
	// meaning such that moving from t to dt skips over
	// at most probability mass dp
	const double dp = 0.01;
	double max_dt = INFINITY;
	for(double p = 1; p - dp > 0; p -= dp) {
		const double dt = psi.survivalquantile(p - dp ) - psi.survivalquantile(p);
		max_dt = std::min(max_dt,dt);
	}
	return max_dt;
}


double simulate_nmga::phi(absolutetime_t t, interval_t tau) {
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

interval_t simulate_nmga::invphi(absolutetime_t t, double u) {
	return inverse_survival_function(u, tau_precision, [&,t] (double tau) { return phi(t, tau); });
}

void simulate_nmga::update_active_edge_lambdas()
{
	double total = 0;
	/* Recompute lambda for every active edge */
	for(active_edges_entry& e: active_edges) {
		/* Translate t into the edge's frame of reference */
		const double te = current_time - e.source_time;
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
			e.lambda = lambda;
			total += lambda;
		} else {
			e.lambda = 0;
		}
	}
	/* Update sum over all lambdas */
	lambda_total = total;
	if (!std::isfinite(lambda_total))
		throw std::overflow_error("sum over all hazardrates is infinite or NAN");
	if (lambda_total <= 0.0)
		throw all_rates_zero();
}

auto simulate_nmga::draw_active_edge(rng_t& engine) -> std::vector<active_edges_entry>::iterator
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

