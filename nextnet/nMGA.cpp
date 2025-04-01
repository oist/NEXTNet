#include "stdafx.h"
#include "nMGA.h"
#include "random.h"
#include "network.h"
#include "types.h"
#include "utility.h"

network &simulate_nmga::get_network() const
{
    return nw;
}

const class transmission_time &simulate_nmga::transmission_time() const
{
    return psi;
}

const class transmission_time *simulate_nmga::reset_time() const
{
    return rho;
}

bool simulate_nmga::is_infected(node_t node) const
{
    return (infected.find(node) != infected.end());
}

void simulate_nmga::add_infections(const std::vector<std::pair<node_t, absolutetime_t>> &v)
{
    for (const auto &ve : v) {
        /* Add infection to sorted list of outside infections */
        outside_infections_entry inf;
        inf.node = ve.first;
        inf.time = ve.second;
        outside_infections.push(inf);
    }
}

absolutetime_t simulate_nmga::next(rng_t &engine)
{
    /* If the next event was already determined, just return it */
    if (next_event)
        return next_event->time;

    /* If the current time wasn't yet set, start at the earliest time at which
     * and edge becomes active. */
    absolutetime_t base_time = current_time;
    if (std::isnan(base_time)) {
        base_time = outside_infections.empty() ? INFINITY : outside_infections.top().time;
        for (const active_edges_entry &e : active_edges)
            base_time = std::min(base_time, e.source_time);
    }

    /* If there are no active edges and no outside infections, there is no next event */
    if (active_edges.empty() && outside_infections.empty())
        base_time = INFINITY;
    if (std::isinf(base_time))
        return INFINITY;

    /* Find the time of the next event */
    const bool use_exact_algorithm = ((p.approximation_threshold < 0) ||
                                      (active_edges.size() <= (unsigned int)p.approximation_threshold));

    /* First, draw the time of the next event */
    double tau = NAN;
    if (active_edges.empty()) {
        /* No active edges, but there are outside infections (see below) */

        tau = INFINITY;
    } else if (use_exact_algorithm) {
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
        for (;; base_time += max_dt) {
            try {
                /* First, update hazard rates lambda and lambda_total */
                update_active_edge_lambdas(base_time);

                /* Then, draw the time of the next event
                 * Note: Here, this draw *does* depend on the hazard rates
                 * Only accept time increments that dont exceed the maximum
                 * allowed time step!
                 */
                tau = next_time_approximation(engine);
                if (tau <= max_dt)
                    break;
            } catch (const all_rates_zero &e) {
                /* All rates were zero. This typically happens if
                 * the age distribution hasn't convereged when we switch
                 * to the approximate algorithm. Since the rates are zero,
                 * we assume it's going to be a while since the next event
                 * occurs, and skip ahead max_dt time units.
                 */
            }
        }
    }

    /* Check if we would jump over any outside infections. If there is an
     * outside infection, we discard the tau we just generated but that's OK.
     */
    while (!outside_infections.empty()) {
        /* Check if the outside infection occurs before the generated tau */
        const outside_infections_entry inf = outside_infections.top();
        if (inf.time > base_time + tau)
            break;

        /* Skip and remove if the node is already infected */
        if (infected.find(inf.node) != infected.end()) {
            outside_infections.pop();
            continue;
        }

        /* Found next event, store & return its time */
        next_event = epidemic_event_t{
            .kind        = epidemic_event_kind::outside_infection,
            .source_node = -1,
            .node        = inf.node,
            .time        = inf.time
        };
        return next_event->time;
    }
    assert((tau >= 0) && (tau < INFINITY));
    const double next_time = base_time + tau;

    /* Now determine which event takes place.
     * Since the exact algorithm for drawing the next time does not update the hazard rates, do so now.
     * But in case of the approximate algorith, re-use the hazard rates that were used to draw tau.
     */
    if (use_exact_algorithm)
        update_active_edge_lambdas(next_time);
    const active_edges_t::iterator edge_i = draw_active_edge(engine);

    /* Found next event, store & return its time */
    assert(edge_i->kind != epidemic_event_kind::outside_infection);
    next_event = epidemic_event_t{
        .kind        = edge_i->kind,
        .source_node = edge_i->source,
        .node        = edge_i->target,
        .time        = next_time
    };
    return next_event->time;
}

std::optional<epidemic_event_t> simulate_nmga::step(rng_t &engine, absolutetime_t maxtime, event_filter_t evf)
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
        const epidemic_event_t ev = *next_event;
        current_time              = next_event->time;
        next_event                = std::nullopt;

        /* Remove corresponding active edge or outside infection */
        switch (ev.kind) {
            case epidemic_event_kind::outside_infection: {
                if (outside_infections.empty() || (outside_infections.top().time != ev.time) || (outside_infections.top().node != ev.node))
                    throw std::logic_error("next outside infection changed between next() and step()");
                outside_infections.pop();
                break;
            }

            case epidemic_event_kind::infection:
            case epidemic_event_kind::reset: {
                /* Re-find active edge.
                 * NOTE: It would be more efficient to pass the iterator from next() to here. However,
                 * that would make it impossible to check whether the edge is actually still part of the
                 * active edges set, which is a reasonable safe-guard against coding mistakes
                 */
                const active_edges_entry ev_edge = { .kind = ev.kind, .source = ev.source_node, .target = ev.node };
                active_edges_t::iterator edge_i  = active_edges.find(ev_edge);
                if (edge_i == active_edges.end())
                    throw std::logic_error("edge became inactive between next() and step()");
                remove_active_edge(edge_i);
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
            case epidemic_event_kind::infection:
            case epidemic_event_kind::outside_infection: {
                /* Infection event */

                /* If the node is already infected, ignore the event */
                if (infected.find(ev.node) != infected.end())
                    continue;

                /* Mark node as infected */
                infected.insert(ev.node);

                /* Make recovery self-loop if there's a reset time distribution ... */
                if (rho) {
                    active_edges_entry e;
                    e.kind        = epidemic_event_kind::reset;
                    e.source      = ev.node;
                    e.source_time = current_time;
                    e.target      = ev.node;
                    add_active_edge(e);
                }

                /* ... and outgoing edges active */
                const int neighbours = nw.outdegree(ev.node);
                for (int j = 0; j < neighbours; ++j) {
                    const node_t neighbour = nw.neighbour(ev.node, j);
                    if (neighbour < 0) {
                        /* This should never happen unless the graph reported the wrong number of outgoing edges */
                        throw std::logic_error(std::string("neighbour ") + std::to_string(j + 1) +
                                               " of node " + std::to_string(ev.node) + " is invalid");
                    }
                    active_edges_entry e;
                    e.kind        = epidemic_event_kind::infection;
                    e.source      = ev.node;
                    e.source_time = current_time;
                    e.target      = neighbour;
                    add_active_edge(e);
                }
                break;
            }

            case epidemic_event_kind::reset: {
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
                for (auto i = active_edges.begin(); i != active_edges.end();) {
                    /* Skip edges originating anywhere else */
                    if (i->source != ev.node) {
                        ++i;
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

void simulate_nmga::notify_infected_node_neighbour_added(network_event_t event, rng_t &engine)
{
    throw std::logic_error("nGMA currently does not support dynamic networks");
}

interval_t simulate_nmga::next_time_exact(rng_t &engine)
{
    /* Determine time of next event by inverting the global survival function phi */
    const double u   = unif01_dist(engine);
    const double tau = invphi(current_time, u);
    return tau;
}

interval_t simulate_nmga::next_time_approximation(rng_t &engine)
{
    /* Determine time of next event using equation (8) of Boguna. */
    return std::exponential_distribution<double>(lambda_total)(engine);
}

double simulate_nmga::find_maximal_dt(const class transmission_time &psi)
{
    // We find dt such that F(t + dt) - F(t) < dp,
    // meaning such that moving from t to dt skips over
    // at most probability mass dp
    const double dp = 0.01;
    double max_dt   = INFINITY;
    for (double p = 1; p - dp > 0; p -= dp) {
        const double dt = psi.survivalquantile(p - dp) - psi.survivalquantile(p);
        max_dt          = std::min(max_dt, dt);
    }
    return max_dt;
}

double simulate_nmga::phi(absolutetime_t t, interval_t tau)
{
    /* This implements our version of equation 4 from the Boguna paper */
    if (std::isinf(tau))
        return 0;
    double r = 1;
    for (const active_edges_entry &e : active_edges) {
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
            case epidemic_event_kind::outside_infection:
            case epidemic_event_kind::infection:
                p = psi.survivalprobability(te + tau - tp, tp, 1);
                break;
            case epidemic_event_kind::reset:
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

interval_t simulate_nmga::invphi(absolutetime_t t, double u)
{
    return inverse_survival_function(u, p.tau_precision, [&, t](double tau) { return phi(t, tau); });
}

void simulate_nmga::update_active_edge_lambdas(double time)
{
    double total = 0;
    /* Recompute lambda for every active edge */
    for (const active_edges_entry &ec : active_edges) {
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
        active_edges_entry &e = const_cast<active_edges_entry &>(ec);
        /* Translate t into the edge's frame of reference */
        const double te = time - e.source_time;
        /* For edges not yet active, lambda is zero */
        if (te >= 0) {
            /* Compute lambda, check that it's valid and update */
            double lambda = NAN;
            switch (e.kind) {
                case epidemic_event_kind::outside_infection:
                case epidemic_event_kind::infection:
                    lambda = psi.hazardrate(te);
                    break;
                case epidemic_event_kind::reset:
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

auto simulate_nmga::draw_active_edge(rng_t &engine) -> active_edges_t::iterator
{
    const double q = unif01_dist(engine) * lambda_total;
    double l       = 0.0;
    for (auto i = active_edges.begin(); i != active_edges.end(); ++i) {
        l += i->lambda;
        if (l >= q)
            return i;
    }
    throw std::logic_error("inconsistent state, failed to draw an edge");
}
