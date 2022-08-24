#include "stdafx.h"
#include "nMGA.h"
#include "random.h"
#include "graph.h"
#include "types.h"
#include "utility.h"

double simulate_nmga::find_maximal_dt(transmission_time& psi) {
    // We find dt such that F(t + dt) - F(t) < dp,
    // meaning such that moving from t to dt skips over
    // at most probability mass dp
    const double dp = 0.1;
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
        const double p = psi.survivalprobability(te + tau - tp, tp, 1);
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
            const double lambda = psi.hazardrate(te);
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
    if (lambda_total <= 0.0)
        throw std::underflow_error("all active edges report a hazardrate of zero (or negative)");
    else if (std::isinf(lambda_total))
        throw std::overflow_error("sum over all hazardrates is infinite");
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

std::pair<node_t, absolutetime_t> simulate_nmga::step(rng_t& engine)
{
    while (true) {
        /* If the current time wasn't yet set, start at the earliest time at which
        * and edge becomes active. */
        if (std::isnan(current_time)) {
            double t = INFINITY;
            for(const active_edges_entry& e: active_edges)
                t = std::min(t, e.source_time);
            current_time = t;
        }
        
        /* If there are no active edges, the next event is at infinity */
        if (active_edges.empty())
            current_time = INFINITY;
        if (std::isinf(current_time))
            return std::make_pair(-1, current_time);
        
        /* Find the time of the next event */
        node_t node;
        std::vector<active_edges_entry>::iterator edge;
        if ((approximation_threshold >= 0) && (active_edges.size() <= (unsigned int)approximation_threshold)) {
            /* Exact version */

            /* First, draw the time of the next event
            * Note: This step does not use the harard rates lambda, we
            * thus do not have to update them before drawing tau
            */
            const double tau = next_time_exact(engine);

            /* Then, update the current time */
            current_time += tau;        

            /* and update lambdas and lambda_total. */
            update_active_edge_lambdas();

            /* Now determine which event takes place */
            edge = draw_active_edge(engine);
            node = edge->target;
        } else {
            /* Approximate version */

            /* The following fixes an issue in the original NMGA algorithm */
            double tau = NAN;
            while (true) {
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
                
                /* Make maximum allowed time step and try again */
                current_time += maximal_dt;
            }
            
            /* Now determine which event takes place.
            * Note that we do not recompute the hazard rates, which amounts
            * to doing the draw *before* updating the current time
            */
            edge = draw_active_edge(engine);
            node = edge->target;

            /* Finally, update current time */
            current_time += tau;    
        }
        
        /* Remove selected edge from active edge list */
        remove_active_edge(edge);

        /* If the node is already infected, ignore the event */
        if (infected.find(node) != infected.end())
            continue;

        /* Mark node as infected and make outgoing edges active */
        infected.insert(node);
        const int neighbours = network.outdegree(node);
        for(int j=0; j < neighbours; ++j) {
            const node_t neighbour = network.neighbour(node, j);
            if (neighbour < 0) {
                /* This should never happen unless the graph reported the wrong number of outgoing edges */
                throw std::logic_error(std::string("neighbour ") + std::to_string(j + 1) +
                                        " of node " + std::to_string(node) + " is invalid");
            }
            active_edges_entry e;
            e.source = node;
            e.source_time = current_time;
            e.target = neighbour;
            add_active_edge(e);
        }
        
        return std::make_pair(node, current_time);
    }
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

void simulate_nmga::add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v)
{
    for(const auto& ve: v) {
        /* Iterator over source node's neighbours */
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
            e.source_time = ve.second;
            e.source = ve.first;
            e.target = neighbour;
            active_edges.push_back(e);
        }
    }
}
