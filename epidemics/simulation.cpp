//
//  simulation.cpp
//  epidemics_update
//
//  Created by Cure Samuel Cyrus on 2022/02/19.
//

#include "stdafx.h"
#include "simulation.h"
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

double simulate_nmga::phi(absolutetime_t t, interval_t tau) {
    /* This implements our version of equation 4 from the Boguna paper */
    if (std::isinf(tau))
        return 0;
    double r = 1;
    for(const active_edges_entry& e: active_edges) {
        const double te = t - e.source_time;
        /* Skip edges not yet active at relative time tau */
        if (te + tau < 0)
            continue;
        /* Multiply with the survival probability of the current edge,
         * taking into account that it might activate between time t
         * and tau, in which case the 'last firing time' is zero
         */
        r *= psi.survivalprobability(te + tau, std::max(te, 0.0), 1);
    }
    return r;
}

interval_t simulate_nmga::invphi(absolutetime_t t, double u) {
    return inverse_survival_function(u, tau_precision, [&,t] (double tau) { return phi(t, tau); });
}

double simulate_nmga::lambda_total(std::vector<double>* lambda_finite_cumulative,
                                   std::vector<std::size_t>* lambda_infinite) {
    double lambda_total = 0;
    double lambda_finite_total = 0;
    for(std::size_t i=0; i < active_edges.size(); ++i) {
        const active_edges_entry& e = active_edges[i];
        const double te = current_time - e.source_time;
        double lambda = 0;
        /* For edges not yet active, lambda is zero */
        if (te >= 0)
            lambda = psi.hazardrate(te);
        lambda_total += lambda;
        /* Collect both the grand total and the total across all finite lambdas */
        if (std::isfinite(lambda))
            lambda_finite_total += lambda;
        else if (lambda_infinite)
            lambda_infinite->push_back(i);
        if (lambda_finite_cumulative)
            lambda_finite_cumulative->push_back(lambda_finite_total);
    }
    return lambda_total;
}

std::pair<node_t, absolutetime_t> simulate_nmga::step(rng_t& engine) {
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
    double tau;
    if ((approximation_threshold >= 0) && (active_edges.size() <= (unsigned int)approximation_threshold))
        tau = next_time_exact(engine);
    else
        tau = next_time_approximation(engine);
    
    /* Advance current time */
    current_time += tau;
    
    /* Determine which event takes place */
    active_edges_cumulative_finite_lambdas.clear();
    active_edges_infinite_lambdas.clear();
    lambda_total(&active_edges_cumulative_finite_lambdas, &active_edges_infinite_lambdas);
    std::size_t edge_i;
    if (active_edges_infinite_lambdas.empty()) {
        const double ltotal = active_edges_cumulative_finite_lambdas.back();
        if (ltotal == 0) {
            /* All lambdas were zero, draw an active edge uniformly */
            std::uniform_int_distribution<std::size_t> idist(0, active_edges.size() - 1);
            edge_i = idist(engine);
        } else {
            /* All lambdas were finite and some non-zero , draw by inverting the CDF */
            std::uniform_real_distribution<double> pdist(0, ltotal);
            const double p = pdist(engine);
            edge_i = 0;
            while ((active_edges_cumulative_finite_lambdas.at(edge_i) <= p))
                ++edge_i;
        }
    } else {
        /* Some lambdas were infinite, draw uniformly from those */
        std::uniform_int_distribution<std::size_t> idist(0, active_edges_infinite_lambdas.size() - 1);
        edge_i = active_edges_infinite_lambdas.at(idist(engine));
    }
    const node_t node = active_edges.at(edge_i).target;

    
    /* Mark node as infected and make edges to non-infected neighbours active */
    infected.insert(node);
    for(int j=0; ; ++j) {
        const node_t neighbour = network.neighbour(node, j);
        if (neighbour < 0)
            break;
        if (infected.find(neighbour) != infected.end())
            continue;
        active_edges_entry e;
        e.source = node;
        e.source_time = current_time;
        e.target = neighbour;
        add_active_edge(e);
    }
    
    /* Remove active edges whose target just because infected */
    for(auto it = active_edges.begin(); it != active_edges.end();) {
        /* Removing an edge actually replaced the edge"s entry with a different once,
         * therefore we don't want to increment the iterator after removing something
         */
        if (it->target == node)
            remove_active_edge(it);
        else
            ++it;
    }
    
    return std::make_pair(node, current_time);
}

interval_t simulate_nmga::next_time_exact(rng_t& engine) {
    /* Determine time of next event by inverting the global survival function phi */
    const double u = unif01_dist(engine);
    const double tau = invphi(current_time, u);
    
    return tau;
}

interval_t simulate_nmga::next_time_approximation(rng_t& engine) {
    /* Determine time of next event using equation (8) of Boguna.
     * If lambda_total is infinite, the next event happens at the current time, and
     * we therefore return tau = 0
     */
    const double ltotal = lambda_total();
    if (std::isinf(ltotal))
        return 0;
    std::exponential_distribution<double> tau_dist(ltotal);
    const double tau = tau_dist(engine);
    
    return tau;
}

void simulate_nmga::add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v) {
    for(const auto& ve: v) {
        /* Iterator over source node's neighbours */
        for(int i=0; ; ++i) {
            /* Query next neighbour, abort if no more neighbours exist */
            const node_t neighbour = network.neighbour(ve.first, i);
            if (neighbour < 0)
                break;
            
            /* Add future active edge */
            active_edges_entry e;
            e.source_time = ve.second;
            e.source = ve.first;
            e.target = neighbour;
            active_edges.push_back(e);
        }
    }
}


#if 0
std::vector<double> simulatePath(std::vector<double>& infection_times, int n_max,const lognormal_beta& infection_distribution, rng_t& engine){
    
    double absolute_time = 0;
    std::vector<double> time_trajectory({});

    /* Start with one infected individual at time t=0*/
    for (int population = 1; population <n_max; population++ ){

        
        /*----Add new infection times----*/
        int new_infections = poissrnd(infection_distribution.r0,engine);

        for (int i =0; i < new_infections; i++)
        infection_times.push_back( infection_distribution.sample(engine) + absolute_time);  // the infection times are not relative to the age of the infected indiv. They are relative to the general time of the system.


        /*----Determine the next infection time in the population----*/
        if (infection_times.size()==0) { // If there are no new infections times, the simulation is over.
            break;
        }
        
        long index = min_element(infection_times.begin(),infection_times.end()) - infection_times.begin();//Find the position of minimum of infection_times
        double next_infection_time = infection_times[index]; //Find minimum
        
        /*----Update trajectory----*/
        absolute_time = next_infection_time;
        time_trajectory.push_back(absolute_time);
        
        /*----Remove that infection time from the list----*/
        infection_times.erase(infection_times.begin()+index); // Remove next_infection_time;

    }
    
    return time_trajectory;
    
}
#endif

void print_matrix(const std::vector<std::vector<double>>& A){
    long rows = A.size();
    long cols = A[0].size();
    
    for (int i = 0; i<rows; i++) {
        for (int j =0; j<cols; j++) {
            std::cout << ceil(A[i][j]*10)/10 << "   ";
        }
        std::cout << std::endl;
    }
}

void simulatePaths_MeanField(double mean, double variance, int degree,int nb_paths,double size, rng_t& engine){
        
    for (int path=1; path<= nb_paths; path++) {
        std::string file_nb = std::to_string(path);
        

        erdos_reyni network(size,degree, engine);
        transmission_time_lognormal psi(mean, variance); 
        simulate_next_reaction simulation(network, psi);
        simulation.add_infections({ std::make_pair(0, 0.0)});
        
        std::vector<double> time_trajectory({});
        std::vector<double> vertex_path({});
        for (int i =0 ; i< size; i++) {
            auto point = simulation.step(engine);
            if (point.second != INFINITY) {
                time_trajectory.push_back(point.second);
                continue;
            }
            break;
        }
        std::cout << path<< std::endl;
        std::string filename("data");
        std::string ext(".dat");
        exportData(time_trajectory,filename+file_nb+ext);
    }
    
 

}



void generatePaths_NMGA(double mean, double variance, int degree,int nb_paths,double size, rng_t& engine,double threshold){
    
    std::string filename ="data";
    std::string ext= ".dat";
    
    
    for (int path=0; path< nb_paths; path++) {
        std::string file_nb = std::to_string(path);

        erdos_reyni network(size, degree, engine);
        transmission_time_lognormal psi(mean, variance); 
        simulate_nmga simulation(network, psi);
        simulation.approximation_threshold = threshold;
    //    simulate_next_reaction simulation(network);
        simulation.add_infections({ std::make_pair(0, 0.0)});

        std::vector<double> time_trajectory({});
        std::vector<double> vertex_path({});
        for (int i =0 ; i< size; i++) {
            auto point = simulation.step(engine);
            if (point.second != INFINITY) {
                vertex_path.push_back(point.first);
                time_trajectory.push_back(point.second);
                continue;
            }
            break;
        }
        std::cout << path<< std::endl;
        exportData(time_trajectory,filename+file_nb+ext);
    }
}

void generatePaths_next_reaction(double mean, double variance, int degree,int nb_paths,double size, rng_t& engine){
    
    std::string filename ="data";
    std::string ext= ".dat";
    
    
    for (int path=0; path< nb_paths; path++) {
        std::string file_nb = std::to_string(path);

        erdos_reyni network(size, degree, engine);
        transmission_time_lognormal psi(mean, variance); 
        simulate_next_reaction simulation(network, psi);
        simulation.add_infections({ std::make_pair(0, 0.0)});

        std::vector<double> time_trajectory({});
        std::vector<double> vertex_path({});
        for (int i =0 ; i< size; i++) {
            auto point = simulation.step(engine);
            if (point.second != INFINITY) {
                time_trajectory.push_back(point.second);
                continue;
            }
            break;
        }
        std::cout << path<< std::endl;
        exportData(time_trajectory,filename+file_nb+ext);
    }
}

