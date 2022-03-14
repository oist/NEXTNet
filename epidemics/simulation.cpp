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

std::pair<node_t, absolutetime_t> simulate_next_reaction::step() {
    while (true) {
        /* If there are no more infection times, stop */
        if (infectiontimes.empty())
            return std::make_pair(-1, INFINITY);

        /* Fetch the next putatively infected node */
        const auto next = infectiontimes.top();
        infectiontimes.pop();
        
        /* Lazily enqueue next sibling of the putatively infected node if necessary */
        if (next.neighbour_index >= 0) {
            const auto sibling = network.neighbour(next.source_node, next.neighbour_index+1);
            if ((sibling.first >= 0) && (std::isfinite(sibling.second))) {
                /* Create sibling's infection times entry and add to queue */
                infectiontimes_entry e;
                e.time = next.source_time + sibling.second;
                e.node = sibling.first;
                e.source_time = next.source_time;
                e.source_node = next.source_node;
                e.neighbour_index = next.neighbour_index + 1;
                infectiontimes.push(e);
            }
        }
        
        /* Check if the putatively infected node is already infected, if so we're done */
        if (infected.find(next.node) != infected.end())
            continue;
        
        /* Mark the node as infected */
        infected.insert(next.node);
        
        /* Add the infecte node's first neigbhour to the infection times queue */
        const auto neighbour = network.neighbour(next.node, 0);
        if ((neighbour.first >= 0) && (std::isfinite(neighbour.second))) {
            infectiontimes_entry e;
            e.time = next.time + neighbour.second;
            e.node = neighbour.first;
            e.source_time = next.time;
            e.source_node = next.node;
            e.neighbour_index = 0;
            infectiontimes.push(e);
        }
        
        return std::make_pair(next.node, next.time);
    }
}

void simulate_next_reaction::add_infections(const std::vector<std::pair<node_t, absolutetime_t>>& v) {
    for(const auto& ve: v) {
        infectiontimes_entry e;
        e.time = ve.second;
        e.node = ve.first;
        infectiontimes.push(e);
    }
}

std::vector<double> simulatePath(std::vector<double>& infection_times, int n_max,const lognormal_beta& infection_distribution, rng_t engine){
    
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

void simulateManyPaths(int nb_paths, rng_t& engine){
    int size = 10000;
    int degree = 1;
    double mean = 1;
    double variance = 10;
    
    
    for (int path=1; path<= nb_paths; path++) {
        std::string file_nb = std::to_string(path);
        
        erdos_reyni network(size,degree,lognormal_beta(mean,variance,degree), engine);
        
        simulate_next_reaction simulation(network);
        simulation.add_infections({ std::make_pair(0, 0.0)});
        
        std::vector<double> time_trajectory({});
        std::vector<double> vertex_path({});
        for (int i =0 ; i< size; i++) {
            auto point = simulation.step();
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
