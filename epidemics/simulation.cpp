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
#include "nMGA.h"
#include "NextReaction.h"
#include "analysis.h"



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
    
    std::string adjalist ="adjaLIST.csv";
    
    for (int path=0; path< nb_paths; path++) {
        std::string file_nb = std::to_string(path);

        erdos_reyni network(size, degree, engine);
        export_adjacency_list(network.neighbours,adjalist);
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
