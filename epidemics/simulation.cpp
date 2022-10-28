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
			std::optional<event_t> ev = simulation.step(engine);
			if (!ev) break;
			time_trajectory.push_back(ev->time);
        }
        std::cout << path<< std::endl;
        std::string filename("data");
        std::string ext(".dat");
        exportData(time_trajectory,filename+file_nb+ext);
    }
    
 

}



void generatePaths_NMGA(double mean, double variance, int degree,int nb_paths,double size, rng_t& engine,double threshold){
    
//    std::string filename ="datanmga";
//    std::string ext= ".dat";
    //exportData(time_trajectory,filename+file_nb+ext);
    
    
    for (int path=0; path< nb_paths; path++) {
        std::string file_nb = std::to_string(path);

        erdos_reyni network(size, degree, engine);
        transmission_time_lognormal psi(mean, variance); 
        simulate_nmga simulation(network, psi);
        simulation.approximation_threshold = threshold;
    //    simulate_next_reaction simulation(network);
        for (int i =0; i<1; i++) {
            simulation.add_infections({ std::make_pair(i, 0.0)});
        }
        


        //std::vector<double> time_trajectory({});
        //std::vector<double> vertex_path({});
        for (int i =0 ; i< size; i++) {
            auto ev = simulation.step(engine);
            if (ev) {
                //vertex_path.push_back(ev->node);
                //time_trajectory.push_back(ev->time);
                continue;
            }
            break;
        }
        //std::cout << path<< std::endl;
        //exportData(time_trajectory,filename+file_nb+ext);
    }
}

void generatePaths_next_reaction(double mean, double variance, int degree,int nb_paths,double size, rng_t& engine){
    
    std::string filename ="data";
    std::string ext= ".dat";
    
    std::string adjalist ="adjaLIST.csv";
    
    for (int path=0; path< nb_paths; path++) {
        std::string file_nb = std::to_string(path);

        erdos_reyni network(size, degree, engine);
        export_adjacency_list(network.adjacencylist,adjalist);
        transmission_time_lognormal psi(mean, variance); 
        simulate_next_reaction simulation(network, psi);
        simulation.add_infections({ std::make_pair(0, 0.0)});

        std::vector<double> time_trajectory({});
        std::vector<double> vertex_path({});
        for (int i =0 ; i< size; i++) {
            std::optional<event_t> ev = simulation.step(engine);
			if (!ev) break;
			time_trajectory.push_back(ev->time);
        }
        std::cout << path<< std::endl;
        exportData(time_trajectory,filename+file_nb+ext);
    }
}
