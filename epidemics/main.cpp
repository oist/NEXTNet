//
//  main.cpp
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "simulation.h"
#include "graph.h"

using namespace std;

/*  MersenneTwister random number generator */
/*static*/ mt19937 engine;


int main(int argc, const char * argv[]) {

//    simulateManyPaths(1000, mersenneTwister);
    
    int size = 3000;
    int degree = 3;
    double mean = 10;
    double variance = 1.2;
    
    lognormal_beta my_log_norm(mean,variance,degree);
    erdos_reyni network(size,degree, my_log_norm, engine);
    
    simulate_nmga simulation(network,my_log_norm);
    simulation.approximation_threshold = 100;
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
    
    exportData(time_trajectory,"data.dat");
    
    
//
//    erdos_reyni network(size,degree,lognormal_beta(mean,variance,degree) , mersenneTwister);
//
//    vector<vector<pair<node_t,interval_t>>> adjacency_list = network.neighbours;
//
//    export_adjacency_list(adjacency_list, "adjacency_list.csv");
//
//
//
//
//    simulator simulation(network);
//    simulation.add_infections({make_pair(0, 0.0)});
//
//    vector<double> time_trajectory({});
//    vector<double> vertex_path({});
//    for (int i =0 ; i< size; i++) {
//        time_trajectory.push_back(simulation.step().second);
//        vertex_path.push_back(simulation.step().first);
//    }
//
//    exportData(time_trajectory,"data2.dat");
//
////    exportData(vertex_path, "path.dat");
    return 0;
}
