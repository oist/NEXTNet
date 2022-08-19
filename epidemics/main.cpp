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
#include "nMGA.h"

using namespace std;

/*  MersenneTwister random number generator */
/*static*/ mt19937 engine;


int main(int argc, const char * argv[]) {

    int size = 1000;
    int degree = 3;
    double mean = 10;
    double variance = 1.0;
    int threshold = 100;

    // User's input
    if(argc == 2){
        size = atoi(argv[1]);
    } else if (argc == 3){
            size = atoi(argv[1]);
            threshold = atoi(argv[2]);
    }
    


    erdos_reyni network(size,degree, engine);
    transmission_time_lognormal psi(mean, variance); 
    simulate_nmga simulation(network, psi,threshold);


    simulation.add_infections({ std::make_pair(0, 0.0)});
    
    std::vector<double> time_trajectory({});
    std::vector<double> vertex_path({});
    for (int i =0 ; i< size; i++) {
        auto point = simulation.step(engine);
        if (point.second != INFINITY) {
            continue;
        }
        break;
    }


    cout << " done \n";

    return 0;
}

