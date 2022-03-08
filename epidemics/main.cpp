//
//  main.cpp
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#include "stdafx.h"
#include "random.h"
#include "savingData.h"
#include "simulation.h"
#include "testing.h"
#include "graph.h"

using namespace std;

/*  MersenneTwister random number generator */
/*static*/ mt19937 mersenneTwister;


int main(int argc, const char * argv[]) {
    
    setDebugEnabled(false);// testing or not.
//
//
//    /*--------- Initialise Parameters ----------*/
//
//    Tau tau(10,1,2); //Mean Variance R0
//
////    int size = 10;
////    double degree = 3;
////    simulatePathNetwork(size,degree, tau, mersenneTwister);
//    vector<double> infection_times({});
//    vector <double> trajectory = simulatePath(infection_times, 10000, tau, mersenneTwister);
    
    
    int size = 10000;
    int degree = 3;
    double mean = 10;
    double variance = 1;
    
    
    erdos_reyni network(size,degree,lognormal_beta(mean,variance,degree) , mersenneTwister);
    
    simulator simulation(network);
    simulation.add_infections({make_pair(0, 0.0)});
    
    vector<double> trajectory({});
    
    for (int i =0 ; i<size; i++) {
        trajectory.push_back(simulation.step().second);
    }
    
    exportData(trajectory, "data.dat");
    
    exportData(trajectory, "data.dat");
    return 0;
}
