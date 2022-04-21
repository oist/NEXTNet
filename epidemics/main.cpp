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

 
    // PARAMETERS
    
    
    int size = 1000; // Size network or max pop for mean field.
    int degree = 5; // Average Degree
    double mean = 10; // Average infection time
    double variance = 3;
    
    int nb_paths = 1;
    //double threshold=100;
    
    
     /*-----nMGA-----*/
    //generatePaths_NMGA( mean, variance, degree,nb_paths, size, engine,threshold);

  

    /*----Next Reaction on ERDOS REYNI network---*/
    
    generatePaths_next_reaction( mean, variance, degree,nb_paths, size, engine);
    
    
    /*---Mean Field---*/
    
    //simulatePaths_MeanField( mean, variance, degree,nb_paths, size, engine);

    return 0;
}

