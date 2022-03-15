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
    
    
    int size = 5000; // Size network or max pop for mean field.
    int degree = 4; // Average Degree
    double mean = 10; // Average infection time
    double variance = 2;
    
    int nb_paths = 1;
    
    
     /*-----nMGA-----*/

    //generatePaths_next_reaction( men, variance, degree,nb_paths, size, engine, int threshold= 10000);
  

    /*----Next Reaction on ERDOS REYNI network---*/
    
    generatePaths_next_reaction( mean, variance, degree,nb_paths, size, engine);
    
    
    /*---Mean Field---*/
    
    //simulatePaths_MeanField( mean, variance, degree,nb_paths, size, engine);

    return 0;
}

