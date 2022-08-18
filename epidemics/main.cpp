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

    int n = 1000;
    if(argc == 2){
        n = atoi(argv[1]);
    }
    

    
    int degree = 3;
    double mean = 10;
    double variance = 1.0;
    
   
    return 0;
}

