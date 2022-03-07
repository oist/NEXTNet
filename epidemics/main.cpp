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
#include "Tau.h"
#include "boguna.h"

using namespace std;

/*  MersenneTwister random number generator */
/*static*/ mt19937 mersenneTwister;


int main(int argc, const char * argv[]) {
    
    setDebugEnabled(false);// testing or not.


    /*--------- Initialise Parameters ----------*/

    Tau tau(1,1,2); //Mean Variance R0
    
    int size = 10;
    double degree = 3;
    simulatePathNetwork(size,degree, tau, mersenneTwister);
    
    return 0;
}
