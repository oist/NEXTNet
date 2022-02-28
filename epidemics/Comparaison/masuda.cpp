/* Laplacian Gillespie Algorithm */


/*-------*\
 
/!\ Only works for distribution that are completely monotone /!\.
- Realistic infective functions do not satisfy this property.
- e.g Log-normal is not completely monotone.
- Gamma distirbution when alpha < 1 is.

\*------*/

#include "stdafx.h"
#include "random.h"
#include "masuda.h"


using namespace std;


vector<double> simulateLGA(int N,double alpha,double beta,mt19937& mersenneTwister){
    vector<double> trajectory({});
    
    
    /*-----Initialisation-----*/
    
    int i; // counter
    double time = 0;
    double dt= 0;
    vector<double> lambda(N,0.0); // rate of Poisson processes
    double sum_lambda = 0.0;
    

    /*1. Initialize each of the N processes by drawing the rate λ_i*/
    for (i=0 ; i<N ; i++) {
    lambda[i] = rgamma(alpha,beta,mersenneTwister);
    sum_lambda += lambda[i];
    }
    
    
    /*-----Simulation------*/
    
    /* 2. Draw the time of the next event*/
    double u= rand(0,1,mersenneTwister);
    dt = -log(u)/ sum_lambda;
    time += dt;
    
    /*3.  Select the process i that has generated the event. */
    double r = rand(0,1,mersenneTwister);
    double cumulant = lambda[0];
    i=0;
    while (r >= cumulant){
        cumulant += lambda[i];
        i++;
    }
    
    /* 4. Draw a new rate λ_i according to p_i[λ_i]*/
    sum_lambda -= lambda[i];
    lambda[i] = rgamma(alpha,1,mersenneTwister);
    sum_lambda += lambda[i];
    
    /* 5. Increase population size */
    double new_rate = rgamma(alpha,beta,mersenneTwister);
    lambda.push_back(new_rate);
    N++;
    sum_lambda += new_rate;
        
    return trajectory;
}


