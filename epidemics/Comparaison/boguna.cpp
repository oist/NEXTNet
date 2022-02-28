/*----------------------------------------------------------------------------------------*/
/* non-Markovian Gillespie Algorithm  */

/* DOI ORIGINAL PAPER: https://doi.org/10.1103/PhysRevE.90.042108  */

/* QUOTE FROM THE PAPER:
 
 " we consider that each active link (connecting a susceptible-infected pair) defines a statistically independent infection process following the distribution psi(tau).
 That is, a susceptible individual connected to a single infected individual will become infected after a random time distributed by psi(tau) from the moment the link became active.
 If the susceptible individual is connected to more than one infected neighbor, each active link is considered as statistically independent so that the infection event will take place at the time of the first firing event of any of the current active links." */


/* CURRENT IMPLEMENTATION:
 - Consider a mean-field degree regime: Each individual has Poi(R0) neighbours and no infected nodes share the same neighbour.
 */


/*----------------------------------------------------------------------------------------*/
#include "stdafx.h"
#include "random.h"
#include "boguna.h"
#include "Tau.h"

using namespace std;



/* Statistically exact version - Not in their paper because of step 2.*/
vector<double> simulateNMGA_exact(int n_max,Tau& tau,mt19937& mersenneTwister){
    vector<double> trajectory({});
//    double absolute_time = 0;
//    vector<double> time_of_activation({}); // t_i the time elapsed since the last event of the ith process.
//
//    /* Start with one infected individual at time t=0*/
//    for (int population = 1; population <n_max; population++ ){
//
//    /*----Add new infection times----*/
//    int new_infections = poissrnd(tau.r0,mersenneTwister);
//    vector<double> new_infection_times = beta_normalised(new_infections,tau.mean,tau.variance,mersenneTwister);
//    for (int i =0; i < new_infections; i++)
//        time_of_activation.push_back(new_infection_times[i]+absolute_time);
//    }
//
//    /* Draw the time until the next event, Δt, by solving Φ(Δt|{tj}) = u*/
//    double u = rand(0,1,mersenneTwister);
    
    
    
    
 


    
    return trajectory;
    
}

/* Only statistically exact version in the limit n-> \infty - version in their paper.*/
vector<double> simulateNMGA(int n_max,Tau& tau,mt19937& mersenneTwister){
    vector<double> trajectory({});
    return trajectory;
}
