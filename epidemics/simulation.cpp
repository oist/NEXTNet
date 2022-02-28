//
//  simulation.cpp
//  epidemics_update
//
//  Created by Cure Samuel Cyrus on 2022/02/19.
//

#include "stdafx.h"
#include "simulation.h"
#include "random.h"
#include "Tau.h"

using namespace std;

vector<double> simulatePath(vector<double>& infection_times, int n_max, Tau& tau, mt19937& mersenneTwister ){
    
    double absolute_time = 0;
    vector<double> time_trajectory({});

    /* Start with one infected individual at time t=0*/
    for (int population = 1; population <n_max; population++ ){

        
        /*----Add new infection times----*/
        int new_infections = poissrnd(tau.r0,mersenneTwister);
        vector<double> new_infection_times = beta_normalised(new_infections,tau.mean,tau.variance,mersenneTwister);
        for (int i =0; i < new_infections; i++)
            infection_times.push_back(new_infection_times[i]+ absolute_time);  // the infection times are not relative to the age of the infected indiv. They are relative to the general time of the system.


        /*----Determine the next infection time in the population----*/
        if (infection_times.size()==0) { // If there are no new infections times, the simulation is over.
            break;
        }
        
        long index = min_element(infection_times.begin(),infection_times.end()) - infection_times.begin();//Find the position of minimum of infection_times
        double next_infection_time = infection_times[index]; //Find minimum
        
        /*----Update trajectory----*/
        absolute_time = next_infection_time;
        time_trajectory.push_back(absolute_time);
        
        /*----Remove that infection time from the list----*/
        infection_times.erase(infection_times.begin()+index); // Remove next_infection_time;

    }
    
    return time_trajectory;
    
}

/* For each infected individual i:
generate the nb of infections Y_i
generate the Y_i infection times.
Append them all in the infectiontime array.
*/
//vector<double> intialiseInfectionTimes(int number_of_infected, double r0,double mean, double variance){
//    
//    vector<double> infection_times({});
//    
//    for (int i=0; i < number_of_infected; i++){
//        int Y = poissrnd(r0,);
//        vector<double> times = beta_normalised(Y,mean,variance);
//        for (int j=0; j < Y; j++)
//            infection_times.push_back(times[j]);
//    }
//    return infection_times;
//}



/*
This is not yet what Boguna et al propose in their algorithm, but is a first step.
 */
vector<double> simulatePathRenewal(vector<double>& infection_times, int n_max, Tau& tau, mt19937& mersenneTwister ){
    
//    double absolute_time = 0;
    vector<double> time_trajectory({});
//
//    /* Start with one infected individual at time t=0*/
//    for (int population = 1; population <n_max; population++ ){
//
//
//
//    }
    
    return time_trajectory;
}
