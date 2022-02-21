//
//  simulation.cpp
//  epidemics_update
//
//  Created by Cure Samuel Cyrus on 2022/02/19.
//

#include "stdafx.h"
#include "simulation.h"
#include "random.h"

using namespace std;

vector<double> simulatePath(vector<double>& infection_times, int n_max, double& absolute_time, double mean, double variance,double R0){
    
    vector<double> time_trajectory({0});
    
    for (int population = 1; population <n_max; population++ ){
        

        /* Determine the next infection time in the population */
        long index = min_element(infection_times.begin(),infection_times.end()) - infection_times.begin(); //Find the position of minimum of infection_times
        double next_infection_time = infection_times[index]; //Find minimum
        
        /* Update trajectory*/
        absolute_time += next_infection_time;
        time_trajectory.push_back(absolute_time);
        
        /* Remove that infection time from the list */
        infection_times.erase(infection_times.begin()+index); // Remove next_infection_time;
        
        /* Advance the time of the system by updating the entire list*/
        for (int i =0; i< (int) infection_times.size() ; i++)
            infection_times[i] -= next_infection_time;
        
        
        /* Add new infection times*/
        int new_infections = poissrnd(R0);
        vector<double> new_times = beta_normalised(new_infections,mean,variance);
        for (int i =0; i < new_infections; i++)
            infection_times.push_back(new_times[i]);
        
    }
    
    return time_trajectory;
    
}

/* For each infected individual i:
generate the nb of infections Y_i
generate the Y_i infection times.
Append them all in the infectiontime array.
*/
vector<double> intialiseInfectionTimes(int number_of_infected, double r0,double mean, double variance){
    
    vector<double> infection_times({});
    
    for (int i=0; i < number_of_infected; i++){
        int Y = poissrnd(r0);
        vector<double> times = beta_normalised(Y,mean,variance);
        for (int j=0; j < Y; j++)
            infection_times.push_back(times[j]);
    }
    return infection_times;
}
