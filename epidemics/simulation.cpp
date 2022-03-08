//
//  simulation.cpp
//  epidemics_update
//
//  Created by Cure Samuel Cyrus on 2022/02/19.
//

#include "stdafx.h"
#include "simulation.h"
#include "random.h"
#include "graph.h"
#include "types.h"

using namespace std;

vector<double> simulatePath(vector<double>& infection_times, int n_max,const lognormal_beta& infection_distribution, rng_t engine){
    
    double absolute_time = 0;
    vector<double> time_trajectory({});

    /* Start with one infected individual at time t=0*/
    for (int population = 1; population <n_max; population++ ){

        
        /*----Add new infection times----*/
        int new_infections = poissrnd(infection_distribution.r0,engine);

        for (int i =0; i < new_infections; i++)
        infection_times.push_back( infection_distribution.sample(engine) + absolute_time);  // the infection times are not relative to the age of the infected indiv. They are relative to the general time of the system.


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

void print_matrix(vector<vector<double>>& A){
    long rows = A.size();
    long cols = A[0].size();
    
    for (int i = 0; i<rows; i++) {
        for (int j =0; j<cols; j++) {
            cout << ceil(A[i][j]*10)/10 << "   ";
        }
        cout << endl;
    }
}
