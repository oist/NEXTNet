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
        vector<double> new_infection_times = beta_normalised(new_infections,tau,mersenneTwister);
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

vector<double> simulatePathNetwork(int network_size,double degree, Tau& tau, mt19937& mersenneTwister ){
    
//    double absolute_time = 0;
    
    vector<vector<double>> A(network_size, vector<double>(network_size)); // Initialise Adjacency Matrix
    
    vector<vector<int>> neighbours(network_size, vector<int>({}) ); // list of first neighbours to vertex i.
    vector<vector<int>> infection_times(network_size, vector<int>({}) );
    
    vector<double> r =rand( network_size * network_size ,mersenneTwister);// pr. of having a link
    //vector<double> infection_times = beta_normalised(network_size * network_size,tau,mersenneTwister); // infection time between nodes
    
    
    /*--------------Initialisation--------------

     Construct erdos-Reyni graph and for each link we add an infection time:
     => A[i][j] is the time node i (or j) takes to infect node j (or i) once it has itself been infected by someone else.
     
     */
    
    double p = degree/network_size; // probability of an edge: if network_size ->infty and degree-> fixed then we get Poisson Graph.
    
    for (int i=0; i<network_size; i++) {
        for (int j=0; j<i; j++) {
            if (r[i*network_size+j]>p) {
                A[i][j]=0;
                A[j][i]=0;
            }
            else{
                A[i][j]=beta_normalised(tau, mersenneTwister);
                A[j][i]=A[i][j];
                neighbours[i].push_back(j);
                neighbours[j].push_back(i);
                
                infection_times[i].push_back(A[i][j]);
                infection_times[j].push_back(A[i][j]);
            }
        }
        
    }
    
    cout << "Graph generated" << endl;
    
    
    int number_of_infected = 1;
    vector<int> infected_nodes({0}); // Start with vertex 0 be the initial infected, without loss of generality.
    
    while (number_of_infected < network_size) {
        
    }
    cout << "Entire network has been infected" << endl;
    return r;
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
