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

using namespace std;

/* Infection times of an individual */
double beta_normalised( double mean, double variance);

int main(int argc, const char * argv[]) {

    /*--------- Initialise Parameters ----------*/


//    double absolute_time= 0;
//    int number_of_infected = 1;
    int number_of_paths = 100;
    

    int n_max =3000;

    double mean = 5; // Mean time of secondary infections.
    double variance = 3; // Variance "" "" ""
    double r0 = 3.0; // Average number of secondary infection.
    cout << argc;
    if (argc > 1) {
       n_max = atoi(argv[1]); // Optional user's choice of n_max
    }
    
    // if (argc == 5) {
    //     mean = 1.0 * atoi(argv[2]);
    //     variance = 1.0* atoi(argv[3]); // User's choice
    //     r0 = 1.0 * atoi(argv[4]);
    // }
    
    string filename("path");
    string csv(".csv");
    
    for (int path = 0; path < number_of_paths ; path++ ){
        
        /*--------- Initialise Infection times of the starting population ----------*/
        vector<double> infection_times ={0}; // First infected occurs at t=0
        
        // Not needed for now:
        //vector<double> infection_times = intialiseInfectionTimes(number_of_infected, r0, mean, variance);

        /*--------- Begin the simulation ----------*/
        vector<double> time_trajectory = simulatePath(infection_times, n_max, mean, variance,r0);
        //cout << path;

        /*--------- Save Data ----------*/
        string s = to_string(path+1);
        exportData(time_trajectory,filename+s+csv);
    }

    return 0;
}
