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


    double absolute_time= 0;
    int number_of_infected = 1;
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
    
    string filename("data.csv");
    rng(filename);
    for (int path = 0; path < number_of_paths ; path++ ){
        
        absolute_time = 0;
        /*--------- Initialise Infection times of the starting population ----------*/

        vector<double> infection_times = intialiseInfectionTimes(number_of_infected, r0, mean, variance);

        /*--------- Begin the simulation ----------*/

        vector<double> time_trajectory = simulatePath(infection_times,n_max,absolute_time,mean,variance,r0);
        cout << "\n";
//        cout << "SIMULATION ENDED: Reached " << n_max << " at time t = " << absolute_time;
//        cout <<" t = " << absolute_time;
        cout << path;
    
        /*--------- Save Data ----------*/
        string s = to_string(path+1);
        exportData(time_trajectory,s+filename);
    }

//        vector<double>  v= poissrnd(7.5,10000);
//        double s = 0;
//        for (int i = 0; i < 1000; i++){
//            s+=v[i];
//        }
//        cout << r0 << "\n" << s/1000<<"\n";
    return 0;
}
