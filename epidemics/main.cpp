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
#include "nMGA.h"
#include "NextReaction.h"

using namespace std;

int main(int argc, const char * argv[]) {
    rng_t engine;
    
    
    
    
    

    
//    if (argc!=8 || argc!=9)
//        throw logic_error("Not enough arguments");
    
    int N = atoi(argv[1]);
    double R0 = atof(argv[2]);
    double MEAN = atof(argv[3]);
    double VARIANCE = atof(argv[4]);
    double MEAN_rho = atof(argv[5]);
    double VARIANCE_rho = atof(argv[6]);
    double TMAX = atof(argv[7]);
    string filename("trajectory.csv");
    
    if (argc==9)
        filename = string(argv[8]);
    
    
    
    transmission_time_gamma psi(MEAN, VARIANCE);
    transmission_time_gamma rho(MEAN_rho, VARIANCE_rho);

    /* Simulate using next reaction once times */
    std::vector<double> times, infected;
    erdos_reyni network(2, 1,engine);
    simulate_next_reaction_mean_field simulate(network, N, R0, psi,&rho);
    
    simulate.add_infections({ std::make_pair(0, 0.0)});
    double current_infected = 0;
    // Run simulation, collect transmission times
    while (true) {

        auto point = simulate.step(engine);
        if (!point || (point -> time > TMAX))
            break;

        switch (point-> kind) {
            case event_kind::infection:
            case event_kind::outside_infection:
                current_infected+=1;
                break;
            case event_kind::reset:
                current_infected-=1;
                break;
            default:
                throw std::logic_error("unexpected event kind");
        }

        times.push_back(point->time);
        infected.push_back(current_infected/N);
    }
    
    exportData(times,infected,filename);
    
    
    return 0;
}

