#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "NextReactionMeanField.h"

using namespace std;

int program_sis_meanfield_gamma_gamma(int argc, const char * argv[])
{
    rng_t engine;

    if (argc < 9)
        throw logic_error("Not enough arguments");
    
    int N = atoi(argv[2]);
    double R0 = atof(argv[3]);
    double MEAN = atof(argv[4]);
    double VARIANCE = atof(argv[5]);
    double MEAN_rho = atof(argv[6]);
    double VARIANCE_rho = atof(argv[7]);
    double TMAX = atof(argv[8]);
    string filename("trajectory.csv");
    
    if (argc == 10)
        filename = string(argv[9]);
    
    transmission_time_gamma psi(MEAN, VARIANCE);
    transmission_time_gamma rho(MEAN_rho, VARIANCE_rho);

    /* Simulate using next reaction once times */
    std::vector<double> times, infected;
    simulate_next_reaction_mean_field simulate(N, R0, psi,&rho);
    
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
