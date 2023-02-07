#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "nMGA.h"
#include "NextReaction.h"

using namespace std;

int program_average(int argc, const char * argv[]) {
    rng_t engine;

    int size = atoi(argv[2]);

    //Erdos renyi graph with mean degree R0. If R0 a Barabasi albert network is generated
    double R0 = atof(argv[3]);

    /* Parameters for the epidemic*/ 
	double MEAN_INFECTION = atof(argv[4]);
    double VARIANCE_INFECTION = atof(argv[5]);

    // Initial percentage of infected (float from 0 to 100);
    // double I0 = atof(argv[6]);
    int nb_simulation = atoi(argv[7]);
    string output_filename = argv[8];
    

    transmission_time_gamma psi(MEAN_INFECTION, VARIANCE_INFECTION);

    vector<double> all_times;
    int X = 0;
    for (int s = 0; s < nb_simulation; s++)
    {
        cout << s << "\n";
        erdos_reyni network(size,R0,engine);
        simulate_next_reaction simulation(network, psi);

        for (int node = 0; node < 10; node++)
        {
            simulation.add_infections({ std::make_pair(node, 0.0)});
        }

        while (true) {
            auto point = simulation.step(engine);
            if (!point )
                break;
            if (X % nb_simulation == 0){
                all_times.push_back(point->time);
            }
            X +=1;
        }

    }
    cout << all_times.size();
    sort(all_times.begin(), all_times.end());

    // vector<double> trim;
    // for (int ax = 0; ax < all_times.size(); ax += nb_simulation){
    //     trim.push_back(all_times[ax]);
    // }

    exportData(all_times, output_filename);

    // simulate_trajectory(engine, [input_filename,psi,R0](rng_t& engine) {
    //         struct {
    //             std::unique_ptr<graph> nw;
    //             std::unique_ptr<simulation_algorithm> simulator;
    //         } env;
    //         if (R0==0.0){
	// 		    env.nw.reset(new scale_free(n, engine));
    //         } else {
    //             switch (model)
    //             {
    //             case 0:
    //                 env.simulator.reset(new simulate_next_reaction(*env.nw, psi));
    //                 break;
    //             case 1:
    //                 env.simulator.reset(new simulate_next_reaction(*env.nw, psi,&rho,true,false,SIR));
    //                 break;
    //             case 3:
    //                 env.simulator.reset(new simulate_next_reaction(*env.nw, psi,&rho,true,false,SIR));
    //                 break;
                
    //             default:  
    //                 break;
    //             }
    //         }
    //         return env;
    //         }, I0,TMAX,output_filename);
    // return 0;
    return 0;
}
