#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "nMGA.h"
#include "NextReaction.h"

using namespace std;
// Generate the average trajectory of an SI epidemic on a network ensemble (erdos renyi or barabasi albert)

int program_average(int argc, const char * argv[]) {
    rng_t engine;

    engine.seed(3);


    int size = atoi(argv[2]);

    //Erdos renyi graph with mean degree R0. If R0==0 then a Barabasi albert network is generated
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

    if (R0 == 0.0){
        std::uniform_int_distribution<> dis(0, size-1);
        for (int s = 0; s < nb_simulation; s++)
            {
                cout << s << "\n";
                scale_free network(size,engine);
                simulate_next_reaction simulation(network, psi);
                // simulation.add_infections({ std::make_pair(5822, 0.0)});
                // simulation.add_infections({ std::make_pair(7712, 0.0)});
                // simulation.add_infections({ std::make_pair(3086, 0.0)});
                // simulation.add_infections({ std::make_pair(3954, 0.0)});
                // simulation.add_infections({ std::make_pair(1604, 0.0)});
                for (int node = 0; node < I0; node++)
                {
                    node_t n = dis(engine);
                    simulation.add_infections({ std::make_pair(n, 0.0)});
                    // cout << n << endl;
                }
                while (true) {
                    auto point = simulation.step(engine);
                    if (!point )
                        break;
                    all_times.push_back(point->time);
                    
                }

            }
            // cout << all_times.size();
            sort(all_times.begin(), all_times.end());

            vector<double> trim;
            for (int ax = 0; ax < all_times.size(); ax += nb_simulation){
                trim.push_back(all_times[ax]);
            }
        exportData(trim, output_filename);
    } else {
        for (int s = 0; s < nb_simulation; s++)
            {
                cout << s << "\n";
                erdos_reyni network(size,R0,engine);
                simulate_next_reaction simulation(network, psi);

                for (int node = 0; node < I0; node++)
                {
                    simulation.add_infections({ std::make_pair(node, 0.0)});
                }

                while (true) {
                    auto point = simulation.step(engine);
                    if (!point )
                        break;
                    all_times.push_back(point->time);
                    
                }

            }
            cout << all_times.size();
            sort(all_times.begin(), all_times.end());

            vector<double> trim;
            for (int ax = 0; ax < all_times.size(); ax += nb_simulation){
                trim.push_back(all_times[ax]);
            }
        exportData(trim, output_filename);
    }

    return 0;
}
