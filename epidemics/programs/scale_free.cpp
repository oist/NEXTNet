#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "nMGA.h"
#include "NextReaction.h"

using namespace std;
// Generate the average trajectory of an SI epidemic on a scale-free network

int program_scale_free(int argc, const char * argv[]) {
    rng_t engine;


    int seed_number = atoi(argv[2]);

    engine.seed(seed_number);

    // Network size, e.g. 1e6;
    int size = atoi(argv[3]);

    /* Parameters for the epidemic*/ 
	double MEAN_INFECTION = atof(argv[4]);
    double VARIANCE_INFECTION = atof(argv[5]);

    int nb_simulation = atoi(argv[6]);
    string output_filename = argv[7];
    

    transmission_time_gamma psi(MEAN_INFECTION, VARIANCE_INFECTION);

    vector<double> all_times;
    
    cout << "generating network...\n";

    // scale_free network(size,engine);
    // imported_network network(string("~/Desktop/adjalist.csv"));
    erdos_reyni network(size,3,engine);    
    add_correlation(0.1,network,engine);
    double r = assortativity(network);

    double k1 = 0;
    double k2 = 0;
    double k3 = 0;
    double k4 = 0;

    for (int i = 0; i < size; i++ ){
        double k = network.outdegree(i) * 1.0;
        k1 += k ;
        k2 +=  pow(k,2) ;
        k3 +=  pow(k,3) ;
        k4 += pow(k,4); 
    }

    // Export parameters
    string parameters ("parameters.dat");
    ofstream out;
    out.open(parameters);
    
    out << "r " << r << "\n";
    out << "k1 " << k1/size << "\n";
    out << "k2 " << k2 /size << "\n";
    out << "k3 " << k3/size << "\n";
    out << "k4 " << k4/size << "\n";
    
    out.close();


    std::uniform_int_distribution<> dis(0, size-1);

    // run simulation
    cout << "running simulation...\n";
    for (int s = 0; s < nb_simulation; s++)
        {

            cout <<  s <<' / ' << nb_simulation  << " %" << "\r";
            std::cout.flush();
            
            simulate_next_reaction simulation(network, psi);

            node_t rand_node = dis(engine);
            simulation.add_infections({ std::make_pair(rand_node, 0.0)});
                            
            while (true) {
                auto point = simulation.step(engine);
                if (!point )
                    break;
                all_times.push_back(point->time);
                
            }
        }

        
    // cout << all_times.size();
    sort(all_times.begin(), all_times.end());
    cout << "export" << endl;
    vector<double> trim;
    for (int ax = 0; ax < all_times.size(); ax += nb_simulation){
            trim.push_back(all_times[ax]);
    }
    exportData(trim, output_filename);

    return 0;
}