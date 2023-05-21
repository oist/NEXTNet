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

    // Network size, e.g. 6 +> 1e6;
    double size = pow(10,stod(argv[3]));

    /* Parameters for the epidemic*/ 
	double MEAN_INFECTION = atof(argv[4]);
    double VARIANCE_INFECTION = atof(argv[5]);

    int nb_simulation = atoi(argv[6]);
  

    transmission_time_gamma psi(MEAN_INFECTION, VARIANCE_INFECTION);

    vector<double> all_times;
    
    cout << "generating network...\n";


    scale_free network(size,engine);
    
    string filename = "trajectory.dat";

    const double r = assortativity(network);

    cout << " r: " << r << "\n";

    double k1 = 0;
    double k2 = 0;
    double k3 = 0;
    double k4 = 0;

    for (int i = 0; i < size; i++ ){
        double k = network.outdegree(i) * 1.0;
        k1 += k ;
        k2 +=  pow(k,2) ;
        k3 +=  pow(k,3) ;
        k4 +=  pow(k,4); 
    }

    const double mu = k2/k1 - 1;
    const double mu_r1 = k2/k1 - 1 + r*(k3/k2 -1);
    const double mu_r2 = k2/k1*(1-r) - 1 + r*(k3/k2 -1);
    const double a = MEAN_INFECTION*MEAN_INFECTION / VARIANCE_INFECTION;
    const double b = VARIANCE_INFECTION / MEAN_INFECTION;
    const double Lambda = (pow(mu,1/a) - 1) / b;
    const double Lambda_r1 = (pow(mu_r1,1/a) - 1) / b;
    const double Lambda_r2 = (pow(mu_r2,1/a) - 1) / b;
    

    string parameters = "parameters.dat";

    parameters += to_string(seed_number) + ".dat";
    
    ofstream out;
    out.open(parameters);
    
    out << "L " << Lambda << "\n";
    out << "L_r1 " << Lambda_r1 << "\n";
    out << "L_r2 " << Lambda_r2 << "\n";
    out << "r " << r << "\n";
    out << "k1 " << k1/size << "\n";
    out << "k2 " << k2/size << "\n";
    out << "k3 " << k3/size << "\n";
    out << "k4 " << k4/size << "\n";
    
    out.close();


    std::uniform_int_distribution<> dis(0, size-1);
    std::uniform_real_distribution<> ran(0.0,1.0);
    const bool SHUFFLE_NEIGHBOURS = false;
    const bool EDGES_CONCURRENT = true;
    // run simulation
    cout << "running simulation...\n";
    for (int s = 0; s < nb_simulation; s++)
        {

            cout <<  s <<" / " << nb_simulation << "\r";
            std::cout.flush();

            simulate_next_reaction simulation(network, psi,nullptr,SHUFFLE_NEIGHBOURS,EDGES_CONCURRENT,false);

            for (node_t i = 0; i < 1; i++)
            {
                node_t rand_node = dis(engine);
                simulation.add_infections({ std::make_pair(rand_node, 0)});
            }
            

                            
            while (true) {
                auto point = simulation.step(engine);
                if (!point )
                    break;
                all_times.push_back(point->time);
                
            }
        }

        
    cout << "sort...." << endl;
    sort(all_times.begin(), all_times.end());
    
    cout << "export..." << endl;

    vector<double> trimmed_trajectory;
    vector<double> trim_Y;
    for (int ax = 0; ax < (int) all_times.size(); ax += nb_simulation){
            trimmed_trajectory.push_back(all_times[ax]);
    }

    exportData(trimmed_trajectory, filename);
    return 0;
}