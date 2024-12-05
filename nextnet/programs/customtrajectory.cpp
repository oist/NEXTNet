#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "nMGA.h"
#include "NextReaction.h"

using namespace std;

int program_customtrajectory(int argc, const char * argv[]) {
    rng_t engine;

    /* input filename of the adjacency list csv (the nth row represents the link to the nth node, e.g. 4,5,10,345) */ 
    string input_filename = argv[2];

    // 0 -> SI, 1 -> SIR, 2 -> SIS
    int model = atoi(argv[3]);

    /* Parameters for the epidemic*/ 
	double MEAN_INFECTION = atof(argv[4]);
    double VARIANCE_INFECTION = atof(argv[5]);
    double MEAN_RECOVERY = atof(argv[6]);
    double VARIANCE_RECOVERY = atof(argv[7]);

    // Initial percentage of infected (number from 1 to 100);
    double I0 = atof(argv[8]);

    // Maximum time in the system (only used for the SIS, else INFINITY.)
    double TMAX = atof(argv[9]);
    string output_filename = argv[10];

    transmission_time_gamma psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_gamma rho(MEAN_RECOVERY, VARIANCE_RECOVERY);

    switch (model)
    {
    case 0:
        simulate_trajectory(engine, [input_filename,psi](rng_t& engine) {
            struct {
                std::unique_ptr<network> nw;
                std::unique_ptr<simulation_algorithm> simulator;
            } env;
            env.nw.reset(new imported_network(input_filename));
            env.simulator.reset(new simulate_next_reaction(*env.nw,psi));
            return env;
            }, I0, INFINITY,output_filename);
        break;
    case 1:
        simulate_trajectory(engine, [input_filename,psi,rho](rng_t& engine) {
			struct {
				std::unique_ptr<network> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new imported_network(input_filename));
			simulate_next_reaction::params p;
			p.shuffle_neighbours = false;
			p.edges_concurrent = false;
			p.SIR = true;
            env.simulator.reset(new simulate_next_reaction(*env.nw,psi,&rho,p));
			return env;
		}, I0, INFINITY, output_filename);
        break;
    case 2:
        simulate_trajectory(engine, [input_filename,psi,rho](rng_t& engine) {
			struct {
				std::unique_ptr<network> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new imported_network(input_filename));
			simulate_next_reaction::params p;
			p.shuffle_neighbours = true;
			p.edges_concurrent = false;
			p.SIR = false;
            env.simulator.reset(new simulate_next_reaction(*env.nw,psi,&rho,p));
			return env;
		}, I0, TMAX, output_filename);
        break;
    default:
        break;
    }

    return 0;
}
