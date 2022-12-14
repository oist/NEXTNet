#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "nMGA.h"
#include "NextReaction.h"

using namespace std;

int program_benchmark(int argc, const char * argv[]) {
    rng_t engine;

    int method = atoi(argv[2]);
    int SIM_MAX = atoi(argv[3]);
	int SIR_ = atoi(argv[4]);
	double TMAX = atof(argv[5]);
    string filename = argv[6];
	
	bool SIR;
	bool shuffle_neighbours;

	if (SIR_){
		SIR = true;
		shuffle_neighbours = false;
		TMAX = INFINITY;
	} else {
		shuffle_neighbours = true;
		SIR = false;
	}


	const double MEAN_INFECTION = 3;
    const double VARIANCE_INFECTION = 10.0;
    const double MEAN_RECOVERY = 8;
    const double VARIANCE_RECOVERY = 5;
	const double R0 = 3;

    transmission_time_lognormal psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_lognormal rho(MEAN_RECOVERY, VARIANCE_RECOVERY);
    switch (method)
    {
    case 0: // Next reaction + ER graph
		measure_runtime(engine, [R0, psi, rho,shuffle_neighbours,SIR,TMAX](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new erdos_reyni(n, R0, engine));
			env.simulator.reset(new simulate_next_reaction(*env.nw, psi,&rho,shuffle_neighbours,false, SIR));
			return env;
		}, SIM_MAX, TMAX,14, filename);
        //measure_running_time_next_reaction_ER(engine,SIM_MAX,filename);
        break;
    case 1: // Next reaction + BA graph
        measure_runtime(engine, [psi, rho,shuffle_neighbours,SIR,TMAX](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new scale_free(n, engine));
			env.simulator.reset(new simulate_next_reaction(*env.nw, psi,&rho,shuffle_neighbours,false,SIR));
			return env;
		}, SIM_MAX, TMAX,24,filename);
        break;
    case 2: // NMGA + ER graph
        measure_runtime(engine, [R0, psi, rho,SIR,shuffle_neighbours,TMAX](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new erdos_reyni(n, R0, engine));
			env.simulator.reset(new simulate_nmga(*env.nw, psi,&rho,SIR,shuffle_neighbours));
			return env;
		}, SIM_MAX, TMAX,14, filename);
        break;
    case 3: // NMGA + BA graph
        measure_runtime(engine, [psi, rho,SIR,shuffle_neighbours,TMAX](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new scale_free(n, engine));
			env.simulator.reset(new simulate_nmga(*env.nw, psi, &rho, SIR, shuffle_neighbours));
			return env;
		}, SIM_MAX, TMAX,14, filename);
        break;
   
    default:
        break;
    }

    return 0;
}
