#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "nMGA.h"
#include "NextReaction.h"

using namespace std;

int program_benchmark(int argc, const char * argv[]) {
    rng_t engine;

	bool edges_concurrent;
    int method = atoi(argv[2]);
    int SIM_MAX = atoi(argv[3]);
	string concurrent = argv[4];
    string filename = argv[5];

	if (concurrent == "0"){
		edges_concurrent = false;
	} else {
		edges_concurrent = true;
	}

	const double MEAN_INFECTION = 10;
    const double VARIANCE_INFECTION = 1.0;
    const double MEAN_RECOVERY = 20;
    const double VARIANCE_RECOVERY = 1;
	const double R0 = 3;
	const bool shuffle_neighbours = false;

    transmission_time_lognormal psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_lognormal rho(MEAN_RECOVERY, VARIANCE_RECOVERY);

    switch (method)
    {
    case 0: // Next reaction + ER graph
		measure_runtime(engine, [R0, psi, rho,shuffle_neighbours,edges_concurrent](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new erdos_reyni(n, R0, engine));
			env.simulator.reset(new simulate_next_reaction(*env.nw, psi,nullptr,shuffle_neighbours,edges_concurrent));
			return env;
		}, SIM_MAX, INFINITY, filename);
        //measure_running_time_next_reaction_ER(engine,SIM_MAX,filename);
        break;
    case 1: // Next reaction + BA graph
        measure_runtime(engine, [psi, rho,shuffle_neighbours,edges_concurrent](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new scale_free(n, engine));
			env.simulator.reset(new simulate_next_reaction(*env.nw, psi,nullptr,shuffle_neighbours,edges_concurrent));
			return env;
		}, SIM_MAX, INFINITY, filename);
        break;
    case 2: // NMGA + ER graph
        measure_runtime(engine, [R0, psi, rho](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new erdos_reyni(n, R0, engine));
			env.simulator.reset(new simulate_nmga(*env.nw, psi));
			return env;
		}, SIM_MAX, INFINITY, filename);
        break;
    case 3: // NMGA + BA graph
        measure_runtime(engine, [psi, rho](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new scale_free(n, engine));
			env.simulator.reset(new simulate_nmga(*env.nw, psi));
			return env;
		}, SIM_MAX, INFINITY, filename);
        break;
   
    default:
        break;
    }

    return 0;
}
