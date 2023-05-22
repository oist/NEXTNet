#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "nMGA.h"
#include "NextReaction.h"

using namespace std;

int program_benchmark(int argc, const char * argv[]) {
    rng_t engine;

    std::string method = argv[2];
	std::string ensemble = argv[3];
	int MAX_POWER = atoi(argv[4]);
    int SIM_MAX = atoi(argv[5]);
	bool SIR = (bool) atoi(argv[6]);
	double TMAX = atof(argv[7]);
	double MEAN_INFECTION = atof(argv[8]);
    double VARIANCE_INFECTION = atof(argv[9]);
    double MEAN_RECOVERY = atof(argv[10]);
    double VARIANCE_RECOVERY = atof(argv[11]);
	double R0 = atof(argv[12]);
	int M = atoi(argv[13]);
    bool CONCURRENT_EDGES = (bool) atoi(argv[14]);
	
	bool SHUFFLE_NEIGHBOURS;
	std::string filename = ensemble;
	if (ensemble == "ER"){
		filename += "_" + std::to_string((int) R0);
	} else if (ensemble == "BA"){
		filename += "_" + std::to_string((int) M);
	} else {
		throw std::logic_error("unrecognised ensemble name");
	}

	filename += "_" + method + "_" + std::to_string(MAX_POWER) + "_" + std::to_string(SIM_MAX) + "_" + std::to_string((int) MEAN_INFECTION) + "_" + std::to_string((int) VARIANCE_INFECTION) + "_" + std::to_string((int) MEAN_RECOVERY) + "_" + std::to_string((int) VARIANCE_RECOVERY) + "_";

	if (CONCURRENT_EDGES){
		filename += "CON";
	} else{
		filename += "SEQ";
	}

	if (SIR || CONCURRENT_EDGES){
		SHUFFLE_NEIGHBOURS = false;
		TMAX = INFINITY;
	} else {
		SHUFFLE_NEIGHBOURS = true;
	}

	// const double MEAN_INFECTION = 3;
    // const double VARIANCE_INFECTION = 10.0;
    // const double MEAN_RECOVERY = 8;
    // const double VARIANCE_RECOVERY = 5;
	// const double R0 = 3;

    transmission_time_lognormal psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_lognormal rho(MEAN_RECOVERY, VARIANCE_RECOVERY);

	if ((method == "NR") && (ensemble == "ER")){
		measure_runtime(engine, [R0, psi, rho,SHUFFLE_NEIGHBOURS,SIR,TMAX,CONCURRENT_EDGES](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new erdos_reyni(n, R0, engine));
			env.simulator.reset(new simulate_next_reaction(*env.nw, psi,&rho,SHUFFLE_NEIGHBOURS,CONCURRENT_EDGES, SIR));
			return env;
		}, SIM_MAX, TMAX,MAX_POWER, filename);

    } else if ((method == "NR") && (ensemble == "BA")){
        measure_runtime(engine, [psi, rho,SHUFFLE_NEIGHBOURS,SIR,TMAX,M,CONCURRENT_EDGES](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new scale_free(n, engine,M));
			env.simulator.reset(new simulate_next_reaction(*env.nw, psi,&rho,SHUFFLE_NEIGHBOURS,CONCURRENT_EDGES,SIR));
			return env;
		}, SIM_MAX, TMAX,MAX_POWER,filename);

	} else if ((method == "NMGA") && (ensemble == "ER")){
        measure_runtime(engine, [R0, psi, rho,SIR,SHUFFLE_NEIGHBOURS,TMAX](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new erdos_reyni(n, R0, engine));
			env.simulator.reset(new simulate_nmga(*env.nw, psi,&rho,SIR,SHUFFLE_NEIGHBOURS));
			return env;
		}, SIM_MAX, TMAX,MAX_POWER, filename);

	} else if ((method == "NMGA") && (ensemble == "BA")){
        measure_runtime(engine, [psi, rho,SIR,SHUFFLE_NEIGHBOURS,TMAX,M](rng_t& engine, int n) {
			struct {
				std::unique_ptr<graph> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new scale_free(n, engine,M));
			env.simulator.reset(new simulate_nmga(*env.nw, psi, &rho, SIR, SHUFFLE_NEIGHBOURS));
			return env;
		}, SIM_MAX, TMAX,MAX_POWER, filename);
	} else {
		return -1;
    }

    return 0;
}
