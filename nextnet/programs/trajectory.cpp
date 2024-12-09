#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "nMGA.h"
#include "NextReaction.h"

using namespace std;

int program_trajectory(int argc, const char * argv[]) {
    rng_t engine;
    
    /* The first argument is the epidemic type
    * "SIS","SI","SIR"
    * */ 
    string method = argv[2];

    double R0 = atof(argv[3]);
    // if (nw != "MF" && nw != "BA")
    //     R0 = atof(argv[3]);
    /* Parameters for the epidemic*/ 
    
	double MEAN_INFECTION = atof(argv[4]);
    double VARIANCE_INFECTION = atof(argv[5]);
    double MEAN_RECOVERY = atof(argv[6]);
    double VARIANCE_RECOVERY = atof(argv[7]);
	int n = atoi(argv[8]);
    double I0 = atof(argv[9]);
    double TMAX = atof(argv[10]);

    string filename = argv[11];

    transmission_time_lognormal psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_lognormal rho(MEAN_RECOVERY, VARIANCE_RECOVERY);

    if (method == "SI"){
        simulate_trajectory(engine, [R0, psi,n](rng_t& engine) {
			struct {
				std::unique_ptr<network> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new erdos_reyni(n, R0, engine));
			simulate_nmga::params p;
            env.simulator.reset(new simulate_nmga(*env.nw, psi, nullptr, p));
			return env;
		}, I0, INFINITY, filename);
    } else if (method == "SIS")
    {
        simulate_trajectory(engine, [R0,psi,rho,n](rng_t& engine) {
			struct {
				std::unique_ptr<network> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new erdos_reyni(n, R0, engine));
			simulate_nmga::params p;
			env.simulator.reset(new simulate_nmga(*env.nw, psi, &rho, p));
			return env;
		}, I0, TMAX, filename);
    } else if (method == "SIR"){
        simulate_trajectory(engine, [R0,psi,rho,n](rng_t& engine) {
			struct {
				std::unique_ptr<network> nw;
				std::unique_ptr<simulation_algorithm> simulator;
			} env;
			env.nw.reset(new erdos_reyni(n, R0, engine));
			simulate_nmga::params p;
			p.SIR = true;
			env.simulator.reset(new simulate_nmga(*env.nw, psi, &rho, p));
            return env;

		}, I0, INFINITY, filename);
    } else {
        throw std::logic_error("invalid simulation type, enter SI,SIS,SIR as an argument.");
    }

    return 0;
}
