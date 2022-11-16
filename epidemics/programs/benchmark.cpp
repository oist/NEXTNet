#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "nMGA.h"

using namespace std;

int program_benchmark(int argc, const char * argv[]) {
    rng_t engine;

    int method = atoi(argv[1]);
    int SIM_MAX = atoi(argv[2]);
    string filename = argv[3];


    const double MEAN_INFECTION = 10;
    const double VARIANCE_INFECTION = 1.0;
    const double MEAN_RECOVERY = 20;
    const double VARIANCE_RECOVERY = 1;
	const double R0 = 3;

    transmission_time_lognormal psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_lognormal rho(MEAN_RECOVERY, VARIANCE_RECOVERY);

    switch (method)
    {
    case 0: // Next reaction + ER graph
		measure_runtime(engine, [R0, psi, rho](rng_t& engine, int n) {
			erdos_reyni network(n, R0, engine);
			return new simulate_nmga(network,psi);
		}, SIM_MAX, INFINITY, filename);
        //measure_running_time_next_reaction_ER(engine,SIM_MAX,filename);
        break;
    case 1:
        measure_running_time_next_reaction_BA(engine,SIM_MAX,filename);
        break;
    case 2:
        measure_running_time_nMGA_ER(engine,SIM_MAX,filename);
        break;
    case 3:
        measure_running_time_nMGA_BA(engine,SIM_MAX,filename);
        break;
   
    default:
        break;
    }

    return 0;
}
