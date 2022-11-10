#include "tests/stdafx.h"
#include "tests/simulate.h"
#include "tests/analytical.h"
#include "tests/plot.h"

#include "random.h"
#include "NextReaction.h"
#include "nMGA.h"

#if ENABLE_PLOTTING
TEST_CASE("Plot large-population SIR mean-field (NextReaction)", "[nextreaction]") {
    using namespace std::string_literals;
    std::mt19937 engine;

    const int M = 1000;
    const int Nfully = 1000;
    const int Nerdos = 1000;
    const double T = 25;
    const int X = 400;
    const double R0 = 3;
    // Every infecteced node has N-1 neighbours, of which in the large-population limit N-2 are susceptible.
    // To trigger subsequent infections amgonst these N-2 susceptible neighbours, we must infect each neighbour
    // with probability R0/(N-2). Note that for the first infected node, this is not strictly speaking correct,
    // since it's infection has no source, and it will thus create R0(N-1)/(N-2) > R0 subsequent infections.
    // This only causes a relative error of |1 - (N-1)/(N-2)| =~= 1/N though.
    const double pfully = R0/(Nfully-2);
    const double MEAN = 10;
    const double VARIANCE = 1;
    transmission_time_gamma psifully(MEAN, VARIANCE, 1.0-pfully);
    transmission_time_gamma psi(MEAN, VARIANCE);

    /* Simulate using next reaction M times */
    auto rfully = simulate_SIR<fully_connected, simulate_next_reaction>(engine, psifully, T, M, 1, Nfully);
    auto racyclic = simulate_SIR<acyclic, simulate_next_reaction>(engine, psi, T, M, 1, R0+1, true);
    auto rerdos = simulate_SIR<erdos_reyni, simulate_next_reaction>(engine, psi, T, M, 1, Nerdos, R0);

    /* Evaluate analytical solution */
    std::pair<std::vector<double>, std::vector<double>> analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for(std::size_t i=0; i < X; ++i) {
        const double t = (double)T * i / (X-1);
        analytical.first.push_back(t);
        analytical.second.push_back(sol.N(t));
    }

    plot("nextreaction.sir.mean.pdf", "Large-population SIR mean-field [NextReaction]", [&](auto& gp, auto& p) {
        p.add_plot1d(rfully, "with lines title 'next reaction fully-connected (N="s + std::to_string(Nfully) + ")'"s);
        p.add_plot1d(racyclic, "with lines title 'next reaction acyclic'"s);
        p.add_plot1d(rerdos, "with lines title 'next reaction Erdös-Reyni (N="s + std::to_string(Nerdos) + ")'"s);
        p.add_plot1d(analytical, "with lines title 'analytical'");
    });
}
#endif

#if ENABLE_PLOTTING
TEST_CASE("Plot SIS single trajectory (NextReaction)", "[nextreaction]") {
       using namespace std::string_literals;
    std::mt19937 engine;

    const int N = 10000;
    const double T = 200;
    const double R0 = 3;

    const double MEAN = 3;
    const double VARIANCE = 1;
    const double MEAN_rho = 10;
    const double VARIANCE_rho = 1;
    transmission_time_gamma psi(MEAN, VARIANCE);
    transmission_time_gamma rho(MEAN_rho, VARIANCE_rho);

    /* Simulate using next reaction once times */
    std::vector<double> t_sim, y_sim_new, y_sim_total;
    simulate_SIS<erdos_reyni, simulate_next_reaction>(engine, psi,rho, t_sim, y_sim_new,y_sim_total, T, N, R0);

    plot("nextreaction.sis.single.pdf", "SIS single trajectory [NextReaction]", [&](auto& gp, auto& p) {
        p.add_plot1d(std::make_pair(t_sim, y_sim_total), "with lines title 'next reaction'"s);
    });
}
#endif
