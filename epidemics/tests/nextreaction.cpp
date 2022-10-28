#include "tests/stdafx.h"
#include "tests/simulate.h"
#include "tests/analytical.h"

#include "random.h"
#include "NextReaction.h"
#include "nMGA.h"

#if ENABLE_PLOTTING
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
#endif

#if ENABLE_PLOTTING
TEST_CASE("Plot large-population SIR mean-field (NextReaction)", "[nextreaction]") {
    using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t M = 1000;
    const std::size_t Nfully = 1000;
    const std::size_t Nerdos = 1000;
    const std::size_t T = 35;
    const std::size_t X = 400;
    const double R0 = 2;
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
    std::vector<double> t_analytical;
    std::vector<double> y_analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for(std::size_t i=0; i < X; ++i) {
        t_analytical.push_back((double)T * i / (X-1));
        y_analytical.push_back(sol.N(t_analytical.back()));
    }

    plt::figure_size(1600, 1200);
    plt::named_plot("next reaction fully-connected (N="s + std::to_string(Nfully) + ")", rfully.first, rfully.second);
    plt::named_plot("next reaction acyclic", racyclic.first, racyclic.second);
    plt::named_plot("next reaction Erdös-Reyni (N="s + std::to_string(Nerdos) + ")", rerdos.first, rerdos.second);
    plt::named_plot("analytical", t_analytical, y_analytical);
    plt::title("Large-population SIR mean-field [NextReaction]");
    plt::legend();
    std::filesystem::create_directory("tests.out");
    plt::save("tests.out/nextreaction.sir.mean.pdf");
    plt::close();
}
#endif

#if ENABLE_PLOTTING
TEST_CASE("Plot SIS single trajectory (NextReaction)", "[nextreaction]") {
       using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t N = 5000;
    const std::size_t T = 400;
    const std::size_t X = 200;
    const double R0 = 5;
    // Every infecteced node has N-1 neighbours, of which in the large-population limit N-2 are susceptible.
    // To trigger subsequent infections amgonst these N-2 susceptible neighbours, we must infect each neighbour
    // with probability R0/(N-2). Note that for the first infected node, this is not strictly speaking correct,
    // since it's infection has no source, and it will thus create R0(N-1)/(N-2) > R0 subsequent infections.
    // This only causes a relative error of |1 - (N-1)/(N-2)| =~= 1/N though.
    const double p = R0/(N-2);

    const double MEAN = 10;
    const double VARIANCE = 50;
    const double MEAN_rho = 70;
    const double VARIANCE_rho = 1;
    transmission_time_gamma psi(MEAN, VARIANCE, 1.0-p);
    transmission_time_gamma rho(MEAN_rho, VARIANCE_rho);

    /* Simulate using next reaction once times */
    std::vector<double> t_sim, y_sim_new, y_sim_total;
    simulate_SIS<fully_connected, simulate_next_reaction>(engine, psi,rho, t_sim, y_sim_new,y_sim_total, T, N);

    /* Evaluate analytical solution */
    std::vector<double> t_analytical;
    std::vector<double> y_analytical;

    for(std::size_t i=0; i < X; ++i) {
        const double t = (double)T * i / (X-1);
        t_analytical.push_back(t);
        y_analytical.push_back(N*(1-1/R0));
    }

    plt::figure_size(800, 600);
    plt::named_plot("next reaction", t_sim, y_sim_total);
    plt::named_plot("analytical", t_analytical, y_analytical);
    plt::title("SIS single trajectory [NextReaction]");
    plt::legend();
    std::filesystem::create_directory("tests.out");
    plt::save("tests.out/nextreaction.sis.single.pdf");
    plt::close();
}
#endif
