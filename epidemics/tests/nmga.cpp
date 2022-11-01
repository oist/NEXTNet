#include "tests/stdafx.h"
#include "tests/simulate.h"
#include "tests/analytical.h"
#include "tests/plot.h"

#include "random.h"
#include "nMGA.h"

#if ENABLE_PLOTTING
TEST_CASE("Plot large-population SIR mean-field (nMGA)", "[nMGA]") {
    using namespace std::string_literals;
    std::mt19937 engine;

    const std::size_t M = 1000;
    const std::size_t T = 35;
    const std::size_t X = 400;
    const double R0 = 2;
    const double MEAN = 10;
    const double VARIANCE = 1;
    transmission_time_gamma psi(MEAN, VARIANCE);

    /* Simulate using next reaction M times */
    auto racyclic = simulate_SIR<acyclic, simulate_nmga>(engine, psi, T, M, 1, R0+1, true);

    /* Evaluate analytical solution */
    std::pair<std::vector<double>, std::vector<double>> analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for(std::size_t i=0; i < X; ++i) {
        const double t = (double)T * i / (X-1);
        analytical.first.push_back(t);
        analytical.second.push_back(sol.N(t));
    }

    plot("nmga.sir.mean.pdf", "Large-population SIR mean-field [nMGA]", [&](auto& gp, auto& p) {
        gp << "set logscale y";
        p.add_plot1d(racyclic, "with lines title 'next reaction acyclic'"s);
        p.add_plot1d(analytical, "with lines title 'analytical'");
    });
}
#endif
