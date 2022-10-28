#include "tests/stdafx.h"
#include "tests/simulate.h"
#include "tests/analytical.h"

#include "random.h"
#include "nMGA.h"

#if ENABLE_PLOTTING
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
#endif

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
    std::vector<double> t_analytical;
    std::vector<double> y_analytical;
    meanfield_infpop_gamma sol = meanfield_infpop_gamma::mean_variance(MEAN, VARIANCE, R0);
    for(std::size_t i=0; i < X; ++i) {
        t_analytical.push_back((double)T * i / (X-1));
        y_analytical.push_back(sol.N(t_analytical.back()));
    }

    plt::figure_size(1600, 1200);
    plt::named_plot("next reaction (acyclic)", racyclic.first, racyclic.second);
    plt::named_plot("analytical", t_analytical, y_analytical);
    plt::title("Large-population SIR mean-field [nMGA]");
    plt::legend();
    std::filesystem::create_directory("tests.out");
    plt::save("tests.out/ngma.sir.mean.pdf");
    plt::close();
}
#endif
