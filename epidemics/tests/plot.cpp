#if ENABLE_PLOTTING

#include "tests/stdafx.h"
#include "tests/gnuplot-iostream.h"
#include "tests/plot.h"

namespace gp = gnuplotio;

void plot(const std::string filename, const std::string title, std::function<void(gp::Gnuplot&, gp::PlotGroup&)> body) {
    std::filesystem::create_directory("tests.out");
    gp::Gnuplot gp;
    gp::PlotGroup group = gp.plotGroup();
    gp << "set terminal pdf\n";
    gp << "set output 'tests.out/" << filename << "'\n";
    gp << "set title '" << title << "'\n";
    body(gp, group);
    gp << group;
    gp.flush();
}

#endif
