#pragma once

#if ENABLE_PLOTTING

#    include "nextnet/tests/stdafx.h"
#    include "nextnet/tests/gnuplot-iostream.h"

namespace gp = gnuplotio;

void plot(const std::string filename, const std::string title, std::function<void(gp::Gnuplot &, gp::PlotGroup &)> body);

#endif
