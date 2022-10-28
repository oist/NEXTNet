#pragma once

#if ENABLE_PLOTTING

#include "tests/stdafx.h"
#include "tests/gnuplot-iostream.h"

namespace gp = gnuplotio;

void plot(const std::string filename, const std::string title, std::function<void(gp::PlotGroup&)> body);

#endif
