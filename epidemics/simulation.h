#pragma once

#include "stdafx.h"
#include "types.h"
#include "graph.h"

/* Simulates path */

std::vector<double> simulatePath(std::vector<double>& infection_times, int n_max, const lognormal_beta& infection_distribution, rng_t engine);





void print_matrix(std::vector<std::vector<double>>& A);

void simulateManyPaths(int nb_paths, rng_t& engine);
