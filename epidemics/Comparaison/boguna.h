#pragma once

#include "stdafx.h"
#include "Tau.h"

std::vector<double> simulateNMGA_exact(int N,double alpha,double beta,std::mt19937& mersenneTwister);

/* Only statistically exact version in the limit n-> \infty - version in their paper.*/
std::vector<double> simulateNMGA(int n_max,Tau& tau,std::mt19937& mersenneTwister);
