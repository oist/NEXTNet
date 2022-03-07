#pragma once

#include "stdafx.h"
#include "Tau.h"

/* Simulates path */
std::vector<double> simulatePath(std::vector<double>& infection_times, int n_max, Tau& tau ,std::mt19937& mersenneTwister);

std::vector<double> simulatePathRenewal(std::vector<double>& infection_times, int n_max, Tau& tau ,std::mt19937& mersenneTwister);

//std::vector<double> intialiseInfectionTimes(int number_of_infected, double r0,double mean, double variance);


std::vector<double> simulatePathNetwork(int network_size,double degree, Tau& tau, std::mt19937& mersenneTwister );

void print_matrix(std::vector<std::vector<double>>& A);
