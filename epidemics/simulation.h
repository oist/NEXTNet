#pragma once

#include "stdafx.h"


/* Simulates path */
std::vector<double> simulatePath(std::vector<double>& infection_times, int n_max, double mean, double variance,double r0);

std::vector<double> intialiseInfectionTimes(int number_of_infected, double r0,double mean, double variance);
