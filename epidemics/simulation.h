#pragma once

#include "stdafx.h"

std::vector<double> simulatePath(std::vector<double>& infection_times, int n_max, double& absolute_time, double mean, double variance,double R0);

std::vector<double> intialiseInfectionTimes(int number_of_infected, double r0,double mean, double variance);
