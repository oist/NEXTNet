#pragma once

#include "stdafx.h"
#include "types.h"
#include "graph.h"

//--------------------------------------
//---------GENERATING PATHS-------------
//--------------------------------------

/* Simulates several paths in the mean field regime */
void simulatePaths_MeanField(double mean, double variance, int degree,int nb_paths,double size, rng_t& engine);

/* Simulates path with the non-Markovian Gillespie Algorithm following Boguna et al */
void generatePaths_NMGA(double mean, double variance, int degree,int nb_paths,double size, rng_t& engine,double threshold);

/* Simulates a next reaction scheme on a Erdos Reyni Network */
void generatePaths_next_reaction(double mean, double variance, int degree,int nb_paths,double size, rng_t& engine);



