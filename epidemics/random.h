//
//  randomGenerators.h
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#pragma once

#include "stdafx.h"
#include "savingData.h"

typedef std::mt19937 rng_t;

/*  Create uniformly distributed random numbers */
std::vector<double> randu(int n,rng_t& mersenneTwister);

/* Create a vector of size N where each element is an IID lognormal number */
std::vector<double> beta_normalised( int n,double mean, double variance,rng_t& mersenneTwister);

//std::vector<double> poissrnd(double lambda,int n);
int poissrnd(double lambda,rng_t& mersenneTwister);

/* Resets the random number generator*/
//To do
//void rng() ;
