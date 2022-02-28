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
std::vector<double> rand(int n,rng_t& mersenneTwister);
double rand(double a, double b, rng_t& mersenneTwister);

/* Create a vector of size N where each element is an IID lognormal number */
std::vector<double> beta_normalised( int n,double mean, double variance,rng_t& mersenneTwister);

//std::vector<double> poissrnd(double lambda,int n);
int poissrnd(double lambda,rng_t& mersenneTwister);

// Gamma random number
std::vector<double> rgamma( int n, double a, double b, rng_t& mersenneTwister);
double rgamma(double a, double b, rng_t& mersenneTwister);


// Uniform random number
/* Resets the random number generator*/
//To do
//void rng() ;

/* CDF Gaussian*/
double normcdf( double x );

/*Inverse CDF Gaussian*/
double norminv( double x );

/* Testing */

void testNormCdf();

void testNormInv();
