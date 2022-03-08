//
//  randomGenerators.h
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#pragma once

#include "stdafx.h"
#include "savingData.h"
#include "Tau.h"

typedef std::mt19937 rng_t;

/*  Create uniformly distributed random numbers */
std::vector<double> rand(int n,rng_t& mersenneTwister);
double rand(double a, double b, rng_t& mersenneTwister);

/* Create a vector of size N where each element is an IID lognormal number */
std::vector<double> beta_normalised( int n,Tau& tau,rng_t& mersenneTwister);

double beta_normalised(Tau& tau, rng_t& mersenneTwister);


//std::vector<double> poissrnd(double lambda,int n);
int poissrnd(double lambda,rng_t& mersenneTwister);

// Gamma random number
std::vector<double> rgamma( int n, double a, double b, rng_t& mersenneTwister);
double rgamma(double a, double b, rng_t& mersenneTwister);

// Erdos-Renyi Graph
//initialise_adjacency_matrix(std::vector<std::vector<int>>& A, std::vector<int>& K,int n, double degree,rng_t& mersenneTwister);
//

// Uniform random number
/* Resets the random number generator*/
//To do
//void rng() ;


double pdf_log_normal(double t,double mean ,double variance);

double cdf_log_normal(double t,double mean,double variance);

void test_pdf_logv();
void test_cdf_logv();


///**
// *  Computes the cumulative
// *  distribution function of the
// *  normal distribution
// */
//double normcdf( double x );
//
///**
// *  Computes the inverse of normcdf
// */
//double norminv( double x );