//
//  randomGenerators.h
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#pragma once

#include "stdafx.h"
#include "types.h"
#include "analysis.h"

//--------------------------------------
//----------INFECTION TIMES------------
//--------------------------------------

class beta {
public:
    /*
     * sample() returns a single (iid) sample from the
     * (normalized versio of) the distribution.
     */
    virtual interval_t sample(rng_t& engine) const = 0;

    /*
     * sample_next() returns the next single sample from the
     * the distribution condtioned on the last infection time of the individual.
     * /!\ returning infinity is possible and represnets no further infection.
     */
    virtual interval_t sample_next(interval_t last, rng_t& engine) const = 0;
    
    /*
     * sample_next_conditional() returns same as above but takes into consideration the
     * current number of healthy neighbours.
     * i.e if there are no more healthy neighbours
     * it raises an error as no more infections are possible.
     */

    virtual interval_t sample_next_conditional(interval_t last, int healthy, rng_t& engine) const = 0;
};


class lognormal_beta : public beta {
public:
    const double mean;
    const double variance;
    const double r0;
    const double mu = 2 * log(mean) - 0.5 * log( pow(mean,2.0)+ variance );
    const double sigma = sqrt( log( 1 + variance/pow(mean,2.0)));

public:
    lognormal_beta(double _mean, double _variance, double _r0)
        :mean(_mean), variance(_variance), r0(_r0)
    {}
    
    virtual interval_t sample(rng_t& engine) const;

    virtual interval_t sample_next(interval_t last, rng_t& engine) const;
    
    virtual interval_t sample_next_conditional(interval_t last, int healthy, rng_t& engine) const;

private:
    mutable std::lognormal_distribution<interval_t> log_distribution = std::lognormal_distribution<interval_t>(mu,sigma);
};

/*  Create uniformly distributed random numbers */
std::vector<double> rand(int n,rng_t& mersenneTwister);
double rand(double a, double b, rng_t& mersenneTwister);

///* Create a vector of size N where each element is an IID lognormal number */
//std::vector<double> beta_normalised( int n,Tau& tau,rng_t& mersenneTwister);
//
//double beta_normalised(Tau& tau, rng_t& mersenneTwister);


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

