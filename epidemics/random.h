//
//  randomGenerators.h
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#pragma once

#include "stdafx.h"
#include "types.h"
// #include "analysis.h"

//--------------------------------------
//----------INFECTION TIMES------------
//--------------------------------------

class transmission_time {
public:
    /**
     * The probability of no infection, i.e. of the infection time being infinite.
     * The total probability represented by the density function *plus*
     * pinfinity() must be one, meaning that the integral over the density
     * must be 1 - pinfinity().
     */
    const double pinfinity = 0.0;

    transmission_time(double pinf = 0.0);
    
    virtual ~transmission_time();

    /**
     * Samples from the distribution with survival function Psi( tau | t,m ) and pdf psi(tau)
     * May return infinity to indicate that no further infections occur if pinfinity() > 0.
     */
    virtual interval_t sample(rng_t&, interval_t t, int m) const;

    /**
     * Probability Density Function psi(tau) of the samples.
     */
    virtual double density(interval_t tau) const = 0;

    /**
     * "hazard rate" of the process, denoted by lambda in Boguna
     * and defined as lambda(tau) = psi(tau) / survivalprobability(tau).
     */
    virtual double hazardrate(interval_t) const;
	
	/**
	 * Upper bound of the "hazard rate" for times up t. Let the hazard
	 * rate be lambda(tau) = psi(tau) / survivalprobability(tau), then this
	 * function returns a value m(tau) such that lambda(t) <= m(t)
	 * for t in [0, tau]. For t = infinity, this function returns a global
	 * upper bound on lambda(tau) if one exists, otherwise infinity.
	 */
	virtual double hazardbound(interval_t) const;

    /**
     * Evaluates the survival function Psi(tau), 
     * i.e. the probability that a single edges does not fire within time interval [t , t + tau).
     * Since infinity is larger than any finite tau, this probability is always >= pinfinity().
     */
    virtual double survivalprobability(interval_t tau) const = 0;

    /**
     * Evaluates the survival function  Psi( tau | t, m), 
     * i.e. the probability that none of m edges fire within the time interval [t, t+tau).
     * Since infinity is larger than any finite tau, this probability is always >= pinfinity().
     */
    virtual double survivalprobability(interval_t tau, interval_t t, int m) const;

    /**
     * Evaluates the inverse of the survival function Psi(tau),
     * i.e returns the time interval given a probability in 
     */
    virtual interval_t survivalquantile(double u) const;

    /**
     * Evaluates the inverse of the survival function Psi(tau | t, m),
     * i.e returns the time interval given a probability in 
     */
    virtual interval_t survivalquantile(double u, interval_t t, int m) const;
};

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:EXPONENTIAL------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

class transmission_time_exponential : public transmission_time {
    const double lambda;
    const double mean = 1/lambda;
    const double variance = 1/pow(lambda,2);

public:
    transmission_time_exponential(double _lambda)
        :lambda(_lambda)
    {}

    virtual double density(interval_t tau) const;
    virtual double hazardrate(interval_t) const;
	virtual double hazardbound(interval_t) const;
    virtual double survivalprobability(interval_t tau) const;
    virtual double survivalprobability(interval_t tau, interval_t t, int m) const;
    virtual interval_t survivalquantile(double u) const;
    virtual interval_t survivalquantile(double u, interval_t t, int m) const;
};

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:GENERIC----------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

template<typename DistributionType>
class transmission_time_generic_boost : public transmission_time {
    typedef DistributionType distribution_t;
    const distribution_t distribution;

public:
    transmission_time_generic_boost(const distribution_t& d, double pinfinity = 0.0)
        :transmission_time(pinfinity)
        ,distribution(d)
    {}

    virtual double density(interval_t tau) const {
        // TODO: Should we scale according to pinfinity here?
        // Currently, the density is always normalized while the survivalfunction
        // takes pfinfinity into account, which is weird.
        return pdf(distribution, tau);
    }

    // use the implementations of the conditional probabilities
    // in terms of the unconditional probabilities from the base class
    using transmission_time::survivalprobability;
    using transmission_time::survivalquantile;
    
    virtual double survivalprobability(interval_t tau) const {
        // distribution has a scaled version of the specified CDF on tau in [0, inf)
        // and point mass pinfinity at tau = infinity
        if (std::isinf(tau))
            return pinfinity;
        return pinfinity + (1.0 - pinfinity) * cdf(complement(distribution, tau));
    }

    virtual interval_t survivalquantile(double u) const {
        // distribution has a scaled version of the specified CDF on tau in [0, inf)
        // and point mass pinfinity at tau = infinity, so pinfinity is the lower bound
        // of survivalprobability() and any u-quantile for any u < pinfinity is infinity.
        if (u < pinfinity)
            return INFINITY;
        return quantile(complement(distribution, (u - pinfinity) / (1.0 - pinfinity)));
    }
};

// Define a specific policy:
typedef bm::policies::policy<
      bm::policies::overflow_error<bm::policies::ignore_error>,
      bm::policies::underflow_error<bm::policies::ignore_error>
> ignore_error_policy;

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:LOG NORMAL-------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

typedef bm::lognormal_distribution<double, ignore_error_policy> boost_lognormal_dist;

class transmission_time_lognormal
    :public transmission_time_generic_boost<boost_lognormal_dist>
{
    static double mu(const double mean, double variance) {
        return 2 * log(mean) - 0.5 * log( pow(mean,2.0)+ variance );
    }
    
    static double sigma(const double mean, double variance) {
        return sqrt( log( 1 + variance/pow(mean,2.0)));
    }
    
public:
    const double mean;
    const double variance;

    transmission_time_lognormal(double m, double v, double pinf = 0.0)
        :transmission_time_generic_boost(boost_lognormal_dist(mu(m, v), sigma(m, v)), pinf)
        ,mean(m), variance(v)
    {}
};

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:GAMMA------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

typedef bm::gamma_distribution<double, ignore_error_policy> boost_gamma_dist;

class transmission_time_gamma
    :public transmission_time_generic_boost<boost_gamma_dist>
{
    static double shape(const double mean, double variance) {
        return std::pow(mean, 2.0) / variance;
    }

    static double scale(const double mean, double variance) {
        return variance / mean;
    }

public:
    const double mean;
    const double variance;

    transmission_time_gamma(double m, double v, double pinf = 0.0)
        :transmission_time_generic_boost(boost_gamma_dist(shape(m, v), scale(m, v)), pinf)
        ,mean(m), variance(v)
    {}
	
	virtual double hazardbound(interval_t) const;
};

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------SUB RNGS---------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/**
 * @brief Creates n sub-RNGs from a root random number generator
 *
 * This can be used to make code that depends on randum numbers
 * parallelizlable. Typically, n is the number of work units,
 * and the work unit i would then use RNG i.
 */
struct sub_rngs {
    sub_rngs(std::size_t n, rng_t& engine);

    rng_t& operator[](std::size_t i) {
        return rngs[i];
    }

private:
    std::vector<rng_t> rngs;
};

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- VARIOUS HELPER FUNCTIONS ---------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/* Uniformly distributed random numbers */
std::vector<double> rand(int n, rng_t& engine);

double rand(double a, double b, rng_t& engine);

int poissrnd(double lambda, rng_t& engine);

/* Gamma random number */
std::vector<double> rgamma( int n, double a, double b, rng_t& engine);

double rgamma(double a, double b, rng_t& engine);

/* Zipfian random numbers */
int zipf(int alpha, int n,rng_t& engine);

/* Log-normal distribution */

double pdf_log_normal(double t,double mean, double variance);

double cdf_log_normal(double t,double mean, double variance);

