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
	 *
	 * m is the weight/multiplicity of the edge which scales the hazardrate,
	 * see also survivalprobability / survivalquantile.
     */
    virtual interval_t sample(rng_t&, interval_t t, double m) const;

    /**
     * Probability Density Function psi(tau) of the samples.
     */
    virtual double density(interval_t tau) const = 0;

	/**
	 * Probability Density Function psi(tau) of the samples.
	 *
	 * m is the weight/multiplicity of the edge which scales the hazardrate,
	 * see also survivalprobability / survivalquantile.
	 */
	virtual double density(interval_t tau, double m) const;

    /**
     * "hazard rate" of the process, denoted by lambda in Boguna
     * and defined as lambda(tau) = psi(tau) / survivalprobability(tau).
	 *
	 * The reported hazardrate is for edges with weight/multiplicity 1.0,
	 * and must be scaled for edges with a different weight/mutiplicity.
     */
    virtual double hazardrate(interval_t) const;

	/**
	 * Upper bound of the "hazard rate" for times up t. Let the hazard
	 * rate be lambda(tau) = psi(tau) / survivalprobability(tau), then this
	 * function returns a value m(tau) such that lambda(t) <= m(t)
	 * for t in [0, tau]. For t = infinity, this function returns a global
	 * upper bound on lambda(tau) if one exists, otherwise infinity.
	 *
	 * The reported hazardbound is for edges with weight/multiplicity 1.0,
	 * and must be scaled for edges with a different weight/mutiplicity.
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
	 *
	 * m is used both to support weighted networks where it is equal to the weight
	 * (i.e. multiplicity) of the edge, and to support the sequential edge mode of the
	 * next reaction schema.
     */
    virtual double survivalprobability(interval_t tau, interval_t t, double m) const;

    /**
     * Evaluates the inverse of the survival function Psi(tau).
     */
    virtual interval_t survivalquantile(double u) const;

    /**
     * Evaluates the inverse of the survival function Psi(tau | t, m).
     */
    virtual interval_t survivalquantile(double u, interval_t t, double m) const;
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

    virtual double density(interval_t tau) const override;
    virtual double hazardrate(interval_t) const override;
	virtual double hazardbound(interval_t) const override;
    virtual double survivalprobability(interval_t tau) const override;
    virtual double survivalprobability(interval_t tau, interval_t t, double m) const override;
    virtual interval_t survivalquantile(double u) const override;
    virtual interval_t survivalquantile(double u, interval_t t, double m) const override;
};

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:GENERIC----------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

template<typename MathDistributionType, typename RandomDistributionType>
class transmission_time_generic_boost : public transmission_time {
    typedef MathDistributionType math_distribution_t;
    typedef RandomDistributionType random_distribution_t;
    const math_distribution_t math_distribution;
    const random_distribution_t random_distribution;

public:
    transmission_time_generic_boost(const math_distribution_t& md, const random_distribution_t& rd, double pinfinity = 0.0)
        :transmission_time(pinfinity)
        ,math_distribution(md)
        ,random_distribution(rd)
    {}

    virtual interval_t sample(rng_t& e, interval_t t, double m) const override {
        if ((t < 0) || !std::isfinite(t))
            throw std::range_error("condition t must be non-negative and finite");
        if (m < 0.0)
            throw std::range_error("m must be positive");
        // For m == 1 and t == 0.0 we sample from the base distribution,
        // and may thus use the random distribution class' faster implementation.
        // Otherwise, we sample from a modified distribution and fall back
        // to CDF inversion.
		if (m == 0.0)
			return INFINITY;
        if ((m == 1.0) && (t == 0.0)) {
            if ((pinfinity > 0) && std::bernoulli_distribution(pinfinity)(e))
                return INFINITY;
            return random_distribution_t(random_distribution.param())(e);
        }
        return transmission_time::sample(e, t, m);
    }

    virtual double density(interval_t tau) const override {
        if (tau < 0)
            return 0;
        return (1.0 - pinfinity) * pdf(math_distribution, tau);
    }

    // use the implementations of the conditional probabilities
    // in terms of the unconditional probabilities from the base class
    using transmission_time::survivalprobability;
    using transmission_time::survivalquantile;
    
    virtual double survivalprobability(interval_t tau) const override {
        // distribution has a scaled version of the specified CDF on tau in [0, inf)
        // and point mass pinfinity at tau = infinity
        if (tau < 0)
            return 1.0;
        if (std::isinf(tau))
            return pinfinity;
        return pinfinity + (1.0 - pinfinity) * cdf(complement(math_distribution, tau));
    }

    virtual interval_t survivalquantile(double u) const override {
        // distribution has a scaled version of the specified CDF on tau in [0, inf)
        // and point mass pinfinity at tau = infinity, so pinfinity is the lower bound
        // of survivalprobability() and any u-quantile for any u < pinfinity is infinity.
        if ((u < 0) || (u > 1) || !std::isfinite(u))
            throw std::range_error("u-quantile undefined for u not in [0,1]");
        if (u < pinfinity)
            return INFINITY;
        return quantile(complement(math_distribution, (u - pinfinity) / (1.0 - pinfinity)));
    }
};

// Define a specific policy:
typedef bm::policies::policy<
      bm::policies::overflow_error<bm::policies::ignore_error>,
      bm::policies::underflow_error<bm::policies::ignore_error>
> ignore_error_policy;

#define TRANSMISSION_TIME_GENERIC_BOOST(name, boost_tpl, std_tpl) \
    class name: public transmission_time_generic_boost<boost::math::boost_tpl<double, ignore_error_policy>, std::std_tpl<double>> { \
        typedef boost::math::boost_tpl<double, ignore_error_policy> bm_dist_t; \
        typedef std::std_tpl<double> std_dist_t;

#define TRANSMISSION_TIME_GENERIC_BOOST_SAMENAME(name) \
    TRANSMISSION_TIME_GENERIC_BOOST(transmission_time_ ## name, name ## _distribution, name ## _distribution)

/* Defines a distribution parameterized in terms of its mean and variance */
#define TRANSMISSION_TIME_GENERIC_MEAN_VARIANCE(name, boost_tpl, std_tpl, p1name, p1formula, p2name, p2formula) \
    TRANSMISSION_TIME_GENERIC_BOOST(name, boost_tpl, std_tpl) \
    static double p1name(const double m, const double v) { return p1formula; } \
    static double p2name(const double m, const double v) { return p2formula; } \
    public: \
    const double mean; \
    const double variance; \
    name(double m, double v, double pinf = 0.0) \
        :transmission_time_generic_boost(bm_dist_t(p1name(m, v), p2name(m, v)), \
                                         std_dist_t(p1name(m, v), p2name(m, v)), \
                                         pinf) \
        ,mean(m), variance(v) \
    {}

/* Defines a distribution parameterized by two double parameters */
#define TRANSMISSION_TIME_GENERIC_2PARAM(name, boost_tpl, std_tpl, p1name, p2name, meanformula, varformula) \
    TRANSMISSION_TIME_GENERIC_BOOST(name, boost_tpl, std_tpl) \
    public: \
    const double mean; \
    const double variance; \
    name(double p1name, double p2name, double pinf = 0.0) \
        :transmission_time_generic_boost(bm_dist_t(p1name, p2name), \
                                         std_dist_t(p1name, p2name), \
                                         pinf) \
        ,mean(meanformula), variance(varformula) \
    {}

#define TRANSMISSION_TIME_GENERIC_1PARAM(name, boost_tpl, std_tpl, p1name, meanformula) \
    TRANSMISSION_TIME_GENERIC_BOOST(name, boost_tpl, std_tpl) \
    public: \
    const double mean; \
    name(double p1name, double p2name, double pinf = 0.0) \
        :transmission_time_generic_boost(bm_dist_t(p1name), std_dist_t(p1name), pinf) \
        ,mean(meanformula) \
    {}

#define TRANSMISSION_TIME_GENERIC_SAMENAME_MEAN_VARIANCE(name, p1name, p1formula, p2name, p2formula) \
    TRANSMISSION_TIME_GENERIC_MEAN_VARIANCE(transmission_time_ ## name, name ## _distribution, name ## _distribution, \
                                            p1name, p1formula, p2name, p2formula)

#define TRANSMISSION_TIME_GENERIC_SAMENAME_2PARAM(name, p1name, p2name, meanformula, varformula) \
    TRANSMISSION_TIME_GENERIC_2PARAM(transmission_time_ ## name, name ## _distribution, name ## _distribution, \
                                     p1name, p2name, meanformula, varformula)

#define TRANSMISSION_TIME_GENERIC_SAMENAME_1PARAM(name, p1name,  meanformula) \
    TRANSMISSION_TIME_GENERIC_1PARAM(transmission_time_ ## name, name ## _distribution, name ## _distribution, \
                                     p1name, meanformula)

/*-----------TRANSMISSION TIMES DEFINED IN TERMS OF BOOST AND C++ STDLIB ---------*/

TRANSMISSION_TIME_GENERIC_SAMENAME_MEAN_VARIANCE(lognormal,
                                                 mu,  2 * log(m) - 0.5 * log(pow(m,2.0) + v),
                                                 sigma, sqrt(log( 1 + v/pow(m,2.0))))
};

TRANSMISSION_TIME_GENERIC_SAMENAME_MEAN_VARIANCE(gamma,
                                                 shape, pow(m, 2.0) / v,
                                                 scale, v / m)
    virtual double hazardbound(interval_t) const;
};

TRANSMISSION_TIME_GENERIC_SAMENAME_2PARAM(weibull, shape, scale,
                                          scale * std::tgamma(1 + 1/shape),
                                          pow(scale,2) * (std::tgamma(1 + 2/shape)- pow(std::tgamma(1 + 1/shape),2)))
};


/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME: POLYNOMIAL RATE-------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

class transmission_time_polynomial_rate : public transmission_time
{
public:
    template<typename ...Args>
    transmission_time_polynomial_rate(Args&& ...args)
        :transmission_time_polynomial_rate(std::vector<double>(std::forward<Args>(args)...), 0)
    {}

private:
    explicit
    transmission_time_polynomial_rate(std::vector<double>&& _coeffs, int dummy);

public:
    const std::vector<double> coeffs;

    virtual double density(interval_t tau) const override;

    virtual double hazardrate(interval_t) const override;

    double totalhazard(interval_t) const;

    virtual double hazardbound(interval_t) const override;

    virtual double survivalprobability(interval_t tau) const override;
};


/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME: DETERMINISTIC---------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

class transmission_time_deterministic : public transmission_time
{
public:

    transmission_time_deterministic(double v) :value(v) {}

    const double value;

    virtual interval_t sample(rng_t&, interval_t t, double m) const override;

    virtual double density(interval_t tau) const override;

    virtual double hazardrate(interval_t) const override;
	
	virtual double hazardbound(interval_t) const override;

    virtual double survivalprobability(interval_t tau) const override;

    virtual double survivalprobability(interval_t tau, interval_t t, double m) const override;

    virtual interval_t survivalquantile(double u) const override;

    virtual interval_t survivalquantile(double u, interval_t t, double m) const override;
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

