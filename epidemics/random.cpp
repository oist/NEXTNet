//
//  randomGenerators.cpp
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#include "stdafx.h"
#include "random.h"
#include "utility.h"


/*----------------------------------------------------*/
/*----------------------CLASS-------------------------*/
/*----------------TRANSMISSION TIME-------------------*/
/*----------------------------------------------------*/

transmission_time::transmission_time(double pinf)
    :pinfinity(pinf)
{
    if ((pinfinity < 0.0) || (pinfinity > 1.0) || !std::isfinite(pinfinity))
        throw std::range_error("pinfinity must lie within [0, 1]");
}

transmission_time::~transmission_time() {
    // Nothing to do
}

/* "this ->"" is also used to ensure that if some of the functions are 
* redefined in herited classes then they are they ones being used and
* not the by-default implementations.*/

interval_t transmission_time::sample(rng_t& rng, interval_t t, int m) const {
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("t must be non-negative and finite");
    if (m < 1)
        throw std::range_error("m must be positive");
    // sample by inverting the survival function
    const double u = std::uniform_real_distribution<double>(0, 1)(rng);
    return this->survivalquantile(u, t, m);
}

double transmission_time::hazardrate(interval_t tau) const {
    return this->density(tau) / this->survivalprobability(tau);
}

double transmission_time::hazardbound(interval_t tau) const {
	throw std::runtime_error("hazardbound is not implemented for this transmission_time");
}

double transmission_time::survivalprobability(interval_t tau, interval_t t, int m) const {
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("t must be non-negative and finite");
    if (m < 1)
        throw std::range_error("m must be positive");
    // By default compute the conditional survival probability
    // using the unconditional survival function Psi(tau) based on
    // Psi(tau | t, m) = (Psi(t + tau) / Psi(t))^m
    return std::pow((this->survivalprobability(t + tau) / this->survivalprobability(t)), m);
}

interval_t transmission_time::survivalquantile(double u) const {
    if ((u < 0) || (u > 1) || !std::isfinite(u))
        throw std::range_error("u-quantile undefined for u not in [0,1]");
    // pinfinity() is the lower bound of survivalprobability(), so the u-quantile
    // for any u < pinfinity is infinity.
    if (u < pinfinity)
        return INFINITY;
    // By default, numerically invert the survival function
    return inverse_survival_function(u, 1e-6, [&] (double tau) { return this->survivalprobability(tau); });
}

interval_t transmission_time::survivalquantile(double u, interval_t t, int m) const {
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("t must be non-negative and finite");
    if (m < 1)
        throw std::range_error("m must be positive");
    /* We consider i.i.d. tau_1, ..., tau_m and ask for the probability
     *
     *   p(tau) = P[ tau_1, ..., tau_m >= t + tau | tau_1, ..., tau_m >= t ].
     *
     * In terms of the (unconditional) CDF F of tau this yields
     *
     *                                   m
     *            (   1 - F(t + tau)   )
     *   p(tau) = ( ------------------ )
     *            (       1 - F(t)     )
     *
     * and solving p(tau) = u gives
     *
     *  t + tau = F^-1[ 1 - u^(1/m) * (1 - F(t)) ].
     *
     * In terms of the (unconditional) survival function G = 1 - F we get
     *
     *  t + tau = G^-1[ u^(1/m) * G(t) ].
     */
    const double up = this->survivalprobability(t) * std::pow(u, 1.0 / double(m));
    const interval_t t_plus_tau = this->survivalquantile(up);
    if (t_plus_tau < t)
        throw std::logic_error("encountered invalid result when inverting the survival function");
    return (t_plus_tau - t);
}


/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:EXPONENTIAL------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

double transmission_time_exponential::density(interval_t tau) const {
    return lambda*exp(-lambda*tau);
}

double transmission_time_exponential::hazardrate(interval_t tau) const {
    return lambda;
}

double transmission_time_exponential::hazardbound(interval_t) const {
	return lambda;
}

double transmission_time_exponential::survivalprobability(interval_t tau) const {
    return exp(-lambda*tau);
}

double transmission_time_exponential::survivalprobability(interval_t tau, interval_t t, int m) const {
    return exp(-lambda * m * tau);
}

double transmission_time_exponential::survivalquantile(double u) const {
    return -mean * log(u);
}

double transmission_time_exponential::survivalquantile(double u, interval_t t, int m) const {
    return -mean/m * log(u);
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:GAMMA------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

double transmission_time_gamma::hazardbound(interval_t) const {
	return mean / variance;
}

double transmission_time_weibull::hazardbound(interval_t) const {
	return INFINITY;
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:POLYNOMIAL RATE--------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

transmission_time_polynomial_rate::transmission_time_polynomial_rate(std::vector<double>&& _coeffs, int dummy)
    :coeffs(std::move(_coeffs))
{
    bool all_zero = true;
    for(double c: coeffs) {
        if ((c < 0.0) || (!std::isfinite(c)))
            throw std::range_error("coefficients must be non-negative");
        all_zero = all_zero && (c == 0.0);
    }
    if (all_zero)
        throw std::range_error("polynomial cannot be identically zero");
}

double transmission_time_polynomial_rate::density(interval_t tau) const {
    return hazardrate(tau) * std::exp(-totalhazard(tau));
}

double transmission_time_polynomial_rate::hazardrate(interval_t tau) const {
    double r=0.0;
    for(std::size_t i=0; i < coeffs.size(); ++i)
        r += std::pow(tau, i) * coeffs[i];
    return r;
}

double transmission_time_polynomial_rate::totalhazard(interval_t tau) const {
    double s=0.0;
    for(std::size_t i=0; i < coeffs.size(); ++i)
        s += std::pow(tau, i+1) * coeffs[i] / (i+1);
    return s;
}

double transmission_time_polynomial_rate::hazardbound(interval_t) const {
    return INFINITY;
}

double transmission_time_polynomial_rate::survivalprobability(interval_t tau) const {
    return std::exp(-totalhazard(tau));
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:DETERMINISTIC------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

double transmission_time_deterministic::sample(rng_t&, interval_t t, int m) const{
    if (t>value)
        throw std::logic_error("present time cannot be larger than sampled time");
    return value;
}

double transmission_time_deterministic::density(interval_t tau) const {
        return (tau==value) ? INFINITY : 0;
}

double transmission_time_deterministic::hazardrate(interval_t) const {
    throw std::runtime_error("not implemented");
    return -1;
}

double transmission_time_deterministic::hazardbound(interval_t) const {
    throw std::runtime_error("not implemented");
    return -1;
}

double transmission_time_deterministic::survivalprobability(interval_t tau) const {
    throw std::runtime_error("not implemented");
    return -1;
}

double transmission_time_deterministic::survivalprobability(interval_t tau, interval_t t, int m) const {
    throw std::runtime_error("not implemented");
    return -1;
}

interval_t transmission_time_deterministic::survivalquantile(double u) const {
    throw std::runtime_error("not implemented");
    return -1;
}

interval_t transmission_time_deterministic::survivalquantile(double u, interval_t t, int m) const {
    return -1;
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------SUB RNGS---------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

struct sub_seed_adapter {
    typedef std::uint_least32_t result_type;

    sub_seed_adapter(rng_t& engine_) :engine(engine_) {}

    template<typename It>
    void generate(It begin, It end) {
        for(It i = begin; i != end; ++i)
            *i = engine();
    }

    rng_t& engine;
};

sub_rngs::sub_rngs(std::size_t n, rng_t& engine) {
    // Generate RNGs
    rngs.reserve(n);
    sub_seed_adapter s(engine);
    for(std::size_t i = 0; i < n; ++i)
        rngs.emplace_back(s);
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- VARIOUS HELPER FUNCTIONS ---------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

std::vector<double> rgamma( int n, double a, double b, rng_t& engine) {
    std::gamma_distribution<double> gam(a,b);
    
    std::vector<double> vec(n,0.0); // Initialise empty vector of size n.
    for (int i = 0; i< n; i++)
        vec[i] = gam(engine);

    return vec;
}

double rgamma(double a, double b, rng_t& engine) {
    std::gamma_distribution<double> gam(a,b);
    return gam(engine);
}

std::vector<double> rand( int n, rng_t& engine) {
    std::vector<double> vec(n,0.0);
    std::uniform_real_distribution<> dis(0,1);
    for (int i =0; i<n; i++) {
        vec[i]=dis(engine);
    }
    return vec;
}

double rand(double a, double b, rng_t& engine) {
    std::uniform_real_distribution<> dis(a,b);
    return dis(engine);
}

int poissrnd(double lambda, rng_t& engine) {
    typedef std::poisson_distribution<int> pois_int_t;
    pois_int_t pois(lambda);
    return(pois(engine));
}

int zipf(int alpha, int n, rng_t& engine)
{
  static int first = 1;// Static first time flag
  static double harmonic= 0;
  static std::vector<double> c(n+1,0);  // cumulant function
  double z;                     // Uniform random number (0 < z < 1)
  int zipf_value;            // Computed exponential value to be returned
  int low = 1;
  int high = n;
  int mid;           // Binary-search bounds
  
  
  // Compute normalization constant on first call only
  if (first == 1)
  {
    for (int i = 1; i <= n; i++) {
      harmonic += 1 / pow(i,alpha);
    }
    
    c[0]=0;
    for (int i=1; i<=n; i++)
      c[i] = c[i-1] + 1.0 / (pow(i,alpha+1)*harmonic);
    
    first = 0;
  }
  
  // Pull a uniform random number (0 < z < 1)
  do
  {
    z = rand(0,1,engine);
  }
  while ((z == 0) || (z == 1));
  
  
  while (true) {
    mid = floor((low+high)/2);
    if (c[mid-1] < z && z <= c[mid] ) {
      zipf_value = mid;
      break;
    } else if (c[mid] >= z) {
      high = mid-1;
    } else {
      low = mid+1;
    }
  }
  
  // Assert that zipf_value is between 1 and N
  assert((zipf_value >=1) && (zipf_value <= n));
  
  return(zipf_value);
}

const double PI = 3.1415926535897932384626433832795028842;

double pdf_log_normal(double t,double mean, double variance) {
    if (t==0) {
        return 0;
    }
    double mu = 2 * log(mean) - 0.5 * log( pow(mean,2.0) + variance );
    double sigma = sqrt( log( 1 + variance/pow(mean,2.0)));
    
    
    return 1/ (t*sigma*sqrt(2*PI)) * exp( -pow(log(t)-mu,2) / (2*sigma*sigma) );
}
double cdf_log_normal(double t,double mean, double variance) {
    if (t==0) {
        return 0;
    }
    double mu = 2 * log(mean) - 0.5 * log( pow(mean,2.0) + variance );
    double sigma = sqrt( log( 1 + variance/pow(mean,2.0)));
    return 0.5 * (1 + erf( (log(t)-mu) / (sqrt(2)*sigma) ));
}
