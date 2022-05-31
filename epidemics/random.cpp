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

/* "this ->"" is also used to ensure that if some of the functions are 
* redefined in herited classes then they are they ones being used and
* not the by-default implementations.*/

interval_t transmission_time::sample(rng_t& rng, interval_t t, int m) {
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("t must be non-negative and finite");
    if (m < 1)
        throw std::range_error("m must be positive");
    // sample by inverting the survival function
    const double u = std::uniform_real_distribution<double>(0, 1)(rng);
    return this->survivalquantile(u, t, m);
}

double transmission_time::hazardrate(interval_t tau) {
    return this->density(tau) / this->survivalprobability(tau);
}

double transmission_time::survivalprobability(interval_t tau, interval_t t, int m) {
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("t must be non-negative and finite");
    if (m < 1)
        throw std::range_error("m must be positive");
    // By default compute the conditional survival probability
    // using the unconditional survival function Psi(tau) based on
    // Psi(tau | t, m) = (Psi(t + tau) / Psi(t))^m
    return std::pow((this->survivalprobability(t + tau) / this->survivalprobability(t)), m);
}

interval_t transmission_time::survivalquantile(double u) {
    if ((u < 0) || (u > 1) || !std::isfinite(u))
        throw std::range_error("u-quantile undefined for u not in [0,1]");
    // pinfinity() is the lower bound of survivalprobability(), so the u-quantile
    // for any u < pinfinity is infinity.
    if (u < pinfinity)
        return INFINITY;
    // By default, numerically invert the survival function
    return inverse_survival_function(u, 1e-6, [&] (double tau) { return this->survivalprobability(tau); });
}

interval_t transmission_time::survivalquantile(double u, interval_t t, int m) {
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("t must be non-negative and finite");
    if (m < 1)
        throw std::range_error("m must be positive");
    /* We consider i.i.d. tau_1, ..., tau_m and ask for the probability
     *
     *   p(tau) = P[ tau_1, ..., tau_m >= t + tau | tau_1, ..., tau_m >= t ].
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
     * In terms of the (unconditional) survival function G = 1 -F we get
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

double transmission_time_exponential::density(interval_t tau) {
    return lambda*exp(-lambda*tau);
}

double transmission_time_exponential::hazardrate(interval_t tau) {
    return lambda;
}

double transmission_time_exponential::survivalprobability(interval_t tau) {
    return exp(-lambda*tau);
}

double transmission_time_exponential::survivalprobability(interval_t tau, interval_t t, int m) {
    return exp(-lambda * m * tau);
}

double transmission_time_exponential::survivalquantile(double u) {
    return -mean * log(u);
}

double transmission_time_exponential::survivalquantile(double u, interval_t t, int m) {
    return -mean/m * log(u);
}

std::vector<double> rgamma( int n, double a, double b, std::mt19937& mersenneTwister){
    
    std::gamma_distribution<double> gam(a,b);
    
    std::vector<double> vec(n,0.0); // Initialise empty vector of size n.

    for (int i = 0; i< n; i++)
        vec[i] = gam(mersenneTwister);

    return vec;
}

double rgamma(double a, double b, std::mt19937& mersenneTwister){
    std::gamma_distribution<double> gam(a,b);
    return gam(mersenneTwister);
}


/*  Create uniformly distributed random numbers using the Mersenne Twister algorithm. */
std::vector<double> rand( int n, std::mt19937& mersenneTwister){
    std::vector<double> vec(n,0.0);
    std::uniform_real_distribution<> dis(0,1);
    for (int i =0; i<n; i++) {
        vec[i]=dis(mersenneTwister);
    }
    return vec;
}

double rand(double a, double b, std::mt19937& mersenneTwister){
    std::uniform_real_distribution<> dis(a,b);
    return dis(mersenneTwister);
}

int poissrnd(double lambda, std::mt19937& mersenneTwister) {
    // poisson_distribution<int>(lambda)(mersenneTwister)
    typedef std::poisson_distribution<int> pois_int_t;
    pois_int_t pois(lambda);
    return(pois(mersenneTwister));
}


//void initialise_adjacency_matrix(vector<vector<int>>& A,vector<int>& K,int n, double degree,mt19937& mersenneTwister){
//    //vector<vector<int>> A(n, vector<int>(n)); //adjacency matrix n by n;
//    
//    vector<double> r =rand(n*n,mersenneTwister);
//    vector<double> infection_times = beta_normalised(n, Tau, mersenneTwister);
//    double p = degree/n;
//    for (int i=0; i<n; i++) {
//        for (int j=0; j<n; j++) {
//            if (r[i*n+j]>p) {
//                A[i][j]=0;
//            }
//            else{
//                A[i][j]=
//                K[i]+=1; // Count Neighbours
//            }
//        }
//    }
//}

const double PI = 3.1415926535897932384626433832795028842;

double pdf_log_normal(double t,double mean,double variance){
    if (t==0) {
        return 0;
    }
    double mu = 2 * log(mean) - 0.5 * log( pow(mean,2.0) + variance );
    double sigma = sqrt( log( 1 + variance/pow(mean,2.0)));
    
    
    return 1/ (t*sigma*sqrt(2*PI)) * exp( -pow(log(t)-mu,2) / (2*sigma*sigma) );
}
double cdf_log_normal(double t,double mean,double variance){
    if (t==0) {
        return 0;
    }
    double mu = 2 * log(mean) - 0.5 * log( pow(mean,2.0) + variance );
    double sigma = sqrt( log( 1 + variance/pow(mean,2.0)));
    return 0.5 * (1 + erf( (log(t)-mu) / (sqrt(2)*sigma) ));
}


