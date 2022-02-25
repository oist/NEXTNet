//
//  randomGenerators.cpp
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#include "random.h"
#include "stdafx.h"

using namespace std;

/*-----RANDOM NUMBER GENERATOR-----*/

/*  MersenneTwister random number generator */
//static mt19937 mersenneTwister(1);


/* Generate infection times of a single individual
 In this example we generate a time tau where LOG(tau) is normally distributed with mean mu and variance sigma^2 */
vector<double> beta_normalised( int n, double mean, double variance, mt19937& mersenneTwister){
    
    double mu = 2 * log(mean) - 0.5 * log( pow(mean,2.0) + variance );
    double sigma = sqrt( log( 1 + variance/pow(mean,2.0)));
    
    lognormal_distribution<double> log_distribution(mu,sigma);
    
    vector<double> vec(n,0.0); // Initialise empty vector of size n.
    
    for (int i = 0; i< n; i++)
        vec[i] = log_distribution(mersenneTwister);

    return vec;
}



/*  Create uniformly distributed random numbers using the Mersenne Twister algorithm. */
vector<double> randu( int n,mt19937& mersenneTwister){
    vector<double> vec(n,0.0);
    uniform_real_distribution<> dis(0,1);
    for (int i =0; i<n; i++) {
        vec[i]=dis(mersenneTwister);
    }
    return vec;
}

int poissrnd(double lambda,mt19937& mersenneTwister) {
    // poisson_distribution<int>(lambda)(mersenneTwister)
    typedef poisson_distribution<int> pois_int_t;
    pois_int_t pois(lambda);
    return(pois(mersenneTwister));
}

#if 0
// My own random number generators
vector<double> poissrnd(double lambda,int n){
    vector<double> vec(n, 0.0);
    vector<double> u = randu(n);
    for (int i =0; i<n; i++) {
        int k=0;
        double bound = u[i] * exp(lambda);
        double invcdf = 1.0;
        double factorial = 1.0;
        while (bound > invcdf) {
            k++;
            factorial *= k;
            invcdf += pow(lambda,k) / factorial;
        }
        vec[i]= k;
    }
    return  vec;
}

int poissrnd(double lambda){
    int k=1;
    vector<double> u = randu(1);
    double bound = u[0] * exp(lambda);
    double invcdf = 1.0;
    double factorial = 1.0;
    while (bound > invcdf) {
        factorial *= k;
        invcdf += pow(lambda,k) / factorial;
        k++;
    }
    return  k;
}

#endif
