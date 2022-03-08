//
//  randomGenerators.cpp
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#include "random.h"
#include "stdafx.h"
#include "testing.h"


using namespace std;

/*-----RANDOM NUMBER GENERATOR-----*/

/*  MersenneTwister random number generator */
//static mt19937 mersenneTwister(1);


/* Generate infection times of a single individual
 In this example we generate a time tau where LOG(tau) is normally distributed with mean mu and variance sigma^2 */
//vector<double> beta_normalised( int n, Tau& tau, mt19937& mersenneTwister){
//    
//    double mu = 2 * log(tau.mean) - 0.5 * log( pow(tau.mean,2.0) + tau.variance );
//    double sigma = sqrt( log( 1 + tau.variance/pow(tau.mean,2.0)));
//    
//    lognormal_distribution<double> log_distribution(mu,sigma);
//    
//    vector<double> vec(n,0.0); // Initialise empty vector of size n.
//    
//    for (int i = 0; i< n; i++)
//        vec[i] = log_distribution(mersenneTwister);
//
//    return vec;
//}

//double beta_normalised(Tau& tau, mt19937& mersenneTwister){
//    
//    double mu = 2 * log(tau.mean) - 0.5 * log( pow(tau.mean,2.0) + tau.variance );
//    double sigma = sqrt( log( 1 + tau.variance/pow(tau.mean,2.0)));
//    lognormal_distribution<double> log_distribution(mu,sigma);
//    
//    return log_distribution(mersenneTwister);
//}


vector<double> rgamma( int n, double a, double b, mt19937& mersenneTwister){
    
    gamma_distribution<double> gam(a,b);
    
    vector<double> vec(n,0.0); // Initialise empty vector of size n.

    for (int i = 0; i< n; i++)
        vec[i] = gam(mersenneTwister);

    return vec;
}

double rgamma(double a, double b, mt19937& mersenneTwister){
    gamma_distribution<double> gam(a,b);
    return gam(mersenneTwister);
}


/*  Create uniformly distributed random numbers using the Mersenne Twister algorithm. */
vector<double> rand( int n,mt19937& mersenneTwister){
    vector<double> vec(n,0.0);
    uniform_real_distribution<> dis(0,1);
    for (int i =0; i<n; i++) {
        vec[i]=dis(mersenneTwister);
    }
    return vec;
}

double rand(double a, double b, mt19937& mersenneTwister){
    uniform_real_distribution<> dis(a,b);
    return dis(mersenneTwister);
}






int poissrnd(double lambda,mt19937& mersenneTwister) {
    // poisson_distribution<int>(lambda)(mersenneTwister)
    typedef poisson_distribution<int> pois_int_t;
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


void test_cdf_logv() {
    ASSERT_APPROX_EQUAL( cdf_log_normal(10000, 1 ,1), 1, 0.01 );
    ASSERT_APPROX_EQUAL( cdf_log_normal(0, 1, 1), 0, 0.01 );
}
void test_pdf_logv() {
    ASSERT_APPROX_EQUAL( pdf_log_normal(0, 2, 3), 0, 0.01 );
    ASSERT_APPROX_EQUAL( pdf_log_normal(100000, 2, 3), 0, 0.01 );
}
