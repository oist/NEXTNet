/*----------------------------------------------------------------------------------------*/
/* non-Markovian Gillespie Algorithm  */

/* DOI ORIGINAL PAPER: https://doi.org/10.1103/PhysRevE.90.042108  */

/* QUOTE FROM THE PAPER:
 
 " we consider that each active link (connecting a susceptible-infected pair) defines a statistically independent infection process following the distribution psi(tau).
 That is, a susceptible individual connected to a single infected individual will become infected after a random time distributed by psi(tau) from the moment the link became active.
 If the susceptible individual is connected to more than one infected neighbor, each active link is considered as statistically independent so that the infection event will take place at the time of the first firing event of any of the current active links." */


/* CURRENT IMPLEMENTATION:
 - Consider a mean-field degree regime: Each individual has Poi(R0) neighbours and no infected nodes share the same neighbour.
 */


/*----------------------------------------------------------------------------------------*/
#include "stdafx.h"
#include "random.h"
#include "boguna.h"
#include "Tau.h"
#include "testing.h"

using namespace std;

/* Only statistically exact version in the limit n-> \infty // version of their paper.*/
// If starting with one individual infection_times needs to be thermalised!.
// Elements of infection_times are the time at which they contracted the disease.
// MEAN FIELD VERSION.
//TO DO
void simulateNMGA(vector<double>& trajectory, vector<double>& infection_times, int n_max,Tau& tau,mt19937& mersenneTwister){
//    /*-----PARAMETERS------*/
//    long n = infection_times.size(); // number of renewal processes running.
//    double absolute_time = 0; // calendar time
//    double dt; // increments
//
//    double sum_lambda = 0;
//    vector<double> lambda({}); // instantaneous rate of the n processes at calendar time t
//
//    /*-----INITIALISATION------*/ // Optional // lazy implementation
//    if (n==0) {
//        // Solving PHI[dt|{tj}] with {tj}={0}; is equiv to draw a random number from the log-normal distribution. (or whatever distribution for the infeection time we are using).
//
//        absolute_time = beta_normalised(tau.mean,tau.variance,mersenneTwister);
//        /*----Add new infection times----*/
//
//        double new_active_links = poissrnd(tau.r0, mersenneTwister); // neighbours of the newly infected.
//        n = new_active_links;
//
//        for (int i=0; i< new_active_links; i++) {
//            infection_times.push_back(absolute_time);
//            lambda.push_back(insta_rate(absolute_time,tau));
//            cout << lambda[i]<<endl;
//        }
//        sum_lambda = insta_rate(absolute_time,tau) * new_active_links;
//        cout <<sum_lambda<<endl;
//    }
//
//
//    /*----SIMULATION------*/
//    for (int population =0; population< n_max; population++) {
//
//        /* Draw the next time until the next event*/
//        double u = rand(0,1,mersenneTwister);
//        dt = -log(u)/sum_lambda;
//
//
//        /* Select the process that has generated the event*/
//        vector<double> cumulant(n,lambda[0]/sum_lambda);
//        double v = rand(0,1,mersenneTwister);
//        int i = 0;
//        while (v > cumulant[i]){
//            i++;
//            cumulant[i]+= cumulant[i-1] + lambda[i]/sum_lambda;
//        }
//        int selected_process = i;
//
//        /*Update the time and the rates since the last event*/
//        absolute_time += dt ;
//
        
//
//    }
}


double insta_rate(double t,Tau& tau){
    if (t>0) {
        double psi = pdf_log_normal(t, tau.mean, tau.variance);
        double PSI = (1-cdf_log_normal(t, tau.mean, tau.variance)) ;
        return psi/PSI;
    }
    return 0;

}


//NETWORK VERSION: TO DO
//vector<double> simulateNMGA(vector<double>& infection_times, int n_max,Tau& tau,mt19937& mersenneTwister){
//
