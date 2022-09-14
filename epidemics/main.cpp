//
//  main.cpp
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "simulation.h"
#include "graph.h"
#include "nMGA.h"
#include "NextReaction.h"
using namespace std;

/*  MersenneTwister random number generator */
/*static*/ mt19937 engine;


int main(int argc, const char * argv[]) {

    engine.seed(1);
//    std::string path( std::filesystem::current_path() );
//    cout << path ;
//    for (int letter=0; letter<std::string("epidemics").size(); letter++) {
//        path.pop_back();
//    }
    cout << "starts"<< endl;
    imported_network nw(string("/Users/curesamuelcyrus/Documents/epidemics/mathematica/adjacency_list_mathematica.csv"));
    cout << "done"<< endl;
    
    cout << " old r : " <<  assortativity(nw) << endl ;
    
    add_correlation(0.4,nw,engine);
    
    cout << " new r : " <<  assortativity(nw) << endl ;
    export_adjacency_list(nw.adjacencylist,"network.txt"); 
    //export_adjacency_list(nwk.adjacencylist,"correlated_network.txt"); 

    // cout << "start...\n";

    // int size = 1000000;

    // // number_of_simulations:
    // int n = 1;

    // // User's input
    // if(argc == 2){
    //     size = atoi(argv[1]);
    // } else if (argc == 3){
    //     size = atoi(argv[1]);
    //     n = atoi(argv[2]);
    // }

    // double mean = 10;
    // double variance = 1.0;
    // transmission_time_lognormal psi(mean, variance);

    // scale_free network(size, engine);

    // cout << "network generated \n";
    // std::vector<long> degree;
    // for (int i=0; i<network.adjacencylist.size(); i++) {
    //     degree.push_back(network.adjacencylist[i].size());
    // }

    // ofstream out;

    // out.open(string("degree.dat"));

    // for (int i =0; i<degree.size(); i++){
    //     out << degree[i] << "\n";
    // }
    // out.close();



    // vector<double> time_trajectory({0.0});
    // vector<int> path({});
    // for (int i = 0; i < n; i++)
    // {
    //     cout << i << "\n" ;
    //     simulate_next_reaction simulation(network, psi);
    //     std::uniform_int_distribution<int> dist(0,size-1);

    //     node_t initial_infected = dist(engine);
    //    // path.push_back(initial_infected);
    //     simulation.add_infections({ std::make_pair(initial_infected, 0.0)});
    //     for (int i =0 ; i< size; i++) {
    //         auto point = simulation.step(engine);
    //         if (point.second != INFINITY) {
    //             time_trajectory.push_back(point.second);
    //             //path.push_back(point.first);
    //             continue;
    //         }
    //         break;
    //     }

    // }
    // cout << "sorting....\n";
    // sort(time_trajectory.begin(),time_trajectory.end());


    // out.open(string("average_trajectory.dat"));

    // for (int i =0; i<time_trajectory.size(); i++){
    //     out << time_trajectory[i] << " " << double(i/double(n)) << "\n";
    // }
    // out.close();
    
    // out.open(string("path.dat"));

    // for (int i =0; i<time_trajectory.size(); i++){
    //     out << path[i] << "\n";
    // }
    // out.close();


    return 0;
}

