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
#include "gnuplot-iostream.h"
using namespace std;

/*  MersenneTwister random number generator */
/*static*/ mt19937 engine;


int main(int argc, const char * argv[]) {

    engine.seed(1);
    cout << "start...\n";

    int size = 1000000;

    double mean = 10;
    double variance = 1.0;
    transmission_time_lognormal psi(mean, variance);
    
    scale_free network(size, engine);

    cout << "network generated \n";   
    std::vector<long> degree;
    for (int i=0; i<network.adjacencylist.size(); i++) {
        degree.push_back(network.adjacencylist[i].size());
    }

    ofstream out;

    out.open(string("degree.dat"));
 
    for (int i =0; i<degree.size(); i++){
        out << degree[i] << "\n";
    }
    out.close();

    // number_of_simulations:
    int n = 1000;

    vector<double> time_trajectory({0.0});
    for (int i = 0; i < n; i++)
    {
        cout << i << "\n" ;
        simulate_next_reaction simulation(network, psi);
        std::uniform_int_distribution<int> dist(0,size-1);

        node_t initial_infected = dist(engine);
        simulation.add_infections({ std::make_pair(initial_infected, 0.0)});
        for (int i =0 ; i< size; i++) {
            auto point = simulation.step(engine);
            if (point.second != INFINITY) {

                auto it = std::upper_bound(time_trajectory.cbegin(),time_trajectory.cend(),point.second);
                time_trajectory.insert(it,point.second);
                continue;
            }
            break;
        }

    }

    out.open(string("average_trajectory.dat"));
 
    for (int i =0; i<time_trajectory.size(); i++){
        out << time_trajectory[i] << " " << double(1/n) << "\n";
    }
    out.close();
    
//     cout << "network generated \n";
    
//     simulate_next_reaction simulation(network, psi);
    
//     vector<double> time_trajectory({});
    
//     simulation.add_infections({ std::make_pair(0, 0.0)});
    
//     for (int i =0 ; i< size; i++) {
//         auto point = simulation.step(engine);
//         if (point.second != INFINITY) {
//             time_trajectory.push_back(point.second);
//             continue;
//         }
//         cout << time_trajectory.size() << "\n";
//         break;
//     }
//     std::string filename ="BA";
//     std::string ext= ".dat";
//     exportData(time_trajectory,path+filename+ext);

    
    
//    std::ifstream fin(path+string("degreeList.dat"));
//
//    std::vector<int> degreeList;
//
//    int element;
//    while (fin >> element)
//    {
//        degreeList.push_back(element);
//    }
//
//    cout << "numbers generated \n";
//    config_model network(degreeList, engine);
//    cout << "network generated \n";
//    //erdos_reyni network(size,degree, engine);
//    transmission_time_lognormal psi(mean, variance);
//    simulate_next_reaction simulation(network, psi);
//
//    vector<double> time_trajectory({});
//    vector<double> trajectory_degree_one({});
//    vector<double> trajectory_degree_four({});
//    vector<double> trajectory_degree_ten({});
//    vector<double> trajectory_degree_seven({});
//
//
//    for (node_t node = 0; node < int(100); node ++) {
//        simulation.add_infections({ std::make_pair(node, 0.0)});
//    }
////    simulation.add_infections({ std::make_pair(0, 0.0)});
//
//
//    for (int i =0 ; i< size; i++) {
//        auto point = simulation.step(engine);
//        if (point.second != INFINITY) {
//            const int degree = (int) network.adjacencylist[point.first].size();
//            if (degree== 1) {
//                trajectory_degree_one.push_back(point.second);
//            } else if (degree == 10) {
//                trajectory_degree_ten.push_back(point.second);
//            } else if (degree == 4) {
//                trajectory_degree_four.push_back(point.second);
//            } else if (degree == 7) {
//                trajectory_degree_seven.push_back(point.second);
//            }
//
//            time_trajectory.push_back(point.second);
//            continue;
//        }
//        cout << time_trajectory.size() << "\n";
//        break;
//    }
//
//    std::string filename ="nmgaAugust";
//    std::string ext= ".dat";
//    exportData(time_trajectory,path+filename+ext);
//    exportData(trajectory_degree_one,path+string("rate_one")+ext);
//    exportData(trajectory_degree_ten,path+string("rate_ten")+ext);
//    exportData(trajectory_degree_four,path+string("rate_four")+ext);
//    exportData(trajectory_degree_seven,path+string("rate_seven")+ext);
//    exportData(degreeList, path+string("degreelist.txt"));
//    export_adjacency_list(network.adjacencylist,path+string("mat.dat"));
//    cout << " done \n";

    return 0;
}

