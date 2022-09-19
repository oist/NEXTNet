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

//int main(int argc, const char * argv[]) {
//
//    engine.seed(1);
//
//    cout << "start...\n";
//
//    int size = 1000;
//
//    // number_of_simulations:
//    int n = 1;
//
//    // User's input
//    if(argc == 2){
//        size = atoi(argv[1]);
//    } else if (argc == 3){
//        size = atoi(argv[1]);
//        n = atoi(argv[2]);
//    }
//
//    double mean = 10;
//    double variance = 1.0;
//    transmission_time_lognormal psi(mean, variance);
//
//    scale_free network(size, engine);
//
//    cout << "network generated \n";
//    std::vector<long> degree;
//    for (int i=0; i<network.adjacencylist.size(); i++) {
//        degree.push_back(network.adjacencylist[i].size());
//    }
//
//    ofstream out;
//
//    out.open(string("degree.dat"));
//
//    for (int i =0; i<degree.size(); i++){
//        out << degree[i] << "\n";
//    }
//    out.close();
//
//
//
//    vector<double> time_trajectory({0.0});
//    vector<int> path({});
//    for (int i = 0; i < n; i++)
//    {
//        cout << i << "\n" ;
//        simulate_next_reaction simulation(network, psi);
//        std::uniform_int_distribution<int> dist(0,size-1);
//
//        node_t initial_infected = dist(engine);
//        path.push_back(initial_infected);
//        simulation.add_infections({ std::make_pair(initial_infected, 0.0)});
//        for (int i =0 ; i< size; i++) {
//            auto point = simulation.step(engine);
//            if (point.second != INFINITY) {
//                time_trajectory.push_back(point.second);
//                path.push_back(point.first);
//                continue;
//            }
//            break;
//        }
//
//    }
//    cout << "sorting....\n";
//    sort(time_trajectory.begin(),time_trajectory.end());
//
//
//    out.open(string("average_")+to_string(size)+string(".traj"));
//
//    for (int i =0; i<time_trajectory.size(); i++){
//        out << time_trajectory[i] << " " << double(i/double(n)) << "\n";
//    }
//    out.close();
//
//    out.open(string("path.dat"));
//
//    for (int i =0; i<time_trajectory.size(); i++){
//        out << path[i] << "\n";
//    }
//    out.close();
//
//
//    return 0;
//}
//



int main(int argc, const char * argv[]) {

    engine.seed(34);
    int size = 50000;
    double r_aim = 0.4;
    double r;
    bool assor = false;

    // number_of_simulations:
    int nb_sim = 50;
    string filename("average_0.traj");
    // User's input
    if(argc == 2){
        size = atoi(argv[1]);
    } else if (argc == 3){
        size = atoi(argv[1]);
        nb_sim = atoi(argv[2]);
    } else if (argc == 4){
        size = atoi(argv[1]);
        nb_sim = atoi(argv[2]);
        r_aim = stod(argv[3]);
        assor = true;
    } else if (argc == 5){
        size = atoi(argv[1]);
        nb_sim = atoi(argv[2]);
        r_aim = stod(argv[3]);
        assor = true;
        filename = string(argv[4]);
    }
//
    //erdos_reyni nw(size,4,engine);

    imported_network nw(string("/Users/curesamuelcyrus/log-normal_2.adj"));

    cout << "old assor : " << assortativity(nw) << "\n";
    if (assor){
        add_correlation(r_aim,nw,engine);
        r = assortativity(nw);
        cout << "new assor : " << r << "\n";
    }
//    if (assor==true)
//        return 0;

    export_adjacency_list(nw.adjacencylist,"/Users/curesamuelcyrus/desktop/simulations/network.adj");

    double mean = 10;
    double variance = 1.0;
    transmission_time_lognormal psi(mean, variance);


    std::vector<long> degree;
    for (long i=0; i< (long) nw.adjacencylist.size(); i++) {
        degree.push_back( (long) nw.adjacencylist[i].size());
    }

    ofstream out;

    out.open(string("degree.dat"));

    for (long i =0; i<degree.size(); i++){
        out << degree[i] << "\n";
    }
    out.close();

    vector<double> time_trajectory({0.0});
    // vector<int> path({});
    for (int i = 0; i < nb_sim; i++)
    {
        cout << i << "\n" ;
        simulate_next_reaction simulation(nw, psi);
        std::uniform_int_distribution<int> dist(0,size-1);

        node_t initial_infected = dist(engine);
       // path.push_back(initial_infected);
        simulation.add_infections({ std::make_pair(initial_infected, 0.0)});
        for (int i =0 ; i< size; i++) {
            auto point = simulation.step(engine);
            if (point.second != INFINITY) {
                time_trajectory.push_back(point.second);
                //path.push_back(point.first);
                continue;
            }
            break;
        }

    }
    cout << "sorting....\n";
    sort(time_trajectory.begin(),time_trajectory.end());


    out.open(filename);

    for (long i =0; i<time_trajectory.size(); i++){
        out << time_trajectory[i] << " " << double(i/double(nb_sim)) << "\n";
    }
    out.close();

    // out.open(string("path.traj"));

    // for (int i =0; i<time_trajectory.size(); i++){
    //     out << path[i] << "\n";
    // }
    // out.close();


    return 0;
}

