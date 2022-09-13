//
//  savingData.cpp
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#include "stdafx.h"
#include "NextReaction.h"
#include "nMGA.h"
#include "random.h"
#include "graph.h"
#include "types.h"
#include "utility.h"
#include "simulation.h"


using namespace std;

void exportData(vector<double>& trajectory,string filename) {
    //create an ofstream
    ofstream out;
    //choose where to write
//    string location("Comparaison/");
    out.open(filename);
    int n = (int) trajectory.size();
   // out << "t" << ", " << "n"<< "\n";
//    out << trajectory[0] <<" " << 1 <<"\n";
    for (int i =0; i<n; i++){
//        if (abs(trajectory[i]-trajectory[i+1])<=0.001) {
//            out << trajectory[i]+0.0111 <<" " << i+1 <<"\n";
//            cout << "OUPS" << endl;
//            continue;
//        }
        out << trajectory[i] <<" " << i+1 <<"\n";
    }
    out.close();
}

void exportData(vector<int>& trajectory,string filename) {
    //create an ofstream
    ofstream out;
    //choose where to write
//    string location("Comparaison/");
    out.open(filename);
    int n = (int) trajectory.size();
   // out << "t" << ", " << "n"<< "\n";
//    out << trajectory[0] <<" " << 1 <<"\n";
    for (int i =0; i<n; i++){
//        if (abs(trajectory[i]-trajectory[i+1])<=0.001) {
//            out << trajectory[i]+0.0111 <<" " << i+1 <<"\n";
//            cout << "OUPS" << endl;
//            continue;
//        }
        out << trajectory[i] <<" " << i+1 <<"\n";
    }
    out.close();
}


void export_adjacency_list(std::vector<std::vector<node_t>>& adjacencyList,string filename) {
    //create an ofstream
    ofstream out;
    
    out.open(filename);
    int n = (int) adjacencyList.size();
    for (int i =0; i<n; i++){
        if ((int) adjacencyList[i].size() > 0) {
            for (int j =0; j < (int) adjacencyList[i].size() -1; j++){
                out << adjacencyList[i][j] << ", ";
            }
            out << adjacencyList[i].back() << "\n" ;
            continue;
        } else {
            out << "\n";
        }
    }
        

    out.close();
}

void export_adjacency_matrix(std::vector<std::vector<node_t>>& adjacencyList,string filename) {
    //create an ofstream
    ofstream out;
    
    

    out.open(filename);
    int n = (int) adjacencyList.size();

    for (int i =0; i<n; i++){
        
        for (int j =0; j < n-1; j++){

            //handle case if k_i = 0
            if ((int) adjacencyList[i].size() == 0){
                for (size_t val = 0; val < n-1; val++)
                    out << 0 << ", ";
                out << 0 << "\n";
                continue;
            } 

            // if i and j are neighbours A[i,j]=1 else 0.
            if(std::find(adjacencyList[i].begin(), adjacencyList[i].end(), j) != adjacencyList[i].end()) {
                out << 1 << ", ";
            } else {
                out << 0 << ", ";
            }
        }

        //handle last element:
        // if i and j are neighbours A[i,j]=1 else 0.
        if(std::find(adjacencyList[i].begin(), adjacencyList[i].end(), n-1) != adjacencyList[i].end()) {
            out << 1 << "\n";
        } else {
            out << 0 << "\n";
        }
    }
        

    out.close();
}



void print_matrix(const std::vector<std::vector<double>>& A){
    long rows = A.size();
    long cols = A[0].size();
    
    for (int i = 0; i<rows; i++) {
        for (int j =0; j<cols; j++) {
            std::cout << ceil(A[i][j]*10)/10 << "   ";
        }
        std::cout << std::endl;
    }
}


double measure_running_time(graph_adjacencylist& network,rng_t& engine){
    //long n0 = network.adjacencylist.size()* 3 / 10 ;
    // initial nb of infected is 30% of the network to avoid epidemics that die out too quickly.
    
    double mean = 10;
    double variance = 1.0;
    
    transmission_time_lognormal psi(mean, variance);

    simulate_nmga simulation(network, psi);
    simulation.approximation_threshold = 1000000;

//    simulate_next_reaction simulation(network, psi);

//    std::uniform_real_distribution<> dis(0,mean * 0.05);
    
    // Infect the first 30% of the network at early times
//    for (int i=0; i<=n0; i++) {
//        const double infection_time = dis(engine);
//        const node_t node = i;
//        simulation.add_infections({ std::make_pair(node, infection_time)});
//    }

    
    // Start measuring performance
    auto start = std::chrono::high_resolution_clock::now();

    // Epidemic spreading
    
    //generatePaths_next_reaction(mean,variance,degree,1,size, engine);

    auto stop = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    
    return (double) duration.count();
}

//simulate_nmga returnNMGA(graph_adjacencylist& network, transmission_time& psi){
//    simulate_nmga simulation(network, psi);
//    return simulation;
//}
//
//simulate_next_reaction returnNextReac(graph_adjacencylist& network, transmission_time& psi){
//    simulate_next_reaction simulation(network, psi);
//    return simulation;
//}

void generate_data_running_time(rng_t& engine,int size,bool isNMGA){

    vector<double> time_trajectory{};
    double mean = 5.0;
    double variance = 5.0;
    int degree = 3;
    
    const auto start = std::chrono::high_resolution_clock::now();
    
    if (isNMGA == true) {
        generatePaths_NMGA(mean,variance,degree,1,size,engine,100);
    } else {
        generatePaths_next_reaction(mean, variance,degree,1,size,engine);
    }
    
    
    auto stop = std::chrono::high_resolution_clock::now();
    
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    const double time = duration.count() / 1000;
    time_trajectory.push_back(time);
    cout <<size << " "<< time << "\n";
}
