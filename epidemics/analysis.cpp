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
#include "algorithm.h"
#include "simulation.h"

#include "analysis.h"

using namespace std;

void exportData(vector<double>& X,vector<double>& Y,string filename) {
    ofstream out;
    out.open(filename);
    int n = (int) X.size();
    if (n != (int) Y.size()) throw std::logic_error("Arrays are not the same size.");
    for (int i =0; i<n; i++){
        out << X[i] <<" " << Y[i] <<"\n";
    }
    out.close();
}

void exportData(vector<int>& X,vector<double>& Y,string filename) {
    ofstream out;
    out.open(filename);
    int n = (int) X.size();
    if (n != (int) Y.size()) throw std::logic_error("Arrays are not the same size.");
    for (int i =0; i<n; i++){
        out << X[i] <<" " << Y[i] <<"\n";
    }
    out.close();
}


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
    std::size_t n = (int) adjacencyList.size();

    for (std::size_t i =0; i<n; i++){
        for (std::size_t j =0; j < n-1; j++){
            //handle case if k_i = 0
            if ((int) adjacencyList[i].size() == 0){
                for (std::size_t val = 0; val < n-1; val++)
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

void measure_runtime(rng_t& engine, simulation_factory_t factory, int SMAX, double TMAX, string filename)
{
	//fill sizes
	vector<int> sizes;
	for (int i = 50; i < 15*1e6;)
	{
		sizes.push_back(i);
		i = (int) floor(pow(i,1.1));
	}
	
	//fill times
	vector<double> times;
	for (int N : sizes){
		cout << N << "\n";
		double average = 0;
		for (int i = 0; i < SMAX; i++)
		{
			simulation_algorithm* simulation = factory(engine, N);

			// Start measuring performance
			auto start = std::chrono::high_resolution_clock::now();
	
			// 10% of the network will be infected at t=0;
			const int N0 = floor(simulation->get_network().nodes() / 10);
			for (node_t node = 0; node < N0; node++)
			{
				simulation->add_infections({ std::make_pair(node, 0.0)});
			}
			
			// Run simulation, collect transmission times
			while (true) {

				auto point = simulation->step(engine);
				if (!point || (point -> time > TMAX))
					break;
			}
		
			auto stop = std::chrono::high_resolution_clock::now();
			
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
			const double time = duration.count() ;
			average += time / SMAX;
		}
		times.push_back(average);
	}
		// export data
	exportData(sizes,times,filename);
}

void measure_running_time_next_reaction_ER(rng_t& engine,int SMAX, string filename){
    const double MEAN_INFECTION = 10;
    const double VARIANCE_INFECTION = 1.0;
    const double MEAN_RECOVERY = 20;
    const double VARIANCE_RECOVERY = 1;
    const double TMAX = INFINITY;
    const double R0 = 3;
    transmission_time_lognormal psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_lognormal rho(MEAN_RECOVERY, VARIANCE_RECOVERY);


    //fill sizes
    vector<int> sizes;
    for (int i = 50; i < 15*1e6;)
    {
        sizes.push_back(i);
        i = (int) floor(pow(i,1.1));
    }
    
    //fill times
    vector<double> times;
    for (int N : sizes){
        cout << N << "\n";
        double average = 0;
        for (int i = 0; i < SMAX; i++)
        {

            erdos_reyni network(N,R0,engine);

            // Start measuring performance
            auto start = std::chrono::high_resolution_clock::now();

            simulate_next_reaction simulate(network,psi);

            const int N0 = floor(network.nodes() / 10); // 10% of the network will be infected at t=0;
            for (node_t node = 0; node < N0; node++)
            {
                simulate.add_infections({ std::make_pair(node, 0.0)});
            }
            
            // Run simulation, collect transmission times
            while (true) {

                auto point = simulate.step(engine);
                if (!point || (point -> time > TMAX))
                    break;
            }
        
            auto stop = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            const double time = duration.count() ;
            average += time / SMAX;
        }
        times.push_back(average);
    }
        // export data
    exportData(sizes,times,filename);
}

void measure_running_time_next_reaction_BA(rng_t& engine,int SMAX, string filename){
    const double MEAN_INFECTION = 10;
    const double VARIANCE_INFECTION = 1.0;
    const double MEAN_RECOVERY = 20;
    const double VARIANCE_RECOVERY = 1;
    const double TMAX = INFINITY;
    transmission_time_lognormal psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_lognormal rho(MEAN_RECOVERY, VARIANCE_RECOVERY);


    //fill sizes
    vector<int> sizes;
    for (int i = 50; i < 15*1e6;)
    {
        sizes.push_back(i);
        i = (int) floor(pow(i,1.1));
    }
    
    //fill times
    vector<double> times;
    for (int N : sizes){
        cout << N << "\n";
        double average = 0;
        for (int i = 0; i < SMAX; i++)
        {

            scale_free network(N,engine);

            // Start measuring performance
            auto start = std::chrono::high_resolution_clock::now();

            simulate_next_reaction simulate(network,psi);

            const int N0 = floor(network.nodes() / 10); // 10% of the network will be infected at t=0;
            for (node_t node = 0; node < N0; node++)
            {
                simulate.add_infections({ std::make_pair(node, 0.0)});
            }
            
            // Run simulation, collect transmission times
            while (true) {

                auto point = simulate.step(engine);
                if (!point || (point -> time > TMAX))
                    break;
            }
        
            auto stop = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            const double time = duration.count() ;
            average += time / SMAX;
        }
        times.push_back(average);
    }
        // export data
    exportData(sizes,times,filename);
}


void measure_running_time_nMGA_ER(rng_t& engine,int SMAX, string filename){
    const double MEAN_INFECTION = 10;
    const double VARIANCE_INFECTION = 1.0;
    const double MEAN_RECOVERY = 20;
    const double VARIANCE_RECOVERY = 1;
    const double TMAX = INFINITY;
    const double R0 = 3;
    transmission_time_lognormal psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_lognormal rho(MEAN_RECOVERY, VARIANCE_RECOVERY);


    //fill sizes
    vector<int> sizes;
    for (int i = 50; i < 2647355+1;)
    {
        sizes.push_back(i);
        i = (int) floor(pow(i,1.1));
    }
    
    //fill times
    vector<double> times;
    for (int N : sizes){
        cout << N << "\n";
        double average = 0;
        for (int i = 0; i < SMAX; i++)
        {

            erdos_reyni network(N,R0,engine);

            // Start measuring performance
            auto start = std::chrono::high_resolution_clock::now();

            simulate_nmga simulate(network,psi);

            const int N0 = floor(network.nodes() / 10); // 10% of the network will be infected at t=0;
            for (node_t node = 0; node < N0; node++)
            {
                simulate.add_infections({ std::make_pair(node, 0.0)});
            }
            
            // Run simulation, collect transmission times
            while (true) {

                auto point = simulate.step(engine);
                if (!point || (point -> time > TMAX))
                    break;
            }
        
            auto stop = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            const double time = duration.count() ;
            average += time / SMAX;
        }
        times.push_back(average);
    }
        // export data
    exportData(sizes,times,filename);
}

void measure_running_time_nMGA_BA(rng_t& engine,int SMAX, string filename){
    const double MEAN_INFECTION = 10;
    const double VARIANCE_INFECTION = 1.0;
    const double MEAN_RECOVERY = 20;
    const double VARIANCE_RECOVERY = 1;
    const double TMAX = INFINITY;
    transmission_time_lognormal psi(MEAN_INFECTION, VARIANCE_INFECTION);
    transmission_time_lognormal rho(MEAN_RECOVERY, VARIANCE_RECOVERY);


    //fill sizes
    vector<int> sizes;
    for (int i = 50; i < 2647355+1;)
    {
        sizes.push_back(i);
        i = (int) floor(pow(i,1.1));
    }
    
    //fill times
    vector<double> times;
    for (int N : sizes){
        cout << N << "\n";
        double average = 0;
        for (int i = 0; i < SMAX; i++)
        {

            scale_free network(N,engine);

            // Start measuring performance
            auto start = std::chrono::high_resolution_clock::now();

            simulate_nmga simulate(network,psi);

            const int N0 = floor(network.nodes() / 10); // 10% of the network will be infected at t=0;
            for (node_t node = 0; node < N0; node++)
            {
                simulate.add_infections({ std::make_pair(node, 0.0)});
            }
            
            // Run simulation, collect transmission times
            while (true) {

                auto point = simulate.step(engine);
                if (!point || (point -> time > TMAX))
                    break;
            }
        
            auto stop = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            const double time = duration.count() ;
            average += time / SMAX;
        }
        times.push_back(average);
    }
        // export data
    exportData(sizes,times,filename);
}

