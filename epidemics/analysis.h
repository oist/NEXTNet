#pragma once

#include "stdafx.h"
#include "algorithm.h"
#include "graph.h"

typedef std::function<simulation_algorithm* (rng_t& engine, int n)> simulation_factory_t;

void exportData( std::vector<double>& trajectory,std::string filename);

void exportData( std::vector<double>& X,std::vector<double>& Y,std::string filename);

void exportData( std::vector<int>& X,std::vector<double>& Y,std::string filename);

void exportData( std::vector<int>& trajectory,std::string filename);

void export_adjacency_list(std::vector<std::vector<node_t>>& adjacencyList, std::string filename);
void export_adjacency_matrix(std::vector<std::vector<node_t>>& adjacencyList, std::string filename);

void print_matrix(std::vector<std::vector<double>>& A);

double measure_running_time(graph_adjacencylist& network,rng_t& engine);

void generate_data_running_time(rng_t& engine,int size, bool isNMGA);

// void average_epidemic(graph_adjacencylist& network,rng_t& engine, )
// grjijffw

void measure_runtime(rng_t& engine, simulation_factory_t factory, int SMAX, double TMAX, std::string filename);

void measure_running_time_next_reaction_ER(rng_t& engine,int SMAX, std::string filename);
void measure_running_time_next_reaction_BA(rng_t& engine, int SMAX, std::string filename);
void measure_running_time_nMGA_ER(rng_t& engine, int SMAX, std::string filename);
void measure_running_time_nMGA_BA(rng_t& engine, int SMAX, std::string filename);

template<typename Factory>
void measure_runtime(rng_t& engine, Factory factory,
					 int SMAX, double TMAX, std::string filename)
{
	using namespace std;
	
    //fill sizes
    vector<int> sizes;
    for (int k = 7; k <= 20; k++)
    {
        int i = (int) floor(pow(2,k));
        sizes.push_back(i);
    }
	
	//fill times
	vector<double> times;
	for (int N : sizes){
		cout << N << "\n";
		double average = 0;
		for (int i = 0; i < SMAX; i++)
		{
			// Create simulation environment
			auto s = factory(engine, N);
			simulation_algorithm& simulation = *s.simulator;

			// Start measuring performance
			auto start = std::chrono::high_resolution_clock::now();
	
			// 10% of the network will be infected at t=0;
			const int N0 = floor(simulation.get_network().nodes() / 10);
			for (node_t node = 0; node < N0; node++)
			{
				simulation.add_infections({ std::make_pair(node, 0.0)});
			}
			
			// Run simulation, collect transmission times
			while (true) {

				auto point = simulation.step(engine);
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
