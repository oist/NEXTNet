#pragma once

#include "stdafx.h"
#include "algorithm.h"
#include "graph.h"

typedef std::function<simulation_algorithm* (rng_t& engine, int n)> simulation_factory_t;

void exportData( std::vector<double>& trajectory,std::string filename);

void exportData( std::vector<double>& X,std::vector<double>& Y,std::string filename);

void exportData( std::vector<double>& X,std::vector<double>& Y,std::string filename);

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
