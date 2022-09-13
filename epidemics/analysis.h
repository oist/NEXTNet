#pragma once

#include "stdafx.h"
#include "graph.h"


void exportData( std::vector<double>& trajectory,std::string filename);

void exportData( std::vector<int>& trajectory,std::string filename);

void export_adjacency_list(std::vector<std::vector<node_t>>& adjacencyList, std::string filename);
void export_adjacency_matrix(std::vector<std::vector<node_t>>& adjacencyList, std::string filename);

void print_matrix(std::vector<std::vector<double>>& A);

double measure_running_time(graph_adjacencylist& network,rng_t& engine);

void generate_data_running_time(rng_t& engine,int size, bool isNMGA);

// void average_epidemic(graph_adjacencylist& network,rng_t& engine, )
// grjijffw



