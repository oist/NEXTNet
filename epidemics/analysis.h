#pragma once

#include "stdafx.h"
#include "types.h"

//ostream& file,

void exportData( std::vector<double>& trajectory,std::string filename);

void export_adjacency_list(std::vector<std::vector<node_t>>& adjacencyList, std::string filename);



void print_matrix(std::vector<std::vector<double>>& A);
