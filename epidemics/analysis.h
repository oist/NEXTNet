#pragma once

#include "stdafx.h"
#include "types.h"

//ostream& file,

void exportData( std::vector<double>& trajectory,std::string filename);

void export_adjacency_list(std::vector<std::vector<std::pair<node_t,interval_t>>>& adjacencyList, std::string filename);

