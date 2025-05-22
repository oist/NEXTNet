//
//  graph.hpp
//  epidemics
//
//  Created by Cure Samuel Cyrus on 2022/03/07.
//

#pragma once

#include "nextnet/stdafx.h"
#include "nextnet/types.h"
#include "nextnet/network.h"

//---------------------------------------------------
//-----Output networks-------------------------------
//---------------------------------------------------

void output_adjacencylist(std::ostream &dst, network &nw, bool include_weights = true, bool include_coords = true,
                          bool include_meta = true, bool include_header = true, char csep = ' ', char dsep = ',', char wsep = ':');

void output_network_meta(std::ostream &dst, network &nw, char dsep = ',');
