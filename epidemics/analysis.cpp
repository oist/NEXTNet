//
//  savingData.cpp
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#include "stdafx.h"
#include "analysis.h"
#include "types.h"


using namespace std;

void exportData(vector<double>& trajectory,string filename) {
    //create an ofstream
    ofstream out;
    //choose where to write
//    string location("Comparaison/");
    out.open(filename);
    int n = (int) trajectory.size();
   // out << "t" << ", " << "n"<< "\n";
    for (int i =0; i<n; i++)
        out << trajectory[i] <<" " << i+1 <<"\n";

    out.close();
}

void export_adjacency_list(vector<vector<pair<node_t,interval_t>>>& adjacencyList,string filename) {
    //create an ofstream
    ofstream out;
    
    out.open(filename);
    int n = (int) adjacencyList.size();
    for (int i =0; i<n; i++){
        if ((int) adjacencyList[i].size() > 0) {
            for (int j =0; j < (int) adjacencyList[i].size() -1; j++){
                out << adjacencyList[i][j].first << ", ";
            }
            out << adjacencyList[i].back().first << "\n" ;
            continue;
        } else {
            out << "\n";
        }
    }
        

    out.close();
}



