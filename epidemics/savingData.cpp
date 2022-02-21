//
//  savingData.cpp
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#include "stdafx.h"
#include "savingData.h"

using namespace std;

void exportData(vector<double>& trajectory,string filename) {
    //create an ofstream
    ofstream out;
    //choose where to write
    out.open(filename);
    // string path = "/Users/samuelcure/Documents/Epidemics-Mean-Field/data/";
    int n = (int) trajectory.size();
    out << "t" << ", " << "n"<< "\n";
    for (int i =0; i<n; i++)
        out << trajectory[i] <<", " << i+1 <<"\n";

    out.close();
}
