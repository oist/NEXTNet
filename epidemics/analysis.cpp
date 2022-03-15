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

//
//vector<double> interpolate(vector<double>& t_vec,vector<double>& n_vec,double dt){
//    
//    double nb_points = t_vec.back() / dt ;
//    
//    vector<double> result(nb_points,0);
//    
//    double t_min = t_vec[0];
//    double t_max = t_vec.back();
////    double n_min = n_vec[0];
////    double n_max = n_vec.back();
//    
//    int pos = 0;
//    
//    for (int i=0; i<nb_points; i++) {
//        
//        double t =  t_min + i * (t_max-t_min) / nb_points ;
//        
//        while (t > t_vec[pos+1])
//            pos++ ;
//        
//        
//        
//        double t0 = t_vec[pos];
//        double n0 = n_vec[pos];
//
//        double t1 = t_vec[pos+1];
//        double n1 = n_vec[pos+1];
//        
//        
//        
//        
//    }
//
//
//
//    return result;
//}


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




