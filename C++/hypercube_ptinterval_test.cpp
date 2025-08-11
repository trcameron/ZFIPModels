#include "ip_models.hpp"
#include <fstream>
#include <chrono>
using namespace std;
using namespace std::chrono;
/*
	main function
*/
int main(int argc,char *argv[]){
    // d (dimension)
	if(argc < 2){
		cout << "Usage: d (dimension)" << endl;
		return 1;
	}
	const int d = atoi(argv[1]);
    // hypercube graph
    Graph *g = hypercube_graph(d);
    // ip data variables
    ptiIPdata ptidata;
    // timing variables
    double pti_time = 0;
    // propagation time interval
	auto start = high_resolution_clock::now();
	pt_interval(g,pow(2,d-1),ptidata);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop-start);
	pti_time += duration.count()*1E-6;
    cout << "PT INTERVAL" << endl;
	for(auto it=ptidata.ptMap->cbegin(); it!=ptidata.ptMap->cend(); ++it){
        cout << it->first << ": ";
        for(auto u=it->second.cbegin(); u!=it->second.cend(); ++u){
            cout << *u << " ";
        }
        cout << endl;
	}
    // cdlete graph and ip data
    delete g;
    delete ptidata.ptMap;
    // print timings
    cout << "pti_time: " << pti_time << endl;
    // return
    return 0;
}