#include "testing.hpp"
#include "ip_models.hpp"
#include <fstream>
#include <chrono>
using namespace std;
using namespace std::chrono;
/*
	main function
*/
int main(int argc,char *argv[]){
	// order
	if(argc < 3){
		cout << "Usage: gentreeg order and option" << endl;
		return 1;
	}
	const int order = atoi(argv[1]);
	const char *opt = argv[2];
	// graph variable 
	Graph *g;
	// file stream variables
	int status;
	FILE* pipe;
	char line[PATH_MAX];
	string str;
	// testing variables
	amfIPdata amfdata;
	int max = 0;
	vector<int> diam;
	vector<int> deg;
	vector<string> sparse6;
	// timing variables
	int count = 0;
	double amf_time = 0;
	// call gentreeg for given order
	pipe = gentreeg_call(order,opt);
	// for each non-isomorphic tree
	while(fgets(line, PATH_MAX, pipe) != NULL){
		// new graph
		g = new Graph;
		// sparse6 string
		for(unsigned int i=0; i<sizeof_char(line); ++i){
			str += line[i];
		}
		g->read_sparse6(str);
		count += 1;
		// all minimal forts
		auto start = high_resolution_clock::now();
		all_minimal_forts(g,amfdata);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop-start);
		amf_time += duration.count()*1E-6;
		// store diameter, degree, and graph string
		if(amfdata.mforts->size() > max){
			max = amfdata.mforts->size();
			diam.clear();
			deg.clear();
			sparse6.clear();
			diam.push_back(g->tree_diameter());
			deg.push_back(g->max_degree());
			sparse6.push_back(str);
			cout << "max: " << max << ", sparse6: " << str << endl;
		}
		else if(amfdata.mforts->size()==max){
			diam.push_back(g->tree_diameter());
			deg.push_back(g->max_degree());
			sparse6.push_back(str);
		}
		// delete amfdata and graph and clear string
		delete amfdata.mforts;
		delete g;
		str.clear();
	}
	cout << "amf_time: " << amf_time/count << endl;
	cout << "maximum number of minimal forts: " << max << endl;
	cout << "diameter: ";
	for(auto i=diam.cbegin(); i!=diam.cend(); ++i){
		cout << *i << " ";
	}
	cout << endl;
	cout << "max degree: ";
	for(auto i=deg.cbegin(); i!=deg.cend(); ++i){
		cout << *i << " ";
	}
	cout << endl;
	cout << "graph6: ";
	for(auto i=sparse6.cbegin(); i!=sparse6.cend(); ++i){
		cout << *i << " ";
	}
	cout << endl;
	return 0;
}