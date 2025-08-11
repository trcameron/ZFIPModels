#include "ip_models.hpp"
#include <fstream>
using namespace std;
/*
	main function
*/
int main(int argc,char *argv[]){
    // file stream variables
	int status;
	FILE* pipe;
	char line[PATH_MAX];
	// testing variables
	fcIPdata fcdata;
	fzfIPdata fzfdata;
	ftIPdata ftdata;
	// output file
	std::cout.precision(3);
	ofstream myfile;
    // for each bhhs graph, find zero forcing, prop time, throttling number,
	// fractional zero forcing, all minimal forts, and fort number
	myfile.open("../csv_files/bhhs_test.csv");
	myfile << "G_{i}, ft, Z*, M, Z" << endl;
    for(int i=1; i<=7; ++i){
        // graph
		Graph *g = new Graph;
        g->read_edge("../test_graphs/BHHS_graphs/E"+to_string(i)+".edg");
		myfile << to_string(i);
		// fort number
		ft_num_ip(g,ftdata);
		myfile << ", " + to_string(ftdata.val);
		delete ftdata.dforts;
		// fractional zero forcing
		fzf_ip(g,fzfdata);
		myfile << ", " + to_string(fzfdata.val);
		delete fzfdata.weights;
        // zero forcing (FCM)
		fort_cover_ip(g,fcdata);
        myfile << ", " + to_string(fcdata.val-1) + ", " + to_string(fcdata.val) + "\n";
		delete fcdata.zf_set;
		// clear graph
		delete g;
    }
    // close file
	myfile.close();
    // return
    return 0;
}