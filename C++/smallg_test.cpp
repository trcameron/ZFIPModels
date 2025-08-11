#include "ip_models.hpp"
#include "testing.hpp"
#include <chrono>
#include <fstream>
using namespace std;
using namespace std::chrono;
/*
	main function
*/
int main(int argc,char *argv[]){
	// lower and upper bound
	if(argc < 3){
		cout << "Usage: order and option" << endl;
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
	fcIPdata fcdata;
	infIPdata infdata;
	fzfIPdata fzfdata;
	tsIPdata tsdata;
	ftIPdata ftdata;
	amfIPdata amfdata;
	ptiIPdata ptidata;
	// timing variables
	double inf_zf_time=0, inf_pt_time=0, inf_th_time=0;
	double tsm_zf_time = 0, tsm_pt_time = 0, tsm_PT_time = 0, tsm_th_time = 0;
	double fcm_zf_time = 0, fzf_time=0, amf_time=0, pti_time=0, ft_time=0;
	int count=0;
	time_point<high_resolution_clock> start, stop;
	microseconds duration;
	// output file
	std::cout.precision(3);
	ofstream myfile;
	// call geng for given order
	pipe = geng_call(order,opt);
	// for each non-isomorphic graph, find zero forcing, prop time, throttling number,
	// fractional zero forcing, all minimal forts, and fort number
	myfile.open("../csv_files/smallg_test"+to_string(order)+opt+".csv");
	myfile << "graph, Z, pt, PT, th, PTI, Z*, # min forts, ft" << endl;
	while(fgets(line, PATH_MAX, pipe) != NULL){
		// new graph
		g = new Graph;
		// graph6 string
		for(unsigned int i=0; i<sizeof_char(line); ++i){
			str += line[i];
		}
		g->read_graph6(str);
		myfile << str;
		count += 1;
		// zero forcing (IM)
		start = high_resolution_clock::now();
		infection_ip(g,g->order-1,infdata,'Z');
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		inf_zf_time += duration.count()*1E-6;
		myfile << ", " << to_string(infdata.val);
		delete infdata.zf_set;
		delete infdata.forcings;
		// zero forcing (FCM)
		start = high_resolution_clock::now();
		fort_cover_ip(g,fcdata);
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		fcm_zf_time += duration.count()*1E-6;
		delete fcdata.zf_set;
		// zero forcing (TSM)
		start = high_resolution_clock::now();
		time_step_ip(g,g->order-1,tsdata,'Z');
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		tsm_zf_time += duration.count()*1E-6;
		delete tsdata.zf_set;
		delete tsdata.forcings;
		// minimum propogation time (IM)
		start = high_resolution_clock::now();
		infection_ip(g,g->order-1,infdata,'p');
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		inf_pt_time += duration.count()*1E-6;
		myfile << ", " << to_string(infdata.val);
		delete infdata.zf_set;
		delete infdata.forcings;
		// minimum propogation time (TSM)
		start = high_resolution_clock::now();
		time_step_ip(g,g->order-1,tsdata,'p');
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		tsm_pt_time += duration.count()*1E-6;
		delete tsdata.zf_set;
		delete tsdata.forcings;
		// maximum propogation time (TSM)
		start = high_resolution_clock::now();
		time_step_ip(g,g->order-1,tsdata,'P');
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		tsm_PT_time += duration.count()*1E-6;
		myfile << ", " << to_string(tsdata.val);
		delete tsdata.zf_set;
		delete tsdata.forcings;
		// throttling number (IM)
		start = high_resolution_clock::now();
		infection_ip(g,g->order-1,infdata,'T');
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		inf_th_time += duration.count()*1E-6;
		myfile << ", " << to_string(infdata.val);
		delete infdata.zf_set;
		delete infdata.forcings;
		// throttling number (TSM)
		start = high_resolution_clock::now();
		time_step_ip(g,g->order-1,tsdata,'T');
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		tsm_th_time += duration.count()*1E-6;
		delete tsdata.zf_set;
		delete tsdata.forcings;
		// propagation time interval
		start = high_resolution_clock::now();
		pt_interval(g,g->order-1,ptidata);
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop-start);
		pti_time += duration.count()*1E-6;
		myfile << ", ";
		for(auto iter=ptidata.ptMap->begin(); iter!=ptidata.ptMap->end(); ++iter){
			myfile << to_string(iter->first) << " ";
		}
		delete ptidata.ptMap;
		// fractional zero forcing
		start = high_resolution_clock::now();
		fzf_ip(g,fzfdata);
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		fzf_time += duration.count()*1E-6;
		myfile << ", " << to_string(fzfdata.val);
		delete fzfdata.weights;
		// all minimal forts
		start = high_resolution_clock::now();
		all_minimal_forts(g,amfdata);
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		amf_time += duration.count()*1E-6;
		myfile << ", " << to_string(amfdata.mforts->size());
		delete amfdata.mforts;
		// fort number
		start = high_resolution_clock::now();
		ft_num_ip(g,ftdata);
		stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
		ft_time += duration.count()*1E-6;
		myfile << ", " << to_string(ftdata.val) << endl;
		delete ftdata.dforts;
		// delete graph and clear graph6 string
		delete g;
		str.clear();
	}
	// report timings
	cout << "inf_zf_time: " << inf_zf_time/count << endl;
	cout << "fcm_zf_time: " << fcm_zf_time/count << endl;
	cout << "tsm_zf_time: " << tsm_zf_time/count << endl;
	cout << "inf_pt_time " << inf_pt_time/count << endl;
	cout << "tsm_pt_time " << tsm_pt_time/count << endl;
	cout << "tsm_PT_time " << tsm_PT_time/count << endl;
	cout << "inf_th_time: " << inf_th_time/count << endl;
	cout << "tsm_th_time: " << tsm_th_time/count << endl;
	cout << "pti_time: " << pti_time/count << endl;
	cout << "z*_time: " << fzf_time/count << endl;
	cout << "amf_time: " << amf_time/count << endl;
	cout << "ft_time: " << ft_time/count << endl;
	// close file
	myfile.close();
	// close pipe
	status = pclose(pipe);
	return 0;
}