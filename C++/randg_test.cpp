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
	// lower and upper bound
	if(argc < 4){
		cout << "Usage: order, num, p (0-10)" << endl;
		return 1;
	}
	const int order = atoi(argv[1]), num = atoi(argv[2]), p = atoi(argv[3]);
	// graph
	Graph *g = new Graph;
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
	double inf_zf_time = 0, inf_pt_time = 0, inf_th_time = 0;
	bool inf_zf_to = false, inf_pt_to = false, inf_th_to = false;
	double tsm_zf_time = 0, tsm_pt_time = 0, tsm_PT_time = 0, tsm_th_time = 0;
	bool tsm_zf_to = false, tsm_pt_to = false, tsm_PT_to = false, tsm_th_to = false;
	double fcm_zf_time = 0, fzf_time=0, amf_time=0, pti_time=0, ft_time=0;
	bool fcm_zf_to = false, fzf_to = false, amf_to = false, pti_to = false, ft_to = false;
	int count=0;
	time_point<high_resolution_clock> start, stop;
	microseconds duration;
	// output file
	std::cout.precision(3);
	ofstream myfile;
	// call geng for given order and edge range
	pipe = genrang_call(order,num,p);
	// for each non-isomorphic graph, find zero forcing, prop time, throttling number,
	// fractional zero forcing, all minimal forts, and fort number
	myfile.open("../csv_files/randg_test"+to_string(order)+to_string(num)+to_string(p)+".csv");
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
		if(!inf_zf_to){
			start = high_resolution_clock::now();
			infection_ip(g,g->order-1,infdata,'Z');
			stop = high_resolution_clock::now();
			duration = duration_cast<microseconds>(stop - start);
			inf_zf_time += duration.count()*1E-6;
			myfile << ", " << to_string(infdata.val);
			delete infdata.zf_set;
			delete infdata.forcings;
			if(infdata.status == GRB_TIME_LIMIT){
				inf_zf_to = true;
			}
		}
		// zero forcing (FCM)
		if(!fcm_zf_to){
			start = high_resolution_clock::now();
			fort_cover_ip(g,fcdata);
			stop = high_resolution_clock::now();
			duration = duration_cast<microseconds>(stop - start);
			fcm_zf_time += duration.count()*1E-6;
			delete fcdata.zf_set;
			if(fcdata.status == GRB_TIME_LIMIT){
				fcm_zf_to = true;
			}
		}
		// zero forcing (TSM)
		if(!tsm_zf_to){
			start = high_resolution_clock::now();
			time_step_ip(g,g->order-1,tsdata,'Z');
			stop = high_resolution_clock::now();
			duration = duration_cast<microseconds>(stop - start);
			tsm_zf_time += duration.count()*1E-6;
			delete tsdata.zf_set;
			delete tsdata.forcings;
			if(tsdata.status == GRB_TIME_LIMIT){
				tsm_zf_to = true;
			}
		}
		// minimum propogation time (IM)
		if(!inf_pt_to){
			start = high_resolution_clock::now();
			infection_ip(g,g->order-1,infdata,'p');
			stop = high_resolution_clock::now();
			duration = duration_cast<microseconds>(stop - start);
			inf_pt_time += duration.count()*1E-6;
			myfile << ", " << to_string(infdata.val);
			delete infdata.zf_set;
			delete infdata.forcings;
			if(infdata.status == GRB_TIME_LIMIT){
				inf_pt_to = true;
			}
		}
		// minimum propogation time (TSM)
		if(!tsm_pt_to){
			start = high_resolution_clock::now();
			time_step_ip(g,g->order-1,tsdata,'p');
			stop = high_resolution_clock::now();
			duration = duration_cast<microseconds>(stop - start);
			tsm_pt_time += duration.count()*1E-6;
			delete tsdata.zf_set;
			delete tsdata.forcings;
			if(tsdata.status == GRB_TIME_LIMIT){
				tsm_pt_to = true;
			}
		}
		// maximum propogation time (TSM)
		if(!tsm_PT_to){
			start = high_resolution_clock::now();
			time_step_ip(g,g->order-1,tsdata,'P');
			stop = high_resolution_clock::now();
			duration = duration_cast<microseconds>(stop - start);
			tsm_PT_time += duration.count()*1E-6;
			myfile << ", " << to_string(tsdata.val);
			delete tsdata.zf_set;
			delete tsdata.forcings;
			if(tsdata.status == GRB_TIME_LIMIT){
				tsm_PT_to = true;
			}
		}
		// throttling number (IM)
		if(!inf_th_to){
			start = high_resolution_clock::now();
			infection_ip(g,g->order-1,infdata,'T');
			stop = high_resolution_clock::now();
			duration = duration_cast<microseconds>(stop - start);
			inf_th_time += duration.count()*1E-6;
			myfile << ", " << to_string(infdata.val);
			delete infdata.zf_set;
			delete infdata.forcings;
			if(infdata.status == GRB_TIME_LIMIT){
				inf_th_to = true;
			}
		}
		// throttling number (TSM)
		if(!tsm_th_to){
			start = high_resolution_clock::now();
			time_step_ip(g,g->order-1,tsdata,'T');
			stop = high_resolution_clock::now();
			duration = duration_cast<microseconds>(stop - start);
			tsm_th_time += duration.count()*1E-6;
			delete tsdata.zf_set;
			delete tsdata.forcings;
			if(tsdata.status == GRB_TIME_LIMIT){
				tsm_th_to = true;
			}
		}
		// propagation time interval
		if(!pti_to){
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
			if(ptidata.status == GRB_TIME_LIMIT){
				pti_to = true;
			}
		}
		// fractional zero forcing
		if(!fzf_to){
			start = high_resolution_clock::now();
			fzf_ip(g,fzfdata);
			stop = high_resolution_clock::now();
			duration = duration_cast<microseconds>(stop - start);
			fzf_time += duration.count()*1E-6;
			myfile << ", " << to_string(fzfdata.val);
			delete fzfdata.weights;
			if(fzfdata.status == GRB_TIME_LIMIT){
				fzf_to = true;
			}
		}
		// all minimal forts
		if(!amf_to){
			start = high_resolution_clock::now();
			all_minimal_forts(g,amfdata);
			stop = high_resolution_clock::now();
			duration = duration_cast<microseconds>(stop - start);
			amf_time += duration.count()*1E-6;
			myfile << ", " << to_string(amfdata.mforts->size());
			delete amfdata.mforts;
			if(amfdata.status == GRB_TIME_LIMIT){
				amf_to = true;
			}
		}
		// fort number
		if(!ft_to){
			start = high_resolution_clock::now();
			ft_num_ip(g,ftdata);
			stop = high_resolution_clock::now();
			duration = duration_cast<microseconds>(stop - start);
			ft_time += duration.count()*1E-6;
			myfile << ", " << to_string(ftdata.val);
			delete ftdata.dforts;
			if(ftdata.status == GRB_TIME_LIMIT){
				ft_to = true;
			}
		}
		// new line in myfile
		myfile << endl;
		// delete graph and clear graph6 string
		delete g;
		str.clear();
	}
	// report timings
	if(inf_zf_to){
		cout << "inf_zf_time: TIME LIMIT" << endl;
	}
	else{
		cout << "inf_zf_time: " << inf_zf_time/count << endl;
	}
	if(fcm_zf_to){
		cout << "fcm_zf_time: TIME LIMIT" << endl;
	}
	else{
		cout << "fcm_zf_time: " << fcm_zf_time/count << endl;
	}
	if(tsm_zf_to){
		cout << "tsm_zf_time: TIME LIMIT" << endl;
	}
	else{
		cout << "tsm_zf_time: " << tsm_zf_time/count << endl;
	}
	if(inf_pt_to){
		cout << "inf_pt_time: TIME LIMIT" << endl;
	}
	else{
		cout << "inf_pt_time " << inf_pt_time/count << endl;
	}
	if(tsm_pt_to){
		cout << "tsm_pt_time: TIME LIMIT" << endl;
	}
	else{
		cout << "tsm_pt_time " << tsm_pt_time/count << endl;
	}
	if(tsm_PT_to){
		cout << "tsm_PT_time: TIME LIMIT" << endl;
	}
	else{
		cout << "tsm_PT_time " << tsm_PT_time/count << endl;
	}
	if(inf_th_to){
		cout << "inf_th_time: TIME LIMIT" << endl;
	}
	else{
		cout << "inf_th_time: " << inf_th_time/count << endl;
	}
	if(tsm_th_to){
		cout << "tsm_th_time: TIME LIMIT" << endl;
	}
	else{
		cout << "tsm_th_time: " << tsm_th_time/count << endl;
	}
	if(pti_to){
		cout << "pti_time: TIME LIMIT" << endl;
	}
	else{
		cout << "pti_time: " << pti_time/count << endl;
	}
	if(fzf_to){
		cout << "fzf_time: TIME LIMIT" << endl;
	}
	else{
		cout << "fzf_time: " << fzf_time/count << endl;
	}
	if(amf_to){
		cout << "amf_time: TIME LIMIT" << endl;
	}
	else{
		cout << "amf_time: " << amf_time/count << endl;
	}
	if(ft_to){
		cout << "ft_time: TIME LIMIT" << endl;
	}
	else{
		cout << "ft_time: " << ft_time/count << endl;
	}
	// close file
	myfile.close();
	// close pipe
	status = pclose(pipe);
	return 0;
}