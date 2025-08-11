#include "ip_models.hpp"
#include <cstdlib>
#include <fstream>
using namespace std;
/*
    rand_int: random integer between a and b
*/
int rand_int(const int a,const int b){
    return rand() % (b - a + 1) + a;
}
/*
    min_rank_vs: minimum rank of vertex sum with graphs of order 8 whose minimum rank is known
*/
int min_rank_vs(const Graph *g,const Graph *h,const int u,const int v,const int mrg,const int mrh){
    // fort cover ip data
    fcIPdata fcdata;
    // g-u, min rank and rank spread
    Graph *vd1 = vert_del(g,u);
    fort_cover_ip(vd1,fcdata);
    int mr1 = vd1->order - fcdata.val;
    int rs1 = mrg - mr1;
    delete vd1;
    delete fcdata.zf_set;
    // h-v, min rank and rank spread
    Graph *vd2 = vert_del(h,v);
    fort_cover_ip(vd2,fcdata);
    int mr2 = vd2->order - fcdata.val;
    int rs2 = mrh - mr2;
    delete vd2;
    delete fcdata.zf_set;
    // sum of rank spreads
    int s = rs1 + rs2;
    if(s > 2){
        s = 2;
    }
    // return minimum rank
    return mr1 + mr2 + s;
}
/*
    max_null_vs: maximum nullity of vertex sum with graphs of order 8 whose maximum nullity is known
*/
int max_null_vs(const Graph *g,const Graph *h,const int u,const int v,const int mng, const int mnh){
    // return maximum nullity
    return (g->order + h->order - 1) - min_rank_vs(g,h,u,v,g->order-mng,h->order-mnh);
}
/*
	main function
*/
int main(int argc,char *argv[]){
    // seed random variable
    srand(time(NULL));
    // file stream variables
	int status;
	FILE* pipe;
	char line[PATH_MAX];
	// testing variables
	fcIPdata fcdata;
	fzfIPdata fzfdata;
	ftIPdata ftdata;
    int count = 0, mng, mnh;
	// output file
	std::cout.precision(3);
	ofstream myfile;
    // for each bhhs graph, find zero forcing, prop time, throttling number,
	// fractional zero forcing, all minimal forts, and fort number
	myfile.open("../csv_files/bhhs_vs_test.csv");
	myfile << "Count, G_{i}, G_{j}, u, v, ft, Z*, M, Z" << endl;
    // graph g_{i}
    for(int i=1; i<=7; ++i){
        Graph *g = new Graph;
        g->read_edge("../test_graphs/BHHS_graphs/E"+to_string(i)+".edg");
        // zero forcing number of g_{i} for max nullity
        fort_cover_ip(g,fcdata);
        mng = fcdata.val - 1;
        delete fcdata.zf_set;
        // graph g_{j}
        for(int j=1; j<=7; ++j){
            Graph *h = new Graph;
            h->read_edge("../test_graphs/BHHS_graphs/E"+to_string(j)+".edg");
            // zero forcing number of g_{j} for max nullity
            fort_cover_ip(h,fcdata);
            mnh = fcdata.val - 1;
            delete fcdata.zf_set;
            // generate random vertices u and v
            int u = rand_int(0,g->order-1);
            int v = rand_int(0,h->order-1);
            count += 1;
            myfile << to_string(count) + ", " + to_string(i) + ", " + to_string(j) + ", " + to_string(u) + ", " + to_string(v);
            // vertex sum graph
            Graph *vs = vert_sum(g,h,u,v);
            // fort number
		    ft_num_ip(vs,ftdata);
		    myfile << ", " + to_string(ftdata.val);
		    delete ftdata.dforts;
		    // fractional zero forcing
		    fzf_ip(vs,fzfdata);
		    myfile << ", " + to_string(fzfdata.val);
		    delete fzfdata.weights;
            // maximum nullity
            myfile << ", " + to_string(max_null_vs(g,h,u,v,mng,mnh));
            // zero forcing (FCM)
		    fort_cover_ip(vs,fcdata);
            myfile << ", " + to_string(fcdata.val) + "\n";
		    delete fcdata.zf_set;
            // delete vertex sum graph
            delete vs;
            // delete graph g_{j}
            delete h;
        }
        // delete graph g_{i}
        delete g;
    }
    // close file
	myfile.close();
    // return
    return 0;
}