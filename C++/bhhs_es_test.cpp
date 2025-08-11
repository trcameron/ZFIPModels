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
    min_rank_es: minimum rank of edge sum with graphs of order 8 whose minimum rank is known
*/
int min_rank_es(const Graph *g,const Graph *h,const int u,const int v,const int mrg,const int mrh){
    // fort cover ip data
    fcIPdata fcdata;
    // g-u, min rank and rank spread
    Graph *vd1 = vert_del(g,u);
    fort_cover_ip(vd1,fcdata);
    int mr1 = vd1->order - fcdata.val;
    //cout << "mr1: " << mr1 << endl;
    int rs1 = mrg - mr1;
    //cout << "rs1: " << rs1 << endl;
    delete vd1;
    delete fcdata.zf_set;
    // h-v, min rank and rank spread
    Graph *vd2 = vert_del(h,v);
    fort_cover_ip(vd2,fcdata);
    int mr2 = vd2->order - fcdata.val;
    //cout << "mr2: " << mr2 << endl;
    int rs2 = mrh - mr2;
    //cout << "rs2: " << rs2 << endl;
    delete vd2;
    delete fcdata.zf_set;
    // return minimum rank
    if(rs1==2 || rs2==2){
        return mrg + mrh;
    }
    else{
        return mrg + mrh + 1;
    }
}
/*
    max_null_es: maximum nullity of edge sum with graphs of order 8 whose maximum nullity is known
*/
int max_null_es(const Graph *g,const Graph *h,const int u,const int v,const int mng, const int mnh){
    // return maximum nullity
    return (g->order + h->order) - min_rank_es(g,h,u,v,g->order-mng,h->order-mnh);
}
/*
	main function
*/
int main(int argc,char *argv[]){
    /*
    fcIPdata fcdata;
    Graph *g = new Graph;
    g->read_edge("../test_graphs/BHHS_graphs/E"+to_string(6)+".edg");
    fort_cover_ip(g,fcdata);
    cout << "zg: " << fcdata.val << endl;
    delete fcdata.zf_set;
    Graph *h = new Graph;
    h->read_edge("../test_graphs/BHHS_graphs/E"+to_string(4)+".edg");
    fort_cover_ip(h,fcdata);
    cout << "zh: " << fcdata.val << endl;
    delete fcdata.zf_set;
    int u = 4, v = 4, mng = 3, mnh = 3;
    int mr = min_rank_es(g,h,u,v,g->order - mng,h->order - mnh);
    cout << "mr: " << mr << endl;
    return 0;
    */
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
	myfile.open("../csv_files/bhhs_es_test.csv");
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
            // edge sum graph
            Graph *es = edge_sum(g,h,u,v);
            // fort number
		    ft_num_ip(es,ftdata);
		    myfile << ", " + to_string(ftdata.val);
		    delete ftdata.dforts;
		    // fractional zero forcing
		    fzf_ip(es,fzfdata);
		    myfile << ", " + to_string(fzfdata.val);
		    delete fzfdata.weights;
            // maximum nullity
            myfile << ", " + to_string(max_null_es(g,h,u,v,mng,mnh));
            // zero forcing (FCM)
		    fort_cover_ip(es,fcdata);
            myfile << ", " + to_string(fcdata.val) + "\n";
		    delete fcdata.zf_set;
            // delete edge sum graph
            delete es;
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