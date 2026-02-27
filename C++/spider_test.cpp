#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_set>
#include "graph.hpp"
#include "ip_models.hpp"
using namespace std;
/*
    gen_part: generate integer partitions of rem into parts in non-increasing order, 
    where each part is in [1,maxPart], no integer is bigger than maxAllowed, and there at least minParts in each parition.
    Accumulate into cur. 
*/
void gen_parts(int rem, int maxPart, int maxAllowed, int minParts, vector<int>& cur, vector<vector<int>>& out){
    if (rem == 0) {
        if ((int)cur.size() >= minParts) out.push_back(cur);
        return;
    }

    // If we still need parts but even all 1's can't reach minParts? (not needed)
    // If rem is positive, we can always add 1's, so no feasibility issue here.

    int upper = rem < maxPart ? rem : maxPart;
    if (upper > maxAllowed) upper = maxAllowed;

    for (int x = upper; x >= 1; --x) {
        cur.push_back(x);
        gen_parts(rem - x, x, maxAllowed, minParts, cur, out);
        cur.pop_back();
    }
}
/*
    spider_leg_parts: returns all leg-length multisets for spiders on n vertices,
    partitions of n-1 with at least 3 parts and max leg length 5
*/
void spider_leg_parts(int n,vector<vector<int>>& legs) {
    legs.clear();
    vector<int> cur;
    gen_parts(n - 1, n - 1, 5, 3, cur, legs);
}
/*
    spider_graph: build spider graph for leg vector
*/
void spider_graph(vector<int> leg,Graph *g){
    int n=1;
    for(auto l=leg.cbegin(); l!=leg.cend(); ++l){
        n += *l;
    }

    g->order = n;
    g->size = 0;
    g->adj = new unordered_set<int>[n];

    const int center = 0;
    int next_id = 1;
    for(auto l=leg.cbegin(); l!=leg.cend(); ++l){
        int prev = center;
        for(int t=0; t<*l; ++t){
            g->addEdge(prev,next_id);
            prev = next_id;
            next_id++;
        }
    }
}
/*
    main file
*/
int main() {
    int count = 0;
    ofstream myfile;
    myfile.open("../csv_files/spider_test.csv");
	myfile << "count, order, minforts" << endl;
    for(int n=4; n<=31; ++n){
        vector<vector<int>> legs;
        spider_leg_parts(n,legs);
        for(auto leg=legs.cbegin(); leg!=legs.cend(); ++leg){
            Graph *g = new Graph();
            spider_graph(*leg,g);
            amfIPdata amfdata;
	        all_minimal_forts(g,amfdata);
            count ++;
            myfile << to_string(count) << ", " << to_string(g->order) << ", " << to_string(amfdata.mforts->size()) << endl;
		    delete amfdata.mforts;
            delete g;
        }
    }
    myfile.close();
    cout << "# of graphs: " << count << endl;
    return 0;
}