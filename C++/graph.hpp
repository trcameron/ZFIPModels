#ifndef GRAPH_H
#define GRAPH_H
#include <iostream>
#include <fstream>
#include <stack>
#include <unordered_set>
#include <vector>
using namespace std;
/*
	Graph: graph class
*/
class Graph{
private:
	int* data_to_n(int *data,int &order,int &num){
		int *edges;
		if(data[0]<=62){
			order = data[0];
			num -= 1;
			edges = &data[1];
			return edges;
		}
		if(data[1]<=62){
			order = (data[1]<<12) + (data[2]<<6) + data[3];
			num -= 4;
			edges = &data[4];
			return edges;
		}
		order = (data[2]<<30) + (data[3]<<24) + (data[4]<<18) + (data[5]<<12) + (data[6]<<6) + data[7];
		num -= 8;
		edges = &data[8];
		return edges;
	}
	bool* bits(const int *edges,const int num){
		bool *b = new bool[num*6];
		for(int i=0; i<num; ++i){
			for(int j=5; j>=0; --j){
				b[i*6+(5-j)] = (edges[i]>>j) & 1;
			}
		}
		return b;
	}
	void parse_data(const int order,const int num,const int *edges,vector<bool> &b,vector<int> &x){
	    int k = 1;
	    while (1 << k < order) {
	        k++;
	    }
		
		vector<int> chunks(edges,edges+num);
	    int d = 0;  // partial data word
	    int dLen = 0;  // how many unparsed bits are left in d
		
		while(true){
			if(dLen < 1){
				if(chunks.empty()){
					return;
				}
				d = chunks.front();
				chunks.erase(chunks.begin());
				dLen = 6;
			}
			dLen --;
			int top_bit = (d >> dLen) & 1;  // grab top remaining bit
			
			int partial_x = d & ((1 << dLen) - 1);  // partially built up value of x
			int xLen = dLen;  // how many bits included so far in x
			while (xLen < k) {  // now grab full chunks until we have enough
				if(chunks.empty()){
					return;
				}
				d = chunks.front();
				chunks.erase(chunks.begin());
				dLen = 6;
	            partial_x = (partial_x << 6) + d;
	            xLen += 6;
			}
			partial_x >>= (xLen - k);  // shift back the extra bits
			dLen = xLen - k;
			b.push_back(top_bit);
			x.push_back(partial_x);
		}
	}
	void dfsUtil(bool* visited,int cur_node,int count,int &maxCount,int &farth_node) const{
		visited[cur_node] = true;
		count++;
		for(auto iter=adj[cur_node].cbegin(); iter!=adj[cur_node].cend(); ++iter){
			if(!visited[*iter]){
				if(count >= maxCount){
					maxCount = count;
					farth_node = *iter;
				}
				dfsUtil(visited,*iter,count,maxCount,farth_node);
			}
		}
	}
	void dfs(int cur_node,int &maxCount,int &farth_node) const{
		bool visited[order];
		int count = 0;
		for(int i=0; i<order; ++i){
			visited[i] = false;
		}
		dfsUtil(visited,cur_node,count,maxCount,farth_node);
	}
public:
	int order, size;
	unordered_set<int>* adj;		// nodes are 0,1,...,order-1; *(adj+i) is an unordered_set<int> that contains all nodes adjacent to i. 
	Graph(){
		order = 0; size = 0;
	}
	void read_graph6(const string &line){
		int num = line.size();		// number of characters in line
		int *data = new int[num];	// data for graph
		int i,j,k;			// indices for loops
		// shift line characters by 63
		for(i=0; i<num; i++){
			data[i] = line[i] - 63;
		}
		// get order and edge data
		int *edges = data_to_n(data,order,num);
		// gets bits representing each edge
		bool *b = bits(edges,num);		
		// add vertices
		adj = new unordered_set<int>[order];
		// add edges
		k = 0;
		for(j=0; j<order; ++j){
			for(i=0; i<j; ++i){
				if(b[k]){
					addEdge(i,j);
				}
				k++;
			}
		}
		// free memory
		delete[] data;
		delete[] b;
	}
	void read_sparse6(const string &line){
		int num = line.size()-1;	// number of characters in line, ignoring : character in front
		int *data = new int[num];	// data for graph
		int i,j,k;			// indices for loops
		// shift line characters by 63
		for(i=0; i<num; i++){
			data[i] = line[i+1] - 63;
		}
		// get order and edge data
		int *edges = data_to_n(data,order,num);
		// parse
		vector<bool> b;
		vector<int> x;
		parse_data(order,num,edges,b,x);
		// add vertices
		adj = new unordered_set<int>[order];
		// add edges
		j = 0;
		for(i = 0; i<b.size(); ++i){
			if(b[i]){
				j++;
			}
			// Padding with ones can cause overlarge number here
			if(x[i] >= order || j >= order){
				break;
			}
			else if(x[i] > j){
				j = x[i];
			}
			else{
				addEdge(x[i],j);
			}
		}
		// free memory
		delete[] data;
	}
	void read_edge(const string file_name){
		int i, j, k;
		fstream file;
		file.open(file_name);
		file >> order;
		file >> size;
		adj = new unordered_set<int>[order];
		for(k=0; k<size; ++k){
			file >> i;
			file >> j;
			(*(adj+i)).insert(j);
			(*(adj+j)).insert(i);
		}
		file.close();
	}
	void addEdge(const int u,const int v){
		if((*(adj+u)).find(v)==(*(adj+u)).cend()){
			(*(adj+u)).insert(v);
			(*(adj+v)).insert(u);
			size += 1;	
		}
	}
	void delEdge(const int u,const int v){
		if((*(adj+u)).find(v)!=(*(adj+u)).cend()){
			(*(adj+u)).erase(v);
			(*(adj+v)).erase(u);
			size -= 1;
		}
	}
	int deg(const int u) const{
		return (*(adj+u)).size();
	}
	bool is_connected(void) const{
		bool visited[order];
		for(int i=0; i<order; ++i){
			visited[i] = false;
		}
		stack<int> active;
		active.push(0);
		int v;
		while(!active.empty()){
			v = active.top();
			active.pop();
			if(!visited[v]){
				visited[v] = true;
				for(auto iter=adj[v].cbegin(); iter!=adj[v].cend(); ++iter){
					active.push(*iter);
				}
			}
		}
		bool check = visited[0];
		for(int i=1; i<order; ++i){
			check = check and visited[i];
		}
		return check;
	}
	int tree_diameter(void) const{
		int maxCount = INT_MIN;
		int cur_node = 0, farth_node = 0;
		dfs(cur_node,maxCount,farth_node);
		dfs(farth_node,maxCount,farth_node);
		return maxCount;
	}
	int max_degree(void) const{
		int max = adj->size();
		for(int i=1; i<order; ++i){
			if((*(adj + i)).size() > max){
				max = (*(adj + i)).size();
			}
		}
		return max;
	}
	void print(void) const{
		cout << "order: " << order << ", " << "size: " << size << endl;
		//nodeIter iter;
		for(int i=0; i<order; ++i){
			cout << i << ": ";
			//nbhd(i,iter);
			for(auto iter=adj[i].cbegin(); iter!=adj[i].cend(); ++iter){
				cout << *iter << " ";
			}
			cout << endl;
		}
	}
	void clear(void){
		order = 0; size = 0;
		delete[] adj;
	}
	~Graph(){
		order = 0; size = 0;
		delete[] adj;
	}
};
/*
	corona_prod: corona product of two graphs
*/
Graph* corona_prod(const Graph *g,const Graph *h);
/*
	cart_prod: cartesian product of two graphs
*/
Graph* cart_prod(const Graph *g,const Graph *h);
/*
	vert_del: vertex deletion from graph
*/
Graph* vert_del(const Graph *g,const int u);
/*
	vert_sum: vertex sum of two graphs
*/
Graph* vert_sum(const Graph *g,const Graph *h,const int u,const int v);
/*
	edge_sum: edge sum of two graphs
*/
Graph* edge_sum(const Graph *g,const Graph *h,const int u,const int v);
/*
	path: path graph of order n
*/
Graph* path_graph(const int n);
/*
	cycle: cycle graph of order n
*/
Graph* cycle_graph(const int n);
/*
 complete: complete graph of order n
*/
Graph* complete_graph(const int n);
/*
	hypercube: hypercube graph of dimension d
*/
Graph* hypercube_graph(const int d);
/*
	kneser: the K(n,k) Kneser graph
*/
Graph* kneser_graph(const int n,const int k);
/*
	petersen: petersen graph
*/
Graph* petersen(void);
/*
	nsun: Cn corona K1
*/
Graph* nsun_graph(const int n);
/*
	sun_link: link of k 5sun graphs
*/
Graph* sun_link_graph(const int k);
/*
	copy: return copy of given graph
*/
Graph* copy_graph(const Graph *g);
/*
	zf_closure: zero-forcing closure of given set on given graph
*/
void zf_closure(const Graph *g,unordered_set<int>* filled,int &pt);
/*
	node_comp: node complement of given set on given graph
*/
//unordered_set<int>* node_comp(const Graph *g,const unordered_set<int> *init);
/*
	wavefront: dynamic programming style improvement of the brute force algorithm to compute the zero-forcing number of a given graph
*/
int wavefront(const Graph *g);
#endif