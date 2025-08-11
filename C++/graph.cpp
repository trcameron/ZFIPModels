#include "graph.hpp"
#include <algorithm>
#include <list>
#include <map>
#include <queue>
using namespace std;
/*
	corona_prod: corona product of two graphs
*/
Graph* corona_prod(const Graph *g,const Graph *h){
	Graph *cp = new Graph;
	cp->order = g->order + (g->order)*(h->order);
	cp->adj = new unordered_set<int>[cp->order];
	for(int u=0; u<g->order; ++u){
		for(auto iter=g->adj[u].cbegin(); iter!=g->adj[u].cend(); ++iter){
			cp->addEdge(u,*iter);
		}
		for(int v=0; v<h->order; ++v){
			cp->addEdge(u,g->order + u*h->order + v);
			for(auto iter=h->adj[v].cbegin(); iter!=h->adj[v].cend(); ++iter){
				cp->addEdge(g->order + u*h->order + v,g->order + u*h->order + *iter);
			}
		}
	}
	return cp;
}
/*
	cart_prod: cartesian product of two graphs
*/
Graph* cart_prod(const Graph *g,const Graph *h){
	Graph *cp = new Graph;
	cp->order = (g->order)*(h->order);
	cp->adj = new unordered_set<int>[cp->order];
	for(int u=0; u<g->order; ++u){
		for(int v=0; v<h->order; ++v){
			for(auto iter=h->adj[v].cbegin(); iter!=h->adj[v].cend(); ++iter){
				cp->addEdge(u*h->order+v,u*h->order+(*iter));
			}
		}
	}
	for(int v=0; v<h->order; ++v){
		for(int u=0; u<g->order; ++u){
			for(auto iter=g->adj[u].cbegin(); iter!=g->adj[u].cend(); ++iter){
				cp->addEdge(u*h->order+v,(*iter)*h->order+v);
			}
		}
	}
	return cp;
}
/*
	vert_del: vertex deletion from graph
*/
Graph* vert_del(const Graph *g,const int u){
	Graph *vd = new Graph;
	vd->order = g->order - 1;
	vd->adj = new unordered_set<int>[vd->order];
	for(int i=0; i<u; ++i){
		for(auto it=g->adj[i].cbegin(); it!=g->adj[i].cend(); ++it){
			if(*it < u){
				vd->addEdge(i,*it);
			}
			else if(*it > u){ 
				vd->addEdge(i,*it-1);
			}
		}
	}
	for(int i=u+1; i<g->order; ++i){
		for(auto it=g->adj[i].cbegin(); it!=g->adj[i].cend(); ++it){
			if(*it < u){
				vd->addEdge(i-1,*it);
			}
			else if(*it > u){
				vd->addEdge(i-1,*it-1);
			}
		}
	}
	return vd;
}
/*
	vert_sum: vertex sum of two graphs
*/
Graph* vert_sum(const Graph *g,const Graph *h,const int u,const int v){
	Graph *vs = new Graph;
	vs->order = (g->order-1)+(h->order-1) + 1;
	vs->adj = new unordered_set<int>[vs->order];
	for(int i=0; i<g->order; ++i){
		for(auto iter=g->adj[i].cbegin(); iter!=g->adj[i].cend(); ++iter){
			vs->addEdge(i,*iter);
		}
	}
	for(int i=0; i<h->order; ++i){
		for(auto iter=h->adj[i].cbegin(); iter!=h->adj[i].cend(); ++iter){
			if(i<v and *iter<v){
				vs->addEdge(g->order + i,g->order + *iter);
			}
			else if(i<v and *iter>v){
				vs->addEdge(g->order + i,g->order + *iter - 1);
			}
			else if(i>v and *iter<v){
				vs->addEdge(g->order + i - 1,g->order + *iter);
			}
			else if(i>v and *iter>v){
				vs->addEdge(g->order + i - 1,g->order + *iter - 1);
			}
			else if(i==v and *iter<v){
				vs->addEdge(u,g->order + *iter);
			}
			else if(i==v and *iter>v){
				vs->addEdge(u,g->order + *iter - 1);
			}
		}
	}
	return vs;
}
/*
	edge_sum: edge sum of two graphs
*/
Graph* edge_sum(const Graph *g,const Graph *h,const int u,const int v){
	Graph *es = new Graph;
	es->order = g->order + h->order;
	es->adj = new unordered_set<int>[es->order];
	for(int i=0; i<g->order; ++i){
		for(auto iter=g->adj[i].cbegin(); iter!=g->adj[i].cend(); ++iter){
			es->addEdge(i,*iter);
		}
	}
	for(int i=0; i<h->order; ++i){
		for(auto iter=h->adj[i].cbegin(); iter!=h->adj[i].cend(); ++iter){
			es->addEdge(g->order + i,g->order + *iter);
		}
	}
	es->addEdge(u,g->order + v);
	return es;
}
/*
	path: path graph of order n
*/
Graph* path_graph(const int n){
	Graph *p = new Graph;
	p->order = n;
	p->adj = new unordered_set<int>[n];
	for(int i=0; i<n-1; ++i){
		p->addEdge(i,i+1);
	}
	return p;
}
/*
	cycle: cycle graph of order n
*/
Graph* cycle_graph(const int n){
	Graph *c = path_graph(n);
	c->addEdge(0,n-1);
	return c;
}
/*
 complete: complete graph of order n
*/
Graph* complete_graph(const int n){
	Graph *k = new Graph;
	k->order = n;
	k->adj = new unordered_set<int>[n];
	for(int i=0; i<n; ++i){
		for(int j=i+1; j<n; ++j){
			k->addEdge(i,j);
		}
	}
	return k;
}
/*
	hypercube: hypercube graph of dimension d
*/
Graph* hypercube_graph(const int d){
	// path graph of order 2
	Graph *p = path_graph(2);
	// initiate old q as path graph of order 2
	Graph *old_q = copy_graph(p);
	// new q pointer
	Graph *new_q = old_q;
	// recursively update new q
	for(int i=2; i<=d; ++i){
		new_q = cart_prod(old_q,p);
		delete old_q;
		old_q = new_q;
	}
	// delete p
	delete p;
	// return new q (old q)
	return new_q;
}
/*
	generate_k_subsets: generates all size k subsets from a set of size n
*/
void generate_k_subsets(vector<unordered_set<int>> *subsets,const int n,const int k){
	int p[n];
	bool v[n];
	for(int i=0; i<n; ++i){
		p[i] = i;
		if(i<(n-k)){
			v[i] = false;
		}
		else{
			v[i] = true;
		}
	}
	unordered_set<int> current_subset;
	do{
		for(int i=0; i<n; ++i){
			if(v[i]){
				current_subset.insert(i);
			}
		}
		subsets->push_back(current_subset);
		current_subset.clear();
	}while(next_permutation(v,v+n));
}
/*
	are_disjoint: checks if to unordered sets are disjoint
*/
bool are_disjoint(const unordered_set<int> *set1,const unordered_set<int> *set2){
	for(auto u=set1->cbegin(); u!=set1->cend(); ++u){
		if(set2->find(*u) != set2->cend()){
			return false;
		}
	}
	return true;
}
/*
	kneser: the K(n,k) Kneser graph
*/
Graph* kneser_graph(const int n,const int k){
	vector<unordered_set<int>> *subsets = new vector<unordered_set<int>>;
	generate_k_subsets(subsets,n,k);
	Graph *g = new Graph;
	g->order = subsets->size();
	g->adj = new unordered_set<int>[g->order];
	for(int i=0; i<g->order; ++i){
		for(int j=i+1; j<g->order; ++j){
			if(are_disjoint(&subsets->at(i),&subsets->at(j))){
				g->addEdge(i,j);
			}
		}
	}
	delete subsets;
	return g;
}
/*
	petersen: petersen graph
*/
Graph* petersen(void){
	Graph *g = new Graph;
	g->order = 10;
	g->adj = new unordered_set<int>[10];
	for(int i=0; i<5; ++i){
		g->addEdge(i,(i+1)%5);
		g->addEdge(i,i+5);
	}
	g->addEdge(5,7);
	g->addEdge(5,8);
	g->addEdge(6,8);
	g->addEdge(6,9);
	g->addEdge(7,9);
	return g;
}
/*
	nsun: Cn corona K1
*/
Graph* nsun_graph(const int n){
	// cycle of order n
	Graph *c = cycle_graph(n);
	// complete graph of order 1
	Graph *k = complete_graph(1);
	// nsun graph
	Graph *g = corona_prod(c,k);
	// delete cycle and complete graph
	delete c;
	delete k;
	// return nsun graph
	return g;
}
/*
	sun_link: link of k 5sun graphs
*/
Graph* sun_link_graph(const int k){
	Graph *g = nsun_graph(5);
	Graph *h = copy_graph(g);
	Graph *vs;

	const int v=9;
	int u = 6;
	for(int i=1; i<k; ++i){
		vs = vert_sum(g,h,u,v);
		u = g->order + 6;
		delete g;
		g = vs;
	}
	delete h;
	return g;
}
/*
	copy: return copy of given graph
*/
Graph* copy_graph(const Graph *g){
	Graph *c = new Graph;
	c->order = g->order;
	c->adj = new unordered_set<int>[c->order];
	for(int u=0; u<c->order; ++u){
		for(auto v=g->adj[u].cbegin(); v!=g->adj[u].cend(); ++v){
			c->addEdge(u,*v);
		}
	}
	return c;
}
/*
	zf_closure: zero-forcing closure of given set on given graph
*/
void zf_closure(const Graph *g,unordered_set<int>* filled,int &pt){
	// for each possible time step
	int count, node, t;
	unordered_set<int> active;
	for(t=0; t<g->order; ++t){
		// fill active set
		for(auto u=filled->cbegin(); u!=filled->cend(); ++u){
			count = 0;
			for(auto v=g->adj[*u].cbegin(); v!=g->adj[*u].cend(); ++v){
				if(filled->find(*v)==filled->cend()){
					count += 1;
					node = *v;
				}
			}
			if(count == 1){
				active.insert(node);
			}
		}
		// if active is empty, break
		if(active.empty()){
			break;
		}
		else{
			// apply all active forcings
			filled->insert(active.cbegin(),active.cend());
			active.clear();
		}
	}
	// set propagation time
	if(filled->size()==g->order){
		pt = t;
	}
	else{
		pt = INT_MAX;
	}
}
/*
	node_comp: node complement of given set on given graph
*/
/*
unordered_set<int>* node_comp(const Graph *g,const unordered_set<int> *init){
	// initialize complement
	unordered_set<int>* comp = new unordered_set<int>;
	for(int i=0; i<g->order; ++i){
		if(init->find(i)==init->cend()){
			comp->insert(i);
		}
	}
	return comp;
}
*/
/*
	wavefront: dynamic programming style improvement of the brute force algorithm to compute the zero-forcing number of a given graph
*/
int wavefront(const Graph *g){
	list<pair<unordered_set<int>,int>> cl_pairs;
	cl_pairs.push_back(pair<unordered_set<int>,int>({},0));
	unordered_set<int> s;
	unordered_set<int> *s_new;
	int r, r_new, card;
	for(int i=1; i<=g->order; ++i){
		for(auto j=cl_pairs.cbegin(); j!=cl_pairs.cend(); ++j){
			s = j->first;
			r = j->second;
			for(int v=0; v<g->order; ++v){
				s_new = new unordered_set<int>(s.cbegin(),s.cend());
				s_new->insert(v);
				s_new->insert(g->adj[v].cbegin(),g->adj[v].cend());
				zf_closure(g,s_new,card);
				r_new = r + (int)(s.find(v)==s.cend());
				card = -1;
				for(auto iter=g->adj[v].cbegin(); iter!=g->adj[v].cend(); ++iter){
					card += (int)(s.find(*iter)==s.cend());
				}
				if(card > 0){
					r_new += card;
				}
				if(r_new <= i){
					for(auto iter=cl_pairs.cbegin(); iter!=cl_pairs.cend(); ++iter){
						if(iter->first==*s_new and iter->second<=i){
							goto DELETE;
						}
					}
					if(s_new->size()==g->order){
						goto RETURN;
					}
					cl_pairs.push_back(pair<unordered_set<int>,int>(*s_new,r_new));
				}
				DELETE:
				delete s_new;
			}
		}
	}
	RETURN:
	delete s_new;
	return r_new;
}
/*
	main function
*/
/*
int main(){
	Graph *g = hypercube_graph(3);
	Graph *h = hypercube_graph(3);
	Graph *vs = vert_sum(g,h,0,0);
	cout << "order: " << vs->order << endl;
	cout << "size: " << vs->size << endl;
	delete g;
	delete h;
	delete vs;
	return 0;
}
*/