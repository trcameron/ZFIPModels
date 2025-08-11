#ifndef IP_MODELS_H
#define IP_MODELS_H
#include "graph.hpp"
#include "gurobi_c++.h"
#include <map>
using namespace std;
/*
	IP_MAX_TIME: maximum time (seconds) allowed for ip solver
*/
const double IP_MAX_TIME = 7200;
/*
	fcIPdata: fort cover IP data
*/
struct fcIPdata{
	int status, val;
	unordered_set<int> *zf_set;
};
/*
	fraczfIPdata: fractional zero-forcing IP data
*/
struct fzfIPdata{
	int status;
	double val;
	double *weights;
};
/*
	infIPdata: infection IP data
*/
struct infIPdata{
	int status, val;
	unordered_set<int> *zf_set;
	map<pair<int,int>,int> *forcings;
};
/*
 tsIPdata: time step IP data
*/
struct tsIPdata{
	int status, val;
	unordered_set<int> *zf_set;
	map<pair<int,int>,int> *forcings;
};
/*
	ptiIPdata: propagation time interval IP data
*/
struct ptiIPdata{
	int status;
	map<int,unordered_set<int>> *ptMap;
};
/*
	amfIPdata: all minimal forts IP data
*/
struct amfIPdata{
	int status;
	vector<unordered_set<int>> *mforts;
};
/*
	ftIPdata: fort number IP data
*/
struct ftIPdata{
	int status, val;
	vector<unordered_set<int>> *dforts;
};
/*
	fort_cover_ip: ip model for the zero-forcing number of a graph using fort cover constraints
*/
void fort_cover_ip(const Graph *g,fcIPdata &data);
/*
	fzf_ip: ip model for the fractional zero-forcing number of given graph
*/
void fzf_ip(const Graph *g,fzfIPdata &data);
/*
	infection_ip: ip model for the zero-forcing parameter of given graph
*/
void infection_ip(const Graph *g,const int T,infIPdata &data,char type);
/*
	time_step_ip: ip model for the zero-forcing parameter of given graph
*/
void time_step_ip(const Graph *g,const int T,tsIPdata &data,char type);
/*
	pt_interval: algorithm for computing pt interval of given graph
*/
void pt_interval(const Graph *g,const int T,ptiIPdata &data);
/*
	all_minimal_forts: recursive model for generating all minimal forts of given graph
*/
void all_minimal_forts(Graph *g,amfIPdata &data);
/*
	minimal_fort_ip: ip model for generating a minimal fort of the given graph not in the given collection of minimal forts
*/
/*
void minimal_fort_ip(const Graph *g,vector<unordered_set<int>> *forts);
*/
/*
	ft_num_ip: ip model for generating a maximum collection of disjoint forts of the given graph
*/
void ft_num_ip(const Graph *g,ftIPdata &data);
/*
	class violated_fort: Gurobi callback class for adding violated fort constraints to fort_cover_ip model
*/
class violated_fort: public GRBCallback{
	public:
		const Graph *g;
		GRBVar *s, *x;
		GRBModel *model_sub;
		violated_fort(const Graph *xg,GRBVar *xs,GRBVar *xx,GRBModel *xmodel_sub){
			g = xg;
			s = xs;
			x = xx;
			model_sub = xmodel_sub;
		}
	protected:
		void callback(){
			try{
				if(where == GRB_CB_MIPSOL){
					// found an integer solution, update sub model to see if it violates a fort constraint
					model_sub->remove(model_sub->getConstrByName("dynamic_constr"));
					GRBLinExpr expr = 0;
					for(int i=0; i<g->order; ++i){
						if(getSolution(s[i]) > 0.5){
							expr += x[i];
						}
					}
					model_sub->addConstr(expr,GRB_LESS_EQUAL,0,"dynamic_constr");
					model_sub->optimize();
					// if sub model is feasible then there is a fort violation, add constraint (best to add minimum fort violation)
					if(model_sub->get(GRB_IntAttr_Status) == GRB_OPTIMAL){
						expr = 0;
						for(int i=0; i<g->order; ++i){
							if(x[i].get(GRB_DoubleAttr_X) > 0.5){
								expr += s[i];
							}
						}
						addLazy(expr >= 1);
					}
				}
			}catch(GRBException e){
				cout << "Error number: " << e.getErrorCode() << endl;
				cout << e.getMessage() << endl;
			}catch(...){
				cout << "Error during violated_fort callback" << endl;
			}
		}
};
/*
	class frac_violated_fort: Gurobi callback class for adding violated fort constraints to frac_zf_ip model
*/
class frac_violated_fort: public GRBCallback{
	public:
		const Graph *g;
		GRBVar *s, *x;
		GRBModel *model_sub;
		frac_violated_fort(const Graph *xg,GRBVar *xs,GRBVar *xx,GRBModel *xmodel_sub){
			g = xg;
			s = xs;
			x = xx;
			model_sub = xmodel_sub;
		}
	protected:
		void callback(){
			try{
				if(where == GRB_CB_MIPSOL){
					// found an RMP solution, update objective function
					GRBLinExpr expr = 0;
					for(int i=0; i<g->order; ++i){
						expr += getSolution(s[i])*x[i];
					}
					model_sub->setObjective(expr,GRB_MINIMIZE);
					model_sub->optimize();
					// check if RMP is violated by fort constraint
					if(model_sub->get(GRB_DoubleAttr_ObjVal) < 1.0){
						expr = 0;
						for(int i=0; i<g->order; ++i){
							if(x[i].get(GRB_DoubleAttr_X) > 0.5){
								expr += s[i];
							}
						}
						addLazy(expr >= 1);
					}
				}
			}catch(GRBException e){
				cout << "Error number: " << e.getErrorCode() << endl;
				cout << e.getMessage() << endl;
			}catch(...){
				cout << "Error during frac_violated_fort callback" << endl;
			}
		}
};
#endif