#include "graph.hpp"
#include "ip_models.hpp"
using namespace std;
/*
	fort_cover_ip: ip model for the zero-forcing number of a graph using fort cover constraints
*/
void fort_cover_ip(const Graph *g,fcIPdata &data){
	// main model
	GRBEnv* env_main = new GRBEnv(true);
	env_main->set(GRB_IntParam_OutputFlag,0);
	env_main->set(GRB_IntParam_LogToConsole,0);
	env_main->set(GRB_IntParam_LazyConstraints,1);	
	env_main->set(GRB_DoubleParam_TimeLimit,IP_MAX_TIME);
	env_main->start();
	GRBModel* model_main = new GRBModel(*env_main);
	// main model variables
	GRBVar* s = new GRBVar[g->order];
	for(int i=0; i<g->order; ++i){
		s[i] = model_main->addVar(0,1,1,GRB_BINARY);
	}
	// main model objective
	model_main->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
	// sub model
	GRBEnv* env_sub = new GRBEnv(true);
	env_sub->set(GRB_IntParam_OutputFlag,0);
	env_sub->set(GRB_IntParam_LogToConsole,0);
	env_sub->start();
	GRBModel* model_sub = new GRBModel(*env_sub);
	// sub model variables
	GRBVar* x = new GRBVar[g->order];
	for(int i=0; i<g->order; ++i){
		x[i] = model_sub->addVar(0,1,1,GRB_BINARY);
	}
	// sub model objective
	model_sub->set(GRB_IntAttr_ModelSense,GRB_MINIMIZE);
	// sub model static constraints
	GRBLinExpr expr = 0;
	for(int i=0; i<g->order; ++i){
		expr += x[i];
	}
	model_sub->addConstr(expr,GRB_GREATER_EQUAL,1);
	for(int i=0; i<g->order; ++i){
		for(auto iter1=g->adj[i].cbegin(); iter1!=g->adj[i].cend(); ++iter1){
			expr = x[*iter1] - 2*x[i];
			for(auto iter2=g->adj[*iter1].cbegin(); iter2!=g->adj[*iter1].cend(); ++iter2){
				expr += x[*iter2];
			}
			model_sub->addConstr(expr,GRB_GREATER_EQUAL,0);
		}
	}
	// sub model dynamic constraint
	expr = 0;
	model_sub->addConstr(expr,GRB_LESS_EQUAL,0,"dynamic_constr");
	// sub model update
	model_sub->update();
	// set callback function in main model
	violated_fort cb = violated_fort(g,s,x,model_sub);
	model_main->setCallback(&cb);
	// gurobi optimize
	model_main->optimize();
	// extract solution
	data.status = model_main->get(GRB_IntAttr_Status);
	data.val = round(model_main->get(GRB_DoubleAttr_ObjVal));
	data.zf_set = new unordered_set<int>;
	for(int i=0; i<g->order; ++i){
		if(s[i].get(GRB_DoubleAttr_X) > 0.5){
			data.zf_set->insert(i);
		}
	}
	// free memory
	delete env_main;
	delete env_sub;
	delete model_main;
	delete model_sub;
	delete[] s;
	delete[] x;
}
/*
	fzf_ip: ip model for the fractional zero-forcing number of given graph
*/
void fzf_ip(const Graph *g,fzfIPdata &data){
	// main model
	GRBEnv* env_main = new GRBEnv(true);
	env_main->set(GRB_IntParam_OutputFlag,0);
	env_main->set(GRB_IntParam_LogToConsole,0);
	env_main->set(GRB_IntParam_LazyConstraints,1);
	env_main->set(GRB_DoubleParam_TimeLimit,IP_MAX_TIME);
	env_main->start();
	GRBModel* model_main = new GRBModel(*env_main);
	// main model variables
	GRBVar* s = new GRBVar[g->order];
	for(int i=0; i<g->order; ++i){
		s[i] = model_main->addVar(0,1,1,GRB_CONTINUOUS);
	}
	GRBVar z = model_main->addVar(0,1,1,GRB_BINARY);		// dummy binary variable so Gurobi will use callback
	// main model objective
	model_main->set(GRB_IntAttr_ModelSense,GRB_MINIMIZE);
	// sub model
	GRBEnv* env_sub = new GRBEnv(true);
	env_sub->set(GRB_IntParam_OutputFlag,0);
	env_sub->set(GRB_IntParam_LogToConsole,0);
	env_sub->start();
	GRBModel* model_sub = new GRBModel(*env_sub);
	// sub model variables
	GRBVar* x = new GRBVar[g->order];
	for(int i=0; i<g->order; ++i){
		x[i] = model_sub->addVar(0,1,0,GRB_BINARY);
	}
	// sub model constraints
	GRBLinExpr expr = 0;
	for(int i=0; i<g->order; ++i){
		expr += x[i];
	}
	model_sub->addConstr(expr,GRB_GREATER_EQUAL,1);
	for(int i=0; i<g->order; ++i){
		for(auto iter1=g->adj[i].cbegin(); iter1!=g->adj[i].cend(); ++iter1){
			expr = x[*iter1] - 2*x[i];
			for(auto iter2=g->adj[*iter1].cbegin(); iter2!=g->adj[*iter1].cend(); ++iter2){
				expr += x[*iter2];
			}
			model_sub->addConstr(expr,GRB_GREATER_EQUAL,0);
		}
	}
	// sub model update
	model_sub->update();
	// set callback function in main model
	frac_violated_fort cb = frac_violated_fort(g,s,x,model_sub);
	model_main->setCallback(&cb);
	// gurobi optimize
	model_main->optimize();
	// extract solution
	data.status = model_main->get(GRB_IntAttr_Status);
	data.val = model_main->get(GRB_DoubleAttr_ObjVal);
	data.weights = new double[g->order];
	for(int i=0; i<g->order; ++i){
		data.weights[i] = s[i].get(GRB_DoubleAttr_X);
	}
	// free memory
	delete env_main;
	delete env_sub;
	delete model_main;
	delete model_sub;
	delete[] s;
	delete[] x;
}
/*
	infection_ip: ip model for the zero-forcing parameter of given graph
*/
void infection_ip(const Graph *g,const int T,infIPdata &data,char type){
	// arc map
	map<pair<int,int>,int> arcMap;
	int id = 0;
	for(int i=0; i<g->order; ++i){
		for(auto iter=g->adj[i].cbegin(); iter!=g->adj[i].cend(); ++iter){
			arcMap[pair<int,int>(i,*iter)] = id;
			id += 1;
		}
	}
	// gurobi model
	GRBEnv* env = new GRBEnv(true);
	env->set(GRB_IntParam_OutputFlag,0);
	env->set(GRB_IntParam_LogToConsole,0);
	env->set(GRB_DoubleParam_TimeLimit,IP_MAX_TIME);
	env->start();
	GRBModel* model = new GRBModel(*env);
	// gurobi variables
	GRBVar* s = new GRBVar[g->order];
	GRBVar* x = new GRBVar[g->order];
	GRBVar* y = new GRBVar[2*g->size];
	for(int i=0; i<g->order; ++i){
		s[i] = model->addVar(0,1,1,GRB_BINARY);
		x[i] = model->addVar(0,T,0,GRB_INTEGER);
	}
	for(int i=0; i<2*g->size; ++i){
		y[i] = model->addVar(0,1,0,GRB_BINARY);
	}
	GRBVar z;
	if(type == 'Z'){
		z = model->addVar(0,T,0,GRB_INTEGER);
	}
	else if(type == 'p'){
		z = model->addVar(0,T,1.0/(2.0*T),GRB_INTEGER);
	}
	else if(type == 'T'){
		z = model->addVar(0,T,1,GRB_INTEGER);
	}
	else{
		cout << "infection_ip: invalid type parameter" << endl;
		return;
	}
	// gurobi set objective
	model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
	// constraint 1
	GRBLinExpr expr;
	for(int i=0; i<g->order; ++i){
		expr = s[i];
		for(auto iter=g->adj[i].cbegin(); iter!=g->adj[i].cend(); ++iter){
			id = arcMap[pair<int,int>(*iter,i)];
			expr += y[id];
		}
		model->addConstr(expr==1);
	}
	// constraint 2
	int u,v;
	for(auto iter=arcMap.cbegin(); iter!=arcMap.cend(); ++iter){
		u = iter->first.first;
		v = iter->first.second;
		id = iter->second;
		expr = x[u] - x[v] + (T+1)*y[id];
		model->addConstr(expr<=T);
	}
	// constraint 3
	for(auto iter1=arcMap.cbegin(); iter1!=arcMap.cend(); ++iter1){
		u = iter1->first.first;
		v = iter1->first.second;
		id = iter1->second;
		for(auto iter2=g->adj[u].cbegin(); iter2!=g->adj[u].cend(); ++iter2){
			if(*iter2 != v){
				expr = x[*iter2] - x[v] + (T+1)*y[id];
				model->addConstr(expr<=T);	
			}
		}
	}
	// constraint 4
	for(int i=0; i<g->order; ++i){
		expr = x[i] - z;
		model->addConstr(expr<=0);
	}
	// gurobi optimize
	model->optimize();
	// extract solution
	data.status = model->get(GRB_IntAttr_Status);
	if(type == 'Z' or type == 'T'){
		data.val = round(model->get(GRB_DoubleAttr_ObjVal));
	}
	else{
		data.val = round(z.get(GRB_DoubleAttr_X));
	}
	data.zf_set = new unordered_set<int>;
	for(int i=0; i<g->order; ++i){
		if(s[i].get(GRB_DoubleAttr_X) > 0.5){
			data.zf_set->insert(i);
		}
	}
	data.forcings = new map<pair<int,int>,int>;
	for(auto iter=arcMap.cbegin(); iter!=arcMap.cend(); ++iter){
		u = iter->first.first;
		v = iter->first.second;
		id = iter->second;
		if(y[id].get(GRB_DoubleAttr_X) > 0.5){
			data.forcings->insert({pair<int,int>(u,v),round(x[v].get(GRB_DoubleAttr_X))});
		}
	}
	// free memory
	delete env;
	delete model;
	delete[] s;
	delete[] x;
	delete[] y;
}
/*
	time_step_ip: ip model for the zero-forcing parameter of given graph
*/
void time_step_ip(const Graph *g,const int T,tsIPdata &data,char type){
	// arc map
	map<pair<int,int>,int> arcMap;
	int id = 0;
	for(int i=0; i<g->order; ++i){
		for(auto iter=g->adj[i].cbegin(); iter!=g->adj[i].cend(); ++iter){
			arcMap[pair<int,int>(i,*iter)] = id;
			id += 1;
		}
	}
	// gurobi model
	GRBEnv* env = new GRBEnv(true);
	env->set(GRB_IntParam_OutputFlag,0);
	env->set(GRB_IntParam_LogToConsole,0);
	env->set(GRB_DoubleParam_TimeLimit,IP_MAX_TIME);
	env->start();
	GRBModel* model = new GRBModel(*env);
	// gurobi variables
	GRBVar** x = new GRBVar*[T+1];
	for(int i=0; i<=T; ++i){
		x[i] = new GRBVar[g->order];
		for(int j=0; j<g->order; ++j){
			x[i][j] = model->addVar(0,1,0,GRB_BINARY);
		}
	}
	GRBVar** y = new GRBVar*[T];
	for(int i=0; i<T; ++i){
		y[i] = new GRBVar[2*g->size];
		for(int j=0; j<2*g->size; ++j){
			y[i][j] = model->addVar(0,1,0,GRB_BINARY);
		}
	}
	GRBVar* z = new GRBVar[T];
	for(int i=0; i<T; ++i){
		z[i] = model->addVar(0,1,0,GRB_BINARY);
	}
	// gurobi set objective
	double eps = 1.0/(2.0*T);
	GRBLinExpr expr = x[0][0];
	for(int j=1; j<g->order; ++j){
		expr += x[0][j];
	}
	if(type == 'Z'){
		model->setObjective(expr,GRB_MINIMIZE);
	}
	else if(type == 'p'){
		for(int i=0; i<T; ++i){
			expr += eps*z[i];
		}
		model->setObjective(expr,GRB_MINIMIZE);
	}
	else if(type == 'P'){
		for(int i=0; i<T; ++i){
			expr -= eps*z[i];
		}
		model->setObjective(expr,GRB_MINIMIZE);
	}
	else if(type == 'T'){
		for(int i=0; i<T; ++i){
			expr += z[i];
		}
		model->setObjective(expr,GRB_MINIMIZE);
	}
	else{
		cout << "time_step_ip: invalid type parameter" << endl;
		return;
	}
	// constraint 1
	for(int v=0; v<g->order; ++v){
		expr =  x[0][v];
		for(int i=0; i<T; ++i){
			for(auto iter=g->adj[v].cbegin(); iter!=g->adj[v].cend(); ++iter){
				int j = arcMap[pair<int,int>(*iter,v)];
				expr += y[i][j];
			}
		}
		model->addConstr(expr==1);
	}
	// constraints 2 and 3
	for(auto a=arcMap.cbegin(); a!=arcMap.cend(); ++a){
		int u = (a->first).first, v = (a->first).second, j = a->second;
		for(int i=0; i<T; ++i){
			expr = y[i][j] - x[i][u];
			model->addConstr(expr<=0);
		}
		for(auto iter=g->adj[u].cbegin(); iter!=g->adj[u].cend(); ++iter){
			if(*iter != v){
				for(int i=0; i<T; ++i){
					expr = y[i][j] - x[i][*iter];
					model->addConstr(expr<=0);
				}
			}
		}
	}
	// constraint 4
	for(int v=0; v<g->order; ++v){
		for(int i=0; i<T; ++i){
			expr = x[i+1][v] - x[i][v];
			for(auto iter=g->adj[v].cbegin(); iter!=g->adj[v].cend(); ++iter){
				int j = arcMap[pair<int,int>(*iter,v)];
				expr -= y[i][j];
			}
			model->addConstr(expr==0);
		}
	}
	// constraint 5
	for(int i=0; i<T; ++i){
		for(auto a=arcMap.cbegin(); a!=arcMap.cend(); ++a){
			int u = (a->first).first, v = (a->first).second;
			expr = x[i][u] - x[i][v];
			for(auto iter=g->adj[u].cbegin(); iter!=g->adj[u].cend(); ++iter){
				if(*iter != v){
					expr += x[i][*iter];	
				}
			}
			for(auto iter=g->adj[v].cbegin(); iter!=g->adj[v].cend(); ++iter){
				int j = arcMap[pair<int,int>(*iter,v)];
				expr -= y[i][j];
			}
			model->addConstr(expr<=g->deg(u)-1);
		}
	}
	// constraint 6
	eps = 1.0/g->order;
	for(int i=1; i<=T; ++i){
		expr = eps*(x[i][0] - x[i-1][0]);
		for(int j=1; j<g->order; ++j){
			expr += eps*(x[i][j] - x[i-1][j]);
		}
		expr -= z[i-1];
		model->addConstr(expr<=0);
	}
	// constraint 7
	for(int i=1; i<=T; ++i){
		expr = z[i-1] - (x[i][0] - x[i-1][0]);
		for(int j=1; j<g->order; ++j){
			expr -= (x[i][j] - x[i-1][j]);
		}
		model->addConstr(expr<=0);
	}
	// gurobi optimize
	model->optimize();
	// extract solution
	data.status = model->get(GRB_IntAttr_Status);
	if(type == 'Z' or type == 'T'){
		data.val = round(model->get(GRB_DoubleAttr_ObjVal));
	}
	else{
		data.val = round(z[0].get(GRB_DoubleAttr_X));
		for(int i=1; i<T; ++i){
			data.val += round(z[i].get(GRB_DoubleAttr_X));
		}
	}
	data.zf_set = new unordered_set<int>;
	for(int j=0; j<g->order; ++j){
		if(x[0][j].get(GRB_DoubleAttr_X) > 0.5){
			data.zf_set->insert(j);
		}
	}
	data.forcings = new map<pair<int,int>,int>;
	for(auto iter=arcMap.cbegin(); iter!=arcMap.cend(); ++iter){
		int u = iter->first.first;
		int v = iter->first.second;
		int j = iter->second;
		for(int i=0; i<T; ++i){
			if(y[i][j].get(GRB_DoubleAttr_X) > 0.5){
				data.forcings->insert({pair<int,int>(u,v),i});
				break;
			}
		}
	}
	// free memory
	delete env;
	delete model;
	for(int i=0; i<=T; ++i){
		delete[] x[i];
	}
	delete[] x;
	for(int i=0; i<T; ++i){
		delete[] y[i];
	}
	delete[] y;
	delete[] z;
}
/*
	pt_interval: algorithm for computing pt interval of given graph
*/
void pt_interval(const Graph *g,const int T,ptiIPdata &data){
	// arc map
	map<pair<int,int>,int> arcMap;
	int id = 0;
	for(int i=0; i<g->order; ++i){
		for(auto iter=g->adj[i].cbegin(); iter!=g->adj[i].cend(); ++iter){
			arcMap[pair<int,int>(i,*iter)] = id;
			id += 1;
		}
	}
	// gurobi model
	GRBEnv* env = new GRBEnv(true);
	env->set(GRB_IntParam_OutputFlag,0);
	env->set(GRB_IntParam_LogToConsole,0);
	env->set(GRB_DoubleParam_TimeLimit,IP_MAX_TIME);
	env->start();
	GRBModel* model = new GRBModel(*env);
	// gurobi variables
	GRBVar** x = new GRBVar*[T+1];
	for(int i=0; i<=T; ++i){
		x[i] = new GRBVar[g->order];
		for(int j=0; j<g->order; ++j){
			x[i][j] = model->addVar(0,1,0,GRB_BINARY);
		}
	}
	GRBVar** y = new GRBVar*[T];
	for(int i=0; i<T; ++i){
		y[i] = new GRBVar[2*g->size];
		for(int j=0; j<2*g->size; ++j){
			y[i][j] = model->addVar(0,1,0,GRB_BINARY);
		}
	}
	GRBVar* z = new GRBVar[T];
	for(int i=0; i<T; ++i){
		z[i] = model->addVar(0,1,0,GRB_BINARY);
	}
	// gurobi set objective
	double eps = 1.0/(2.0*T);
	GRBLinExpr expr = x[0][0];
	for(int j=1; j<g->order; ++j){
		expr += x[0][j];
	}
	for(int i=0; i<T; ++i){
		expr += eps*z[i];
	}
	model->setObjective(expr,GRB_MINIMIZE);
	// constraint 1
	for(int v=0; v<g->order; ++v){
		expr =  x[0][v];
		for(int i=0; i<T; ++i){
			for(auto iter=g->adj[v].cbegin(); iter!=g->adj[v].cend(); ++iter){
				int j = arcMap[pair<int,int>(*iter,v)];
				expr += y[i][j];
			}
		}
		model->addConstr(expr==1);
	}
	// constraints 2 and 3
	for(auto a=arcMap.cbegin(); a!=arcMap.cend(); ++a){
		int u = (a->first).first, v = (a->first).second, j = a->second;
		for(int i=0; i<T; ++i){
			expr = y[i][j] - x[i][u];
			model->addConstr(expr<=0);
		}
		for(auto iter=g->adj[u].cbegin(); iter!=g->adj[u].cend(); ++iter){
			if(*iter != v){
				for(int i=0; i<T; ++i){
					expr = y[i][j] - x[i][*iter];
					model->addConstr(expr<=0);
				}
			}
		}
	}
	// constraint 4
	for(int v=0; v<g->order; ++v){
		for(int i=0; i<T; ++i){
			expr = x[i+1][v] - x[i][v];
			for(auto iter=g->adj[v].cbegin(); iter!=g->adj[v].cend(); ++iter){
				int j = arcMap[pair<int,int>(*iter,v)];
				expr -= y[i][j];
			}
			model->addConstr(expr==0);
		}
	}
	// constraint 5
	for(int i=0; i<T; ++i){
		for(auto a=arcMap.cbegin(); a!=arcMap.cend(); ++a){
			int u = (a->first).first, v = (a->first).second;
			expr = x[i][u] - x[i][v];
			for(auto iter=g->adj[u].cbegin(); iter!=g->adj[u].cend(); ++iter){
				if(*iter != v){
					expr += x[i][*iter];	
				}
			}
			for(auto iter=g->adj[v].cbegin(); iter!=g->adj[v].cend(); ++iter){
				int j = arcMap[pair<int,int>(*iter,v)];
				expr -= y[i][j];
			}
			model->addConstr(expr<=g->deg(u)-1);
		}
	}
	// constraint 6
	eps = 1.0/g->order;
	for(int i=1; i<=T; ++i){
		expr = eps*(x[i][0] - x[i-1][0]);
		for(int j=1; j<g->order; ++j){
			expr += eps*(x[i][j] - x[i-1][j]);
		}
		expr -= z[i-1];
		model->addConstr(expr<=0);
	}
	// constraint 7
	for(int i=1; i<=T; ++i){
		expr = z[i-1] - (x[i][0] - x[i-1][0]);
		for(int j=1; j<g->order; ++j){
			expr -= (x[i][j] - x[i-1][j]);
		}
		model->addConstr(expr<=0);
	}
	// constraint 8
	int k = 0;
	expr = z[0];
	for(int i=1; i<T; ++i){
		expr += z[i];
	}
	model->addConstr(expr>=k);
	// generate all integers in prop time interval
	model->optimize();
	data.ptMap = new map<int,unordered_set<int>>;
	unordered_set<int> zf_set;
	while(model->get(GRB_IntAttr_Status)==GRB_OPTIMAL){
		// store propagation time and minimum zero forcing set
		k = 0;
		for(int i=0; i<T; ++i){
			if(z[i].get(GRB_DoubleAttr_X)>0.5){
				k += 1;
			}
		}
 		for(int i=0; i<g->order; ++i){
			if(x[0][i].get(GRB_DoubleAttr_X) > 0.5){
				zf_set.insert(i);
			}
		}
		data.ptMap->insert(pair<int,unordered_set<int>>(k,zf_set));
		zf_set.clear();
		// update model
		k += 1;
		model->addConstr(expr>=k);
		// re-optimize
		model->optimize();
	}
	data.status = model->get(GRB_IntAttr_Status);
	// free memory
	delete env;
	delete model;
	for(int i=0; i<=T; ++i){
		delete[] x[i];
	}
	delete[] x;
	for(int i=0; i<T; ++i){
		delete[] y[i];
	}
	delete[] y;
	delete[] z;
}
/*
	all_minimal_forts: recursive model for generating all minimal forts of given graph
*/
void all_minimal_forts(Graph *g,amfIPdata &data){
	// gurobi model
	GRBEnv* env = new GRBEnv(true);
	env->set(GRB_IntParam_OutputFlag,0);
	env->set(GRB_IntParam_LogToConsole,0);
	env->set(GRB_DoubleParam_TimeLimit,IP_MAX_TIME);
	env->start();
	GRBModel* model = new GRBModel(*env);
	// gurobi variables
	GRBVar* x = new GRBVar[g->order];
	for(int i=0; i<g->order; ++i){
		x[i] = model->addVar(0,1,1,GRB_BINARY);
	}
	// gurobi set objective
	model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
	// gurobi static constraints
	GRBLinExpr expr = 0;
	for(int i=0; i<g->order; ++i){
		expr += x[i];
	}
	model->addConstr(expr,GRB_GREATER_EQUAL,1);
	for(int i=0; i<g->order; ++i){
		for(auto iter1=g->adj[i].cbegin(); iter1!=g->adj[i].cend(); ++iter1){
			expr = x[*iter1] - x[i];
			for(auto iter2=g->adj[*iter1].cbegin(); iter2!=g->adj[*iter1].cend(); ++iter2){
				if(*iter2 != i){
					expr += x[*iter2];
				}
			}
			model->addConstr(expr,GRB_GREATER_EQUAL,0);
		}
	}
	// generate all minimal forts
	model->optimize();
	data.mforts = new vector<unordered_set<int>>;
	unordered_set<int> fort;
	while(model->get(GRB_IntAttr_Status)==GRB_OPTIMAL){
		// store minimal fort
		for(int i=0; i<g->order; ++i){
			if(x[i].get(GRB_DoubleAttr_X) > 0.5){
				fort.insert(i);
			}
		}
		// add to collection
		data.mforts->push_back(fort);
		// add constraint
		expr = 0;
		for(auto iter=fort.cbegin(); iter!=fort.cend(); ++iter){
			expr += x[*iter];
		}
		model->addConstr(expr,GRB_LESS_EQUAL,fort.size()-1);
		// delete fort
		fort.clear();
		// re-optimize
		model->optimize();
	}
	data.status = model->get(GRB_IntAttr_Status);
	// free memory
	delete env;
	delete model;
	delete[] x;
}
/*
	ft_num_ip: ip model for generating a maximum collection of disjoint forts of the given graph
*/
void ft_num_ip(const Graph *g,ftIPdata &mydata){
	// gurobi model
	GRBEnv* env = new GRBEnv(true);
	env->set(GRB_IntParam_OutputFlag,0);
	env->set(GRB_IntParam_LogToConsole,0);
	env->set(GRB_DoubleParam_TimeLimit,IP_MAX_TIME);
	env->start();
	GRBModel* model = new GRBModel(*env);
	// gurobi variables
	GRBVar** x = new GRBVar*[g->order];
	GRBVar* z = new GRBVar[g->order];
	for(int i=0; i<g->order; ++i){
		x[i] = new GRBVar[g->order];
		z[i] = model->addVar(0,1,0,GRB_BINARY);
		for(int j=0; j<g->order; ++j){
			x[i][j] = model->addVar(0,1,0,GRB_BINARY);
		}
	}
	// gurobi objective
	GRBLinExpr expr = z[0];
	for(int i=1; i<g->order; ++i){
		expr += z[i];
	}
	model->setObjective(expr,GRB_MAXIMIZE);
	//  constraint 1: row i of x is non-empty if z[i] is turned on
	for(int i=0; i<g->order; ++i){
		expr = z[i];
		for(int u=0; u<g->order; ++u){
			expr -= x[i][u];
		}
		model->addConstr(expr,GRB_LESS_EQUAL,0);
	}
	// constraint 2: each row of x must represent a fort or an empty set
	for(int i=0; i<g->order; ++i){
		for(int v=0; v<g->order; ++v){
			for(auto iter1=g->adj[v].cbegin(); iter1!=g->adj[v].cend(); ++iter1){
				expr = x[i][*iter1] - x[i][v];
				for(auto iter2=g->adj[*iter1].cbegin(); iter2!=g->adj[*iter1].cend(); ++iter2){
					if(*iter2 != v){
						expr += x[i][*iter2];
					}
				}
				model->addConstr(expr,GRB_GREATER_EQUAL,0);
			}
		}
	}
	// constraint 3: no vertx is in more than 1 fort
	for(int u=0; u<g->order; ++u){
		expr = x[0][u];
		for(int i=1; i<g->order; ++i){
			expr += x[i][u];
		}
		model->addConstr(expr,GRB_LESS_EQUAL,1);
	}
	// gurobi optimize
	model->optimize();
	// extract solution
	mydata.status = model->get(GRB_IntAttr_Status);
	mydata.val = round(model->get(GRB_DoubleAttr_ObjVal));
	mydata.dforts = new vector<unordered_set<int>>;
	unordered_set<int> fort;
	for(int i=0; i<g->order; ++i){
		if(z[i].get(GRB_DoubleAttr_X) > 0.5){
			for(int j=0; j<g->order; ++j){
				if(x[i][j].get(GRB_DoubleAttr_X) > 0.5){
					fort.insert(j);
				}
			}
			mydata.dforts->push_back(fort);
			fort.clear();
		}
	}
	// free memory
	delete env;
	delete model;
	delete[] z;
	for(int i=0; i<g->order; ++i){
		delete[] x[i];
	}
	delete[] x;
}
/*
	main function
*/
//int main(){
	/* hypercube graph
	const int d = 4;
	Graph *g = hypercube_graph(d);

	amfIPdata amfdata;
	all_minimal_forts(g,amfdata);
	cout << "# minimal forts: " << amfdata.mforts->size() << endl;
	for(int n = d; n<=12; ++n){
		int count = 0;
		for(auto it=amfdata.mforts->cbegin(); it!=amfdata.mforts->cend(); ++it){
			if(it->size() == n){
				count ++;
			}
		}
		cout << "# minimal forts of size " << n << " : " << count << endl;
	}
	return 0;
	*/
	/* Tracy Graph
	const int n = 12;
	Graph *g = new Graph;
	g->order = n;
	g->adj = new unordered_set<int>[n];

	g->addEdge(0,4);
	g->addEdge(0,5);
	g->addEdge(1,4);
	g->addEdge(1,5);
	g->addEdge(1,6);
	g->addEdge(2,6);
	g->addEdge(2,7);
	g->addEdge(3,6);
	g->addEdge(3,7);

	g->addEdge(4,8);
	g->addEdge(4,9);
	g->addEdge(5,8);
	g->addEdge(5,9);
	g->addEdge(5,10);
	g->addEdge(6,10);
	g->addEdge(6,11);
	g->addEdge(7,10);
	g->addEdge(7,11);

	ftIPdata ftdata;
	fzfIPdata fzfdata;
	zfIPdata zfdata;

	fort_cover_ip(g,zfdata);
	cout << "Z = " << zfdata.val << endl;

	fzf_ip(g,fzfdata);
	cout << "Z* = " << fzfdata.val << endl;

	ft_num_ip(g,ftdata);
	cout << "ft = " << ftdata.val << endl;

	for(auto it = ftdata.dforts->cbegin(); it!= ftdata.dforts->cend(); ++it){
		for(auto u = it->cbegin(); u!= it->cend(); ++u){
			cout << *u << " ";
		}
		cout << endl;
	}
	return 0;
	*/
//}