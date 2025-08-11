#####################################################
#   Zero-Forcing IP Models Make File 				#
#####################################################
include make.inc

graph: graph.hpp graph.cpp
	$(CPP) $(CARGS) -o graph graph.cpp $(INC) $(LIB)

ip_models: graph.hpp graph.cpp ip_models.hpp ip_models.cpp
	$(CPP) $(CARGS) -o ip_models graph.cpp ip_models.cpp $(INC) $(LIB)

testing: testing.hpp testing.cpp
	$(CPP) $(CARGS) -o testing testing.cpp $(INC) $(LIB)

smallg_test: graph.hpp graph.cpp ip_models.hpp ip_models.cpp testing.hpp testing.cpp smallg_test.cpp
	$(CPP) $(CARGS) -o smallg_test graph.cpp ip_models.cpp testing.cpp smallg_test.cpp $(INC) $(LIB)

randg_test: graph.hpp graph.cpp ip_models.hpp ip_models.cpp testing.hpp testing.cpp randg_test.cpp
	$(CPP) $(CARGS) -o randg_test graph.cpp ip_models.cpp testing.cpp randg_test.cpp $(INC) $(LIB)

bhhs_test: graph.hpp graph.cpp ip_models.hpp ip_models.cpp bhhs_test.cpp
	$(CPP) $(CARGS) -o bhhs_test graph.cpp ip_models.cpp bhhs_test.cpp $(INC) $(LIB)

bhhs_vs_test: graph.hpp graph.cpp ip_models.hpp ip_models.cpp bhhs_vs_test.cpp
	$(CPP) $(CARGS) -o bhhs_vs_test graph.cpp ip_models.cpp bhhs_vs_test.cpp $(INC) $(LIB)

bhhs_es_test: graph.hpp graph.cpp ip_models.hpp ip_models.cpp bhhs_es_test.cpp
	$(CPP) $(CARGS) -o bhhs_es_test graph.cpp ip_models.cpp bhhs_es_test.cpp $(INC) $(LIB)

tree_fort_test: graph.hpp graph.cpp ip_models.hpp ip_models.cpp testing.hpp testing.cpp tree_fort_test.cpp
	$(CPP) $(CARGS) -o tree_fort_test graph.cpp ip_models.cpp testing.cpp tree_fort_test.cpp $(INC) $(LIB)

hypercube_ptinterval_test: graph.hpp graph.cpp ip_models.hpp ip_models.cpp hypercube_ptinterval_test.cpp
	$(CPP) $(CARGS) -o hypercube_ptinterval_test graph.cpp ip_models.cpp hypercube_ptinterval_test.cpp $(INC) $(LIB)
	
uninstall:
	@rm -f graph
	@rm -f ip_models
	@rm -f testing
	@rm -f smallg_test
	@rm -f randg_test
	@rm -f bhhs_test
	@rm -f bhhs_vs_test
	@rm -f bhhs_es_test
	@rm -f tree_fort_test
	@rm -f hypercube_ptinterval_test