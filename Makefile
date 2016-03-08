# LIBS = -lm -lpthread  -DIL_STD

# FILES = vne_reallocation.cc cplex_solver.cc util.cc vne_solution_builder.cc
FILES = util.cc vnr.cc
# LIBS = -lboost
LIBS = -lpthread
INCLUDES = -I/usr/local/include 

all:
	g++ -O3  $(LIB_PATHS) $(INCLUDE_PATHS) $(FILES) $(LIBS) -o vne_reallocation

dbg:
	g++ -DDBG -g  $(LIB_PATHS) $(INCLUDE_PATHS) $(FILES) $(LIBS) -o vne_reallocation

debug:
	g++ -g $(LIB_PATHS) $(INCLUDE_PATHS) $(FILES) $(LIBS) -o vne_reallocation
