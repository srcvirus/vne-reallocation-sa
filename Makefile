FILES = util.cc vnr.cc
LIBS = -lpthread
INCLUDE_PATHS =

all:
	g++ -O3  $(LIB_PATHS) $(INCLUDE_PATHS) $(FILES) $(LIBS) -o vne_reallocation

dbg:
	g++ -DDBG -g  $(LIB_PATHS) $(INCLUDE_PATHS) $(FILES) $(LIBS) -o vne_reallocation

debug:
	g++ -g $(LIB_PATHS) $(INCLUDE_PATHS) $(FILES) $(LIBS) -o vne_reallocation
