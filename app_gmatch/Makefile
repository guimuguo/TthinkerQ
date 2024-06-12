CCOMPILE=g++ 
CPPFLAGS=-lpthread -march=native -O3

all: run

run: run.cpp
	$(CCOMPILE) -std=c++14 -fopenmp -lpthread run.cpp intersection/computesetintersection.cpp $(CPPFLAGS)  -o run

clean:
	-rm run
