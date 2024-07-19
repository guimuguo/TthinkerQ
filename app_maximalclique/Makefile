CCOMPILE=g++ 
CPPFLAGS=-lpthread -march=native -g

all: run

run: run.cpp
	$(CCOMPILE) -std=c++14 -fopenmp -lpthread run.cpp $(CPPFLAGS)  -o run

clean:
	rm run
