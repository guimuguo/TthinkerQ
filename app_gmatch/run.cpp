#include "app_gmatch.h"
#include "FormatGraph.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>

// #define GEN_GRAPH

using namespace std;

// #include <chrono>
// using namespace std::chrono;

int main(int argc, char *argv[])
{

#ifdef GEN_GRAPH
    char *file_path = "/home/lyuan/graph_data/GSE10158_q";
    FormatGraph fg;
    fg.loadGraphFromFile(file_path);
    fg.graphDump();
#else
    char *file_path = "../data/amazon.graph";
    int num_compers = 32;

    GMWorker worker(num_compers);
    worker.load_data(file_path);

    worker.run();
#endif

    return 0;
}

