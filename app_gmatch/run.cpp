#include "app_gmatch.h"
#include "FormatGraph.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>

#define GEN_GRAPH

using namespace std;

// #include <chrono>
// using namespace std::chrono;

int main(int argc, char *argv[])
{

#ifdef GEN_GRAPH
    // char *file_path = "/home/lyuan/graph_data/GSE10158_q";
    // FormatGraph fg;
    // fg.loadGraphFromFile(file_path);
    // fg.graphDump();

    std::string file_path = "xxx";
    std::string query_path = "xxx";
    FormatGraph g(file_path);
    g.Preprocess();
    g.writeGraphFile(file_path);
    g.sampleQueryGraph(query_path);
#else
    char *file_path = "/home/lyuan/graph_data/gmatch_data/GSE1730.graph";
    int num_compers = 32;

    GMWorker worker(num_compers);
    worker.load_data(file_path);

    worker.run();
#endif

    return 0;
}

