#include "mc.h"

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

    // std::string file_path = "/home/lyuheng/TKDE-revise/TthinkerQ/data/data_graphs/amazon/com-amazon.ungraph.graph";
    char * file_path_c    = "/home/lyuheng/TKDE-revise/TthinkerQ/app_maximalclique/test_graph.graph";

    int num_compers = 32;

    MCWorker worker(num_compers);
    worker.load_data(file_path_c);

    worker.run();

    return 0;
}

