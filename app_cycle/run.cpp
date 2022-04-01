
#include "cycle_app.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>

using namespace std;

// #include <chrono>
// using namespace std::chrono;

int main(int argc, char *argv[])
{
    char *file_path = "../data/GSE1730_q";
    int num_compers = 32;

    CEWorker worker(num_compers);
    worker.load_data(file_path);

    worker.run();

    return 0;
}

