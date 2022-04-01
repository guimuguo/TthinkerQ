#include "hpcycle_app.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>


using namespace std;

#include <chrono>
using namespace std::chrono;

int main(int argc, char *argv[])
{
    const char *file_path = "../data/GSE1730_q";
    int num_compers = 32;

    HCEWorker worker(num_compers);

    auto start_t = steady_clock::now();

    worker.load_data(file_path, num_compers);

    auto end_t = steady_clock::now();

    float duration = (float)duration_cast<milliseconds>(end_t - start_t).count() / 1000;

    cout<<"Total build Index Graph time: "<<duration<<" seconds"<<endl;

    worker.run();

    return 0;
}

