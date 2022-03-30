# include <iostream>
# include "format_graph.h"

using namespace std;

int main(int argc, char **argv){
    /*
        argv[1] = filepath of original graph file to be read
        argv[2] = filepath where the output graph is to be saved
    */
    if(argc < 2){
        cout << "arg1 = graph to be read, arg2 = output file path for new format graph" << endl;
        return -1;
    }

    FormatGraph fg;
    fg.loadGraphFromFile(argv[1]);
    fg.graphDump(argv[2]);

    return 0;
}