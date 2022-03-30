/*this code is used to output the graph in required format
the required format for SCS Application is:
vertex_id neighbor_id

//example
0 1
0 2
0 3
1 2
1 3
*/

#pragma once

#include "data.h"
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <stdlib.h> 
#include <fstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>

#define LABELSIZE 80
#define QUERY_GRAPH_SIZE 32
#define QUERY_GRAPH_NUM 10
#define ADD_PROB 1

using namespace std;

int comp_int(const void *e1, const void *e2)
{
    int n1, n2;
    n1 = *(int *) e1;
    n2 = *(int *) e2;
    if (n1>n2)
        return 1;
    else if (n1<n2)
        return -1;
    else
        return 0;
}


class FormatGraph
{ 
public: 

    int mnum_of_vertices, mnum_of_edges;
    int **mppadj_lists, *mpadj_list_buf;

    int loadGraphFromFile(char *file_path)
    {  
        Data *pdata;
        Transaction *ptransaction;
        int num_of_vertices, nbuf_size, nvertex_no, nbuf_pos;
        int *padj_lens, nsize, *ptemp, nlist_len, i, nmax_deg;

        pdata = new Data(file_path);

        nsize = 5000; // capacity of padj_lens
        padj_lens = new int[nsize]; // buffer for keeping each transaction for parsing
        memset(padj_lens, 0, sizeof(int)*nsize);

        num_of_vertices = 0;
        nbuf_size = 0; // ---> all transactions are saved in one buffer.... (concat adj-lists)
        nmax_deg = 0; // track longest transaction, just for printing at Line 1855...
        ptransaction = pdata->getNextTransaction();
        while(ptransaction)
        {
            if(nmax_deg<ptransaction->length)
                nmax_deg = ptransaction->length;
            if(num_of_vertices>=nsize)
            {
                ptemp = new int[2*nsize]; // ptemp is used to double the capacity of padj_lens, (vector)
                memcpy(ptemp, padj_lens, sizeof(int)*nsize);
                memset(&ptemp[nsize], 0, sizeof(int)*nsize);
                delete []padj_lens;
                padj_lens = ptemp;
                nsize *= 2;
            }
            padj_lens[num_of_vertices] = ptransaction->length;
            num_of_vertices++;
            nbuf_size += ptransaction->length; // ---> all transactions are saved in one buffer....
            ptransaction = pdata->getNextTransaction(); // line by line transaction reading
        }

        mppadj_lists = new int*[num_of_vertices]; // mppadj_lists[] keeps the start position (pointer) of each vertex's adj-list in mpadj_list_buf
        mpadj_list_buf = new int[num_of_vertices+nbuf_size];
        nvertex_no = 0;
        nbuf_pos = 0;

        ptransaction = pdata->getNextTransaction(); // note the rewind operation inside getNextTransaction()
        while(ptransaction)
        {
            nlist_len = 1;
            mppadj_lists[nvertex_no] = &mpadj_list_buf[nbuf_pos];
            for(i=0; i<ptransaction->length; i++)
            {
                if(ptransaction->t[i]!=nvertex_no) // remove self-loop
                    mppadj_lists[nvertex_no][nlist_len++] = ptransaction->t[i];
            }
            nbuf_pos += nlist_len;
            mppadj_lists[nvertex_no][0] = nlist_len-1; // from position 1 onwards, neighbors are kept; position 0 keeps the number of neighbors (adj-list length)
            // temporarially don't need to sort.
            qsort(&mppadj_lists[nvertex_no][1], mppadj_lists[nvertex_no][0], sizeof(int), comp_int); // adj-lists are sorted !!!

            nvertex_no++;
            ptransaction = pdata->getNextTransaction();
        }
        delete pdata;
        delete []padj_lens;

        printf("maximum vertex degree: %d\n", nmax_deg); // max_deg just for printing...

        assert(nbuf_size%2 == 0);
        mnum_of_vertices = num_of_vertices;
        mnum_of_edges = nbuf_size/2;

        cout<<"# of vertices : "<<mnum_of_vertices<<endl;
        
        return num_of_vertices;
    }

    void graphDump(char * outfile_path)
    {   
        ofstream fout;
        char file[200];
        strcpy(file, outfile_path);
        fout.open(file);
        int *labels = new int[mnum_of_vertices];

        for(int i=0; i<mnum_of_vertices; i++) {
            int degree = mppadj_lists[i][0];
            if(degree > 0){
                for(int j=1; j<=degree; j++){
                    fout << i << " " << mppadj_lists[i][j] << "\n";
                }
            }
        }

        fout.close();

        delete[] labels;
    }
};