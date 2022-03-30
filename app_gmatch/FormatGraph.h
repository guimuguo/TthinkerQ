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

#define LABELSIZE 4
#define QUERY_GRAPH_SIZE 8
#define QUERY_GRAPH_NUM 20
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

    void graphDump()
    {   
        ofstream fout;
        char file[200];
        strcpy(file, "/home/lyuan/graph_data/gmatch_data/hyves.graph");
        fout.open(file);
        int *labels = new int[mnum_of_vertices];
        
        fout<<"t "<<mnum_of_vertices<<" "<<mnum_of_edges<<endl;

        for(int i=0; i<mnum_of_vertices; i++) {
            int label = rand() % LABELSIZE; // [0,LS-1]
            labels[i] = label;
            int degree = mppadj_lists[i][0];
            fout<<"v "<<i<<" "<<label<<" "<<degree<<endl;
        }

        for(int i=0; i<mnum_of_vertices; i++) {
            for(int j=1; j<=mppadj_lists[i][0]; j++) {
                if(mppadj_lists[i][j]>i) {
                    fout<<"e "<<i<<" "<<mppadj_lists[i][j]<<endl; 
                }
            }
        }
        fout.close();

        // Generate Query Graph

        unordered_set<int> contained_vertices;
        unordered_map<int, vector<int> > query_adj;
        unordered_map<int, int> vid2int;
        unordered_map<int, int> int2vid;
        int newId = 0, num_q_edges=0;
        int qcnt = 0;
        
        while(qcnt < QUERY_GRAPH_NUM) {

            sprintf(file, "/home/lyuan/TthinkerQ/client/console_version/hyves/query_graph_%d_%d.graph", QUERY_GRAPH_SIZE, qcnt);
            fout.open(file);

            int rand_vtx = rand()%mnum_of_vertices;
            contained_vertices.insert(rand_vtx);
            int nxt_rand_vtx;
            bool try_flag = true;

            while(contained_vertices.size() < QUERY_GRAPH_SIZE) 
            {
                int deg = mppadj_lists[rand_vtx][0];
                if(deg == 0) {
                    try_flag = false;
                    break;
                }
                nxt_rand_vtx = mppadj_lists[rand_vtx][rand()%deg+1];
                if(contained_vertices.find(nxt_rand_vtx) != contained_vertices.end()) {
                    try_flag = false;
                    break;
                }
                contained_vertices.insert(nxt_rand_vtx);
                rand_vtx = nxt_rand_vtx;
            }
            if(!try_flag) {
                contained_vertices.clear();
                fout.close();
                continue;
            }
            
            for(auto& ve: contained_vertices)
            {
                vid2int[ve] = newId;
                int2vid[newId] = ve;
                newId++;

                for(int j=1; j<=mppadj_lists[ve][0]; j++) {
                    int ne = mppadj_lists[ve][j];
                    if(ne > ve && contained_vertices.find(ne) != contained_vertices.end()) {
                        double r = (rand()%100)/100;
                        if(r < ADD_PROB) {
                            if(query_adj.find(ve) == query_adj.end())
                                query_adj[ve] = vector<int>();
                            if(query_adj.find(ne) == query_adj.end())
                                query_adj[ne] = vector<int>();
                            query_adj[ve].push_back(ne);
                            query_adj[ne].push_back(ve);
                            
                            num_q_edges++;
                        }
                    }
                }
            }
            for(auto& ve: contained_vertices)
                sort(query_adj[ve].begin(), query_adj[ve].end(), [&vid2int](const int& x, const int& y) {
                    return vid2int[x] < vid2int[y];
                });

            if(num_q_edges < 10) 
            {
                contained_vertices.clear();
                fout.close();
                query_adj.clear();
                int2vid.clear();
                vid2int.clear();
                num_q_edges = 0;
                newId = 0;
                continue;
            }

            fout<<"t "<<QUERY_GRAPH_SIZE<<" "<<num_q_edges<<endl;
            for(int i=0; i<QUERY_GRAPH_SIZE; i++) 
            {
                int ve = int2vid[i];
                int deg = query_adj[ve].size();
                fout<<"v "<<i<< " "<< labels[ve] << " " <<deg<<endl;
            }
            for(int i=0; i<QUERY_GRAPH_SIZE; i++) 
            {
                int ve = int2vid[i];
                for(auto &ne: query_adj[ve]) {
                    if(vid2int[ne] > i) {
                        fout<<"e "<<i<<" "<<vid2int[ne]<<endl;
                    }
                }
            }
            fout.close();
            contained_vertices.clear();
            query_adj.clear();
            int2vid.clear();
            vid2int.clear();
            num_q_edges = 0;
            newId = 0;
            qcnt++;
        } 

        delete[] labels;
    }
};