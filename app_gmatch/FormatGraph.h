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

#include <sys/time.h>

#include <cstdlib>
#include <iostream>

class Timer {
public:
    Timer() {}
    ~Timer() {}

    void StartTimer() { start_timestamp_ = wtime(); }
    void EndTimer() { end_timestamp_ = wtime(); }

    double GetElapsedMicroSeconds() const { return end_timestamp_ - start_timestamp_; }

    inline void PrintElapsedMicroSeconds(const std::string& time_tag) const {
        std::cout << std::fixed << "finish " << time_tag << ", elapsed_time=" << (end_timestamp_ - start_timestamp_) / 1000.0 << "ms" << std::endl;
    }

    double wtime() {
        double time[2];
        struct timeval time1;
        gettimeofday(&time1, NULL);

        time[0] = time1.tv_sec;
        time[1] = time1.tv_usec;

        return time[0] * (1.0e6) + time[1];
    }

private:
    double start_timestamp_;
    double end_timestamp_;
};



#define LABELSIZE 4
#define QUERY_GRAPH_SIZE 8
#define QUERY_GRAPH_NUM 20
#define ADD_PROB 1

/************************

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
*/

#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <queue>
#include <utility>
#include <fstream>
#include <string.h>

typedef unsigned int ui;
typedef unsigned long long int ull;

typedef unsigned int uintV;
typedef unsigned long long int uintE;


class FormatGraph
{
public:
    FormatGraph() {}
    FormatGraph(const std::string &filename);

    virtual ~FormatGraph() {}

    size_t GetEdgeCount() const { return edge_count_; }
    size_t GetVertexCount() const { return vertex_count_; }
    void SetVertexCount(size_t vertex_count) { vertex_count_ = vertex_count; }
    void SetEdgeCount(size_t edge_count) { edge_count_ = edge_count; }

    uintE *GetRowPtrs() const { return row_ptrs_; }
    uintV *GetCols() const { return cols_; }

    void SetRowPtrs(uintE *row_ptrs) { row_ptrs_ = row_ptrs; }
    void SetCols(uintV *cols) { cols_ = cols; }

    void GetMaxDegree();

    void readBinFile(const std::string &filename);
    void readGraphFile(const std::string &filename);
    void readSnapFile(const std::string &filename);
    void readLVIDFile(const std::string &filename);

    void writeBinFile(const std::string &filename);
    void writeGraphFile(const std::string &filename);

    void writePeregrineFile(const std::string &filename);

    void writeQuasiCliqueFile(const std::string &filename);

    void Preprocess();

    void ReMapVertexId();

    void sampleQueryGraph(const std::string &filename);


private:
    uintE *row_ptrs_;
    uintV *cols_;
    size_t vertex_count_;
    size_t edge_count_;

    std::vector<int> labels;
};

FormatGraph::FormatGraph(const std::string &filename)
{
    std::string suffix = filename.substr(filename.rfind(".") + 1);
    if (suffix == "bin")
        readBinFile(filename);
    else if (suffix == "graph")
        readGraphFile(filename);
    else if (suffix == "txt")
        readSnapFile(filename);
    else if (suffix == "lvid")
        readLVIDFile(filename);
    else {
        std::cout << "Cannot read graph file based on its suffix ..." << std::endl;
        assert(false);
    }
}

void FormatGraph::GetMaxDegree()
{
    ull max_deg = 0;
    for (size_t i = 0; i < vertex_count_; ++i)
    {
        max_deg = std::max(max_deg, row_ptrs_[i + 1] - row_ptrs_[i]);
    }
    std::cout << "max degree=" << max_deg << std::endl;
}

void FormatGraph::ReMapVertexId()
{
    // 1. find the root vertex with largest degree
    size_t max_degree = 0;
    uintV root = 0;
    for (uintV i = 0; i < vertex_count_; i++) {
        if (row_ptrs_[i + 1] - row_ptrs_[i] > max_degree) {
            max_degree = row_ptrs_[i + 1] - row_ptrs_[i];
            root = i;
        }
    }
    // 2. bfs from the root vertex, make sure connected
    //    order: higher degree, more connections to the visited vertices
    std::queue<uintV> queue;
    std::vector<bool> visited(vertex_count_, false);
    queue.push(root);
    visited[root] = true;
    uintV new_vid = 0;
    std::vector<uintV> old_to_new(vertex_count_);
    std::vector<uintV> new_to_old(vertex_count_);
    uintE *new_row_ptrs_ = new uintE[vertex_count_ + 1];
    uintV *new_cols_ = new uintV[edge_count_];

    while (!queue.empty()) {
        size_t size = queue.size();
        std::vector<uintV> same_level_vertices;
        for (size_t i = 0; i < size; i++) {
            uintV front = queue.front();
            same_level_vertices.push_back(front);
            queue.pop();
            for (size_t j = row_ptrs_[front]; j < row_ptrs_[front + 1]; ++j)
            {  
                uintV ne = cols_[j];
                if (!visited[ne]) {
                    visited[ne] = true;
                    queue.push(ne);
                }
            }
        }
        std::vector<std::tuple<size_t, size_t, uintV>> weights;  // degree, connections, vid
        for (size_t i = 0; i < size; i++) {
            uintV v = same_level_vertices[i];
            size_t connections = 0;
            for (size_t j = row_ptrs_[v]; j < row_ptrs_[v + 1]; ++j)
            {
                uintV ne = cols_[j];
                if (visited[ne])
                    connections++;
            }
            weights.emplace_back(row_ptrs_[v + 1] - row_ptrs_[v], connections, v);
        }
        std::sort(weights.begin(), weights.end(), [](const auto& a, const auto& b) {
            if (std::get<0>(a) != std::get<0>(b))
                return std::get<0>(a) > std::get<0>(b);
            else if (std::get<1>(a) != std::get<1>(b))
                return std::get<1>(a) > std::get<1>(b);
            else if (std::get<2>(a) != std::get<2>(b))
                return std::get<2>(a) < std::get<2>(b);
            return false;
        });
        for (const auto& w : weights) {
            old_to_new[std::get<2>(w)] = new_vid;
            new_to_old[new_vid] = std::get<2>(w);
            new_vid++;
        }
    }
    auto offsets = new uintE[vertex_count_ + 1];
    memset(offsets, 0, sizeof(uintE) * (vertex_count_ + 1));

    for (size_t i = 0; i < vertex_count_; ++i)
    {
        offsets[i] = row_ptrs_[new_to_old[i] + 1] - row_ptrs_[new_to_old[i]];
    }
    uintE prefix = 0;
    for (size_t i = 0; i < vertex_count_ + 1; ++i) {
        new_row_ptrs_[i] = prefix;
        prefix += offsets[i];
        offsets[i] = new_row_ptrs_[i];
    }
    for (size_t i = 0; i < vertex_count_; ++i)
    {
        for (size_t j = row_ptrs_[i]; j < row_ptrs_[i + 1] ; ++j)
        {
            uintV ne = cols_[j];
            new_cols_[offsets[old_to_new[i]]++] = old_to_new[ne];
        }
    }
    for (uintV u = 0; u < vertex_count_; ++u) {
        std::sort(new_cols_ + new_row_ptrs_[u], new_cols_ + new_row_ptrs_[u + 1]);
    }
    delete[] offsets;
    delete[] row_ptrs_;
    delete[] cols_;
    row_ptrs_ = new_row_ptrs_;
    cols_ = new_cols_;
}

void FormatGraph::readSnapFile(const std::string &filename)
{
    Timer timer;
    timer.StartTimer();

    vertex_count_ = 0;
    edge_count_ = 0;
    row_ptrs_ = NULL;
    cols_ = NULL;

    uintV min_vertex_id = std::numeric_limits<uintV>::max();
    uintV max_vertex_id = std::numeric_limits<uintV>::min();
    std::ifstream file(filename.c_str(), std::fstream::in);
    std::string line;
    uintV vids[2];
    while (getline(file, line)) {
        if (line.length() == 0 || !std::isdigit(line[0]))
            continue;
        std::istringstream iss(line);
        for (int i = 0; i < 2; ++i) {
            iss >> vids[i];
            min_vertex_id = std::min(min_vertex_id, vids[i]);
            max_vertex_id = std::max(max_vertex_id, vids[i]);
        }
        edge_count_++;
    }
    file.close();

    vertex_count_ = max_vertex_id - min_vertex_id + 1;
    edge_count_ *= 2;
    std::cout << "vertex_count=" << vertex_count_ << ", edge_count=" << edge_count_ << std::endl;

    row_ptrs_ = new uintE[vertex_count_ + 1];
    cols_ = new uintV[edge_count_];
    auto offsets = new uintE[vertex_count_ + 1];
    memset(offsets, 0, sizeof(uintE) * (vertex_count_ + 1));
    
    {
        std::ifstream file(filename.c_str(), std::fstream::in);
        std::string line;
        uintV vids[2];
        while (getline(file, line)) {
            if (line.length() == 0 || !std::isdigit(line[0]))
                continue;
            std::istringstream iss(line);
            for (int i = 0; i < 2; ++i)
            {
                iss >> vids[i];
                vids[i] -= min_vertex_id;
            }
            offsets[vids[0]]++;
            offsets[vids[1]]++;
        }
        file.close();
    }
    
    uintE prefix = 0;
    for (size_t i = 0; i < vertex_count_ + 1; ++i) {
        row_ptrs_[i] = prefix;
        prefix += offsets[i];
        offsets[i] = row_ptrs_[i];
    }

    {
        std::ifstream file(filename.c_str(), std::fstream::in);
        std::string line;
        uintV vids[2];
        while (getline(file, line)) {
            if (line.length() == 0 || !std::isdigit(line[0]))
                continue;
            std::istringstream iss(line);
            for (int i = 0; i < 2; ++i)
            {
                iss >> vids[i];
                vids[i] -= min_vertex_id;
            }
            cols_[offsets[vids[0]]++] = vids[1];
            cols_[offsets[vids[1]]++] = vids[0];
        }
        file.close();
    }
    delete[] offsets;
    offsets = NULL;

    for (uintV u = 0; u < vertex_count_; ++u) {
        std::sort(cols_ + row_ptrs_[u], cols_ + row_ptrs_[u + 1]);
    }

    timer.EndTimer();
    timer.PrintElapsedMicroSeconds("reading CSR Snap file");
}

void FormatGraph::readBinFile(const std::string &filename)
{
    vertex_count_ = 0;
    edge_count_ = 0;
    row_ptrs_ = NULL;
    cols_ = NULL;

    Timer timer;
    timer.StartTimer();
    std::cout << "start reading CSR bin file...." << std::endl;
    FILE* file_in = fopen(filename.c_str(), "rb");
    assert(file_in != NULL);
    fseek(file_in, 0, SEEK_SET);
    size_t res = 0;
    size_t uintV_size = 0, uintE_size = 0;
    res += fread(&uintV_size, sizeof(size_t), 1, file_in);
    res += fread(&uintE_size, sizeof(size_t), 1, file_in);
    res += fread(&vertex_count_, sizeof(size_t), 1, file_in);
    res += fread(&edge_count_, sizeof(size_t), 1, file_in);
    assert(uintV_size == sizeof(uintV));
    assert(uintE_size == sizeof(uintE));
    std::cout << "vertex_count=" << vertex_count_ << ", edge_count=" << edge_count_ << std::endl;

    row_ptrs_ = new uintE[vertex_count_ + 1];
    cols_ = new uintV[edge_count_];
    res += fread(row_ptrs_, sizeof(uintE), vertex_count_ + 1, file_in);
    res += fread(cols_, sizeof(uintV), edge_count_, file_in);

    assert(res == 4 + (vertex_count_ + 1) + edge_count_);

    GetMaxDegree();

    fgetc(file_in);
    assert(feof(file_in));
    fclose(file_in);

    timer.EndTimer();
    timer.PrintElapsedMicroSeconds("reading CSR bin file");  
}

void FormatGraph::readGraphFile(const std::string &filename)
{
    vertex_count_ = 0;
    edge_count_ = 0;
    row_ptrs_ = NULL;
    cols_ = NULL;

    std::ifstream file_in(filename);
    if (!file_in.is_open()) 
    {
        std::cout << "Unable to read the graph " << filename << std::endl;
    }

    char type;
    file_in >> type >> vertex_count_ >> edge_count_;
    edge_count_ *= 2;
    std::cout << "vertex_count=" << vertex_count_ << ", edge_count=" << edge_count_ << std::endl;

    row_ptrs_ = new uintE[vertex_count_ + 1];
    cols_ = new uintV[edge_count_];
    row_ptrs_[0] = 0;

    std::vector<uintV> neighbors_offset(vertex_count_, 0);
    
    while (file_in >> type) {
        if (type == 'v') { // Read vertex.
            ui id;
            ui label;
            ui degree;
            file_in >> id >> label >> degree;
            row_ptrs_[id + 1] = row_ptrs_[id] + degree;
        }
        else if (type == 'e') { // Read edge.
            uintV begin;
            uintV end;
            file_in >> begin >> end;

            uintV offset = row_ptrs_[begin] + neighbors_offset[begin];
            cols_[offset] = end;
            offset = row_ptrs_[end] + neighbors_offset[end];
            cols_[offset] = begin;
            neighbors_offset[begin] += 1;
            neighbors_offset[end] += 1;
        }
    }

    file_in.close();

    for (ui i = 0; i < vertex_count_; ++i) {
        std::sort(cols_ + row_ptrs_[i], cols_ + row_ptrs_[i + 1]);
    }
}

void FormatGraph::readLVIDFile(const std::string &filename)
{
    vertex_count_ = 0;
    edge_count_ = 0;
    row_ptrs_ = NULL;
    cols_ = NULL;

    std::cout << "start build csr..." << std::endl;
    // const char* kDelimiters = " ,;\t";
    const char* kDelimiters = "0123456789";
    std::unordered_map<std::string, uintV> ids;
    std::vector<uintV> edge_pairs;
    {
        std::ifstream file(filename.c_str(), std::fstream::in);
        std::string line;
        while (getline(file, line)) {
            if (line.length() == 0 || !std::isdigit(line[0]))
                continue;

            std::vector<std::string> num_strs;
            size_t cur_pos = 0;
            while (cur_pos < line.length()) {
                cur_pos = line.find_first_of(kDelimiters, cur_pos, strlen(kDelimiters));
                if (cur_pos < line.length()) {
                    size_t next_pos = line.find_first_not_of(kDelimiters, cur_pos, strlen(kDelimiters));
                    num_strs.push_back(line.substr(cur_pos, next_pos - cur_pos));
                    assert(next_pos > cur_pos);
                    cur_pos = next_pos;
                }
            }

            for (auto& str : num_strs) {
                assert(str.length());
                for (auto ch : str) {
                    assert(std::isdigit(ch));
                }
            }

            for (auto& str : num_strs) {
                if (ids.find(str) == ids.end()) {
                    ids.insert(std::make_pair(str, vertex_count_++));
                }
                edge_pairs.push_back(ids[str]);
            }
        }
        file.close();
    }
    ids.clear();

    std::cout << "edge pairs size=" << edge_pairs.size() << std::endl;
    assert(edge_pairs.size() % 2 == 0);
    edge_count_ = edge_pairs.size(); // / 2;

    std::vector<uintE> offsets(vertex_count_ + 1, 0);
    for (size_t i = 0; i < edge_pairs.size(); i += 2) {
        offsets[edge_pairs[i]]++;
        offsets[edge_pairs[i + 1]]++;
    }

    row_ptrs_ = new uintE[vertex_count_ + 1];
    cols_ = new uintV[edge_count_];

    uintE prefix = 0;
    for (uintV i = 0; i <= vertex_count_; ++i) {
        row_ptrs_[i] = prefix;
        prefix += offsets[i];
        offsets[i] = row_ptrs_[i];
    }

    for (size_t i = 0; i < edge_pairs.size(); i += 2) 
    {
        cols_[offsets[edge_pairs[i]]++] = edge_pairs[i + 1];
        cols_[offsets[edge_pairs[i + 1]]++] = edge_pairs[i];
    }

    offsets.clear();
    edge_pairs.clear();

#pragma omp parallel for schedule(dynamic)
    for (uintV u = 0; u < vertex_count_; ++u) 
    {
        std::sort(cols_ + row_ptrs_[u], cols_ + row_ptrs_[u + 1]);
    }
    std::cout << "finish building CSR" << std::endl;
}

void FormatGraph::writeBinFile(const std::string &filename)
{
    std::string prefix = filename.substr(0, filename.rfind("."));

    Timer timer;
    timer.StartTimer();
    std::cout << "start write CSR bin file...." << std::endl;
    std::string output_filename = prefix + ".bin";
    FILE* file_out = fopen(output_filename.c_str(), "wb");
    assert(file_out != NULL);
    size_t res = 0;
    size_t uintV_size = sizeof(uintV), uintE_size = sizeof(uintE);
    res += fwrite(&uintV_size, sizeof(size_t), 1, file_out);
    res += fwrite(&uintE_size, sizeof(size_t), 1, file_out);
    res += fwrite(&vertex_count_, sizeof(size_t), 1, file_out);
    res += fwrite(&edge_count_, sizeof(size_t), 1, file_out);
    res += fwrite(row_ptrs_, sizeof(uintE), vertex_count_ + 1, file_out);
    res += fwrite(cols_, sizeof(uintV), edge_count_, file_out);

    assert(res == 4 + (vertex_count_ + 1) + edge_count_);
    fclose(file_out);
    timer.EndTimer();
    timer.PrintElapsedMicroSeconds("writing CSR bin file");
}

void FormatGraph::writeGraphFile(const std::string &filename)
{   
    std::string prefix = filename.substr(0, filename.rfind("."));

    Timer timer;
    timer.StartTimer();
    std::cout << "start write Graph file...." << std::endl;
    std::string output_filename = prefix + ".graph";
    // FILE* file_out = fopen(output_filename.c_str(), "w");
    // assert(file_out != NULL);

    std::ofstream file_out(output_filename);
    file_out << "t " << vertex_count_ << " " << edge_count_ / 2 << '\n';

    // std::vector<int> labels(vertex_count_);
    labels.resize(vertex_count_);

    for (uintV i = 0; i < vertex_count_; ++i)
    {
        int label = rand() % LABELSIZE; // [0,LS-1]
        labels[i] = label;
        file_out << "v " << i << " " << label << " " << row_ptrs_[i + 1] - row_ptrs_[i] << '\n';
    }

    for (uintV i = 0; i < vertex_count_; ++i)
    {
        for (uintE j = row_ptrs_[i]; j < row_ptrs_[i + 1]; ++j)
        {
            if (i < cols_[j])
                file_out << "e " << i << " " << cols_[j] << '\n';
        }
    }

    // fclose(file_out);
    file_out.close();

    timer.EndTimer();
    timer.PrintElapsedMicroSeconds("writing Graph file");
}

void FormatGraph::writePeregrineFile(const std::string &filename) 
{
    std::string prefix = filename.substr(0, filename.rfind("."));

    Timer timer;
    timer.StartTimer();
    std::cout << "start write Graph file...." << std::endl;
    std::string output_filename = prefix + ".peregrine.edges.txt";

    {
        std::ofstream file_out(output_filename);
        for (uintV i = 0; i < vertex_count_; ++i)
        {
            for (uintE j = row_ptrs_[i]; j < row_ptrs_[i + 1]; ++j)
            {
                if (i < cols_[j])
                    file_out << i << " " << cols_[j] << '\n';
            }
        }
        file_out.close();
    }

    
    output_filename = prefix + ".peregrine.labels.txt";
    {
        std::ofstream file_out(output_filename);
        for (uintV i = 0; i < vertex_count_; ++i)
        {
            file_out << i << " " << labels[i] << '\n';
        }
        file_out.close();
    }

    timer.EndTimer();
    timer.PrintElapsedMicroSeconds("writing Peregrine file");
}

void FormatGraph::writeQuasiCliqueFile(const std::string &filename)
{
    std::string prefix = filename.substr(0, filename.rfind("."));
    Timer timer;
    timer.StartTimer();
    std::cout << "start write QuasiClique file...." << std::endl;
    std::string output_filename = prefix + "_q";

    std::ofstream file_out(output_filename);
    for (uintV i = 0; i < vertex_count_; ++i)
    {
        for (uintE j = row_ptrs_[i]; j < row_ptrs_[i + 1]; ++j)
        {
            file_out << cols_[j] << ' ';
        }
        file_out << '\n';
    }
    file_out.close();
    
    timer.EndTimer();
    timer.PrintElapsedMicroSeconds("writing QuasiClique file");
}

void FormatGraph::sampleQueryGraph(const std::string &filename) {

    std::unordered_set<uintV> contained_vertices;
    std::unordered_map<uintV, std::vector<uintV> > query_adj;
    std::unordered_map<uintV, uintV> vid2int;
    std::unordered_map<uintV, uintV> int2vid;
    uintV newId = 0, num_q_edges=0;
    uintV qcnt = 0;

    std::ofstream fout;
    char file[100];
    
    while(qcnt < QUERY_GRAPH_NUM) {

        sprintf(file, "/query_graph_%d_%d.graph", QUERY_GRAPH_SIZE, qcnt);
        fout.open(filename + std::string(file));

        uintV rand_vtx = rand() % vertex_count_;
        contained_vertices.insert(rand_vtx);
        uintV nxt_rand_vtx;
        bool try_flag = true;

        while(contained_vertices.size() < QUERY_GRAPH_SIZE) 
        {
            uintV deg = row_ptrs_[rand_vtx + 1] - row_ptrs_[rand_vtx];
            if (deg == 0) {
                try_flag = false;
                break;
            }
            nxt_rand_vtx = cols_[row_ptrs_[rand_vtx] + rand() % deg];
            if (contained_vertices.find(nxt_rand_vtx) != contained_vertices.end()) {
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

            for(uintE j = row_ptrs_[ve]; j < row_ptrs_[ve+1]; j++) {
                uintV ne = cols_[j];
                if(ne > ve && contained_vertices.find(ne) != contained_vertices.end()) {
                    double r = (rand()%100)/100;
                    if(r < ADD_PROB) {
                        if (query_adj.find(ve) == query_adj.end())
                            query_adj[ve] = std::vector<uintV>();
                        if (query_adj.find(ne) == query_adj.end())
                            query_adj[ne] = std::vector<uintV>();
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

        if(num_q_edges < 15) 
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

        fout << "t "<< QUERY_GRAPH_SIZE <<" "<< num_q_edges << std::endl;
        for(int i=0; i<QUERY_GRAPH_SIZE; i++) 
        {
            int ve = int2vid[i];
            int deg = query_adj[ve].size();
            fout<<"v "<<i<< " "<< labels[ve] << " " <<deg<< std::endl;
        }
        for(int i=0; i<QUERY_GRAPH_SIZE; i++) 
        {
            int ve = int2vid[i];
            for(auto &ne: query_adj[ve]) {
                if(vid2int[ne] > i) {
                    fout<<"e "<<i<<" "<<vid2int[ne]<< std::endl;
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
}

void FormatGraph::Preprocess() {
    // remove dangling nodes
    // self loops
    // parallel edges
    Timer timer;
    timer.StartTimer();
    std::cout << "start preprocess..." << std::endl;
    size_t vertex_count = GetVertexCount();
    size_t edge_count = GetEdgeCount();
    auto row_ptrs = GetRowPtrs();
    auto cols = GetCols();

    auto vertex_ecnt = new uintE[vertex_count + 1];
    memset(vertex_ecnt, 0, sizeof(uintE) * (vertex_count + 1));
    for (uintV u = 0; u < vertex_count; ++u) 
    {
        for (auto j = row_ptrs[u]; j < row_ptrs[u + 1]; ++j) 
        {
            auto v = cols[j];
            bool parallel_edge = (j > row_ptrs[u] && v == cols[j - 1]);
            bool self_loop = u == v;
            if (!parallel_edge && !self_loop) 
            {
                vertex_ecnt[u]++;
            }
        }
    }
    auto nrow_ptrs = new uintE[vertex_count + 1];
    uintE nedge_count = 0;
    for (uintV u = 0; u < vertex_count; ++u) 
    {
        nrow_ptrs[u] = nedge_count;
        nedge_count += vertex_ecnt[u];
    }
    nrow_ptrs[vertex_count] = nedge_count;
    delete[] vertex_ecnt;
    vertex_ecnt = NULL;

    auto ncols = new uintV[nedge_count];
    for (uintV u = 0; u < vertex_count; ++u) 
    {
        auto uoff = nrow_ptrs[u];
        for (uintE j = row_ptrs[u]; j < row_ptrs[u + 1]; ++j) 
        {
            auto v = cols[j];
            bool parallel_edge = j > row_ptrs[u] && v == cols[j - 1];
            bool self_loop = u == v;
            if (!parallel_edge && !self_loop) 
            {
                ncols[uoff++] = v;
            }
        }
    }
    edge_count = nedge_count;
    std::swap(row_ptrs, nrow_ptrs);
    std::swap(cols, ncols);
    delete[] nrow_ptrs;
    nrow_ptrs = NULL;
    delete[] ncols;
    ncols = NULL;

    auto new_vertex_ids = new uintV[vertex_count];
    uintV max_vertex_id = 0;
    for (uintV u = 0; u < vertex_count; ++u) 
    {
        if (row_ptrs[u] == row_ptrs[u + 1]) {
            new_vertex_ids[u] = vertex_count;
        } 
        else 
        {
            new_vertex_ids[u] = max_vertex_id++;
            row_ptrs[new_vertex_ids[u]] = row_ptrs[u];
        }
    }
    for (uintE j = 0; j < edge_count; ++j) 
    {
        cols[j] = new_vertex_ids[cols[j]];
    }
    delete[] new_vertex_ids;
    new_vertex_ids = NULL;
    vertex_count = max_vertex_id;
    row_ptrs[vertex_count] = edge_count;

    timer.EndTimer();
    std::cout << "finish preprocess, time=" << timer.GetElapsedMicroSeconds() / 1000.0 << "ms"
              << ", now vertex_count=" << vertex_count << ",edge_count=" << edge_count << std::endl;

    SetVertexCount(vertex_count);
    SetEdgeCount(edge_count);
    SetRowPtrs(row_ptrs);
    SetCols(cols);
}