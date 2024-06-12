#pragma once

#define MAX_STAGE_NUM 3

#define STAGE_1 1
#define STAGE_2 2
#define STAGE_3 3

// #define ENABLE_FAILING_SET
#define MAXIMUM_QUERY_GRAPH_SIZE 64

#define TIME_THRESHOLD 10000000

#include "../system/workerOL.h"
#include "../system/task.h"

#include "intersection/computesetintersection.h"
// #include "intersection/avx2.hpp"
// #include "intersection/galloping.hpp"

#include <atomic>
#include <bitset>
#include <fstream>
#include <unordered_map>
#include <assert.h>
#include <algorithm>

#include "graph.h"
#include "FilterVertices.h"
#include "GenerateQueryPlan.h"
#include "BuildTable.h"
#include "leapfrogjoin.h"

Graph data_graph;

typedef unsigned long long int ULL;

struct ContextValue
{   
    ui query_vertices_num;
    ui cur_depth;

    // ui valid_candidate_cnt;
    // ui valid_candidate_max_cnt;
    // ui *valid_candidate_idx;

    ui *embedding, *idx_embedding;
    // hash_set<ui> visited;

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set;
    std::unordered_map<ui, ui> reverse_embedding;
#endif

    ContextValue()
    {
        embedding = NULL;
        idx_embedding = NULL;
    }
    ~ContextValue()
    {
        if (embedding != NULL)
            delete[] embedding;
        if (idx_embedding != NULL)
            delete[] idx_embedding;
        // delete[] valid_candidate_idx;
    }
};

#ifdef ENABLE_FAILING_SET
ofbinstream & operator>>(ofbinstream & m, std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> & c)
{
    std::vector<char> buf;
    m >> buf;
    for (ui j=0; j<MAXIMUM_QUERY_GRAPH_SIZE; j++)
        c[j] = ((buf[j>>3] >> (j & 7)) & 1);
    return m;
}
ifbinstream & operator<<(ifbinstream & m, const std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> & c)
{
    std::vector<char> buf((MAXIMUM_QUERY_GRAPH_SIZE + 7) >> 3);
    for (ui j=0; j<MAXIMUM_QUERY_GRAPH_SIZE; j++)
        buf[j>>3] |= (c[j] << (j & 7));
    m << buf;
    return m;
}
#endif

ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{
    m >> c.cur_depth;
    // m >> c.valid_candidate_max_cnt;
    // c.valid_candidate_idx = new ui[c.valid_candidate_max_cnt];
    // c.valid_candidate_cnt = 0;

    m >> c.query_vertices_num;
    c.embedding = new ui[c.query_vertices_num];
    for(ui i = 0; i < c.query_vertices_num; i++)
		m >> c.embedding[i];
    c.idx_embedding = new ui[c.query_vertices_num];
    for(ui i = 0; i < c.query_vertices_num; i++)
        m >> c.idx_embedding[i];

    // m >> c.visited;

#ifdef ENABLE_FAILING_SET
    // m >> c.vec_failing_set;
    size_t size;
    m >> size;
    for (size_t i = 0; i < size; i++) {
        ui key;
        m >> key;
        m >> c.reverse_embedding[key];
    }
#endif
    return m;
}
ifbinstream & operator<<(ifbinstream & m, const ContextValue & c) 
{

    m << c.cur_depth;
    // m << c.valid_candidate_max_cnt;
    m << c.query_vertices_num;

    for(ui i = 0; i < c.query_vertices_num; i++) 
        m << c.embedding[i];
    for(ui i = 0; i < c.query_vertices_num; i++) 
        m << c.idx_embedding[i];
    
    // m << c.visited;
    
#ifdef ENABLE_FAILING_SET
    // m << c.vec_failing_set;
    m << c.reverse_embedding.size();
    for (auto it = c.reverse_embedding.begin(); it != c.reverse_embedding.end(); ++it) {
        m << it->first;
        m << it->second;
    }
#endif
    return m;
}

struct GMQuery
{   
    // given by user
    Graph query_graph;

    //==== resources held for each query ===

    vector<ULL> counters;

    struct timeb start_t, end_t;

    ui **candidates;
    ui *candidates_count;
    ui *bfs_order, *matching_order, *pivot;
    TreeNode *tree;

    Edges ***edge_matrix;

    ui **bn;
    ui *bn_count;

    ui max_candidate_cnt;

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
#endif

    ui cur_stage_num;

    // =======================================

    GMQuery() 
    {
        cur_stage_num = STAGE_1;
        counters.assign(32, 0);

        max_candidate_cnt = 0;

    }

    ~GMQuery()
    {
        for(ui i=0; i<query_graph.getVerticesCount(); ++i) {
            delete[] bn[i]; 
            delete[] candidates[i];
            for (ui j=0; j<query_graph.getVerticesCount(); ++j) {
                delete edge_matrix[i][j];
            }
            delete[] edge_matrix[i];
        }
        delete[] bn;
        delete[] bn_count;
        delete[] candidates_count;
        delete[] candidates;
        delete[] edge_matrix;
        delete[] bfs_order;
        delete[] matching_order;
        delete[] pivot;
        delete[] tree;
    }
};

typedef Task<ContextValue> GMTask;

class GMComper: public Comper<GMTask, GMQuery>
{
public:

    ULL counter;

    ui *temp_buffer;
    bool *visited_arr;
    ui *idx;
    ui *idx_count;
    ui **valid_candidate_idx;


    GMComper(): temp_buffer(NULL), visited_arr(NULL), idx(NULL), idx_count(NULL), valid_candidate_idx(NULL) {}


#ifdef ENABLE_FAILING_SET
    void computeAncestor(Graph &query_graph, ui **bn, ui *bn_cnt, ui *order,
                        std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors) 
    {
        ui query_vertices_num = query_graph.getVerticesCount();
        ancestors.resize(query_vertices_num);

        // Compute the ancestor in the top-down order.
        for (ui i = 0; i < query_vertices_num; ++i) {
            ui u = order[i];
            ancestors[u].set(u);
            for (ui j = 0; j < bn_cnt[i]; ++j) {
                ui u_bn = bn[i][j];
                ancestors[u] |= ancestors[u_bn];
            }
        }
    }
#endif

    void generateBN(Graph &query_graph, ui *order, ui **&bn, ui *&bn_count) 
    {
        ui query_vertices_num = query_graph.getVerticesCount();
        bn_count = new ui[query_vertices_num];
        std::fill(bn_count, bn_count + query_vertices_num, 0);
        bn = new ui *[query_vertices_num];
        for (ui i = 0; i < query_vertices_num; ++i) {
            bn[i] = new ui[query_vertices_num];
        }

        std::vector<bool> visited_vertices(query_vertices_num, false);
        visited_vertices[order[0]] = true;
        for (ui i = 1; i < query_vertices_num; ++i) {
            ui vertex = order[i];

            ui nbrs_cnt;
            const ui *nbrs = query_graph.getVertexNeighbors(vertex, nbrs_cnt);
            for (ui j = 0; j < nbrs_cnt; ++j) {
                ui nbr = nbrs[j];

                if (visited_vertices[nbr]) {
                    bn[i][bn_count[i]++] = nbr;
                }
            }
            visited_vertices[vertex] = true;
        }

        cout << "======= BN ========" << endl;
        for (int i = 1; i < query_vertices_num; ++i)
        {
            for (int j = 0; j < bn_count[i]; ++j)
            {
                cout << bn[i][j] << " ";
            }
            cout << endl;
        }
        cout << "==================" << endl;
    }


    void generateValidCandidateIndex(ui depth, ui *idx_embedding, ui &valid_candidate_cnt, ui *valid_candidate_index,
                                    Edges ***edge_matrix, ui **bn, ui *bn_cnt, ui *order, ui *temp_buffer_)
    {
        
        ui u = order[depth];
        ui previous_bn = bn[depth][0];
        ui previous_index_id = idx_embedding[previous_bn];

        ui valid_candidate_count = 0;

        Edges& previous_edge = *edge_matrix[previous_bn][u];

        valid_candidate_count = previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
        ui* previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

        memcpy(valid_candidate_index, previous_candidates, valid_candidate_count * sizeof(ui));
        
        ui temp_count;
        for (ui i = 1; i < bn_cnt[depth]; ++i) {
            
            VertexID current_bn = bn[depth][i];

            Edges& current_edge = *edge_matrix[current_bn][u];
            ui current_index_id = idx_embedding[current_bn];


            ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];

            ui* current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];


            if (current_candidates_count < valid_candidate_cnt)
                ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index, valid_candidate_count,
                            temp_buffer_, temp_count);
            else
                ComputeSetIntersection::ComputeCandidates(valid_candidate_index, valid_candidate_count, current_candidates, current_candidates_count,
                            temp_buffer_, temp_count);

            for(int i = 0; i < temp_count; ++i)
            {
                valid_candidate_index[i] = temp_buffer_[i];
            }
            valid_candidate_count = temp_count;
        }

        valid_candidate_cnt = valid_candidate_count;
    }


    void generateValidCandidateIndex(ui depth, ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                    Edges ***edge_matrix, ui **bn, ui *bn_cnt, ui *order, ui *temp_buffer_)
    {   

        ui u = order[depth];
        ui previous_bn = bn[depth][0];
        ui previous_index_id = idx_embedding[previous_bn];
        ui valid_candidates_count = 0;


        Edges& previous_edge = *edge_matrix[previous_bn][u];

        valid_candidates_count = previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
        ui* previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

        memcpy(valid_candidate_index[depth], previous_candidates, valid_candidates_count * sizeof(ui));

        
        ui temp_count;
        for (ui i = 1; i < bn_cnt[depth]; ++i) {
            
            VertexID current_bn = bn[depth][i];

            Edges& current_edge = *edge_matrix[current_bn][u];
            ui current_index_id = idx_embedding[current_bn];


            ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];

            ui* current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];


            if (current_candidates_count < valid_candidates_count)
                ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index[depth], valid_candidates_count,
                            temp_buffer_, temp_count);
            else
                ComputeSetIntersection::ComputeCandidates(valid_candidate_index[depth], valid_candidates_count, current_candidates, current_candidates_count,
                            temp_buffer_, temp_count);
          

            // std::swap(temp_buffer, valid_candidate_index[depth]); // all elements are swapped

            for(int i = 0; i < temp_count; ++i)
            {
                valid_candidate_index[depth][i] = temp_buffer_[i];
            }
            valid_candidates_count = temp_count;
        }

        idx_count[depth] = valid_candidates_count;
    }

    double countElaspedTime()
    {
        struct timeb cur_time;
        ftime(&cur_time);
        return cur_time.time-data_graph.gtime_start[thread_id].time+(double)(cur_time.millitm-data_graph.gtime_start[thread_id].millitm)/1000;
    }


    void LFTJ(int enter_depth, Graph &query_graph, Edges ***edge_matrix, ui **candidates,
                ui *candidates_count, ui *order, ui *embedding, ui *idx_embedding,
                ui **bn, ui *bn_count, GMQuery &q)

    {
        // ui *idx = new ui[query_graph.getVerticesCount()];
        // ui *idx_count = new ui[query_graph.getVerticesCount()];
        // ui *embedding = new ui[query_graph.getVerticesCount()];
        // ui *idx_embedding = new ui[query_graph.getVerticesCount()];

        // ui max_candidates_num = candidates_count[0];

        // for (ui i = 1; i < query_graph.getVerticesCount(); ++i) {
        //     ui cur_vertex = i;
        //     ui cur_candidate_num = candidates_count[cur_vertex];

        //     if (cur_candidate_num > max_candidates_num) {
        //         max_candidates_num = cur_candidate_num;
        //     }
        // }

        // ui *temp_buffer = new ui[max_candidates_num];

        // ui *valid_candidate_idx[query_graph.getVerticesCount()];
        // for (ui i = 0; i < query_graph.getVerticesCount(); ++i) {
        //     valid_candidate_idx[i] = new ui[max_candidates_num];
        // }

        //=============================================

        int cur_depth = enter_depth;
        int max_depth = query_graph.getVerticesCount();

        if (cur_depth == 0)
        {

            ui start_vertex = order[0];

            idx[cur_depth] = 0;

            idx_count[cur_depth] = candidates_count[start_vertex];

            for (ui i = 0; i < idx_count[cur_depth]; ++i) {
                valid_candidate_idx[cur_depth][i] = i;
            }
        }
        else
        {  
            idx[cur_depth] = 0;
        
            // compute set intersection
            generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn, bn_count, order, temp_buffer);
  
            
            // initialize visited_arr array 
            for (ui i = 0; i < enter_depth; ++i)
            {
                visited_arr[embedding[order[i]]] = true;
            }
        }

        while (true) {
            while (idx[cur_depth] < idx_count[cur_depth]) {
                ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];

                ui u = order[cur_depth];
                
                ui v = candidates[u][valid_idx];


                if (visited_arr[v]) {
                    idx[cur_depth] += 1;
                    continue;
                }

                embedding[u] = v;
                idx_embedding[u] = valid_idx;

                visited_arr[v] = true;

                idx[cur_depth] += 1;

                if (cur_depth == max_depth - 1) {
                        
                    counter += 1;

                    // print first 10000 results
                    // if (counter < 10000)
                    // {
                    //     for(ui i = 0; i < max_depth; ++i)
                    //     {
                    //         cout << embedding[i] << " ";
                    //     }
                    //     cout << endl;
                    // }

                    visited_arr[v] = false;
                    
                    if(counter % 1000000000 == 0) cout<<counter<<endl;
                    
                    continue;
                }

                // if not timeout, continue search 
                if(countElaspedTime() < TIME_THRESHOLD) 
                {
                    cur_depth += 1;
                    idx[cur_depth] = 0;
                    generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn, bn_count, order, temp_buffer);
                }
                else  // if timeout, start task splitting
                {
                    ui query_vertices_num = query_graph.getVerticesCount();
                    GMTask *t = new GMTask();

                    t->context.query_vertices_num = query_vertices_num;
                    t->context.cur_depth = cur_depth+1;

                    t->context.embedding = new ui[query_vertices_num];
                    memcpy(t->context.embedding, embedding, sizeof(ui)*query_vertices_num);
                    t->context.idx_embedding = new ui[query_vertices_num];
                    memcpy(t->context.idx_embedding, idx_embedding, sizeof(ui)*query_vertices_num);
                      
                    add_task(t);

                    visited_arr[v] = false;
                }

            }
            cur_depth -= 1;
            if (cur_depth < enter_depth)
                break;
            else
            {
                visited_arr[embedding[order[cur_depth]]] = false;
            }
        }

        // ####
        for (ui i = 0; i < enter_depth; ++i)
        {
            visited_arr[embedding[order[i]]] = false;
        }
    }

/**
    void backtrack(ui cur_depth, Graph &query_graph, ui **candidates, ui *candidates_count,
                    Edges ***edge_matrix, ui *order, ui valid_candidate_cnt, 
                    ui *valid_candidate_idx, ui *embedding, ui *idx_embedding, hash_set<ui> &visited, 
                    ui **bn, ui *bn_count, FILE *gfpout, GMQuery &q)
    {   
        struct timeb cur_time;
		double drun_time;

        if(cur_depth == query_graph.getVerticesCount()) {
            // for(ui i=0; i<cur_depth; i++)
            //     fpruif(gfpout, "%d ", embedding[i]);
            // fpruif(gfpout, "\n");

            counter++;
            
            if(counter % 1000000000 == 0) {
                ftime(&q.end_t);
                double totaltime = q.end_t.time-q.start_t.time+(double)(q.end_t.millitm-q.start_t.millitm)/1000;
                cout<<"temp time: "<<totaltime<<endl;
            }

        } else {
            ui u = order[cur_depth];
            if(cur_depth == 0) {
                valid_candidate_cnt = candidates_count[u];
                for (ui i = 0; i < valid_candidate_cnt; ++i) {
                    valid_candidate_idx[i] = i;
                }
            } else {
                generateValidCandidateIndex(cur_depth, idx_embedding, valid_candidate_cnt, valid_candidate_idx, edge_matrix, 
                                            bn, bn_count, order, temp_buffer);
            }


            for(ui i=0; i<valid_candidate_cnt; i++) {

                ui valid_idx = valid_candidate_idx[i];
                ui v = candidates[u][valid_idx];
                // if(visited.find(v) != visited.end()) {
                //     continue;
                // }

                if (visited_arr[v]) continue;

               
                embedding[u] = v;

                idx_embedding[u] = valid_idx;

                // visited.insert(v);
                visited_arr[v] = true;

                ui next_valid_candidate_cnt = 0;

                ui *next_valid_candidate_idx = NULL;
                
                if(cur_depth+1 < query_graph.getVerticesCount())
                {
                    next_valid_candidate_idx = new ui[candidates_count[order[cur_depth+1]]];
                }

                ftime(&cur_time);
                drun_time = cur_time.time-data_graph.gtime_start[thread_id].time+(double)(cur_time.millitm-data_graph.gtime_start[thread_id].millitm)/1000;
                if(drun_time < TIME_THRESHOLD) 
                {
                    backtrack(cur_depth+1, query_graph, candidates, candidates_count, 
                            edge_matrix, order, next_valid_candidate_cnt, next_valid_candidate_idx, 
                            embedding, idx_embedding, visited, bn, bn_count, gfpout, q);
                }
                else {
                    ui query_vertices_num = query_graph.getVerticesCount();
                    GMTask *t = new GMTask();
                    t->context.cur_depth = cur_depth+1;
                    t->context.query_vertices_num = query_vertices_num;

                    t->context.embedding = new ui[query_vertices_num];
                    memcpy(t->context.embedding, embedding, sizeof(ui)*query_vertices_num);
                    t->context.idx_embedding = new ui[query_vertices_num];
                    memcpy(t->context.idx_embedding, idx_embedding, sizeof(ui)*query_vertices_num);
                    t->context.visited = visited;

                    t->context.valid_candidate_cnt = 0;
                    t->context.valid_candidate_max_cnt = 0;
                    t->context.valid_candidate_idx = NULL;
                    if(cur_depth+1 < query_vertices_num) {
                        t->context.valid_candidate_idx = new ui[candidates_count[order[cur_depth+1]]];
                        t->context.valid_candidate_max_cnt = candidates_count[order[cur_depth+1]];
                    }
                    add_task(t);
                }

                // visited.erase(v);
                visited_arr[v] = false;
                if(next_valid_candidate_idx != NULL)
                    delete []next_valid_candidate_idx;
            }
        }
    }

**/

#ifdef ENABLE_FAILING_SET
    void backtrack_fs(ui cur_depth, Graph &query_graph, ui **candidates, ui *candidates_count,
                    Edges ***edge_matrix, ui *order, ui valid_candidate_cnt, 
                    ui *valid_candidate_idx, ui *embedding, ui *idx_embedding, hash_set<ui> &visited, 
                    ui **bn, ui *bn_count, std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &vec_failing_set, 
                    std::unordered_map<ui, ui> &reverse_embedding, 
                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors, FILE *gfpout, GMQuery &q,
                    bool &enter_split_phase)
    {   
        struct timeb cur_time;
		double drun_time;

        if(cur_depth == query_graph.getVerticesCount()) {
            // for(ui i=0; i<cur_depth; i++)
            //     fpruif(gfpout, "%d ", embedding[i]);
            // fpruif(gfpout, "\n");
            counter++;
            if(counter == 100000) {
                ftime(&q.end_t);
                double totaltime = q.end_t.time-q.start_t.time+(double)(q.end_t.millitm-q.start_t.millitm)/1000;
                cout<<"Query "<<get_queryID()<<" first 10^5 total time: "<<totaltime<<endl;
            }

            vec_failing_set[cur_depth - 1].set();
            vec_failing_set[cur_depth - 2] |= vec_failing_set[cur_depth - 1];
        } else {
            ui u = order[cur_depth];
            if(cur_depth == 0) {
                valid_candidate_cnt = candidates_count[u];
                for (ui i = 0; i < valid_candidate_cnt; ++i) {
                    valid_candidate_idx[i] = i;
                }
            } else {
                generateValidCandidateIndex(u, idx_embedding, valid_candidate_cnt, valid_candidate_idx,
                                            edge_matrix, bn[cur_depth], bn_count[cur_depth]);
            }

            if(cur_depth != 0) {
                if(valid_candidate_cnt == 0) {
                    vec_failing_set[cur_depth - 1] = ancestors[u];
                } else {
                    vec_failing_set[cur_depth - 1].reset(); // to 0
                }
            }

            for(ui i=0; i<valid_candidate_cnt; i++) {
                ui valid_idx = valid_candidate_idx[i];
                ui v = candidates[u][valid_idx];
                if(visited.find(v) != visited.end()) {
                    vec_failing_set[cur_depth] = ancestors[u];
                    vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                    continue;
                }
                embedding[u] = v;
                reverse_embedding[v] = u;
                idx_embedding[u] = valid_idx;
                visited.insert(v);

                ui next_valid_candidate_cnt = 0;
                ui *next_valid_candidate_idx = NULL;
                if(cur_depth+1 < query_graph.getVerticesCount())
                    next_valid_candidate_idx = new ui[candidates_count[order[cur_depth+1]]];

                ftime(&cur_time);
                drun_time = cur_time.time-data_graph.gtime_start[thread_id].time+(double)(cur_time.millitm-data_graph.gtime_start[thread_id].millitm)/1000;
                if(drun_time < TIME_THRESHOLD) {
                    backtrack_fs(cur_depth+1, query_graph, candidates, candidates_count, 
                        edge_matrix, order, next_valid_candidate_cnt, next_valid_candidate_idx, 
                        embedding, idx_embedding, visited, bn, bn_count, vec_failing_set, 
                        reverse_embedding, ancestors, gfpout, q, enter_split_phase);
                } else {
                    enter_split_phase = true;
                    ui query_vertices_num = query_graph.getVerticesCount();
                    GMTask *t = new GMTask();
                    t->context.cur_depth = cur_depth+1;
                    t->context.query_vertices_num = query_vertices_num;

                    t->context.embedding = new ui[query_vertices_num];
                    memcpy(t->context.embedding, embedding, sizeof(ui)*query_vertices_num);
                    t->context.idx_embedding = new ui[query_vertices_num];
                    memcpy(t->context.idx_embedding, idx_embedding, sizeof(ui)*query_vertices_num);
                    t->context.visited = visited;

                    t->context.valid_candidate_cnt = 0;
                    t->context.valid_candidate_max_cnt = 0;
                    t->context.valid_candidate_idx = NULL;
                    if(cur_depth+1 < query_vertices_num) {
                        t->context.valid_candidate_idx = new ui[candidates_count[order[cur_depth+1]]];
                        t->context.valid_candidate_max_cnt = candidates_count[order[cur_depth+1]];
                    }
                    t->context.reverse_embedding = reverse_embedding;
                    // t->context.vec_failing_set = vec_failing_set;
                    add_task(t);
                }

                visited.erase(v);
                reverse_embedding.erase(v);

                if(!enter_split_phase) {
                    if(cur_depth != 0) {
                        if (!vec_failing_set[cur_depth].test(u)) { //if 0
                            // case 2.1
                            vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                            i = valid_candidate_cnt; // Pruned by Failing Set
                        } else {
                            // case 2.2
                            vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                        }
                    }
                }
                if(next_valid_candidate_idx != NULL)
                    delete []next_valid_candidate_idx;
            }
        }
    }
#endif


    virtual bool toQuery(string& line, GMQuery& q)
    {   
        ftime(&q.start_t);
        q.query_graph.loadGraphFromFile(line);
        return true;
    }

    virtual bool task_spawn(GMQuery &q)
    {
        if(q.cur_stage_num == STAGE_1) {
            GMTask *t = new GMTask();
            add_task(t);
        } else if(q.cur_stage_num == STAGE_2) {
            GMTask *t = new GMTask();
            add_task(t);
        } else if(q.cur_stage_num == MAX_STAGE_NUM) {
            GMTask *t = new GMTask();
            t->context.cur_depth = 0;
            // t->context.valid_candidate_cnt = 0;
            // t->context.valid_candidate_max_cnt = q.candidates_count[q.matching_order[0]];
            t->context.query_vertices_num = q.query_graph.getVerticesCount();
            // t->context.valid_candidate_idx = new ui[t->context.valid_candidate_max_cnt];
            t->context.embedding = new ui[q.query_graph.getVerticesCount()];
            t->context.idx_embedding = new ui[q.query_graph.getVerticesCount()];
            
#ifdef ENABLE_FAILING_SET
            t->context.vec_failing_set.resize(q.query_graph.getVerticesCount());
            t->context.reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif
            add_task(t);
        }
        return true;
    }

    virtual void compute(ContextT &context, GMQuery &q)
    {   
        if(q.cur_stage_num == STAGE_1) 
        {

            FilterVertices::DPisoFilter(data_graph, q.query_graph, q.candidates, q.candidates_count, 
                                        q.bfs_order, q.tree);      
            FilterVertices::sortCandidates(q.candidates, q.candidates_count, q.query_graph.getVerticesCount());

            for (ui i = 0; i < q.query_graph.getVerticesCount(); ++i)
            {
                q.max_candidate_cnt = std::max(q.max_candidate_cnt, q.candidates_count[i]);
            }

            std::cout << " MAX CANDS : " << q.max_candidate_cnt << std::endl;

        } 
        else if(q.cur_stage_num == STAGE_2) 
        {

            GenerateQueryPlan::generateGQLQueryPlan(data_graph, q.query_graph, q.candidates_count, 
                                                    q.matching_order, q.pivot);

            std::cout<<"======= print matching order =========="<<std::endl;
            for(ui i=0; i<q.query_graph.getVerticesCount(); i++) {
                std::cout<<q.matching_order[i]<<" ";
            }
            std::cout<<std::endl;

            generateBN(q.query_graph, q.matching_order, q.bn, q.bn_count);

            q.edge_matrix = new Edges **[q.query_graph.getVerticesCount()];
            for (ui i = 0; i < q.query_graph.getVerticesCount(); ++i) {
                q.edge_matrix[i] = new Edges *[q.query_graph.getVerticesCount()];
            }
            
            BuildTable::buildTable(data_graph, q.query_graph, q.candidates, q.candidates_count, q.edge_matrix);

#ifdef ENABLE_FAILING_SET
            computeAncestor(q.query_graph, q.bn, q.bn_count, q.matching_order, q.ancestors);
#endif
            
        } 
        else 
        {
            if (temp_buffer == NULL)
            {
                // allocate space for temp_buffer array
                temp_buffer = new ui[q.max_candidate_cnt];

                // allocate space for visited_arr array
                visited_arr = new bool[data_graph.getVerticesCount()];
                memset(visited_arr, false, sizeof(bool)*data_graph.getVerticesCount());

                // allocate space for idx and idx_count array
                idx = new ui[q.query_graph.getVerticesCount()];
                idx_count = new ui[q.query_graph.getVerticesCount()];

                // allocate space for valid candidate 2-dimensional array
                valid_candidate_idx = new ui*[q.query_graph.getVerticesCount()];
                for (ui i = 0; i < q.query_graph.getVerticesCount(); ++i) {
                    valid_candidate_idx[i] = new ui[q.max_candidate_cnt];
                }
            }            

#ifdef ENABLE_FAILING_SET
            ftime(&data_graph.gtime_start[thread_id]);
            counter = q.counters[thread_id];
            bool enter_split_phase = false;
            context.vec_failing_set.resize(q.query_graph.getVerticesCount());
            backtrack_fs(context.cur_depth, q.query_graph, q.candidates, q.candidates_count, q.edge_matrix,
            q.matching_order, context.valid_candidate_cnt, context.valid_candidate_idx, 
            context.embedding, context.idx_embedding, context.visited, q.bn, q.bn_count, context.vec_failing_set,
            context.reverse_embedding, q.ancestors, gfpout, q, enter_split_phase);
            q.counters[thread_id] = counter;
#else
            ftime(&data_graph.gtime_start[thread_id]);
            counter = q.counters[thread_id];


            // backtrack(context.cur_depth, q.query_graph, q.candidates, q.candidates_count, q.edge_matrix,
            // q.matching_order, context.valid_candidate_cnt, context.valid_candidate_idx, 
            // context.embedding, context.idx_embedding, context.visited, q.bn, q.bn_count, gfpout, q);

            // std::cout<<"backtrack"<<std::endl;q

            LFTJ(context.cur_depth, q.query_graph, q.edge_matrix, q.candidates, q.candidates_count, q.matching_order, context.embedding, 
                    context.idx_embedding, q.bn, q.bn_count, q);

            // std::cout<<"LFTJ"<<std::endl;

            q.counters[thread_id] = counter;
#endif
        }
    }

    virtual bool postprocess(GMQuery &q)
    {
        cout<<"[INFO] Query "<<get_queryID()<<", Stage "<<q.cur_stage_num<<"/"<<MAX_STAGE_NUM<<" is done."<<endl;
        if(q.cur_stage_num < MAX_STAGE_NUM) {
            q.cur_stage_num++;
            task_spawn(q);
            return true;
        } else {
            ftime(&q.end_t);
            double totaltime = q.end_t.time-q.start_t.time+(double)(q.end_t.millitm-q.start_t.millitm)/1000;
            cout<<"Query "<<get_queryID()<<" total time: "<<totaltime<<endl;
            ULL total_results = 0;
            for(ui i=0; i<32; i++)
            {
                total_results += q.counters[i];
            }
            cout<<"Query "<<get_queryID()<<" total count: "<< total_results <<endl;

            return false;
        }
    }

    virtual bool is_bigTask(ContextValue &context)
    {
        return false;
    }
};


class GMWorker : public Worker<GMComper>
{
public:
    GMWorker(ui num_compers) : Worker(num_compers)
    {
    }

    ~GMWorker()
    {
    }

    void load_data(char* file_path)
    {
        std::string fp = std::string(file_path);
        data_graph.loadGraphFromFile(fp);
    }
};






