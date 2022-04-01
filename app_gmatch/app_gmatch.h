#pragma once

#define MAX_STAGE_NUM 3

#define STAGE_1 1
#define STAGE_2 2
#define STAGE_3 3

// #define ENABLE_FAILING_SET
#define MAXIMUM_QUERY_GRAPH_SIZE 64

#define TIME_THRESHOLD 0.1

#include "../system/workerOL.h"
#include "../system/task.h"

#include <atomic>
#include <bitset>
#include <fstream>
#include <unordered_map>

#include "graph.h"
#include "FilterVertices.h"
#include "GenerateQueryPlan.h"
#include "BuildTable.h"
#include "leapfrogjoin.h"

Graph data_graph;

typedef unsigned long long int ULL;

struct ContextValue
{   
    int query_vertices_num;
    int cur_depth;
    int valid_candidate_cnt;
    int valid_candidate_max_cnt;

    int *valid_candidate_idx;
    int *embedding, *idx_embedding;
    hash_set<int> visited;

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set;
    std::unordered_map<int, int> reverse_embedding;
#endif

    ContextValue()
    {

    }
    ~ContextValue()
    {
        delete[] embedding;
        delete[] idx_embedding;
        delete[] valid_candidate_idx;
    }
};

#ifdef ENABLE_FAILING_SET
ofbinstream & operator>>(ofbinstream & m, std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> & c)
{
    std::vector<char> buf;
    m >> buf;
    for (int j=0; j<MAXIMUM_QUERY_GRAPH_SIZE; j++)
        c[j] = ((buf[j>>3] >> (j & 7)) & 1);
    return m;
}
ifbinstream & operator<<(ifbinstream & m, const std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> & c)
{
    std::vector<char> buf((MAXIMUM_QUERY_GRAPH_SIZE + 7) >> 3);
    for (int j=0; j<MAXIMUM_QUERY_GRAPH_SIZE; j++)
        buf[j>>3] |= (c[j] << (j & 7));
    m << buf;
    return m;
}
#endif

ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{
    m >> c.cur_depth;
    m >> c.valid_candidate_max_cnt;
    c.valid_candidate_idx = new int[c.valid_candidate_max_cnt];
    c.valid_candidate_cnt = 0;
    m >> c.query_vertices_num;

    c.embedding = new int[c.query_vertices_num];
    for(int i = 0; i < c.query_vertices_num; i++)
		m >> c.embedding[i];
    c.idx_embedding = new int[c.query_vertices_num];
    for(int i = 0; i < c.query_vertices_num; i++)
        m >> c.idx_embedding[i];

    m >> c.visited;

#ifdef ENABLE_FAILING_SET
    // m >> c.vec_failing_set;
    size_t size;
    m >> size;
    for (size_t i = 0; i < size; i++) {
        int key;
        m >> key;
        m >> c.reverse_embedding[key];
    }
#endif
    return m;
}
ifbinstream & operator<<(ifbinstream & m, const ContextValue & c) 
{
    m << c.cur_depth;
    m << c.valid_candidate_max_cnt;
    m << c.query_vertices_num;

    for(int i = 0; i < c.query_vertices_num; i++) 
        m << c.embedding[i];
    for(int i = 0; i < c.query_vertices_num; i++) 
        m << c.idx_embedding[i];
    
    m << c.visited;
    
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

    int **candidates;
    int *candidates_count;
    int *bfs_order, *matching_order, *pivot;
    TreeNode *tree;

    Edges ***edge_matrix;

    int **bn;
    int *bn_count;

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
#endif

    int cur_stage_num;

    // =======================================

    GMQuery() 
    {
        cur_stage_num = STAGE_1;
        counters.assign(32, 0);
    }

    ~GMQuery()
    {
        for(int i=0; i<query_graph.getVerticesCount(); ++i) {
            delete[] bn[i]; 
            delete[] candidates[i];
            for (int j=0; j<query_graph.getVerticesCount(); ++j) {
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


#ifdef ENABLE_FAILING_SET
    void computeAncestor(Graph &query_graph, int **bn, int *bn_cnt, int *order,
                        std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors) 
    {
        int query_vertices_num = query_graph.getVerticesCount();
        ancestors.resize(query_vertices_num);

        // Compute the ancestor in the top-down order.
        for (int i = 0; i < query_vertices_num; ++i) {
            int u = order[i];
            ancestors[u].set(u);
            for (int j = 0; j < bn_cnt[i]; ++j) {
                int u_bn = bn[i][j];
                ancestors[u] |= ancestors[u_bn];
            }
        }
    }
#endif

    void generateBN(Graph &query_graph, int *order, int **&bn, int *&bn_count) 
    {
        int query_vertices_num = query_graph.getVerticesCount();
        bn_count = new int[query_vertices_num];
        std::fill(bn_count, bn_count + query_vertices_num, 0);
        bn = new int *[query_vertices_num];
        for (int i = 0; i < query_vertices_num; ++i) {
            bn[i] = new int[query_vertices_num];
        }

        std::vector<bool> visited_vertices(query_vertices_num, false);
        visited_vertices[order[0]] = true;
        for (int i = 1; i < query_vertices_num; ++i) {
            int vertex = order[i];

            int nbrs_cnt;
            const int *nbrs = query_graph.getVertexNeighbors(vertex, nbrs_cnt);
            for (int j = 0; j < nbrs_cnt; ++j) {
                int nbr = nbrs[j];

                if (visited_vertices[nbr]) {
                    bn[i][bn_count[i]++] = nbr;
                }
            }
            visited_vertices[vertex] = true;
        }
    }

    void generateValidCandidateIndex(int u, int *idx_embedding, int &valid_candidate_cnt, int *valid_candidate_index,
                                    Edges ***edge_matrix, int *bn, int bn_cnt)
    {
        if(bn_cnt == 1) {
            int current_bn = bn[0];
            Edges &current_edge = *edge_matrix[current_bn][u];
            int current_index_id = idx_embedding[current_bn];

            int current_candidates_count =
                current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
            int *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];
            valid_candidate_cnt = current_candidates_count;
            for(int i=0; i<current_candidates_count; i++) {
                valid_candidate_index[i] = current_candidates[i];
            }
        }
        else {
            vector<vector<int> > vecs;

            // vector<int> res;

            for(int i=0; i<bn_cnt; i++) {
                int current_bn = bn[i];

                Edges &current_edge = *edge_matrix[current_bn][u];
                int current_index_id = idx_embedding[current_bn];

                int current_candidates_count =
                    current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
                int *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

                vector<int> vec(current_candidates, current_candidates+current_candidates_count);

                vecs.push_back(vec);
                // if(i == 0)
                // {
                //     res.insert(res.begin(), vec.begin(), vec.end());
                // }
                // else
                // {
                //     vector<int> new_res;
                //     std::set_intersection(res.begin(), res.end(),
                //             vec.begin(), vec.end(),
                //             back_inserter(new_res));
                //     res.swap(new_res);
                // }
            }
            auto intersection = leapfrogJoin(vecs);

            valid_candidate_cnt = intersection.size();
            for(int i=0; i<valid_candidate_cnt; i++) {
                valid_candidate_index[i] = intersection[i];
            }

            // valid_candidate_cnt = res.size();
            // for(int i=0; i<valid_candidate_cnt; i++) {
            //     valid_candidate_index[i] = res[i];
            // }
        }
    }

    void generateValidCandidateIndex(int cur_depth, int *idx_embedding, int *idx_count, int **valid_candidate_idx,
                                    Edges ***edge_matrix, int *bn, int bn_cnt, int *order)
    {   
        int u = order[cur_depth];
        if(bn_cnt == 1) {
            int current_bn = bn[0];
            Edges &current_edge = *edge_matrix[current_bn][u];
            int current_index_id = idx_embedding[current_bn];

            int current_candidates_count =
                current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
            int *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];
            idx_count[cur_depth] = current_candidates_count;
            for(int i=0; i<current_candidates_count; i++) {
                valid_candidate_idx[cur_depth][i] = current_candidates[i];
            }
        }
        else
        {
            vector<int> res;

            for(int i=0; i<bn_cnt; i++) {
                int current_bn = bn[i];

                Edges &current_edge = *edge_matrix[current_bn][u];
                int current_index_id = idx_embedding[current_bn];

                int current_candidates_count =
                    current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
                int *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

                vector<int> vec;
                vec.insert(vec.begin(), current_candidates, current_candidates+current_candidates_count);

                // vecs.push_back(vec);
                if(i == 0)
                {
                    res.insert(res.begin(), vec.begin(), vec.end());
                }
                else
                {
                    vector<int> new_res;
                    std::sort(res.begin(), res.end());
                    std::sort(vec.begin(), vec.end());
                    std::set_intersection(res.begin(), res.end(),
                            vec.begin(), vec.end(),
                            back_inserter(new_res));
                    res.swap(new_res);
                }
            }
            // auto intersection = leapfrogJoin(vecs);

            // valid_candidate_cnt = intersection.size();
            // for(int i=0; i<valid_candidate_cnt; i++) {
            //     valid_candidate_index[i] = intersection[i];
            // }

            idx_count[cur_depth] = res.size();
            for(int i=0; i<idx_count[cur_depth]; i++) {
                valid_candidate_idx[cur_depth][i] = res[i];
            }
        }
    }


    void LFTJ(Graph &query_graph, Edges ***edge_matrix, int **candidates,
                int *candidates_count, int *order, int *embedding, int *idx_embedding,
                int **bn, int *bn_count, GMQuery &q)

    {
        int *idx = new int[query_graph.getVerticesCount()];
        int *idx_count = new int[query_graph.getVerticesCount()];

        bool *visited_vertices = new bool[data_graph.getVerticesCount()];
        std::fill(visited_vertices, visited_vertices + data_graph.getVerticesCount(), false);

        int max_candidates_num = candidates_count[0];

        for (int i = 1; i < query_graph.getVerticesCount(); ++i) {
            int cur_vertex = i;
            int cur_candidate_num = candidates_count[cur_vertex];

            if (cur_candidate_num > max_candidates_num) {
                max_candidates_num = cur_candidate_num;
            }
        }

        int **valid_candidate_idx = new int *[query_graph.getVerticesCount()];
        for (int i = 0; i < query_graph.getVerticesCount(); ++i) {
            valid_candidate_idx[i] = new int[max_candidates_num];
        }

        int cur_depth = 0;
        int max_depth = query_graph.getVerticesCount();
        int start_vertex = order[0];

        idx[cur_depth] = 0;
        idx_count[cur_depth] = candidates_count[start_vertex];

        for (int i = 0; i < idx_count[cur_depth]; ++i) {
            valid_candidate_idx[cur_depth][i] = i;
        }

        while (true) {
            while (idx[cur_depth] < idx_count[cur_depth]) {
                int valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
                int u = order[cur_depth];
                int v = candidates[u][valid_idx];

                embedding[u] = v;
                idx_embedding[u] = valid_idx;
                visited_vertices[v] = true;
                idx[cur_depth] += 1;

                if (cur_depth == max_depth - 1) {
                    counter ++;
                    visited_vertices[v] = false;

                    if(counter % 100000000 == 0) cout<<counter<<endl;
                
                } else {
                    cur_depth += 1;
                    idx[cur_depth] = 0;
                    generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn[cur_depth], bn_count[cur_depth], order);
                }
            }
            cur_depth -= 1;
            if (cur_depth < 0)
                break;
            else
                visited_vertices[embedding[order[cur_depth]]] = false;
        }
    }

    void backtrack(int cur_depth, Graph &query_graph, int **candidates, int *candidates_count,
                    Edges ***edge_matrix, int *order, int valid_candidate_cnt, 
                    int *valid_candidate_idx, int *embedding, int *idx_embedding, hash_set<int> &visited, 
                    int **bn, int *bn_count, FILE *gfpout, GMQuery &q)
    {   
        struct timeb cur_time;
		double drun_time;

        if(cur_depth == query_graph.getVerticesCount()) {
            for(int i=0; i<cur_depth; i++)
                fprintf(gfpout, "%d ", embedding[i]);
            fprintf(gfpout, "\n");

            counter++;
            
            // if(counter == 100000) {
            //     ftime(&q.end_t);
            //     double totaltime = q.end_t.time-q.start_t.time+(double)(q.end_t.millitm-q.start_t.millitm)/1000;
            //     cout<<"Query "<<get_queryID()<<" first 10^5 total time: "<<totaltime<<endl;
            // }

        } else {
            int u = order[cur_depth];
            if(cur_depth == 0) {
                valid_candidate_cnt = candidates_count[u];
                for (int i = 0; i < valid_candidate_cnt; ++i) {
                    valid_candidate_idx[i] = i;
                }
            } else {
                generateValidCandidateIndex(u, idx_embedding, valid_candidate_cnt, valid_candidate_idx,
                                            edge_matrix, bn[cur_depth], bn_count[cur_depth]);
            }

            for(int i=0; i<valid_candidate_cnt; i++) {
                int valid_idx = valid_candidate_idx[i];
                int v = candidates[u][valid_idx];
                if(visited.find(v) != visited.end()) {
                    continue;
                }
                embedding[u] = v;
                idx_embedding[u] = valid_idx;
                visited.insert(v);

                int next_valid_candidate_cnt = 0;

                int *next_valid_candidate_idx = NULL;
                if(cur_depth+1 < query_graph.getVerticesCount())
                    next_valid_candidate_idx = new int[candidates_count[order[cur_depth+1]]];

                ftime(&cur_time);
                drun_time = cur_time.time-data_graph.gtime_start[thread_id].time+(double)(cur_time.millitm-data_graph.gtime_start[thread_id].millitm)/1000;
                if(drun_time < TIME_THRESHOLD) 
                {
                    backtrack(cur_depth+1, query_graph, candidates, candidates_count, 
                            edge_matrix, order, next_valid_candidate_cnt, next_valid_candidate_idx, 
                            embedding, idx_embedding, visited, bn, bn_count, gfpout, q);
                }
                else {
                    int query_vertices_num = query_graph.getVerticesCount();
                    GMTask *t = new GMTask();
                    t->context.cur_depth = cur_depth+1;
                    t->context.query_vertices_num = query_vertices_num;

                    t->context.embedding = new int[query_vertices_num];
                    memcpy(t->context.embedding, embedding, sizeof(int)*query_vertices_num);
                    t->context.idx_embedding = new int[query_vertices_num];
                    memcpy(t->context.idx_embedding, idx_embedding, sizeof(int)*query_vertices_num);
                    t->context.visited = visited;

                    t->context.valid_candidate_cnt = 0;
                    t->context.valid_candidate_max_cnt = 0;
                    t->context.valid_candidate_idx = NULL;
                    if(cur_depth+1 < query_vertices_num) {
                        t->context.valid_candidate_idx = new int[candidates_count[order[cur_depth+1]]];
                        t->context.valid_candidate_max_cnt = candidates_count[order[cur_depth+1]];
                    }
                    add_task(t);
                }

                visited.erase(v);
                if(next_valid_candidate_idx != NULL)
                    delete []next_valid_candidate_idx;
            }
        }
    }

#ifdef ENABLE_FAILING_SET
    void backtrack_fs(int cur_depth, Graph &query_graph, int **candidates, int *candidates_count,
                    Edges ***edge_matrix, int *order, int valid_candidate_cnt, 
                    int *valid_candidate_idx, int *embedding, int *idx_embedding, hash_set<int> &visited, 
                    int **bn, int *bn_count, std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &vec_failing_set, 
                    std::unordered_map<int, int> &reverse_embedding, 
                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors, FILE *gfpout, GMQuery &q,
                    bool &enter_split_phase)
    {   
        struct timeb cur_time;
		double drun_time;

        if(cur_depth == query_graph.getVerticesCount()) {
            // for(int i=0; i<cur_depth; i++)
            //     fprintf(gfpout, "%d ", embedding[i]);
            // fprintf(gfpout, "\n");
            counter++;
            if(counter == 100000) {
                ftime(&q.end_t);
                double totaltime = q.end_t.time-q.start_t.time+(double)(q.end_t.millitm-q.start_t.millitm)/1000;
                cout<<"Query "<<get_queryID()<<" first 10^5 total time: "<<totaltime<<endl;
            }

            vec_failing_set[cur_depth - 1].set();
            vec_failing_set[cur_depth - 2] |= vec_failing_set[cur_depth - 1];
        } else {
            int u = order[cur_depth];
            if(cur_depth == 0) {
                valid_candidate_cnt = candidates_count[u];
                for (int i = 0; i < valid_candidate_cnt; ++i) {
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

            for(int i=0; i<valid_candidate_cnt; i++) {
                int valid_idx = valid_candidate_idx[i];
                int v = candidates[u][valid_idx];
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

                int next_valid_candidate_cnt = 0;
                int *next_valid_candidate_idx = NULL;
                if(cur_depth+1 < query_graph.getVerticesCount())
                    next_valid_candidate_idx = new int[candidates_count[order[cur_depth+1]]];

                ftime(&cur_time);
                drun_time = cur_time.time-data_graph.gtime_start[thread_id].time+(double)(cur_time.millitm-data_graph.gtime_start[thread_id].millitm)/1000;
                if(drun_time < TIME_THRESHOLD) {
                    backtrack_fs(cur_depth+1, query_graph, candidates, candidates_count, 
                        edge_matrix, order, next_valid_candidate_cnt, next_valid_candidate_idx, 
                        embedding, idx_embedding, visited, bn, bn_count, vec_failing_set, 
                        reverse_embedding, ancestors, gfpout, q, enter_split_phase);
                } else {
                    enter_split_phase = true;
                    int query_vertices_num = query_graph.getVerticesCount();
                    GMTask *t = new GMTask();
                    t->context.cur_depth = cur_depth+1;
                    t->context.query_vertices_num = query_vertices_num;

                    t->context.embedding = new int[query_vertices_num];
                    memcpy(t->context.embedding, embedding, sizeof(int)*query_vertices_num);
                    t->context.idx_embedding = new int[query_vertices_num];
                    memcpy(t->context.idx_embedding, idx_embedding, sizeof(int)*query_vertices_num);
                    t->context.visited = visited;

                    t->context.valid_candidate_cnt = 0;
                    t->context.valid_candidate_max_cnt = 0;
                    t->context.valid_candidate_idx = NULL;
                    if(cur_depth+1 < query_vertices_num) {
                        t->context.valid_candidate_idx = new int[candidates_count[order[cur_depth+1]]];
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
            t->context.valid_candidate_cnt = 0;
            t->context.valid_candidate_max_cnt = q.candidates_count[q.matching_order[0]];
            t->context.query_vertices_num = q.query_graph.getVerticesCount();
            t->context.valid_candidate_idx = new int[t->context.valid_candidate_max_cnt];
            t->context.embedding = new int[q.query_graph.getVerticesCount()];
            t->context.idx_embedding = new int[q.query_graph.getVerticesCount()];
            
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
        if(q.cur_stage_num == STAGE_1) {
            FilterVertices::DPisoFilter(data_graph, q.query_graph, q.candidates, q.candidates_count, 
                                        q.bfs_order, q.tree);

        } else if(q.cur_stage_num == STAGE_2) {
            GenerateQueryPlan::generateGQLQueryPlan(data_graph, q.query_graph, q.candidates_count, 
                                                    q.matching_order, q.pivot);

            std::cout<<"======= print matching order =========="<<std::endl;
            for(int i=0; i<q.query_graph.getVerticesCount(); i++) {
                std::cout<<q.matching_order[i]<<" ";
            }
            std::cout<<std::endl;

            generateBN(q.query_graph, q.matching_order, q.bn, q.bn_count);

            q.edge_matrix = new Edges **[q.query_graph.getVerticesCount()];
            for (int i = 0; i < q.query_graph.getVerticesCount(); ++i) {
                q.edge_matrix[i] = new Edges *[q.query_graph.getVerticesCount()];
            }
            
            BuildTable::buildTable(data_graph, q.query_graph, q.candidates, q.candidates_count, q.edge_matrix);

#ifdef ENABLE_FAILING_SET
            computeAncestor(q.query_graph, q.bn, q.bn_count, q.matching_order, q.ancestors);
#endif
            
        } else {
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
            backtrack(context.cur_depth, q.query_graph, q.candidates, q.candidates_count, q.edge_matrix,
            q.matching_order, context.valid_candidate_cnt, context.valid_candidate_idx, 
            context.embedding, context.idx_embedding, context.visited, q.bn, q.bn_count, gfpout, q);

            // LFTJ(q.query_graph, q.edge_matrix, q.candidates, q.candidates_count, q.bfs_order, context.embedding, 
            //         context.idx_embedding, q.bn, q.bn_count, q);
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
            for(int i=0; i<32; i++)
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
    GMWorker(int num_compers) : Worker(num_compers)
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






