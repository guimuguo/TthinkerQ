#pragma once

#include "graph.h"
#include <queue>

#define INVALID_VERTEX_ID 100000000

class TreeNode {
public:
    int id_;
    int parent_;
    int level_;
    int under_level_count_;
    int children_count_;
    int bn_count_;
    int fn_count_;
    int* under_level_;
    int* children_;
    int* bn_;
    int* fn_;
    size_t estimated_embeddings_num_;

public:
    TreeNode() {
        id_ = 0;
        under_level_ = NULL;
        bn_ = NULL;
        fn_ = NULL;
        children_ = NULL;
        parent_ = 0;
        level_ = 0;
        under_level_count_ = 0;
        children_count_ = 0;
        bn_count_ = 0;
        fn_count_ = 0;
        estimated_embeddings_num_ = 0;
    }

    ~TreeNode() {
        delete[] under_level_;
        delete[] bn_;
        delete[] fn_;
        delete[] children_;
    }

    void initialize(const int size) {
        under_level_ = new int[size];
        bn_ = new int[size];
        fn_ = new int[size];
        children_ = new int[size];
    }
};

class FilterVertices
{
public:

    static void allocateBuffer(Graph &data_graph, Graph &query_graph, int** &candidates, int* &candidates_count)
    {
        int query_vertex_num = query_graph.getVerticesCount();
        int candidates_max_num = data_graph.getGraphMaxLabelFrequency();

        candidates_count = new int[query_vertex_num];
        memset(candidates_count, 0, sizeof(int) * query_vertex_num);

        candidates = new int*[query_vertex_num];

        for (int i = 0; i < query_vertex_num; ++i) {
            candidates[i] = new int[candidates_max_num];
        }
    }

    static bool LDFFilter(Graph &data_graph, Graph &query_graph, int** &candidates, int* &candidates_count)
    {   
        allocateBuffer(data_graph, query_graph, candidates, candidates_count);

        for (int i = 0; i < query_graph.getVerticesCount(); ++i) {
            int label = query_graph.getVertexLabel(i);
            int degree = query_graph.getVertexDegree(i);

            int data_vertex_num;
            
            const int* data_vertices = data_graph.getVerticesByLabel(label, data_vertex_num);

            for (int j = 0; j < data_vertex_num; ++j) {
                int data_vertex = data_vertices[j];
            
                if (data_graph.getVertexDegree(data_vertex) >= degree) 
                { // deg_G >= deg_Q
                    
                    candidates[i][candidates_count[i]++] = data_vertex;
                }
            }

            if (candidates_count[i] == 0) {
                return false;
            }
        }

        return true;
    }

    static void computeCandidateWithLDF(Graph &data_graph, Graph &query_graph, int query_vertex,
                                        int &count) 
    {
        int label = query_graph.getVertexLabel(query_vertex);
        int degree = query_graph.getVertexDegree(query_vertex);
        count = 0;
        int data_vertex_num;
        const int* data_vertices = data_graph.getVerticesByLabel(label, data_vertex_num);

        for (int i = 0; i < data_vertex_num; ++i) {
            int v = data_vertices[i];
            if (data_graph.getVertexDegree(v) >= degree) {
                count += 1;
            }
        }
    }

    static int selectDPisoStartVertex(Graph &data_graph, Graph &query_graph, int *&candidates_count) 
    {
        double min_score = data_graph.getVerticesCount();
        int start_vertex = 0;

        for (int i = 0; i < query_graph.getVerticesCount(); ++i) {
            int degree = query_graph.getVertexDegree(i);
            if (degree <= 1)
                continue;

            double cur_score = candidates_count[i] / (double)degree;
            if (cur_score < min_score) {
                min_score = cur_score;
                start_vertex = i;
            }
        }
        return start_vertex;
    }


    static void bfsTraversal(Graph &graph, int root_vertex, TreeNode * &tree, int * &bfs_order) 
    {
        int vertex_num = graph.getVerticesCount();

        std::queue<int> bfs_queue;
        std::vector<bool> visited(vertex_num, false);

        tree = new TreeNode[vertex_num];
        for (int i = 0; i < vertex_num; ++i) {
            tree[i].initialize(vertex_num);
        }
        bfs_order = new int[vertex_num];

        int visited_vertex_count = 0;
        bfs_queue.push(root_vertex);
        visited[root_vertex] = true;
        tree[root_vertex].level_ = 0;
        tree[root_vertex].id_ = root_vertex;

        while(!bfs_queue.empty()) {
            const int u = bfs_queue.front();
            bfs_queue.pop();
            bfs_order[visited_vertex_count++] = u;

            int u_nbrs_count;
            const int* u_nbrs = graph.getVertexNeighbors(u, u_nbrs_count);
            for (int i = 0; i < u_nbrs_count; ++i) {
                int u_nbr = u_nbrs[i];

                if (!visited[u_nbr]) {
                    bfs_queue.push(u_nbr);
                    visited[u_nbr] = true;
                    tree[u_nbr].id_ = u_nbr;
                    tree[u_nbr].parent_ = u;
                    tree[u_nbr].level_ = tree[u].level_ + 1; // parent.level+1
                    tree[u].children_[tree[u].children_count_++] = u_nbr;
                }
            }
        }
    }

    static void generateDPisoFilterPlan(Graph &data_graph, Graph &query_graph, TreeNode *&tree,
                                                    int *&order, int *&candidates_count) 
    {
        int start_vertex = selectDPisoStartVertex(data_graph, query_graph, candidates_count);
        bfsTraversal(query_graph, start_vertex, tree, order);

        int query_vertices_num = query_graph.getVerticesCount();
        std::vector<int> order_index(query_vertices_num);
        for (int i = 0; i < query_vertices_num; ++i) {
            int query_vertex = order[i];
            order_index[query_vertex] = i;
        }

        for (int i = 0; i < query_vertices_num; ++i) {
            int u = order[i];
            tree[u].under_level_count_ = 0;
            tree[u].bn_count_ = 0;
            tree[u].fn_count_ = 0;

            int u_nbrs_count;
            const int* u_nbrs = query_graph.getVertexNeighbors(u, u_nbrs_count);
            for (int j = 0; j < u_nbrs_count; ++j) {
                int u_nbr = u_nbrs[j];
                if (order_index[u_nbr] < order_index[u]) {
                    tree[u].bn_[tree[u].bn_count_++] = u_nbr;
                }
                else {
                    tree[u].fn_[tree[u].fn_count_++] = u_nbr;
                }
            }
        }
    }

    static void pruneCandidates(Graph &data_graph, Graph &query_graph, int query_vertex, int *pivot_vertices, 
                        int pivot_vertices_count, int **candidates, int *candidates_count, 
                        int *flag, int *updated_flag) 
    {
        int query_vertex_label = query_graph.getVertexLabel(query_vertex);
        int query_vertex_degree = query_graph.getVertexDegree(query_vertex);

        int count = 0;
        int updated_flag_count = 0;
        for (int i = 0; i < pivot_vertices_count; ++i) {
            int pivot_vertex = pivot_vertices[i];

            for (int j = 0; j < candidates_count[pivot_vertex]; ++j) {
                int v = candidates[pivot_vertex][j];

                if (v == INVALID_VERTEX_ID)
                    continue;
                int v_nbrs_count;
                const int* v_nbrs = data_graph.getVertexNeighbors(v, v_nbrs_count);

                for (int k = 0; k < v_nbrs_count; ++k) {
                    int v_nbr = v_nbrs[k];
                    int v_nbr_label = data_graph.getVertexLabel(v_nbr);
                    int v_nbr_degree = data_graph.getVertexDegree(v_nbr);

                    if (flag[v_nbr] == count && v_nbr_label == query_vertex_label && v_nbr_degree >= query_vertex_degree) {
                        flag[v_nbr] += 1;

                        if (count == 0) {
                            updated_flag[updated_flag_count++] = v_nbr;
                        }
                    }
                }
            }
            count += 1;
        }

        for (int i = 0; i < candidates_count[query_vertex]; ++i) {
            int v = candidates[query_vertex][i];
            if (v == INVALID_VERTEX_ID)
                continue;

            if (flag[v] != count) {
                candidates[query_vertex][i] = INVALID_VERTEX_ID;
            }
        }

        for (int i = 0; i < updated_flag_count; ++i) {
            int v = updated_flag[i];
            flag[v] = 0;
        }
    }

    static void compactCandidates(int **&candidates, int *&candidates_count, int query_vertex_num) {
        for (int i = 0; i < query_vertex_num; ++i) {
            int query_vertex = i;
            int next_position = 0;
            for (int j = 0; j < candidates_count[query_vertex]; ++j) {
                int data_vertex = candidates[query_vertex][j];

                if (data_vertex != INVALID_VERTEX_ID) {
                    candidates[query_vertex][next_position++] = data_vertex;
                }
            }
            candidates_count[query_vertex] = next_position;
        }
    }


    static bool isCandidateSetValid(int **&candidates, int *&candidates_count, int query_vertex_num) {
        for (int i = 0; i < query_vertex_num; ++i) {
            if (candidates_count[i] == 0)
                return false;
        }
        return true;
    }

    static bool DPisoFilter(Graph &data_graph, Graph &query_graph, int **&candidates, int *&candidates_count, 
                            int *&order, TreeNode *&tree)
    {   
        if (!LDFFilter(data_graph, query_graph, candidates, candidates_count))
            return false;

        generateDPisoFilterPlan(data_graph, query_graph, tree, order, candidates_count);

        int query_vertices_num = query_graph.getVerticesCount();
        int* updated_flag = new int[data_graph.getVerticesCount()];
        int* flag = new int[data_graph.getVerticesCount()];
        std::fill(flag, flag + data_graph.getVerticesCount(), 0);

        // The number of refinement is k. According to the original paper, we set k as 3.
        for (int k = 0; k < 3; ++k) {
            if (k % 2 == 0) {
                for (int i = 1; i < query_vertices_num; ++i) {
                    int query_vertex = order[i];
                    TreeNode& node = tree[query_vertex];
                    pruneCandidates(data_graph, query_graph, query_vertex, node.bn_, node.bn_count_, candidates, candidates_count, flag, updated_flag);
                }
            }
            else {
                for (int i = query_vertices_num - 2; i >= 0; --i) {
                    int query_vertex = order[i];
                    TreeNode& node = tree[query_vertex];
                    pruneCandidates(data_graph, query_graph, query_vertex, node.fn_, node.fn_count_, candidates, candidates_count, flag, updated_flag);
                }
            }
        }

        compactCandidates(candidates, candidates_count, query_graph.getVerticesCount());

        delete[] updated_flag;
        delete[] flag;
        return isCandidateSetValid(candidates, candidates_count, query_graph.getVerticesCount());
    }
};