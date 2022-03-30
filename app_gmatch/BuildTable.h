#pragma once

#include <stdlib.h>
#include "graph.h"
#include <algorithm>

class Edges {
public:
    int* offset_;
    int* edge_;
    int vertex_count_;
    int edge_count_;
    int max_degree_;
public:
    Edges() {
        offset_ = NULL;
        edge_ = NULL;
        vertex_count_ = 0;
        edge_count_ = 0;
        max_degree_ = 0;
    }

    ~Edges() {
        delete[] offset_;
        delete[] edge_;
    }
};

class BuildTable {
public:

    static void buildTable(Graph &data_graph, Graph &query_graph, int **candidates, 
                            int *candidates_count, Edges ***edge_matrix) 
    {   
        int query_vertices_num = query_graph.getVerticesCount();
        int* flag = new int[data_graph.getVerticesCount()];
        int* updated_flag = new int[data_graph.getVerticesCount()];
        std::fill(flag, flag + data_graph.getVerticesCount(), 0);


        for (int i = 0; i < query_vertices_num; ++i) {
            for (int j = 0; j < query_vertices_num; ++j) {
                edge_matrix[i][j] = NULL;
            }
        }

        std::vector<int> build_table_order(query_vertices_num);
        for (int i = 0; i < query_vertices_num; ++i) {
            build_table_order[i] = i;
        }

        std::sort(build_table_order.begin(), build_table_order.end(), [&query_graph](int l, int r) {
            if (query_graph.getVertexDegree(l) == query_graph.getVertexDegree(r)) {
                return l < r;
            }
            return query_graph.getVertexDegree(l) > query_graph.getVertexDegree(r);
        });

        std::vector<int> temp_edges(data_graph.getEdgesCount() * 2);

        for (auto u : build_table_order) {
            int u_nbrs_count;
            const int* u_nbrs = query_graph.getVertexNeighbors(u, u_nbrs_count);

            int updated_flag_count = 0;

            for (int i = 0; i < u_nbrs_count; ++i) {
                int u_nbr = u_nbrs[i];

                if (edge_matrix[u][u_nbr] != NULL)
                    continue;

                if (updated_flag_count == 0) { // record u's CS only once 
                    for (int j = 0; j < candidates_count[u]; ++j) {
                        int v = candidates[u][j];
                        flag[v] = j + 1; // v's ID -> index in u's CS
                        updated_flag[updated_flag_count++] = v;
                    }
                }

                edge_matrix[u_nbr][u] = new Edges;
                edge_matrix[u_nbr][u]->vertex_count_ = candidates_count[u_nbr];
                edge_matrix[u_nbr][u]->offset_ = new int[candidates_count[u_nbr] + 1];

                edge_matrix[u][u_nbr] = new Edges;
                edge_matrix[u][u_nbr]->vertex_count_ = candidates_count[u];
                edge_matrix[u][u_nbr]->offset_ = new int[candidates_count[u] + 1];
                std::fill(edge_matrix[u][u_nbr]->offset_, edge_matrix[u][u_nbr]->offset_ + candidates_count[u] + 1, 0);

                int local_edge_count = 0;
                int local_max_degree = 0;

                for (int j = 0; j < candidates_count[u_nbr]; ++j) {
                    int v = candidates[u_nbr][j];
                    edge_matrix[u_nbr][u]->offset_[j] = local_edge_count;

                    int v_nbrs_count;
                    const int* v_nbrs = data_graph.getVertexNeighbors(v, v_nbrs_count);

                    int local_degree = 0;

                    for (int k = 0; k < v_nbrs_count; ++k) {
                        int v_nbr = v_nbrs[k]; // u's neighors's candidate's neighbor (in G)

                        if (flag[v_nbr] != 0) { // v_nbr in u's CS
                            int position = flag[v_nbr] - 1; // get index in u's CS
                            temp_edges[local_edge_count++] = position;
                            edge_matrix[u][u_nbr]->offset_[position + 1] += 1;
                            local_degree += 1;
                        }
                    }

                    if (local_degree > local_max_degree) {
                        local_max_degree = local_degree;
                    }
                }

                edge_matrix[u_nbr][u]->offset_[candidates_count[u_nbr]] = local_edge_count;
                edge_matrix[u_nbr][u]->max_degree_ = local_max_degree;
                edge_matrix[u_nbr][u]->edge_count_ = local_edge_count;
                edge_matrix[u_nbr][u]->edge_ = new int[local_edge_count];
                std::copy(temp_edges.begin(), temp_edges.begin() + local_edge_count, edge_matrix[u_nbr][u]->edge_);

                edge_matrix[u][u_nbr]->edge_count_ = local_edge_count;
                edge_matrix[u][u_nbr]->edge_ = new int[local_edge_count];

                local_max_degree = 0;
                for (int j = 1; j <= candidates_count[u]; ++j) {
                    if (edge_matrix[u][u_nbr]->offset_[j] > local_max_degree) {
                        local_max_degree = edge_matrix[u][u_nbr]->offset_[j];
                    }
                    edge_matrix[u][u_nbr]->offset_[j] += edge_matrix[u][u_nbr]->offset_[j - 1];
                }

                edge_matrix[u][u_nbr]->max_degree_ = local_max_degree;

                for (int j = 0; j < candidates_count[u_nbr]; ++j) {
                    int begin = j; // index in u_nbr
                    for (int k = edge_matrix[u_nbr][u]->offset_[begin]; k < edge_matrix[u_nbr][u]->offset_[begin + 1]; ++k) {
                        int end = edge_matrix[u_nbr][u]->edge_[k]; // index in u

                        edge_matrix[u][u_nbr]->edge_[edge_matrix[u][u_nbr]->offset_[end]++] = begin;
                    }
                }

                for (int j = candidates_count[u]; j >= 1; --j) { // reset 
                    edge_matrix[u][u_nbr]->offset_[j] = edge_matrix[u][u_nbr]->offset_[j - 1];
                }
                edge_matrix[u][u_nbr]->offset_[0] = 0;
            }

            for (int i = 0; i < updated_flag_count; ++i) {
                int v = updated_flag[i];
                flag[v] = 0;
            }
        }
    }
};
