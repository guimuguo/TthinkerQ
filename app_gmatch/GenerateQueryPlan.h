#pragma once

#include "graph.h"
#include <vector>

class GenerateQueryPlan
{
public:

    static int selectGQLStartVertex(Graph &query_graph, int *candidates_count) 
    {
        /**
         * Select the vertex with the minimum number of candidates as the start vertex.
         * Tie Handling:
         *  1. degree
         *  2. label id
         */

        int start_vertex = 0;

        for (int i = 1; i < query_graph.getVerticesCount(); ++i) {
            int cur_vertex = i;

            if (candidates_count[cur_vertex] < candidates_count[start_vertex]) {
                start_vertex = cur_vertex;
            }
            else if (candidates_count[cur_vertex] == candidates_count[start_vertex]
                    && query_graph.getVertexDegree(cur_vertex) > query_graph.getVertexDegree(start_vertex)) {
                start_vertex = cur_vertex;
            }
        }
        return start_vertex;
    }


    static void updateValidVertices(Graph &query_graph, int query_vertex, std::vector<bool> &visited,
                                    std::vector<bool> &adjacent) 
    {
        visited[query_vertex] = true;
        int nbr_cnt;
        const int* nbrs = query_graph.getVertexNeighbors(query_vertex, nbr_cnt);

        for (int i = 0; i < nbr_cnt; ++i) {
            int nbr = nbrs[i];
            adjacent[nbr] = true;
        }
    }


    static void generateGQLQueryPlan(Graph &data_graph, Graph &query_graph, int *candidates_count,
                                    int *&order, int *&pivot) 
    {
        /**
         * Select the vertex v such that 
         * (1) v is adjacent to the selected vertices;  
         * (2) v has the minimum number of candidates.
         */
        std::vector<bool> visited_vertices(query_graph.getVerticesCount(), false);
        std::vector<bool> adjacent_vertices(query_graph.getVerticesCount(), false);
        order = new int[query_graph.getVerticesCount()];
        pivot = new int[query_graph.getVerticesCount()];

        int start_vertex = selectGQLStartVertex(query_graph, candidates_count);
        order[0] = start_vertex;
        updateValidVertices(query_graph, start_vertex, visited_vertices, adjacent_vertices);

        for (int i = 1; i < query_graph.getVerticesCount(); ++i) {
            int next_vertex;
            int min_value = data_graph.getVerticesCount() + 1;
            for (int j = 0; j < query_graph.getVerticesCount(); ++j) {
                int cur_vertex = j;

                if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
                    if (candidates_count[cur_vertex] < min_value) {
                        min_value = candidates_count[cur_vertex];
                        next_vertex = cur_vertex;
                    }
                    else if (candidates_count[cur_vertex] == min_value && query_graph.getVertexDegree(cur_vertex) > query_graph.getVertexDegree(next_vertex)) {
                        next_vertex = cur_vertex;
                    }
                }
            }
            updateValidVertices(query_graph, next_vertex, visited_vertices, adjacent_vertices);
            order[i] = next_vertex;
        }

        // // Pick a pivot randomly.
        // for (int i = 1; i < query_graph.getVerticesCount(); ++i) {
        //     int u = order[i];
        //     for (int j = 0; j < i; ++j) {
        //         int cur_vertex = order[j];
        //         if (query_graph.checkEdgeExistence(u, cur_vertex, query_graph.getVertexLabel(u))) { // TODO, fix
        //             pivot[i] = cur_vertex;
        //             break;
        //         }
        //     }
        // }
    }
    
};