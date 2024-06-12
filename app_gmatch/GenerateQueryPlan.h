#pragma once

#include "graph.h"
#include <vector>

class GenerateQueryPlan
{
public:

    static ui selectGQLStartVertex(Graph &query_graph, ui *candidates_count) 
    {
        /**
         * Select the vertex with the minimum number of candidates as the start vertex.
         * Tie Handling:
         *  1. degree
         *  2. label id
         */

        ui start_vertex = 0;

        for (ui i = 1; i < query_graph.getVerticesCount(); ++i) {
            ui cur_vertex = i;

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


    static void updateValidVertices(Graph &query_graph, ui query_vertex, std::vector<bool> &visited,
                                    std::vector<bool> &adjacent) 
    {
        visited[query_vertex] = true;
        ui nbr_cnt;
        const ui* nbrs = query_graph.getVertexNeighbors(query_vertex, nbr_cnt);

        for (ui i = 0; i < nbr_cnt; ++i) {
            ui nbr = nbrs[i];
            adjacent[nbr] = true;
        }
    }


    static void generateGQLQueryPlan(Graph &data_graph, Graph &query_graph, ui *candidates_count,
                                    ui *&order, ui *&pivot) 
    {
        /**
         * Select the vertex v such that 
         * (1) v is adjacent to the selected vertices;  
         * (2) v has the minimum number of candidates.
         */
        std::vector<bool> visited_vertices(query_graph.getVerticesCount(), false);
        std::vector<bool> adjacent_vertices(query_graph.getVerticesCount(), false);
        order = new ui[query_graph.getVerticesCount()];
        pivot = new ui[query_graph.getVerticesCount()];

        ui start_vertex = selectGQLStartVertex(query_graph, candidates_count);
        order[0] = start_vertex;
        updateValidVertices(query_graph, start_vertex, visited_vertices, adjacent_vertices);

        for (ui i = 1; i < query_graph.getVerticesCount(); ++i) {
            ui next_vertex;
            ui min_value = data_graph.getVerticesCount() + 1;
            for (ui j = 0; j < query_graph.getVerticesCount(); ++j) {
                ui cur_vertex = j;

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
        // for (ui i = 1; i < query_graph.getVerticesCount(); ++i) {
        //     ui u = order[i];
        //     for (ui j = 0; j < i; ++j) {
        //         ui cur_vertex = order[j];
        //         if (query_graph.checkEdgeExistence(u, cur_vertex, query_graph.getVertexLabel(u))) { // TODO, fix
        //             pivot[i] = cur_vertex;
        //             break;
        //         }
        //     }
        // }
    }
    
};