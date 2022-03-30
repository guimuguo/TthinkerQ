#pragma once

#include <unordered_map>
#include <iostream>
#include <fstream>
#include <vector>
#include <sys/timeb.h>

class Graph {
private:

    int vertices_count_;
    int edges_count_;
    int labels_count_;
    int max_degree_;
    int max_label_frequency_;

    int* offsets_;
    int * neighbors_;
    int* labels_;
    int* reverse_index_offsets_;
    int* reverse_index_;

    int* core_table_;
    int core_length_;

    std::unordered_map<int, int> labels_frequency_;

    int* labels_offsets_;

private:
    void BuildReverseIndex();
    void BuildLabelOffset();

public:
    timeb gtime_start[40];

    Graph() {
        vertices_count_ = 0;
        edges_count_ = 0;
        labels_count_ = 0;
        max_degree_ = 0;
        max_label_frequency_ = 0;
        core_length_ = 0;

        offsets_ = NULL;
        neighbors_ = NULL;
        labels_ = NULL;
        reverse_index_offsets_ = NULL;
        reverse_index_ = NULL;
        core_table_ = NULL;
        labels_frequency_.clear();
    }

    ~Graph() {
        delete[] offsets_;
        delete[] neighbors_;
        delete[] labels_;
        delete[] reverse_index_offsets_;
        delete[] reverse_index_;
        delete[] core_table_;
    }

public:
    void loadGraphFromFile(std::string& file_path);
    void printGraphMetaData();

public:
    const int getLabelsCount() const {
        return labels_count_;
    }

    const int getVerticesCount() const {
        return vertices_count_;
    }

    const int getEdgesCount() const {
        return edges_count_;
    }

    const int getGraphMaxDegree() const {
        return max_degree_;
    }

    const int getGraphMaxLabelFrequency() const {
        return max_label_frequency_;
    }

    const int getVertexDegree(const int id) const {
        return offsets_[id + 1] - offsets_[id];
    }

    const int getLabelsFrequency(const int label) const {
        return labels_frequency_.find(label) == labels_frequency_.end() ? 0 : labels_frequency_.at(label);
    }

    const int getCoreValue(const int id) const {
        return core_table_[id];
    }

    const int get2CoreSize() const {
        return core_length_;
    }

    const int getVertexLabel(const int id) const {
        return labels_[id];
    }

    const int * getVertexNeighbors(const int id, int& count) const {
        count = offsets_[id + 1] - offsets_[id];
        return neighbors_ + offsets_[id];
    }

    const int * getVerticesByLabel(const int id, int& count) const {
        count = reverse_index_offsets_[id + 1] - reverse_index_offsets_[id];
        return reverse_index_ + reverse_index_offsets_[id];
    }

    const int * getNeighborsByLabel(const int id, const int label, int& count) const {
        int offset = id * labels_count_ + label;
        count = labels_offsets_[offset + 1] - labels_offsets_[offset];
        return neighbors_ + labels_offsets_[offset];
    }

    bool checkEdgeExistence(const int u, const int v, const int u_label) const {
        int count = 0;
        const int* neighbors = getNeighborsByLabel(v, u_label, count);
        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }
        return false;
    }

    bool checkEdgeExistence(int u, int v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        int count = 0;
        const int* neighbors = getVertexNeighbors(v, count);

        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

    void buildCoreTable();
};

void getKCore(const Graph *graph, int *core_table) {
    int vertices_count = graph->getVerticesCount();
    int max_degree = graph->getGraphMaxDegree();

    int* vertices = new int[vertices_count];          // Vertices sorted by degree.
    int* position = new int[vertices_count];          // The position of vertices in vertices array.
    int* degree_bin = new int[max_degree + 1];      // Degree from 0 to max_degree.
    int* offset = new int[max_degree + 1];          // The offset in vertices array according to degree.

    std::fill(degree_bin, degree_bin + (max_degree + 1), 0);

    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexDegree(i);
        core_table[i] = degree;
        degree_bin[degree] += 1;
    }

    int start = 0;
    for (int i = 0; i < max_degree + 1; ++i) {
        offset[i] = start;
        start += degree_bin[i];
    }

    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexDegree(i);
        position[i] = offset[degree];
        vertices[position[i]] = i;
        offset[degree] += 1;
    }

    for (int i = max_degree; i > 0; --i) {
        offset[i] = offset[i - 1];
    }
    offset[0] = 0;

    for (int i = 0; i < vertices_count; ++i) {
        int v = vertices[i];

        int count;
        const int * neighbors = graph->getVertexNeighbors(v, count);

        for(int j = 0; j < count; ++j) {
            int u = neighbors[j];

            if (core_table[u] > core_table[v]) {

                // Get the position and vertex which is with the same degree
                // and at the start position of vertices array.
                int cur_degree_u = core_table[u];
                int position_u = position[u];
                int position_w = offset[cur_degree_u];
                int w = vertices[position_w];

                if (u != w) {
                    // Swap u and w.
                    position[u] = position_w;
                    position[w] = position_u;
                    vertices[position_u] = w;
                    vertices[position_w] = u;
                }

                offset[cur_degree_u] += 1;
                core_table[u] -= 1;
            }
        }
    }

    delete[] vertices;
    delete[] position;
    delete[] degree_bin;
    delete[] offset;
}

void Graph::buildCoreTable() {
    core_table_ = new int[vertices_count_];
    getKCore(this, core_table_);
    for (int i = 0; i < vertices_count_; ++i) {
        if (core_table_[i] > 1) {
            core_length_ += 1;
        }
    }
}

void Graph::BuildReverseIndex() {
    reverse_index_ = new int[vertices_count_];
    reverse_index_offsets_= new int[labels_count_ + 1];
    reverse_index_offsets_[0] = 0;

    int total = 0;
    for (int i = 0; i < labels_count_; ++i) {
        reverse_index_offsets_[i + 1] = total;
        total += labels_frequency_[i];
    }

    for (int i = 0; i < vertices_count_; ++i) {
        int label = labels_[i];
        reverse_index_[reverse_index_offsets_[label + 1]++] = i;
    }
}

void Graph::BuildLabelOffset() {
    size_t labels_offset_size = (size_t)vertices_count_ * labels_count_ + 1;
    labels_offsets_ = new int[labels_offset_size];
    std::fill(labels_offsets_, labels_offsets_ + labels_offset_size, 0);


    // sort by label, then by ID
    for (int i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1],
            [this](const int u, const int v) -> bool {
                return labels_[u] == labels_[v] ? u < v : labels_[u] < labels_[v];
            });
    }

    for (int i = 0; i < vertices_count_; ++i) {
        int previous_label = 0;
        int current_label = 0;

        labels_offset_size = i * labels_count_;
        labels_offsets_[labels_offset_size] = offsets_[i];

        for (int j = offsets_[i]; j < offsets_[i + 1]; ++j) {
            current_label = labels_[neighbors_[j]];

            if (current_label != previous_label) {
                for (int k = previous_label + 1; k <= current_label; ++k) {
                    labels_offsets_[labels_offset_size + k] = j;
                }
                previous_label = current_label;
            }
        }

        for (int l = current_label + 1; l <= labels_count_; ++l) {
            labels_offsets_[labels_offset_size + l] = offsets_[i + 1];
        }
    }
}

/**
std::pair<int,int> graph_shrink(std::string &file_path, std::vector<int> &degrees)
{   
    int vertices_count, edges_count;
    std::ifstream infile1(file_path);
    char type;
    infile1 >> type >> vertices_count >> edges_count;
    degrees.assign(vertices_count, 0);
    vertices_count = 0;
    edges_count = 0;

    while (infile1 >> type) {
        if(type == 'v') {
            int id;
            int label;
            int degree;
            infile1 >> id >> label >> degree;
            vertices_count++;
        }
        else if (type == 'e') {  
            int begin;
            int end;
            infile1 >> begin >> end;
            if(begin >= 10000 || end >= 10000) continue;
            degrees[begin]++;
            degrees[end]++;

            edges_count++;
        }
    }
    infile1.close();
    return std::make_pair(vertices_count, edges_count);
}
**/

void Graph::loadGraphFromFile(std::string &file_path) {

    std::ifstream infile(file_path);

    if (!infile.is_open()) {
        std::cout << "Can not open the graph file " << file_path << " ." << std::endl;
        exit(-1);
    }

    char type;
    infile >> type >> vertices_count_ >> edges_count_;

    std::cout << vertices_count_ << " " << edges_count_ << std::endl;

    offsets_ = new int[vertices_count_ +  1];
    offsets_[0] = 0;

    neighbors_ = new int[edges_count_ * 2];
    labels_ = new int[vertices_count_];
    labels_count_ = 0;
    max_degree_ = 0;

    int max_label_id = 0;
    std::vector<int> neighbors_offset(vertices_count_, 0);

    while (infile >> type) {
        if (type == 'v') { // Read vertex.
            int id;
            int label;
            int degree;
            infile >> id >> label >> degree;

            labels_[id] = label;
            offsets_[id + 1] = offsets_[id] + degree;

            if (degree > max_degree_) {
                max_degree_ = degree;
            }

            if (labels_frequency_.find(label) == labels_frequency_.end()) {
                labels_frequency_[label] = 0;
                if (label > max_label_id)
                    max_label_id = label;
            }

            labels_frequency_[label] += 1;
        }
        else if (type == 'e') { // Read edge.
            int begin;
            int end;
            infile >> begin >> end;

            int offset = offsets_[begin] + neighbors_offset[begin];
            neighbors_[offset] = end;

            offset = offsets_[end] + neighbors_offset[end];
            neighbors_[offset] = begin;

            neighbors_offset[begin] += 1;
            neighbors_offset[end] += 1;
        }
    }


    infile.close();
    labels_count_ = (int)labels_frequency_.size() > (max_label_id + 1) ? (int)labels_frequency_.size() : max_label_id + 1;

    for (auto element: labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    for (int i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
    }

    
    BuildReverseIndex();
    
    BuildLabelOffset();
}