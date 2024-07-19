#pragma once

#include <unordered_map>
#include <iostream>
#include <fstream>
#include <vector>
#include <sys/timeb.h>

typedef unsigned int ui;

class Graph {
private:

    ui vertices_count_;
    ui edges_count_;
    ui labels_count_;
    ui max_degree_;
    ui max_label_frequency_;

    ui* offsets_;
    ui * neighbors_;
    ui* labels_;
    ui* reverse_index_offsets_;
    ui* reverse_index_;

    ui* core_table_;
    ui core_length_;

    std::unordered_map<ui, ui> labels_frequency_;

    ui* labels_offsets_;

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
    void pruiGraphMetaData();

public:
    const ui getLabelsCount() const {
        return labels_count_;
    }

    const ui getVerticesCount() const {
        return vertices_count_;
    }

    const ui getEdgesCount() const {
        return edges_count_;
    }

    const ui getGraphMaxDegree() const {
        return max_degree_;
    }

    const ui getGraphMaxLabelFrequency() const {
        return max_label_frequency_;
    }

    const ui getVertexDegree(const ui id) const {
        return offsets_[id + 1] - offsets_[id];
    }

    const ui getLabelsFrequency(const ui label) const {
        return labels_frequency_.find(label) == labels_frequency_.end() ? 0 : labels_frequency_.at(label);
    }

    const ui getCoreValue(const ui id) const {
        return core_table_[id];
    }

    const ui get2CoreSize() const {
        return core_length_;
    }

    const ui getVertexLabel(const ui id) const {
        return labels_[id];
    }

    const ui * getVertexNeighbors(const ui id, ui& count) const {
        count = offsets_[id + 1] - offsets_[id];
        return neighbors_ + offsets_[id];
    }

    const ui * getVerticesByLabel(const ui id, ui& count) const {
        count = reverse_index_offsets_[id + 1] - reverse_index_offsets_[id];
        return reverse_index_ + reverse_index_offsets_[id];
    }

    const ui * getNeighborsByLabel(const ui id, const ui label, ui& count) const {
        ui offset = id * labels_count_ + label;
        count = labels_offsets_[offset + 1] - labels_offsets_[offset];
        return neighbors_ + labels_offsets_[offset];
    }

    bool checkEdgeExistence(const ui u, const ui v, const ui u_label) const {
        ui count = 0;
        const ui* neighbors = getNeighborsByLabel(v, u_label, count);
        ui begin = 0;
        ui end = count - 1;
        while (begin <= end) {
            ui mid = begin + ((end - begin) >> 1);
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

    bool checkEdgeExistence(ui u, ui v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        ui count = 0;
        const ui* neighbors = getVertexNeighbors(v, count);

        ui begin = 0;
        ui end = count - 1;
        while (begin <= end) {
            ui mid = begin + ((end - begin) >> 1);
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

void getKCore(const Graph *graph, ui *core_table) {
    ui vertices_count = graph->getVerticesCount();
    ui max_degree = graph->getGraphMaxDegree();

    ui* vertices = new ui[vertices_count];          // Vertices sorted by degree.
    ui* position = new ui[vertices_count];          // The position of vertices in vertices array.
    ui* degree_bin = new ui[max_degree + 1];      // Degree from 0 to max_degree.
    ui* offset = new ui[max_degree + 1];          // The offset in vertices array according to degree.

    std::fill(degree_bin, degree_bin + (max_degree + 1), 0);

    for (ui i = 0; i < vertices_count; ++i) {
        ui degree = graph->getVertexDegree(i);
        core_table[i] = degree;
        degree_bin[degree] += 1;
    }

    ui start = 0;
    for (ui i = 0; i < max_degree + 1; ++i) {
        offset[i] = start;
        start += degree_bin[i];
    }

    for (ui i = 0; i < vertices_count; ++i) {
        ui degree = graph->getVertexDegree(i);
        position[i] = offset[degree];
        vertices[position[i]] = i;
        offset[degree] += 1;
    }

    for (ui i = max_degree; i > 0; --i) {
        offset[i] = offset[i - 1];
    }
    offset[0] = 0;

    for (ui i = 0; i < vertices_count; ++i) {
        ui v = vertices[i];

        ui count;
        const ui * neighbors = graph->getVertexNeighbors(v, count);

        for(ui j = 0; j < count; ++j) {
            ui u = neighbors[j];

            if (core_table[u] > core_table[v]) {

                // Get the position and vertex which is with the same degree
                // and at the start position of vertices array.
                ui cur_degree_u = core_table[u];
                ui position_u = position[u];
                ui position_w = offset[cur_degree_u];
                ui w = vertices[position_w];

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
    core_table_ = new ui[vertices_count_];
    getKCore(this, core_table_);
    for (ui i = 0; i < vertices_count_; ++i) {
        if (core_table_[i] > 1) {
            core_length_ += 1;
        }
    }
}

void Graph::BuildReverseIndex() {
    reverse_index_ = new ui[vertices_count_];
    reverse_index_offsets_= new ui[labels_count_ + 1];
    reverse_index_offsets_[0] = 0;

    ui total = 0;
    for (ui i = 0; i < labels_count_; ++i) {
        reverse_index_offsets_[i + 1] = total;
        total += labels_frequency_[i];
    }

    for (ui i = 0; i < vertices_count_; ++i) {
        ui label = labels_[i];
        reverse_index_[reverse_index_offsets_[label + 1]++] = i;
    }
}

void Graph::BuildLabelOffset() {
    size_t labels_offset_size = (size_t)vertices_count_ * labels_count_ + 1;
    labels_offsets_ = new ui[labels_offset_size];
    std::fill(labels_offsets_, labels_offsets_ + labels_offset_size, 0);


    // sort by label, then by ID
    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1],
            [this](const ui u, const ui v) -> bool {
                return labels_[u] == labels_[v] ? u < v : labels_[u] < labels_[v];
            });
    }

    for (ui i = 0; i < vertices_count_; ++i) {
        ui previous_label = 0;
        ui current_label = 0;

        labels_offset_size = i * labels_count_;
        labels_offsets_[labels_offset_size] = offsets_[i];

        for (ui j = offsets_[i]; j < offsets_[i + 1]; ++j) {
            current_label = labels_[neighbors_[j]];

            if (current_label != previous_label) {
                for (ui k = previous_label + 1; k <= current_label; ++k) {
                    labels_offsets_[labels_offset_size + k] = j;
                }
                previous_label = current_label;
            }
        }

        for (ui l = current_label + 1; l <= labels_count_; ++l) {
            labels_offsets_[labels_offset_size + l] = offsets_[i + 1];
        }
    }
}

/**
std::pair<ui,ui> graph_shrink(std::string &file_path, std::vector<ui> &degrees)
{   
    ui vertices_count, edges_count;
    std::ifstream infile1(file_path);
    char type;
    infile1 >> type >> vertices_count >> edges_count;
    degrees.assign(vertices_count, 0);
    vertices_count = 0;
    edges_count = 0;

    while (infile1 >> type) {
        if(type == 'v') {
            ui id;
            ui label;
            ui degree;
            infile1 >> id >> label >> degree;
            vertices_count++;
        }
        else if (type == 'e') {  
            ui begin;
            ui end;
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

    offsets_ = new ui[vertices_count_ +  1];
    offsets_[0] = 0;

    neighbors_ = new ui[edges_count_ * 2];
    labels_ = new ui[vertices_count_];
    labels_count_ = 0;
    max_degree_ = 0;

    ui max_label_id = 0;
    std::vector<ui> neighbors_offset(vertices_count_, 0);

    while (infile >> type) {
        if (type == 'v') { // Read vertex.
            ui id;
            ui label;
            ui degree;
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
            ui begin;
            ui end;
            infile >> begin >> end;

            ui offset = offsets_[begin] + neighbors_offset[begin];
            neighbors_[offset] = end;

            offset = offsets_[end] + neighbors_offset[end];
            neighbors_[offset] = begin;

            neighbors_offset[begin] += 1;
            neighbors_offset[end] += 1;
        }
    }


    infile.close();
    labels_count_ = (ui)labels_frequency_.size() > (max_label_id + 1) ? (ui)labels_frequency_.size() : max_label_id + 1;

    for (auto element: labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
    }

    
    BuildReverseIndex();

    // BuildLabelOffset();
}