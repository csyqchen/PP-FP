#pragma once
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "Utils.h"

#ifndef GRAPH_H
#define GRAPH_H

class Graph {
public:
    typedef std::vector<Edge> Neighbors;
    typedef std::unordered_set<std::string> Attributes;

    Graph();
    Graph(const Graph& graph);
    Graph(const Graph* p_graph);
    ~Graph();

    //set data
    void setNodes(std::unordered_set<int> nodes);
    void setEdges(std::unordered_set<Edge, MyEdgeHash> edges);
    void setAttrs(std::unordered_map<int, Attributes> attr);
    void setNbs(std::unordered_map<int, Neighbors> ng);
    void setData(Graph* p_graph);

    // Insert / Delete
    void addEdge(int x, int y);
    void addEdge(Edge e);
    bool addNode(int v);
    void deleteEdge(Edge e);
    void deleteNode(int v);
    void unionGraph(const Graph& g);
    void add_attributes(int v, Attributes atts);
    void remove_attribute(int v, Attributes atts);
    void add_edge(int x, int y);

    const Neighbors& neighbors(int v) const;
    std::vector<int> getNeighborIDs(int v) const;

    // Existence check
    bool hasEdge(int x, int y) const;
    bool hasEdge(Edge e) const;
    bool hasNode(int v) const;

    // Properties
    int getMaxN() const;
    int getM() const;
    int getN() const;
    int deg(int v) const;
    const std::unordered_set<int>& getNodes() const;
    const std::unordered_set<Edge, MyEdgeHash>& getEdges() const;
    std::unordered_map<int, Neighbors> getNeighbors();
    std::unordered_map<int, Attributes> getAttrs();
    //const Neighbors& neighbors(int v);
    const Attributes& attributes(int v);

    //functions
    void load_graph(std::string data_path);
    void load_attribute(std::string data_path);
    void clear();


private:
    std::unordered_map<int, Neighbors> ng;
    std::unordered_set<int> nodes;
    std::unordered_set<Edge, MyEdgeHash> edges;
    std::unordered_map<int, Attributes> attr;
    int maxN;
};

#endif //GRAPH_H
