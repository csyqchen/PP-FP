#ifndef KCORE_H
#define KCORE_H

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <set> // Include for std::set
#include <ctime>
#include "Graph.h"

class Kcore {
public:
    // Move Node struct to the top of the class definition
    std::set<int> built; //查询过conditional tree但无法构成kcore的结点集

    struct Node {
        int id;
        int count;
        int totalCount;
        int level;
        int step;
        int coreness;
        std::set<int> branches; // 记录结点所属的分支
    };

    // Constructors
    Kcore(Graph* p_graph, int v, int k);
    Kcore(Graph* p_graph, int k);

    // Utility functions
    void deleteNode(int v);
    std::unordered_set<int> getNodes();
    Graph::Attributes attrs(int v);
    int deg(int v);

    template<typename T>
    std::unordered_set<T> intersect(std::unordered_set<T>, std::unordered_set<T>);
    template<typename T>
    void intersect_filter(std::unordered_set<T>& set1, const std::unordered_set<T>& set2);

    bool check_kcore(Graph* graph, int k, std::vector<int>& kcore_nodes);
    void search_kcore(Graph* graph);

    template<typename T>
    void backtrack(const std::vector<T>& arr, int k, int start, std::vector<T>& path, std::vector<std::vector<T>>& result);
    template<typename T>
    std::vector<std::vector<T>> combinations(const std::vector<T>& arr, int k);
    
    // Baselines
    double baseline1(std::ofstream& output_file, size_t& success_d, size_t& community);
    double baseline2(std::ofstream& output_file, size_t& success_d, size_t& community);
    void baseline3();
    
    // Our method
    bool isKCore(int k, std::vector<int>& nodes);
    std::tuple<double, int, int> fptree(int k, const std::string& fptreeFile, const std::string& corenessFile, const std::string& outputFileName);
    std::vector<std::vector<Node>> getConditionalPatternBase(int targetNodeID, int rootNode, const std::vector<std::vector<Node>>& branches);
    std::pair<std::vector<int>, int> buildConditionalFPTree(int targetNodeID, int rootNode, const std::vector<std::vector<Node>>& branches);

    // Display results
    void display();

private:
    Graph* graph;
    int v;
    int k;
};

#endif // KCORE_H