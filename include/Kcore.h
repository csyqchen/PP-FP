#ifndef KCORE_H
#define KCORE_H

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <queue>
#include <ctime>
#include "Graph.h"
#include "Utils.h"

struct TreeNode {
    int coreness;
    std::vector<int> nodes;
    std::vector<TreeNode*> children;

    TreeNode(int core) : coreness(core) {}
};

class Kcore {
public:
    std::set<int> built;

    std::vector<std::pair<std::set<std::string>, std::set<int>>> attrIntersect;

    struct Node {
        int id;
        int count;
        int totalCount;
        int level;
        std::set<std::string> attributes;
        // std::unordered_set<std::string> attributes;
        std::set<int> branches;
    };

    // Constructors
    Kcore(Graph* p_graph, int v, int k);
    Kcore(Graph* p_graph, int k);
    Kcore(Graph* p_graph);

    void set_k(int k);
    void set_query_node(int node);

    // Utility functions
    void deleteNode(int v);
    std::unordered_set<int> getNodes();
    Graph::Attributes attrs(int v);
    int getdeg(int v);

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
    double baseline1(std::ofstream& output_file, size_t& success_d, size_t& community, size_t& iteration_count);
    double baseline2(std::ofstream& output_file, size_t& success_d, size_t& community, size_t& iteration_count);
    void baseline3();

    // Index
    std::vector<int> kCoreDecomposeOptimized(const std::vector<int>& component);

    void printTree(TreeNode* node, int depth = 0) const;

    void printTreeToFile(TreeNode* node, std::ofstream& outFile, int depth = 0) const;

    int getNodeCoreness(int node) const;

    // Tree and Inverted Index Access Functions
    TreeNode* getTreeIndexRoot() const;

    // TNode* buildAttributeIndexAndTree(const std::string& attributeOutputFilePath, const std::string& treeOutputFilePath);
    void buildAttributeIndex();
    void handleALevel(int startIdx, int endIdx, int curCoreNum);
    int* getCore();
    TNode** getInvert();
    void traverse(TNode* root, std::ofstream& outFile);
    void buildInvertedIndex(TNode* currentNode);

    std::unordered_set<int> getNodesWithAttribute(const std::string& attribute) const;

    //ablation study1
    std::tuple<double, int, int, int> ablation1(const std::string& outputFileName);
    std::tuple<double, int, int, int> ablation2(const std::string& fptreeFile, const std::string& outputFileName);

    // Our method
    bool isKCore(int k, std::vector<int>& nodes);
    std::tuple<double, int, int, int> fptree(const std::string& fptreeFile, const std::string& outputFileName);
    // std::vector<std::vector<Node>> getConditionalPatternBase(int targetNodeID, int rootNode, const std::vector<std::vector<Node>>& branches);

    // typedef std::unordered_map<std::set<int>, std::set<std::string>> PatternData;

    std::vector<std::pair<std::set<int>, std::set<std::string>>> getConditionalPatternBase(int targetNodeID, int rootNode, const std::vector<std::vector<Node>>& branches, int totalMaxAttributes);

    // std::tuple<std::vector<int>, int, std::set<std::string>> buildConditionalFPTree(int targetNodeID, int rootNode, const std::vector<std::vector<Node>>& branches, int totalMaxAttributes);
    void read_attrmap(const std::string& file_path);
    void read_treeIndex(const std::string& file_path);
    std::pair<bool, int> validateCandidateSet(
        const std::vector<int>& selectedNodes,
        const std::set<std::string>& commattr,
        Graph* graph,
        int v,
        int k,
        const std::string& outputFileName
    );

    std::pair<bool, int> validateCandidateSetAblation(
        const std::vector<int>& selectedNodes,
        const std::set<std::string>& commattr,
        Graph* graph,
        int v,
        int k,
        const std::string& outputFileName
    );

    const std::unordered_set<int>& getSubtreeNodeSet(int targetID, int k) const;
    std::set<int> collectSubtreeNodeSet(TNode* node);
    void precomputeAncestorCoreness(TNode* node, std::unordered_map<int, TNode*> ancestorMap);
    void aggregateNodeSetBottomUp(TNode* node);

    void collectTreeNodes(TNode* root, std::vector<TNode *>& subtreeNodes,  std::unordered_set<TNode*>& visited);
    // Display results
    void display();

    int* decompose();
    int obtainMaxCore();
    int* obtainReverseCoreArr();

    //private:
    Graph* graph;
    int v; //node
    int k;

    int n = 0;
    int m = 0;
    int* coreReverseFang = nullptr;
    int* deg = nullptr;
    //index
    int* core = nullptr;
    
    UNode** unodeArr;
    TNode** invert;
    std::set<TNode*> restNodeSet;
    UnionFind* uf = nullptr;

    TNode* root;
    std::unordered_map<int, int> nodeToCoreIndex;
    std::unordered_map<std::string, std::unordered_set<int>> attributeToNodes;

    std::vector<int> sortedNodes;

    std::unordered_map<std::string, std::unordered_set<int>> attr_map;
    std::unordered_map<int, TNode*> nodeIndex;

};

#endif // KCORE_H
