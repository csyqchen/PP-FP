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

// TreeNode ��Ķ���
struct TreeNode {
    int coreness;  // ��ǰ�ڵ�� k-core ֵ
    std::vector<int> nodes;  // �����Ľڵ㼯��
    std::vector<TreeNode*> children;  // �ӽڵ��б�

    TreeNode(int core) : coreness(core) {}
};

class Kcore {
public:
    std::set<int> built; // ��ѯ�� conditional tree ���޷����� k-core �Ľ�㼯

    struct Node {
        int id;
        int count;
        int totalCount;
        int level;
        std::set<std::string> attributes;
        std::set<int> branches; // ��¼��������ķ�֧
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
    double baseline1(std::ofstream& output_file, size_t& success_d, size_t& community);
    double baseline2(std::ofstream& output_file, size_t& success_d, size_t& community);
    void baseline3();

    // Index
    std::vector<int> kCoreDecomposeOptimized(const std::vector<int>& component);

    void printTree(TreeNode* node, int depth = 0) const;

    // �����ĺ������������ڽ�����ӡ���ļ�
    void printTreeToFile(TreeNode* node, std::ofstream& outFile, int depth = 0) const;

    int getNodeCoreness(int node) const;

    // Tree and Inverted Index Access Functions
    TreeNode* getTreeIndexRoot() const;

    // �������µĺ��������������Ե����������������������ܽ��
    // TNode* buildAttributeIndexAndTree(const std::string& attributeOutputFilePath, const std::string& treeOutputFilePath);  
    void buildAttributeIndex();
    void handleALevel(int startIdx, int endIdx, int curCoreNum);
    int* getCore();
    TNode** getInvert();
    void traverse(TNode* root, std::ofstream& outFile);

    std::unordered_set<int> getNodesWithAttribute(const std::string& attribute) const;


    // Our method
    bool isKCore(int k, std::vector<int>& nodes);
    std::tuple<double, int, int> fptree(const std::string& fptreeFile, const std::string& outputFileName, std::vector<std::string>& result);
    std::vector<std::vector<Node>> getConditionalPatternBase(int targetNodeID, int rootNode, const std::vector<std::vector<Node>>& branches);
    std::tuple<std::vector<int>, int, std::set<std::string>> buildConditionalFPTree(int targetNodeID, int rootNode, const std::vector<std::vector<Node>>& branches, int totalMaxAttributes);
    void read_attrmap(const std::string& file_path);
    void read_treeIndex(const std::string& file_path);
    bool validateCandidateSet(
        const std::vector<int>& selectedNodes,
        const std::set<std::string>& commattr,
        Graph* graph,
        int v,
        int k,
        const std::string& outputFileName
    );
    const std::set<int>& getSubtreeNodeSet(int targetID) const;
    std::set<int> collectSubtreeNodeSet(TNode* node);

    // Display results
    void display();

    //icde decompose
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
    int* core = nullptr; // Equivalent to int[] core in Java
    //int n = -1;
    //int* coreReverseFang = nullptr; // Equivalent to int[] coreReverseFang in Java
    UNode** unodeArr; // Array of UNode objects
    TNode** invert; // Array of TNode objects
    std::set<TNode*> restNodeSet; // Equivalent to Set<TNode> in Java
    UnionFind* uf = nullptr; // Pointer to UnionFind object

    // Tree and Inverted Index Member Variables
    TNode* root;  // ��״�������ڵ�
    std::unordered_map<int, int> nodeToCoreIndex;  // ��������
    std::unordered_map<std::string, std::unordered_set<int>> attributeToNodes;  // ���Ե��ڵ�ĵ�������

    // ������Ա���������ڴ洢�����Ľڵ�
    std::vector<int> sortedNodes;  // ���������Ľڵ�˳�������Ե����ϵ�������

    std::unordered_map<std::string, std::set<int>> attr_map;
    std::unordered_map<int, TNode*> nodeIndex; // �ڵ�������

};

#endif // KCORE_H
