#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "../include/Kcore.h"
#include <queue>
#include "../include/Utils.h"
#include <memory>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <set>
#include <stack>
#include <chrono>


using namespace std;

Kcore::Kcore(Graph* p_graph, int query_node, int k_coreness) :v(query_node), k(k_coreness) {
    graph = p_graph;
}

Kcore::Kcore(Graph* p_graph, int k_coreness) : k(k_coreness) {
    graph = p_graph;
}

Kcore::Kcore(Graph* p_graph) {
    graph = p_graph;
}

void Kcore::set_query_node(int node) {
    v = node;
    TNode* targetNode = nodeIndex[node];
    if (!targetNode) {
        std::cerr << "Error: Target ID " << node << " not found in the tree." << std::endl;
        return;
    }
    // 使用栈来迭代收集子树的所有节点，避免递归
    std::set<int> result;
    std::stack<TNode*> nodeStack;
    nodeStack.push(targetNode);

    while (!nodeStack.empty()) {
        TNode* currentNode = nodeStack.top();
        nodeStack.pop();

        // 合并当前节点的 nodeSet 到结果集中
        result.insert(currentNode->getNodeSet().begin(), currentNode->getNodeSet().end());

        // 将当前节点的子节点压入栈
        for (TNode* child : currentNode->getChildList()) {
            nodeStack.push(child);
        }
    }
    targetNode->setCoreSet(result);
}

void Kcore::set_k(int query_k) {
    k = query_k;
}

void Kcore::deleteNode(int v) {
    graph->deleteNode(v);
}

std::unordered_set<int> Kcore::getNodes() {
    return graph->getNodes();
}

Graph::Attributes Kcore::attrs(int v) {
    return graph->attributes(v);
}
// int Kcore::deg(int v) {
//     return graph->deg(v);
// }

template<typename T>
std::unordered_set<T> Kcore::intersect(std::unordered_set<T> set1, std::unordered_set<T> set2) {
    std::unordered_set<T> ret_set;
    if (set1.size() < set2.size()) {
        for (const auto& elem : set1) {
            if (set2.count(elem)) {
                ret_set.insert(elem);
            }
        }
    }
    else
    {
        for (const auto& elem : set2) {
            if (set1.count(elem)) {
                ret_set.insert(elem);
            }
        }
    }
    return ret_set;
}
template<typename T>
void Kcore::intersect_filter(std::unordered_set<T>& set1, const std::unordered_set<T>& set2) {
    for (auto it = set1.begin(); it != set1.end(); ) {
        if (!set2.count(*it)) {
            it = set1.erase(it);  // 使用迭代器删除元素
        }
        else {
            ++it;
        }
    }
}

template<typename T>
void Kcore::backtrack(const std::vector<T>& arr, int k, int start, std::vector<T>& path, std::vector<std::vector<T>>& result) {
    if (path.size() == k) {
        result.push_back(path);
        return;
    }

    for (int i = start; i < arr.size(); ++i) {
        path.push_back(arr[i]);
        backtrack(arr, k, i + 1, path, result);
        path.pop_back();
    }
}

template<typename T>
std::vector<std::vector<T>> Kcore::combinations(const std::vector<T>& arr, int k) {
    std::vector<std::vector<T>> result;
    std::vector<T> path;
    backtrack(arr, k, 0, path, result);
    return result;
}

void Kcore::search_kcore(Graph* graph) {
    //printf("start searching kcore \n");
    bool changes = true;

    while (changes) {
        changes = false;
        std::queue<int> q;

        // check degree for each node
        for (auto node : graph->getNodes()) {
            if (graph->deg(node) < k) {
                q.push(node);
            }
        }
        if (!q.empty()) changes = true;
        // delete nodes
        while (!q.empty()) {
            int v = q.front();
            q.pop();

            graph->deleteNode(v);
        }
    }

}

bool Kcore::check_kcore(Graph* graph, int k, std::vector<int>& kcore_nodes) {

    if (std::find(graph->getNodes().begin(), graph->getNodes().end(), v) == graph->getNodes().end()) {
        std::cerr << "Query node " << v << " is not present in the graph." << std::endl;
        return false;
    }

    int validNodeCount = 0;  // the number of nodes with degree >= k
    std::vector<int> temp_kcore_nodes;

    for (auto node : graph->getNodes()) {
        int degree = graph->deg(node);
        std::cerr << "Node " << node << " has degree " << degree << " (k=" << k << ")" << std::endl;

        if (degree >= k) {
            temp_kcore_nodes.push_back(node);
            validNodeCount++;
        }
    }

    std::cerr << "Number of nodes with degree >= k: " << validNodeCount << std::endl;

    kcore_nodes = std::move(temp_kcore_nodes);
    return true;
}

double Kcore::baseline1(std::ofstream& output_file, size_t& success_d, size_t& community) {
    clock_t start = clock();
    //find kcore to filter graph
    search_kcore(graph);
    if (!getNodes().count(v)) {
        //printf("query node %d is not in a %d-core\n", v, k);
        output_file << "query node " << v << " is not in a " << k << "-core\n";
        clock_t end = clock();
        double duration = (end - start) / CLOCKS_PER_SEC;
        //printf("baseline1 cost:%2f \n", duration);
        output_file << "baseline1 cost: " << duration << " \n";
        return duration;
    }

    Graph::Attributes query_attributes = attrs(v);
    std::vector<std::string>query_attributes_list(query_attributes.begin(), query_attributes.end());
    std::unordered_map<std::string, std::unordered_set<int>> attr_map;
    //the key is the attribute, and the value is the set of nodes with the attribute

    //fulfill the attr_map
    printf("making the attr map......");
    for (auto item : query_attributes_list) {
        std::unordered_set<int> node_list;
        for (auto node : getNodes()) {
            if (attrs(node).count(item)) node_list.insert(node);

        }
        attr_map[item] = node_list;
    }
    printf(" attr map finished\n");
    //std::cout << "size of graphs:" << attr_map["graphs"].size() <<std::endl;
    //std::vector<std::string>query_attributes_list(query_attributes.begin(), query_attributes.end());
    //std::unique_ptr<Graph> graph_copy = std::make_unique<Graph>();

    success_d = 0;
    community = 0;

    for (size_t d = 1; d <= query_attributes_list.size(); d++) {
        printf("start searching the community with attrs: %zu\n", d);
        output_file << "start searching the community with attrs: " << d << "\n";
        bool find = false;
        std::vector<std::vector<std::string>> attr_combinations = combinations(query_attributes_list, d);
        for (auto attr_candidate : attr_combinations) {
            //merge corresponding nodes
            std::unordered_set<int> nodes_candidate = attr_map[attr_candidate[0]];

            for (size_t i = 1; i < attr_candidate.size(); i++) {
                intersect_filter(nodes_candidate, attr_map[attr_candidate[i]]);
            }
            //std::cout << "size of nodes_candidate:" << nodes_candidate.size() <<std::endl;
            if (nodes_candidate.size() < k + 1) continue;

            std::unique_ptr<Graph> graph_copy = std::make_unique<Graph>();
            try {

                graph_copy->setData(graph);
                // if (graph_copy == nullptr) {
                //     std::cerr << "graph_copy is null" << std::endl;
                //     return;
                // }
                auto nodes = graph_copy->getNodes();
                // if (nodes.empty()) {
                //     std::cerr << "No nodes in the graph" << std::endl;
                //     return;
                // }

                for (auto node : nodes) {
                    if (!nodes_candidate.count(node)) graph_copy->deleteNode(node);
                }

            }
            catch (const std::bad_alloc& e) {
                std::cerr << "Memory allocation failed: " << e.what() << std::endl;
                continue;
            }


            // for (auto node : graph_copy->getNodes()) {
            //     if (!nodes_candidate.count(node)) graph_copy->deleteNode(node);
            // }

            search_kcore(graph_copy.get());
            if (graph_copy->getNodes().count(v)) {
                printf("find the community with max attrs: %zu\n", d);
                output_file << "find the community with max attrs: " << d << "\n";
                for (const auto& str : attr_candidate) {
                    //std::cout << str << " ";
                    output_file << str << " ";
                }
                success_d = d;
                community = graph_copy->getNodes().size();  // 记录社区的节点数
                printf("\n");
                std::cout << "size of nodes:" << community << std::endl;
                output_file << "\n";
                output_file << "size of nodes: " << community << "\n";
                if (community < 20) {
                    for (const auto& node : graph_copy->getNodes()) {
                        std::cout << node << " ";
                        output_file << node << " ";
                    }
                    printf("\n");
                    output_file << "\n";
                }
                //printf("\n");
                find = true;
                break;
            }

        }
        if (!find) {
            output_file << "no " << d << " attributes community, max attrs: " << success_d << "\n";
            break;
        }
    }

    clock_t end = clock();
    //double duration = (end - start) / CLOCKS_PER_SEC;
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    printf("baseline1 cost:%2f \n", duration);
    output_file << "\nbaseline1 cost: " << duration << " \n";
    return duration;
}

double Kcore::baseline2(std::ofstream& output_file, size_t& success_d, size_t& community) {
    //void Kcore::baseline2() {

    clock_t start = clock();
    search_kcore(graph);
    if (!getNodes().count(v)) {
        //printf("query node %d is not in a %d-core\n", v, k);
        output_file << "query node " << v << " is not in a " << k << "-core\n";
        clock_t end = clock();
        double duration = (end - start) / CLOCKS_PER_SEC;
        //printf("baseline2 cost:%2f \n", duration);
        output_file << "baseline2 cost: " << duration << " \n";
        return duration;
    }

    Graph::Attributes query_attributes = attrs(v);
    std::vector<std::string> query_attributes_list(query_attributes.begin(), query_attributes.end());
    std::unordered_map<std::string, std::unordered_set<int>> attr_map;

    printf("making the attr map......");
    for (auto& item : query_attributes_list) {
        std::unordered_set<int> node_list;
        for (auto node : getNodes()) {
            if (attrs(node).count(item)) node_list.insert(node);
        }
        attr_map[item] = node_list;
    }
    printf(" attr map finished\n");

    size_t min_d = 1;
    size_t max_d = query_attributes_list.size();
    success_d = 0;
    community = 0;
    std::vector<std::vector<std::string>> combinations_candidate;

    while (min_d <= max_d) {
        size_t mid_d = (min_d + max_d) / 2;
        //printf("start searching the community with attrs: %d\n", mid_d);
        output_file << "start searching the community with attrs: " << mid_d << "\n";
        bool find = false;
        std::vector<std::vector<std::string>> attr_combinations = combinations(query_attributes_list, mid_d);

        for (auto& attr_candidate : attr_combinations) {
            std::unordered_set<int> nodes_candidate = attr_map[attr_candidate[0]];

            for (size_t i = 1; i < attr_candidate.size(); i++) {
                intersect_filter(nodes_candidate, attr_map[attr_candidate[i]]);
            }

            if (nodes_candidate.size() < k + 1) continue;
            std::unique_ptr<Graph> graph_copy = std::make_unique<Graph>();

            try {
                graph_copy->setData(graph);
                //new
                auto nodes = graph_copy->getNodes();

                // for (auto node : graph_copy->getNodes()) {
                //     if (!nodes_candidate.count(node)) graph_copy->deleteNode(node);
                // }
                for (auto node : nodes) {
                    if (!nodes_candidate.count(node)) graph_copy->deleteNode(node);
                }
            }
            catch (const std::bad_alloc& e) {
                std::cerr << "Memory allocation failed: " << e.what() << std::endl;
                continue;
            }

            search_kcore(graph_copy.get());
            if (graph_copy->getNodes().count(v)) {
                // combinations_candidate.push_back(attr_candidate); // Storing all successful combinations
                //printf("find the community with max attrs: %zu\n", mid_d);
                output_file << "find the community with max attrs: " << mid_d << "\n";
                for (const auto& str : attr_candidate) {
                    //std::cout << str << " ";
                    output_file << str << " ";
                }
                //printf("\n");
                output_file << "\n";
                success_d = mid_d;
                community = graph_copy->getNodes().size();  // 记录满足条件的社区节点数

                //std::cout << "size of nodes:" << graph_copy->getNodes().size() << std::endl;
                output_file << "size of nodes: " << community << "\n";
                if (community < 20) {
                    for (const auto& node : graph_copy->getNodes()) {
                        //std::cout << node << " ";
                        output_file << node << " ";
                    }
                }
                printf("\n");
                find = true;
                break;
            }
        }
        if (find) {
            min_d = mid_d + 1;
            success_d = mid_d;
        }
        else if (!find && mid_d > success_d + 1) {
            max_d = mid_d - 1;
        }
        else if (!find && mid_d <= success_d + 1) {
            //printf("no %zu attributes community, max attrs: %zu\n", mid_d, success_d);
            output_file << "no " << mid_d << " attributes community, max attrs: " << success_d << "\n";
            break;
        }
    }
    clock_t end = clock();
    //double duration = (end - start) / CLOCKS_PER_SEC;
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    //printf("baseline2 cost:%2f \n", duration);
    output_file << "\nbaseline2 cost: " << duration << " \n";
    return duration;
    // Output or use combinations_candidate as needed
    // for (const auto& comb : combinations_candidate) {
    //     for (const auto& attr : comb) {
    //         std::cout << attr << " ";
    //     }
    //     std::cout << std::endl;
    // }
}

std::tuple<double, int, int> Kcore::fptree(const std::string& fptreeFile, const std::string& outputFileName, std::vector<std::string>& result) {

    clock_t start = clock();

    std::ifstream fpFile(fptreeFile);
    if (!fpFile.is_open()) {
        std::cerr << "Cannot open the file: " << fptreeFile << std::endl;
        return std::make_tuple(-1.0, 0, 0);
    }

    std::vector<std::vector<Node>> branches; //store different branches
    std::map<int, std::vector<Node>> levelToBranchMap; //store each level and its corresponding branches
    int rootNode = v;
    std::cerr << "rootNode: " << v << std::endl;


    std::string line;
    std::vector<Node> currentBranch;
    int previousLevel = -1;

    int branchCounter = 0;

    while (std::getline(fpFile, line)) {
        int indentLevel = 0;
        while (indentLevel < line.length() && line[indentLevel] == ' ') {
            indentLevel++;
        }
        line = line.substr(indentLevel);

        std::istringstream iss(line);
        int nodeID, count, totalCount;
        char tempChar;

        iss >> nodeID;
        iss.ignore(20, '(');
        iss >> count >> tempChar >> totalCount;
        iss.ignore(20, '{');

        std::set<std::string> attributes;
        std::string attribute;

        while (std::getline(iss, attribute, ',')) {
            attribute.erase(0, attribute.find_first_not_of(" \t"));
            attribute.erase(attribute.find_last_not_of(" \t") + 1);

            if (!attribute.empty() && attribute.back() == '}') {
                attribute.pop_back();
            }
            if (!attribute.empty()) {
                attributes.insert(attribute);
            }
        }

        std::cerr << "FPNodeID: " << nodeID << ", Count: " << count << ", TotalCount: " << totalCount << "\n";
        std::cerr << "Attributes: ";
        for (const auto& attr : attributes) {
            std::cerr << "[" << attr << "] ";
        }
        std::cerr << "\n";


        int currentLevel = indentLevel / 2;

        Node currentNode = { nodeID, count, totalCount, currentLevel, attributes };

        // if (currentLevel == 1) {
        if (currentLevel <= previousLevel) {
            // new branch 
            if (!currentBranch.empty()) {
                branches.push_back(currentBranch);
                currentBranch.clear();
                branchCounter++;
            }
        }

        // update currentBranch to include prefixed paths
        if (currentLevel > 0 && levelToBranchMap.find(currentLevel - 1) != levelToBranchMap.end()) {
            currentBranch = levelToBranchMap[currentLevel - 1];
            // Inherit the parent path
        }
        currentNode.branches.insert(branchCounter);
        currentBranch.push_back(currentNode);
        levelToBranchMap[currentLevel] = currentBranch;

        previousLevel = currentLevel;
    }

    if (!currentBranch.empty()) {
        branches.push_back(currentBranch);
    }

    fpFile.close();

    bool foundDegk = false;
    int maxAttributes = 0;
    int totalMaxAttributes = branches[1][1].totalCount; //maximum possible attribute number
    std::cout << "totalMaxAttributes" << totalMaxAttributes;
    std::vector<int> selectedNodes;
    std::vector<int> tempSelectedNodes;
    int tempMaxAttributes = 0;
    std::set<std::string> commattr, tempcommattr; //store the current public attributes

    for (const auto& branch : branches) {
        std::vector<int> candidateNodes = { rootNode };
        auto it = branch.cbegin();
        while (it != branch.cend()) {
            const auto& node = *it;

            if (node.count == totalMaxAttributes) {
                candidateNodes.push_back(node.id);
                auto next_it = std::next(it);
                bool nextNodeCondition = (next_it == branch.cend() || next_it->count < totalMaxAttributes);

                if (candidateNodes.size() >= k && nextNodeCondition) {
                    commattr = std::set<std::string>(node.attributes.begin(), node.attributes.end());
                    if (validateCandidateSet(candidateNodes, commattr, graph, v, k, outputFileName)) {
                        foundDegk = true;
                        selectedNodes = candidateNodes;
                        maxAttributes = totalMaxAttributes;
                        std::cerr << "Validate success...\n";
                        break;
                    }
                } //find the current deg_q >= k and the next node's count is lowered
            }
            ++it;
        }
        if (foundDegk) {
            break;
        }
    }

    while (!foundDegk && totalMaxAttributes > 1) {
        totalMaxAttributes--;
        std::cout << "totalMaxAttributes" << totalMaxAttributes;
        std::cerr << "Searching with totalMaxAttributes = " << totalMaxAttributes << std::endl;

        for (const auto& branch : branches) {
            if (foundDegk) break;
            std::vector<int> candidateNodes = { rootNode };
            std::set<std::string> branchCommAttr; //store the common attributes of the current branch

            auto it = branch.cbegin();
            while (it != branch.cend()) {
                const auto& node = *it;

                if (node.count >= totalMaxAttributes) {
                    candidateNodes.push_back(node.id);
                    // branchCommAttr.insert(node.attributes.begin(), node.attributes.end());

                    auto next_it = std::next(it);
                    bool nextNodeCondition = (next_it == branch.cend() || next_it->count < totalMaxAttributes);

                    if (candidateNodes.size() >= k && nextNodeCondition) {
                        branchCommAttr.insert(node.attributes.begin(), node.attributes.end());

                        std::cerr << "Validating candidateNodes directly without conditional tree...\n";

                        // validate candidate solution
                        if (validateCandidateSet(candidateNodes, branchCommAttr, graph, v, k, outputFileName)) {
                            foundDegk = true;
                            selectedNodes = candidateNodes;
                            commattr = branchCommAttr;
                            maxAttributes = totalMaxAttributes;
                            break;
                        }
                    }
                }

                //satisfy the conditions for constructing a conditional tree
                if (node.count < totalMaxAttributes && node.totalCount >= totalMaxAttributes && built.find(node.id) == built.end()) {
                    std::cerr << "Building conditional tree for Node ID: " << node.id << std::endl;

                    auto [tempNodes, threshold, tempCommAttr] = buildConditionalFPTree(node.id, rootNode, branches, totalMaxAttributes);

                    // Check if the size of candidateNodes satisfies k
                    if (tempNodes.size() >= k) {
                        if (threshold >= totalMaxAttributes) {
                            selectedNodes = tempNodes;
                            maxAttributes = threshold;
                            commattr = tempCommAttr;

                            // validate candidate solution
                            if (validateCandidateSet(selectedNodes, commattr, graph, v, k, outputFileName)) {
                                foundDegk = true;
                                break;
                            }
                        } else if (threshold >= tempMaxAttributes) {
                            // preserve better alternative answers
                            tempSelectedNodes = tempNodes;
                            tempMaxAttributes = threshold;
                            tempcommattr = tempCommAttr;
                            std::cerr << "Saving temporary candidate with threshold = " << threshold << std::endl;
                        }
                    } else {
                        std::cerr << "Skipping CandidateNodes with size < k: " << tempNodes.size() << std::endl;
                    }
                }

                if (foundDegk) break;
                ++it;
            }
            if (foundDegk) break;
        }

        if (!foundDegk && tempMaxAttributes >= totalMaxAttributes - 1 && tempMaxAttributes >= 1) {
            selectedNodes = tempSelectedNodes;
            commattr = tempcommattr;
            std::cerr << "Validating temporary candidate with threshold = " << tempMaxAttributes << std::endl;

            std::cerr << "Selected Nodes: ";
            for (const auto& node : selectedNodes) {
                std::cerr << node << " ";
            }
            std::cerr << std::endl;

            std::cerr << "Common Attributes: ";
            for (const auto& attr : commattr) {
                std::cerr << attr << " ";
            }
            std::cerr << std::endl;

            // validate candidate solution
            if (validateCandidateSet(selectedNodes, commattr, graph, v, k, outputFileName)) {
                foundDegk = true;
            } else {
                selectedNodes.clear();
                commattr.clear();
                tempSelectedNodes.clear();
                tempcommattr.clear();
                tempMaxAttributes = 0;
                std::cerr << "Validation failed. Cleared selectedNodes and commattr.\n";
            }
        }
        built.clear();
    }

    clock_t end = clock();
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "\nfptree cost: " << duration << " \n";
    std::ostringstream oss;
    oss<<"query node:"<<v<<"; community: ";
    for (size_t i = 0; i < selectedNodes.size(); ++i) {
        oss << selectedNodes[i]; 
        if (i != selectedNodes.size() - 1) {
            oss << ", ";
        }
    }

    if (selectedNodes.empty()) {
        result.push_back("Community: [No nodes selected]");
    } else {
        result.push_back(" [" + oss.str() + "]");
    }

    return std::make_tuple(duration, maxAttributes, selectedNodes.size());
}

bool Kcore::validateCandidateSet(const std::vector<int>& selectedNodes, const std::set<std::string>& commattr, Graph* graph, int v, int k, const std::string& outputFileName) {
    // update attrSet
    std::set<int> attrSet;

    // output commattr
    std::cerr << "CommAttr: ";
    for (const auto& attr : commattr) {
        std::cerr << attr << " ";
    }
    std::cerr << std::endl;

    std::set<int> tempIntersection; // store the current intersection result
    bool firstAttribute = true;

    for (const auto& attr : commattr) {
        auto it = attr_map.find(attr);
        if (it != attr_map.end()) {
            const std::set<int>& attrNodes = it->second;

            if (firstAttribute) {
                tempIntersection = attrNodes;
                firstAttribute = false;
            } else {
                std::set<int> newIntersection;
                std::set_intersection(
                    tempIntersection.begin(), tempIntersection.end(),
                    attrNodes.begin(), attrNodes.end(),
                    std::inserter(newIntersection, newIntersection.begin())
                );

                tempIntersection = std::move(newIntersection);
            }
        } else {
            std::cerr << "Attribute " << attr << " not found in attr_map.\n";
            tempIntersection.clear();
            break;
        }
    }

    attrSet = std::move(tempIntersection);

    attrSet.insert(selectedNodes.begin(), selectedNodes.end());
    
    std::cout << "attrSet: ";
    for (const auto& node : attrSet) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    // const std::set<int>& coreSet = attrSet;    
    
    const std::set<int>& coreSet = getSubtreeNodeSet(v);    

    std::set<int> candidateSet;
    
    // auto startTime = std::chrono::high_resolution_clock::now();

    std::set_intersection(
        attrSet.begin(), attrSet.end(),
        coreSet.begin(), coreSet.end(),
        std::inserter(candidateSet, candidateSet.begin())
    );

    // auto endTime = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsedTime = endTime - startTime;

    std::cout << "candidateSet: ";
    for (const auto& node : candidateSet) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    // std::cerr << "Time taken by set_intersection: " << elapsedTime.count() << " seconds" << std::endl;

    if (candidateSet.size() >= k + 1) {
        clock_t startTinyGraphTime = clock();

        std::vector<int> sortedCandidates(candidateSet.begin(), candidateSet.end());
        std::sort(sortedCandidates.begin(), sortedCandidates.end());

        Graph* tiny_graph = new Graph();
        for (size_t i = 0; i < sortedCandidates.size(); ++i) {
            for (size_t j = i + 1; j < sortedCandidates.size(); ++j) {
                if (graph->hasEdge(sortedCandidates[i], sortedCandidates[j])) {
                    tiny_graph->addEdge(sortedCandidates[i], sortedCandidates[j]);
                }
            }
        }

        Kcore* kcore = new Kcore(tiny_graph, v, k);
        kcore->search_kcore(tiny_graph);

        if (tiny_graph->hasNode(v)) {
            std::ofstream outputFile(outputFileName);
            std::cout << "Tiny_graph: true" << std::endl;
            clock_t endTinyGraphTime = clock();
            double tinyGraphElapsedTime = static_cast<double>(endTinyGraphTime - startTinyGraphTime) / CLOCKS_PER_SEC;
            std::cerr << "Tiny_graph construction time: " << tinyGraphElapsedTime << " seconds" << std::endl;
            outputFile << "Community Node ID: ";
            for (int node : tiny_graph->getNodes()) {
                outputFile << node << " ";
            }
            outputFile << "\n";
            outputFile << "Community Size: " << tiny_graph->getNodes().size() << "\n";
            outputFile << "Max attribute: " << commattr.size() << "\n";
            outputFile.close();
            delete tiny_graph;
            delete kcore;
            return true;
        } else {
            delete tiny_graph;
            delete kcore;
        }
    } else {
        std::cerr << "Skipping candidateSet with size < k+1: " << candidateSet.size() << std::endl;
    }

    return false;
}

std::vector<std::vector<Kcore::Node>> Kcore::getConditionalPatternBase(int targetNodeID, int rootNode, const std::vector<std::vector<Kcore::Node>>& branches) { 
    std::vector<std::vector<Kcore::Node>> conditionalPatternBase;
    std::set<std::vector<int>> uniquePaths;


    std::cout << "rootNode in getConditionalPatternBase: " << rootNode << std::endl;

    for (const auto& branch : branches) {
        std::vector<Node> prefixPath;
        bool targetFound = false;

        for (auto it = branch.rbegin(); it != branch.rend(); ++it) {
            if (it->id == targetNodeID) {
                targetFound = true;
            }
            if (targetFound) {
                prefixPath.push_back(*it);
            }
        }

        if (!prefixPath.empty()) {
            if (prefixPath.back().id == 0) {
                prefixPath.back().id = rootNode;
            }
            std::reverse(prefixPath.begin(), prefixPath.end());
            // conditionalPatternBase.push_back(prefixPath);
            // new
            std::vector<int> pathIds;
            for (const auto& node : prefixPath) {
                pathIds.push_back(node.id);
            }
            if (uniquePaths.find(pathIds) == uniquePaths.end()) {
                uniquePaths.insert(pathIds);
                conditionalPatternBase.push_back(prefixPath);
            }
        }
    }
    return conditionalPatternBase;
}

std::tuple<std::vector<int>, int, std::set<std::string>> Kcore::buildConditionalFPTree(int targetNodeID, int rootNode, const std::vector<std::vector<Kcore::Node>>& branches, int totalMaxAttributes) {

    std::vector<std::vector<Kcore::Node>> conditionalPatternBase = getConditionalPatternBase(targetNodeID, rootNode, branches);

    std::cout << "rootNode in buildConditionalFPTree: " << rootNode << std::endl;

    std::map<int, int> itemFrequency;
    std::set<std::string> commattr; //store the concatenation of the last node attribute
    commattr.clear();

    //select branches satisfying nodes >= k 
    std::vector<std::vector<Kcore::Node>> validBranches;
    for (const auto& pattern : conditionalPatternBase) {
        if (pattern.size() >= k) {
            validBranches.push_back(pattern);
        }
    }

    //no valid branch found
    if (validBranches.empty()) {
        std::cout << "No valid branches with node count >= " << k << " found for targetNodeID " << targetNodeID << "\n";
        return { {}, 0, {} };
    }

    for (const auto& pattern : validBranches) {
        int lastNodeCount = pattern.back().count;
        const auto& lastNodeAttributes = pattern.back().attributes;

        std::cerr << "Last Node Attributes: ";
        for (const auto& attr : lastNodeAttributes) {
            std::cerr << "[" << attr << "] ";
        }
        std::cerr << "\n";

        commattr.insert(lastNodeAttributes.begin(), lastNodeAttributes.end());

        std::cerr << "CommAttr After Insert: ";
        for (const auto& attr : commattr) {
            std::cerr << "[" << attr << "] ";
        }
        std::cerr << "\n";

        for (const auto& node : pattern) {
            std::cout << "Node ID: " << node.id << " Count: " << node.count << "\n";
            itemFrequency[node.id] += lastNodeCount;
        }
    }

    std::vector<int> frequencies;
    for (const auto& item : itemFrequency) {
        frequencies.push_back(item.second);
    }
    std::sort(frequencies.begin(), frequencies.end(), std::greater<int>());

    int threshold = 0;
    if (frequencies.size() >= k + 1) {
        threshold = frequencies[k];
    } else {
        threshold = 0;
    }

    if (threshold > 0) {
        built.insert(targetNodeID);
        std::cout << "Nodes in conditional FP-Tree for node " << targetNodeID << " with frequency >= " << threshold << ":\n";

        std::vector<int> candidateNodes;
        for (const auto& item : itemFrequency) {
            if (item.second >= threshold) {
                std::cout << "Node ID: " << item.first << " Frequency: " << item.second << "\n";
                candidateNodes.push_back(item.first);
            }
        }

        return { candidateNodes, threshold, commattr };
    } else {
        built.insert(targetNodeID);
        std::cout << "Current built set contains the following nodes: " << std::endl;
        if (built.empty()) {
            std::cout << "Built set is currently empty." << std::endl;
        } else {
            for (const int& node : built) {
                std::cout << "Node ID: " << node << std::endl;
            }
        }
        return { {}, 0, {} };
    }
}

bool Kcore::isKCore(int k, std::vector<int>& nodes) { 
    std::unique_ptr<Graph> subgraph = std::make_unique<Graph>();

    bool queryNodeInSubgraph = false;

    for (size_t i = 0; i < nodes.size(); ++i) {
        int node = nodes[i];
        subgraph->addNode(node);

        if (node == v) {
            queryNodeInSubgraph = true;
        }

        std::vector<int> neighbors = graph->getNeighborIDs(node);

        for (size_t j = i + 1; j < nodes.size(); ++j) {
            int nextNode = nodes[j];

            if (std::find(neighbors.begin(), neighbors.end(), nextNode) != neighbors.end()) {
                subgraph->addEdge(node, nextNode);
                std::cerr << "nextNode " << nextNode << " is a neighbor of " << node << std::endl;
            }
            else {
                std::cerr << "Warning: nextNode " << nextNode << " is not a neighbor of " << node << std::endl;
            }
        }
    }

    search_kcore(subgraph.get());
    return check_kcore(subgraph.get(), k, nodes);
}

std::unordered_map<std::string, std::unordered_set<int>> attributeToNodes;

void Kcore::buildAttributeIndex() {
    for (const auto& node : graph->getNodes()) {
        const auto& attributes = graph->attributes(node);
        for (const auto& attribute : attributes) {
            attributeToNodes[attribute].insert(node);
        }
    }
    printf("Attribute index built.\n");
}

std::unordered_set<int> Kcore::getNodesWithAttribute(const std::string& attribute) const {
    auto it = attributeToNodes.find(attribute);
    if (it != attributeToNodes.end()) {
        return it->second;
    }
    else {
        return {};
    }
}

// read attr map
void Kcore::read_attrmap(const std::string& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open the file: " << file_path << std::endl;
        return;
    }

    // attr_map.clear();
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        size_t colonPos = line.find(':');
        if (colonPos == std::string::npos) continue;

        std::string key = line.substr(0, colonPos);
        std::string values_str = line.substr(colonPos + 1);

        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        values_str.erase(0, values_str.find_first_not_of(" \t"));
        values_str.erase(values_str.find_last_not_of(" \t") + 1);

        std::set<int> values;
        std::istringstream values_iss(values_str);
        std::string value;
        while (values_iss >> value) {
            try {
                values.insert(std::stoi(value));
            } catch (const std::exception& e) {
                std::cerr << "Warning: Invalid value '" << value << "' in file." << std::endl;
                continue;
            }
        }

        if (key == "retrieval" || key == "shape") {
            std::cerr << "Debug - Key: " << key << std::endl;
            std::cerr << "Debug - Raw Line: " << line << std::endl;
            std::cerr << "Debug - Values: ";
            for (const auto& val : values) {
                std::cerr << val << " ";
            }
            std::cerr << std::endl;
        }

        attr_map[key] = values;
    }

    file.close();
}

void Kcore::read_treeIndex(const std::string& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << file_path << std::endl;
        return;
    }

    std::string line;
    std::stack<std::pair<int, TNode*>> nodeStack;
    root = nullptr;

    const int INDENT_UNIT = 2;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        int indentLevel = 0;
        while (indentLevel < line.size() && line[indentLevel] == ' ') {
            indentLevel++;
        }
        indentLevel /= INDENT_UNIT;

        line = line.substr(indentLevel * INDENT_UNIT);

        std::istringstream iss(line);
        std::string key;
        int core = -1;
        int size = 0;
        std::set<int> nodeSet;

        if (std::getline(iss, key, '=') && key == "k") {
            iss >> core;
            iss.ignore(20, ' ');

            if (std::getline(iss, key, '=') && key == "size") {
                iss >> size;
                iss.ignore(20, ':'); // 跳过 " nodes:"

                // 提取节点集合
                int node;
                while (iss >> node) {
                    nodeSet.insert(node);          // 将节点加入集合
                    nodeIndex[node] = nullptr;     // 初始化为 nullptr，稍后更新
                }
            }
        }

        // 创建当前节点
        TNode* currentNode = new TNode(core);
        currentNode->setNodeSet(nodeSet);

        // 更新索引表
        for (int id : nodeSet) {
            nodeIndex[id] = currentNode; // 每个 nodeID 映射到当前 TNode
        }

        // 如果是根节点
        if (nodeStack.empty()) {
            root = currentNode; // 设置根节点
            nodeStack.push({ indentLevel, currentNode });
        }
        else {
            // 检查当前节点的缩进级别，找到正确的父节点
            while (!nodeStack.empty() && nodeStack.top().first >= indentLevel) {
                nodeStack.pop();
            }

            // 添加当前节点为父节点的子节点
            if (!nodeStack.empty()) {
                nodeStack.top().second->getChildList().push_back(currentNode);
            }

            // 将当前节点加入栈
            nodeStack.push({ indentLevel, currentNode });
        }
    }

    file.close();

}

const std::set<int>& Kcore::getSubtreeNodeSet(int targetID) const {
    TNode* targetNode = nodeIndex.at(targetID);
    const std::set<int>& ret= targetNode->getCoreSet();
    return ret;
}

void Kcore::display() {
    // display函数的具体实现
}

