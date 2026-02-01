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
#include <any>
#include <set>
#include <stack>
#include <chrono>
// #include <omp.h>


using namespace std;

//std::unordered_map<std::string, std::set<int>> attr_map;

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

    std::unordered_set<int> result;
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
    // targetNode->setCoreSet(result);
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

    int validNodeCount = 0;  // 用于统计度数大于或等于 k 的节点个数
    std::vector<int> temp_kcore_nodes;

    for (auto node : graph->getNodes()) {
        int degree = graph->deg(node);
        std::cerr << "Node " << node << " has degree " << degree << " (k=" << k << ")" << std::endl;

        if (degree >= k) {
            temp_kcore_nodes.push_back(node);
            validNodeCount++;
        }
    }

    // 输出有效节点的数量
    std::cerr << "Number of nodes with degree >= k: " << validNodeCount << std::endl;

    // 如果有效节点数量小于 k，则不满足 k-core
    // if (validNodeCount < k + 1) {
    //     return false;
    // }
    //mi
    kcore_nodes = std::move(temp_kcore_nodes);
    return true;
}

std::tuple<double, int, int, int> Kcore::fptree(const std::string& fptreeFile, const std::string& outputFileName) {

    // std::map<int, std::pair<int,int>> corenessMap; //存储coreness和step

    clock_t start = clock();
    
    int iteration = 0;

    std::ifstream fpFile(fptreeFile);
    if (!fpFile.is_open()) {
        std::cerr << "Cannot open the file: " << fptreeFile << std::endl;
        return std::make_tuple(-1.0, 0, 0, 0);
    }

    std::vector<std::vector<Node>> branches; // 存储不同的支链
    std::map<int, std::vector<Node>> levelToBranchMap; //存储每个层次对应的分支
    int rootNode = v; //根结点为query node
    //mi std::cerr << "rootNode: " << v << std::endl;


    std::string line;
    std::vector<Node> currentBranch;
    int previousLevel = -1; //用于跟踪上一个结点的level

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
            // attribute.erase(0, attribute.find_first_not_of(" \t"));

            attribute.erase(0, attribute.find_first_not_of(" \t\r\n"));

            attribute.erase(attribute.find_last_not_of(" \t\r\n") + 1);

            // attribute.erase(attribute.find_last_not_of(" \t") + 1);

            if (!attribute.empty() && attribute.back() == '}') {
                attribute.pop_back();
            }

            attribute.erase(attribute.find_last_not_of(" \t\r\n") + 1);

            if (!attribute.empty()) {
                attributes.insert(attribute);
            }
        }

        // std::cerr << "FPNodeID: " << nodeID << ", Count: " << count << ", TotalCount: " << totalCount << "\n";
        // std::cerr << "Attributes: ";
        // for (const auto& attr : attributes) {
        //     std::cerr << "[" << attr << "] ";
        // }
        // std::cerr << "\n";


        int currentLevel = indentLevel / 2;

        Node currentNode = { nodeID, count, totalCount, currentLevel, attributes };

        // if (currentLevel == 1) {
        if (currentLevel <= previousLevel) {
            // 新的支链
            if (!currentBranch.empty()) {
                branches.push_back(currentBranch);
                currentBranch.clear();
                branchCounter++;
            }
        }

        // 更新currentBranch，使其包含前缀路径
        if (currentLevel > 0 && levelToBranchMap.find(currentLevel - 1) != levelToBranchMap.end()) {
            currentBranch = levelToBranchMap[currentLevel - 1];
            // 继承父路径
        }
        currentNode.branches.insert(branchCounter); //将当前分支编号插入结点
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
    int totalMaxAttributes = branches[0][1].totalCount + 1; //最大的可能attr
    //mi std::cout << "totalMaxAttributes:" << totalMaxAttributes << endl;
    std::vector<int> selectedNodes;
    std::vector<int> tempSelectedNodes; //候选社区结点集合
    int tempMaxAttributes = 0; //候选社区最大属性值
    std::set<std::string> commattr, tempcommattr; //存储当前公共属性
    int community_size = 0, temp_size = 0;

    while (!foundDegk && totalMaxAttributes > 1) {

        // totalMaxAttributes--;
        totalMaxAttributes--;
        built.clear();

        std::cout << "totalMaxAttributes:" << totalMaxAttributes << endl;
        std::cerr << "Searching with totalMaxAttributes = " << totalMaxAttributes << std::endl;

        // 遍历分支
        for (const auto& branch : branches) {

            if (foundDegk) break; // 如果已经找到，跳出循环
            std::vector<int> candidateNodes = {};
            std::set<std::string> branchCommAttr; // 存储当前分支的公共属性

            auto it = branch.cbegin();
            while (it != branch.cend()) {
                const auto& node = *it;

                // std::cerr << "Visiting node ID: " << node.id
                // << " (count: " << node.count << ", totalCount: " << node.totalCount << ")" << std::endl;

                if (node.totalCount < totalMaxAttributes) {
                    break;
                } //退出当前支链

                // 直接验证逻辑
                if (node.count >= totalMaxAttributes) {
                    candidateNodes.push_back(node.id);
                    // branchCommAttr.insert(node.attributes.begin(), node.attributes.end());
                    branchCommAttr = std::set<std::string>(node.attributes.begin(), node.attributes.end());


                    auto next_it = std::next(it);
                    bool nextNodeCondition = (next_it == branch.cend() || next_it->count < totalMaxAttributes);

                    if (candidateNodes.size() >= (k + 1) && nextNodeCondition) {
                        //mi std::cerr << "candidate node size: " << candidateNodes.size() << "\n";
                        //mi std::cout << "Candidate Nodes: ";
                        for (int id : candidateNodes) {
                            std::cout << id << " ";
                        }
                        std::cout << std::endl;
                        //mi std::cerr << "Validating candidateNodes directly without conditional tree...\n";

                        auto validateResultA = validateCandidateSet(candidateNodes, branchCommAttr, graph, v, k, outputFileName);
                        iteration ++;

                        // 直接验证，找到了最优解
                        if (validateResultA.first) {
                            community_size = validateResultA.second;
                            foundDegk = true;
                            selectedNodes = candidateNodes;
                            commattr = branchCommAttr;
                            maxAttributes = commattr.size();
                            break; // 退出内层循环
                        }  //应该改为返回验证后的结点的set
                    }
                }

                // 构造条件树逻辑
                if (node.count < totalMaxAttributes && node.totalCount >= totalMaxAttributes && built.find(node.id) == built.end() && foundDegk != true) {
                    //mi std::cerr << "Checking node ID: " << node.id << " with node.count = " << node.count << " and node.totalCount = " << node.totalCount << std::endl;
                    //mi std::cerr << "Building conditional tree for Node ID: " << node.id << std::endl;

                    // 构造条件树
                    // auto [tempNodes, threshold, tempCommAttr] = buildConditionalFPTree(node.id, rootNode, branches, totalMaxAttributes);
                    // auto [tNodes, tempCommAttr] = getConditionalPatternBase(node.id, rootNode, branches, totalMaxAttributes);

                    //返回条件树中所有可能的组合
                    std::vector<std::pair<std::set<int>, std::set<std::string>>> validCandidates = getConditionalPatternBase(node.id, rootNode, branches, totalMaxAttributes);

                    //mi for (const auto& [nodes, attrs] : validCandidates) {
                    //     std::cerr << "{ ";
                    //     for (int n : nodes) std::cerr << n << " ";
                    //     std::cerr << "} with attrs: { ";
                    //     for (const auto& a : attrs) std::cerr << a << " ";
                    //     std::cerr << "}" << std::endl;
                    // }

                    for (const auto& candidate : validCandidates) {

                        std::vector<int> tempNodes(candidate.first.begin(), candidate.first.end());
                        const std::set<std::string>& tempCommAttr = candidate.second;

                        auto validateResult = validateCandidateSet(tempNodes, tempCommAttr, graph, v, k, outputFileName);
                        iteration ++;

                        if (validateResult.first) {
                            selectedNodes = tempNodes;
                            maxAttributes = tempCommAttr.size();
                            commattr = tempCommAttr;
                            community_size = validateResult.second;
                            foundDegk = true;
                            break;
                        } else {
                            //mi std::cerr << "Validation failed for set size =" << tempNodes.size() << "\n";
                        }
                    }

                    if (!foundDegk) {
                        built.insert(node.id);
                        //mi std::cerr << "Build conditional tree failed." << "\n";
                    }
                    if (foundDegk) break; // 外层循环再次检查
                }
                if (!foundDegk) ++it;
                if (foundDegk) break; // 外层循环再次检查
                // totalMaxAttributes--; //属性-1
            }
        }
    }

    clock_t end = clock();
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "\nfptree cost: " << duration << " \n";
    return std::make_tuple(duration, maxAttributes, community_size, iteration);
}

std::pair<bool,int> Kcore::validateCandidateSet(const std::vector<int>& selectedNodes, const std::set<std::string>& commattr, Graph* graph, int v, int k, const std::string& outputFileName) {

    auto startTime = std::chrono::high_resolution_clock::now();

    std::vector<TNode*> subtreeNodes;
    std::unordered_set<TNode*> visited;

    collectTreeNodes(nodeIndex[v], subtreeNodes, visited);

    std::unordered_set<int> candidateSet;  //最终的候选结点集合

    // std::cout << "commattr: ";
    // for (const auto& attr : commattr) {
    //     std::cout << attr << " ";
    // }
    // std::cout << std::endl;

    for (TNode* treenode : subtreeNodes){
        const auto& invertedIndex = treenode->getInvertedIndex();  //treenode对应的map

        std::vector<std::pair<std::string, std::unordered_set<int>>> setToIntersect;

        std::cout << "[COMMATTR] ";

        for (const std::string& attr: commattr) {
            auto it = invertedIndex.find(attr);
            if (it != invertedIndex.end()) {
                std::unordered_set<int> tempIntersect(it->second.begin(), it->second.end());
                setToIntersect.emplace_back(attr, tempIntersect);
            }
            std::cout << attr << " ";
        }
        std::cout << std::endl;

        if (setToIntersect.size() < commattr.size()) {  //属性匹配不完整则跳过
            continue;
        }

        if (!setToIntersect.empty()) {

            std::sort(setToIntersect.begin(), setToIntersect.end(), [](const auto& a, const auto& b) {
                return a.second.size() < b.second.size();
            });

            //根据集合大小从小到大排序

            std::unordered_set<int> intersectSet = setToIntersect[0].second;

            std::cout << "[INIT] attr=" << setToIntersect[0].first
                      << " intersectSet size=" << intersectSet.size()
                      << " nodes: ";

            std::cout << "[TREENODE] ptr=" << treenode
                      << " nodeSet size=" << treenode->getNodeSet().size()
                      << std::endl;


            int shown = 0;
            for (int x : intersectSet) {
                std::cout << x << " ";
                if (++shown >= 5) break;   // 只看前 5 个
            }
            if (intersectSet.size() > 5) std::cout << "...";
            std::cout << std::endl;

            for (int i = 1; i < setToIntersect.size(); i++) {
                std::unordered_set<int> tempSet;
                for (int node : intersectSet) {
                    if (setToIntersect[i].second.find(node) != setToIntersect[i].second.end()) {
                        tempSet.insert(node);
                    }
                }

                std::cout << "  [TEMP] tempSet size=" << tempSet.size()
                          << " nodes: ";

                int shown = 0;
                for (int x : tempSet) {
                    std::cout << x << " ";
                    if (++shown >= 5) break;
                }

                if (tempSet.size() > 5) {
                    std::cout << "...";
                }

                std::cout << std::endl;

                if (tempSet.empty()) {
                    intersectSet.clear();
                    break;
                } //交集如果已为空，提前终止

                intersectSet = tempSet;
            }

            // for (const auto& pair : setToIntersect) {
            //     std::cout << "smallSetzhongtu:";
            //     // 输出每个 unordered_set<int>
            //     for (int node : pair.second) {
            //         std::cout << node << " ";
            //     }
            //     std::cout << std::endl;
            // }
            //
            // std::cout << "smallSet1:";
            // for (const int& node : intersectSet) {
            //     std::cout << node << " ";
            // }
            // std::cout << std::endl;

            // if (intersectSet.empty()){
            //     return {false, 0};
            // }
            candidateSet.insert(intersectSet.begin(), intersectSet.end());
            if (!intersectSet.empty()){
                std::cout << " intersectSet size= " << intersectSet.size() 
                          << " nodes:" << std::endl;
                
                int shown = 0;
                for (int x : intersectSet) {
                    std::cout << x << " ";
                    if (++shown >= 5) break;  // 只打印前 5 个
                }

                if (intersectSet.size() > 5) {
                    std::cout << "...";
                }
                std::cout << std::endl;
            }
        }
    }

    // 输出 selectedNodes 的内容
    //mi std::cout << "selectedNodes: ";
    // for (const int& node : selectedNodes) {
    //     std::cout << node << " ";
    // }
    // std::cout << std::endl;

    candidateSet.insert(selectedNodes.begin(), selectedNodes.end());

    // auto endTime = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsedTime = endTime - startTime;

    //mi std::cout << "candidateSet: ";
    // for (const auto& node : candidateSet) {
    //     std::cout << node << " ";
    // }
    // std::cout << std::endl;

    // 输出时间
    // std::cerr << "Time taken by the last step: " << elapsedTime.count() << " seconds" << std::endl;

    // 判断 candidateSet 的大小是否满足 k + 1
    if (candidateSet.size() >= k + 1) {

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

        // 检查 tiny_graph 是否包含目标节点 v
        clock_t startTinyGraphTime = clock();
        Kcore* kcore = new Kcore(tiny_graph, v, k);
        kcore->search_kcore(tiny_graph);

        if (tiny_graph->hasNode(v)) {
            std::ofstream outputFile(outputFileName);
            //mi std::cout << "Tiny_graph: true" << std::endl;
            clock_t endTinyGraphTime = clock();
            double tinyGraphElapsedTime = static_cast<double>(endTinyGraphTime - startTinyGraphTime) / CLOCKS_PER_SEC;
            //mi std::cerr << "Tiny_graph construction time: " << tinyGraphElapsedTime << " seconds" << std::endl;
            for (int node : tiny_graph->getNodes()) {
                outputFile << "Node ID: " << node << "\n";
            }
            outputFile << "Attributes:" << endl;
            if (commattr.empty()) {
                outputFile << "No attributes\n";
            }
            for (const auto& attr: commattr) {
                outputFile << attr << " ";
            }
            outputFile << "\n";

            int nodeSize = tiny_graph->getNodes().size();
            //mi cout << "VALIDATE:" << nodeSize << std::endl;

            outputFile.close();
            delete tiny_graph;
            delete kcore;

            return {true, nodeSize};
            //返回kcore decomposition后的结点个数
        } else {
            delete tiny_graph;
            delete kcore;
        }
    } else {
        //mi std::cerr << "Skipping candidateSet with size < k+1: " << candidateSet.size() << std::endl;
    }

    return {false, 0};
}

std::vector<std::pair<std::set<int>, std::set<std::string>>> Kcore::getConditionalPatternBase(int targetNodeID, int rootNode, const std::vector<std::vector<Node>>& branches, int totalMaxAttributes) {

    std::vector<std::vector<int>> nCount;   //用于存储结点与共享属性的个数的数组
    std::unordered_map<int, std::set<std::string>> nAttr;   //用于存储结点与共享的属性的哈希表

    for (const auto& branch : branches) {
        bool targetFound = false;
        std::set<std::string> newSet;  //当前的共享属性

        for (auto it = branch.rbegin(); it != branch.rend(); ++it) {
            if (it->id == targetNodeID) {
                targetFound = true;
                newSet = it->attributes;
            }
            if (targetFound) {
                if (nAttr.find(it->id) == nAttr.end()) {
                    nAttr[it->id] = newSet;  //新加入的id
                    nCount.push_back({it->id, static_cast<int>(newSet.size())});
                } else {
                    nAttr[it->id].insert(newSet.begin(), newSet.end());
                    for (auto& row : nCount) {
                        if (row[0] == it->id) {
                            // row[1] += static_cast<int>(newSet.size());
                            row[1] = nAttr[it->id].size();
                            break;
                        }
                    }
                }
            }
        }
    }

    //记录所有prefix path中出现的结点、共享属性个数、属性名

    std::vector<int> keysToRemove;   //因为频率低而移除的结点

    nCount.erase(
        std::remove_if(nCount.begin(), nCount.end(),
            [&keysToRemove, totalMaxAttributes](const std::vector<int> & row) {
                if (row[1] < totalMaxAttributes) {
                    keysToRemove.push_back(row[0]);  //记录被删除的结点ID
                    return true;
                }
                return false;
            }),
            nCount.end());

    for (const int key : keysToRemove) {
        nAttr.erase(key);    //nAttr中也删除对应的结点记录
    }

    std::sort(nCount.begin(), nCount.end(),
            [](const std::vector<int> & a, const std::vector<int> & b) {
                return a[1] > b[1];
            });

    // std::set<int> condMaxSet;
    // std::set<std::string> condMaxAttr;
    // std::vector<std::pair<std::set<int>, std::set<std::string>>> validSets;
    std::set<std::pair<std::set<int>, std::set<std::string>>> uniqueSets;

    for (auto it = nCount.rbegin(); it != nCount.rend(); ++it) {
        int nID = (*it)[0]; //倒序获取nCount中的结点ID，低频到高频

        auto attrIt = nAttr.find(nID);
        if (attrIt != nAttr.end()) {
            const std::set<std::string>& condAttr = attrIt->second;  //nID与targetNode共享的Attr集合
            std::set<int> condSet;   //当前满足的node集合

            for (const auto& candidate: nCount) {    //看其他结点是否包含这个Attr集合
                int candidateID = candidate[0];
                if (auto candidateIt = nAttr.find(candidateID); candidateIt != nAttr.end()) {
                    const std::set<std::string>& tempAttr = candidateIt->second;

                    if (std::includes(tempAttr.begin(), tempAttr.end(), condAttr.begin(), condAttr.end())) {
                        condSet.insert(candidateID);
                    }
                }
            }

            // if (condSet.size() > condMaxSet.size()) {
            //     condMaxSet = condSet;
            //     condMaxAttr = condAttr;
            // }
            if (condSet.size() >= (k + 1)) {
                uniqueSets.emplace(condSet, condAttr);
            }
        }

    }
    // if (condMaxSet.size() >= (k + 1)) {   //返回共享满足至少totalMaxAttr的最大的结点集合，以及对应的属性集合
    //     return {condMaxSet, condMaxAttr};
    // }
    // else {
    //     return {{}, {}};
    // }

    std::vector<std::pair<std::set<int>, std::set<std::string>>> validSets(uniqueSets.begin(), uniqueSets.end());

    std::sort(validSets.begin(), validSets.end(),
        [](const auto& a, const auto& b) {
            return a.first.size() > b.first.size();
        });

    return validSets;
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

            // 检查 nextNode 是否存在于 neighbors 中
            if (std::find(neighbors.begin(), neighbors.end(), nextNode) != neighbors.end()) {
                subgraph->addEdge(node, nextNode);
                //mi std::cerr << "nextNode " << nextNode << " is a neighbor of " << node << std::endl;
            }
            else {
                //mi std::cerr << "Warning: nextNode " << nextNode << " is not a neighbor of " << node << std::endl;
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

void Kcore::read_attrmap(const std::string& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open the file: " << file_path << std::endl;
        return;
    }

    // attr_map.clear();
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) continue; // 跳过空行

        // 去掉行首尾空格
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        // 检查是否包含 ':'，确保格式正确
        size_t colonPos = line.find(':');
        if (colonPos == std::string::npos) continue; // 如果没有 ':'，跳过该行

        std::string key = line.substr(0, colonPos); // 提取 key
        std::string values_str = line.substr(colonPos + 1); // 提取 values 部分

        // 去掉 key 和 values_str 的首尾空格
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        values_str.erase(0, values_str.find_first_not_of(" \t"));
        values_str.erase(values_str.find_last_not_of(" \t") + 1);

        // 分割 values_str 并解析为整数（以空格分隔）
        std::unordered_set<int> values;
        std::istringstream values_iss(values_str);
        std::string value;
        while (values_iss >> value) {
            try {
                values.insert(std::stoi(value)); // 插入整数
            } catch (const std::exception& e) {
                std::cerr << "Warning: Invalid value '" << value << "' in file." << std::endl;
                continue;
            }
        }

        // // 调试输出针对 "retrieval" 和 "shape"
        // if (key == "retrieval" || key == "shape") {
        //     std::cerr << "Debug - Key: " << key << std::endl;
        //     std::cerr << "Debug - Raw Line: " << line << std::endl;
        //     std::cerr << "Debug - Values: ";
        //     for (const auto& val : values) {
        //         std::cerr << val << " ";
        //     }
        //     std::cerr << std::endl;
        // }

        // 插入 attr_map
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
    std::stack<std::pair<int, TNode*>> nodeStack; // 用于维护节点层级关系，存储缩进和节点指针
    root = nullptr; // 树的根节点

    const int INDENT_UNIT = 2; // 定义每一级缩进的空格数

    while (std::getline(file, line)) {
        if (line.empty()) continue; // 跳过空行

        // 计算当前行的缩进级别
        int indentLevel = 0;
        while (indentLevel < line.size() && line[indentLevel] == ' ') {
            indentLevel++;
        }
        indentLevel /= INDENT_UNIT; // 将缩进级别转换为逻辑层级

        // 去除行中的缩进部分
        line = line.substr(indentLevel * INDENT_UNIT);

        // 解析行内容
        std::istringstream iss(line);
        std::string key;
        int core = -1;
        int size = 0;
        std::unordered_set<int> nodeSet;

        // 提取 k 和 size 的值，以及节点集合
        if (std::getline(iss, key, '=') && key == "k") {
            iss >> core; // 提取 k 值
            iss.ignore(20, ' '); // 跳过空格

            if (std::getline(iss, key, '=') && key == "size") {
                iss >> size; // 提取 size 值
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
        TNode* currentNode = new TNode(core, nullptr);
        currentNode->setNodeSet(nodeSet);

        // 更新索引表
        for (int id : nodeSet) {
            nodeIndex[id] = currentNode;
        }

        buildInvertedIndex(currentNode);

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
                // nodeStack.top().second->getChildList().push_back(currentNode);
                TNode* parentNode = nodeStack.top().second;
                parentNode->getChildList().push_back(currentNode);
                currentNode->setParent(parentNode);
            }

            // 将当前节点加入栈
            nodeStack.push({ indentLevel, currentNode });
        }
    }

    // std::unordered_map<int, TNode*> ancestorCorenessMap;
    // precomputeAncestorCoreness(root, ancestorCorenessMap);
    // aggregateNodeSetBottomUp(root);

    file.close();
}

void Kcore::buildInvertedIndex(TNode* currentNode){
    if (!currentNode) return;

    const std::unordered_set<int>& nodeSet = currentNode->getNodeSet();
    auto& invertedIndex = currentNode->getInvertedIndex();

    for (int nodeId : nodeSet) {
        std::vector<std::string> nodeAttributes(attrs(nodeId).begin(), attrs(nodeId).end());

        for (const std::string& attribute : nodeAttributes) {
            invertedIndex[attribute].push_back(nodeId);
        }
    }
}

void Kcore::collectTreeNodes(TNode* root, std::vector<TNode *>& subtreeNodes, std::unordered_set<TNode*>& visited) {
    // if (root == nullptr) {
    //     return;
    // }
    
    if (!root || visited.count(root)) {
        return;
    }
    visited.insert(root);

    subtreeNodes.push_back(root);

    for (TNode* child : root->getChildList()){
        collectTreeNodes(child, subtreeNodes, visited);
    }
}

void Kcore::display() {
    
}
