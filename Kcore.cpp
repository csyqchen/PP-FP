#include "Kcore.h"
#include <queue>
#include "Utils.h"
#include <memory>
#include <map>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>

Kcore::Kcore(Graph* p_graph, int query_node, int k_coreness) :v(query_node), k(k_coreness) {
    graph = p_graph;
}

Kcore::Kcore(Graph* p_graph, int k_coreness) : k(k_coreness) {
    graph = p_graph;
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
int Kcore::deg(int v) {
    return graph->deg(v);
}

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
        } else {
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
    //键是属性，值是具有该属性的结点集合

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

            }catch (const std::bad_alloc& e) {
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
                std::cout << "size of nodes:" << community <<std::endl;
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
                for (auto node : nodes){
                    if (!nodes_candidate.count(node)) graph_copy->deleteNode(node);
                }
            } catch (const std::bad_alloc& e) {
                std::cerr << "Memory allocation failed: " << e.what() << std::endl;
                continue;
            }

            search_kcore(graph_copy.get());
            if (graph_copy->getNodes().count(v)) {
                // combinations_candidate.push_back(attr_candidate); // Storing all successful combinations
                //printf("find the community with max attrs: %zu\n", mid_d);
                output_file << "find the community with max attrs: " << mid_d << "\n";
                for (const auto& str :attr_candidate){
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
        } else if (!find && mid_d > success_d + 1) {
            max_d = mid_d - 1;
        } else if (!find && mid_d <= success_d + 1) {
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

void Kcore::baseline3() {
    clock_t start = clock();
    //find kcore to filter graph
    search_kcore(graph);
    if (!getNodes().count(v)) {
        printf("query node %d is not in a %d-core\n", v, k);
        return;
    }

    Graph::Attributes query_attributes = attrs(v);
    std::vector<std::string>query_attributes_list(query_attributes.begin(), query_attributes.end());
    std::unordered_map<std::string, std::unordered_set<int>> attr_map;

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
 
    size_t min_d = 1;
    size_t max_d = query_attributes_list.size();
    std::vector<std::vector<std::string>> combinations_candidate;

    bool is_searching_up = false;
    size_t previous_mid_d = 0;

    while (min_d <= max_d){
        size_t mid_d = (min_d + max_d) / 2;
        printf("start searching the community with attrs: %zu\n", mid_d);

        //std::vector<std::vector<std::string>> attr_combinations = combinations(query_attributes_list, mid_d);
        std::vector<std::vector<std::string>> attr_combinations;
        if (is_searching_up && !combinations_candidate.empty()){
            //向上搜索且至少存储一个成功组合
            for (auto& prev_comb : combinations_candidate){
                std::vector<std::string> reduced_attr_list = query_attributes_list;
                for (auto& used_attr : prev_comb){
                    reduced_attr_list.erase(std::remove(reduced_attr_list.begin(), reduced_attr_list.end(), used_attr), reduced_attr_list.end());
                }
                //apriori式计算新组合ß
                auto new_combinations = combinations(reduced_attr_list, mid_d - previous_mid_d);
                attr_combinations.insert(attr_combinations.end(), new_combinations.begin(), new_combinations.end());
            }
        }
        else {
            //向下搜索
            attr_combinations = combinations(query_attributes_list, mid_d);
        }

        bool find = false;

        for (auto attr_candidate : attr_combinations){
            std::unordered_set<int> nodes_candidate = attr_map[attr_candidate[0]];

            for (size_t i = 1; i < attr_candidate.size(); i++)
            {
                intersect_filter(nodes_candidate, attr_map[attr_candidate[i]]);
            }

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

            }catch (const std::bad_alloc& e) {
                std::cerr << "Memory allocation failed: " << e.what() << std::endl;
                continue;
            }
            

            // for (auto node : graph_copy->getNodes()) {
            //     if (!nodes_candidate.count(node)) graph_copy->deleteNode(node);
            // }

            search_kcore(graph_copy.get());
            if (graph_copy->getNodes().count(v)) {
                printf("find the community with max attrs: %zu\n", mid_d);
                combinations_candidate.push_back(attr_candidate);
                for (const auto& str : attr_candidate) {
                    std::cout << str << " ";
                }
                printf("\n");
    
                std::cout << "size of nodes:" << graph_copy->getNodes().size() <<std::endl;
                if (graph_copy->getNodes().size() < 20) {
                    for (const auto& node : graph_copy->getNodes()) {
                        std::cout << node << " ";
                    }
                    printf("\n");
                }
                printf("\n");
                find = true;
                //break;
            }
            
        }
        if (find) {
            min_d = mid_d + 1;
            is_searching_up = true;
            previous_mid_d = mid_d;
        } else {
            max_d = mid_d - 1;
            is_searching_up = false;
        }
    }
    clock_t end = clock();
    double duration = (end - start) / CLOCKS_PER_SEC;
    printf("baseline2 cost:%2f \n", duration);

    for (const auto& comb : combinations_candidate){
        for (const auto& attr : comb){
            std::cout << attr << " ";
        }
        std::cout << std::endl;                             
    }
}

std::tuple<double, int, int> Kcore::fptree(int k, const std::string& fptreeFile, const std::string& corenessFile, const std::string& outputFileName) {

    std::map<int, std::pair<int,int>> corenessMap; //存储coreness和step

    clock_t start = clock();

    std::ifstream coreFile (corenessFile);
    if (!coreFile.is_open()) {
        std::cerr << "Cannot open the file: " << corenessFile << std::endl;
        return std::make_tuple(-1.0, 0, 0);
    }

    std::string line;
    getline(coreFile, line); //跳过首行

    while (getline(coreFile, line)) {
        if (line.empty() || line[0] == 'Q') continue;

        size_t last_space = line.rfind(' ');
        if (last_space == std::string::npos) continue;  // No space found, skip the line

        size_t second_last_space = line.rfind(' ', last_space - 1);
        if (second_last_space == std::string::npos) continue;  // No second space found, skip the line

        int coreness = std::stoi(line.substr(last_space + 1));
        int step = std::stoi(line.substr(second_last_space + 1, last_space - second_last_space - 1));
        int nodeID = std::stoi(line.substr(0, second_last_space));
        
        corenessMap[nodeID] = {step, coreness};

    }
    coreFile.close();

    std::ifstream fpFile(fptreeFile);
    if (!fpFile.is_open()) {
        std::cerr << "Cannot open the file: " << fptreeFile << std::endl;
        return std::make_tuple(-1.0, 0, 0);
    }

    std::vector<std::vector<Node>> branches; // 存储不同的支链
    std::map<int, std::vector<Node>> levelToBranchMap; //存储每个层次对应的分支
    int rootNode = v; //根结点为query node
    std::cerr << "rootNode: " << v << std::endl;


    // std::string line;
    std::vector<Node> currentBranch;
    int previousLevel = -1; //用于跟踪上一个结点的level

    int branchCounter = 0; //分支计数器

    while (std::getline(fpFile, line)) {
        int indentLevel = 0;
        while (indentLevel < line.length() && line[indentLevel] == ' ') {
            indentLevel++;
        }
        line = line.substr(indentLevel); // 去除缩进
        
        std::istringstream iss(line);
        int nodeID, count, totalCount;
        char tempChar;

        iss >> nodeID;
        iss.ignore(20, '(');
        iss >> count >> tempChar >> totalCount;

        int currentLevel = indentLevel / 2;

        auto coreData = corenessMap[nodeID];
        Node currentNode = { nodeID, count, totalCount, currentLevel, coreData.first, coreData.second};

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

    //input结束

    // int targetNodeID = 274774;
    // std::cout << "rootNode in fptree: " << rootNode << std::endl;
    // buildConditionalFPTree(targetNodeID, rootNode, branches);

    bool foundKCore = false;
    int maxAttributes = 0;
    int totalMaxAttributes = branches[1][1].totalCount; //最大的可能attr
    std::cout << "totalMaxAttributes" << totalMaxAttributes;
    std::vector<int> selectedNodes;
    std::vector<int> tempSelectedNodes;
    int tempMaxAttributes = 0;

    // for (const auto& branch : branches) {
    //     //std::vector<int> candidateNodes = { rootNode };
    //     int currentCount = branch[0].count;
    //     int branchMaxAttributes = 0; //当前分支的第k个count
    //     bool foundCurrentKCore = false;
    //     std::vector<int> currentSelectedNodes;

        //new
    // for (const auto& branch : branches) {
    //     std::vector<int> candidateNodes = { rootNode };
    //     for (const auto& node : branch) {
    //         if (node.count == totalMaxAttributes && node.coreness >= k) {
    //             candidateNodes.push_back(node.id);
    //             if (candidateNodes.size() >= k + 1) {
    //                 if (isKCore(k, candidateNodes)) {
    //                     foundKCore = true;
    //                     selectedNodes = candidateNodes;
    //                     maxAttributes = totalMaxAttributes;
    //                     break;
    //                 }
    //             }
    //         }
    //     }
    //     if (foundKCore) break;
    // } //mi
    for (const auto& branch : branches) {
        std::vector<int> candidateNodes = { rootNode };
        for (const auto& node : branch) {
            // 只要找到一个符合条件的节点就终止
            if (node.count == totalMaxAttributes && node.coreness >= k) {
                candidateNodes.push_back(node.id);

                // 检查当前节点集合是否满足 k-core 的条件
                if (candidateNodes.size() >= k + 1 && isKCore(k, candidateNodes)) {
                    foundKCore = true;
                    selectedNodes = candidateNodes;
                    maxAttributes = totalMaxAttributes;
                    break;  // 结束内层循环（遍历 branch 中的 node）
                }
            }
        }

        if (foundKCore) {
            break;  // 结束外层循环（遍历 branches）
        }
    }


    while (!foundKCore && totalMaxAttributes > 0) {
        totalMaxAttributes--;
        std::cout << "totalMaxAttributes" << totalMaxAttributes;
        std::cerr << "Searching with totalMaxAttributes = " << totalMaxAttributes << std::endl;
        
        for (const auto& branch : branches) {
            for (const auto& node : branch) {
                if (node.totalCount >= totalMaxAttributes && node.coreness >= k && built.find(node.id)== built.end()) {
                    // 改动过
                    std::cerr << "Building conditional tree for Node ID: " << node.id << std::endl;
                    auto [candidateNodes, threshold] = buildConditionalFPTree(node.id, rootNode, branches);
                    if (!candidateNodes.empty()) {
                        if (threshold == totalMaxAttributes) {
                            foundKCore = true;
                            selectedNodes = candidateNodes;
                            maxAttributes = threshold;
                            break;
                        } else if (threshold > tempMaxAttributes) {
                            tempSelectedNodes = candidateNodes;
                            tempMaxAttributes = threshold;
                        }
                    }
                }
            }
            if (foundKCore) break;
        }

        if (!foundKCore && tempMaxAttributes >= totalMaxAttributes - 1) {
            foundKCore = true;
            selectedNodes = tempSelectedNodes;
            totalMaxAttributes = tempMaxAttributes; //???
            maxAttributes = tempMaxAttributes;
        } else if (!foundKCore) {
            // totalMaxAttributes --;
        }
    }

    if (foundKCore) {
        std::ofstream outputFile(outputFileName);
        if (outputFile.is_open()) {
            outputFile << "Find " << k << "-core community with " << totalMaxAttributes << " attributes\n";
            for (int node : selectedNodes) {
                outputFile << "Node ID: " << node << "\n";
            }
            outputFile.close();
        } else {
            std::cerr << "Cannot open file: " << outputFileName << std::endl;
        }
    } else {
        std::cerr << "No k-core community with at least 1 attribute found." << std::endl;
    }

    clock_t end = clock();
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "\nfptree cost: " << duration << " \n";
    return std::make_tuple(duration, maxAttributes, selectedNodes.size());
}

std::vector<std::vector<Kcore::Node>> Kcore::getConditionalPatternBase(int targetNodeID, int rootNode, const std::vector<std::vector<Kcore::Node>>& branches) {
    std::vector<std::vector<Kcore::Node>> conditionalPatternBase;
    std::set<std::vector<int>> uniquePaths; //存储唯一的路径标识


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

std::pair<std::vector<int>, int> Kcore::buildConditionalFPTree(int targetNodeID, int rootNode, const std::vector<std::vector<Kcore::Node>>& branches) {
    std::vector<std::vector<Kcore::Node>> conditionalPatternBase = getConditionalPatternBase(targetNodeID, rootNode, branches);

    std::cout << "rootNode in buildConditionalFPTree: " << rootNode << std::endl;

    std::map<int, int> itemFrequency;

    for (const auto& pattern : conditionalPatternBase) {
        int lastNodeCount = pattern.back().count;
        for (const auto& node : pattern) {
            std::cout << "Node ID: " << node.id << " Count: " << node.count << "\n";
            // itemFrequency[node.id] ++;
            itemFrequency[node.id] += lastNodeCount;
        }
    }

    std::vector<int> frequencies;
    for (const auto& item : itemFrequency) {
        frequencies.push_back(item.second);
    }

    std::sort(frequencies.begin(), frequencies.end(), std::greater<int>());

    // int threshold = frequencies.size() >= k ? frequencies[k] : 0;
    // int threshold = frequencies.size() >= k ? frequencies[k + 1] : 0;
    int threshold = 0;
    if (frequencies.size() >= k + 1) {
        threshold = frequencies[k];  
        std::cout << "threshold:" << threshold <<  "\n";
    } else {
        threshold = 0;
    }
    if (threshold > 0) {
        std::cout << "Nodes in conditional FP-Tree for node " << targetNodeID << " with frequency >= " << threshold << ":\n";

        std::vector<int> candidateNodes;
        for (const auto& item : itemFrequency) {
            if (item.second >= threshold) {
                std::cout << "Node ID: " << item.first << " Frequency: " << item.second << "\n";
                candidateNodes.push_back(item.first);
            }
        }

        if (isKCore(k, candidateNodes)) {
            std::cout << "The nodes form a " << k << "-core subgraph.\n";
            return {candidateNodes, threshold};
        } else {
            std::cout << "The nodes do not form a " << k << "-core subgraph.\n";
            built.insert(targetNodeID);
            std::cout << "Current built set contains the following nodes: " << std::endl;
            if (built.empty()) {
                std::cout << "Built set is currently empty." << std::endl;
            } else {
            for (const int& node : built) {
                std::cout << "Node ID: " << node << std::endl;
                }
            }
            return {{}, 0};
        }
    }
    else {
        std::cout << "threshold:" << threshold << std::endl;
        built.insert(targetNodeID);
        std::cout << "Current built set contains the following nodes: " << std::endl;
        if (built.empty()) {
            std::cout << "Built set is currently empty." << std::endl;
        } else {
            for (const int& node : built) {
            std::cout << "Node ID: " << node << std::endl;
            }
        }
        return {{}, 0};
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

            // 检查 nextNode 是否存在于 neighbors 中
            if (std::find(neighbors.begin(), neighbors.end(), nextNode) != neighbors.end()) {
                subgraph->addEdge(node, nextNode);
                std::cerr << "nextNode " << nextNode << " is a neighbor of " << node << std::endl;
            } else {
                std::cerr << "Warning: nextNode " << nextNode << " is not a neighbor of " << node << std::endl;
            }
        }
    }

    search_kcore(subgraph.get());
    return check_kcore(subgraph.get(), k, nodes);
}

void Kcore::display(){
    // display函数的具体实现
}

