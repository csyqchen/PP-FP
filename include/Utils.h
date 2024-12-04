
#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <set>
#include <stack>
#include <queue>

//************************************************************ data type **************************************

struct Edge {
    int x, y;

    Edge() : x(-1), y(-1) {}
    Edge(int x, int y) { setNodes(x, y); }

    void updateEdge(int x, int y) { setNodes(x, y); }

    int getX() const { return x; }
    int getY() const { return y; }

    int other(int z) const {
        if (z != x && z != y)
            throw std::invalid_argument("Invalid node value");
        return x + y - z;
    }

    bool operator<(const Edge& e) const {
        return (this->getX() < e.getX() || (this->getX() == e.getX() && this->getY() < e.getY()));
    }

    bool operator==(const Edge& e) const {
        return this->getX() == e.getX() && this->getY() == e.getY();
    }


private:
    void setNodes(int x, int y) {
        if (x > y)
            std::swap(x, y);
        this->x = x;
        this->y = y;
    }
};

struct MyEdgeHash {
    std::size_t operator()(const Edge& s) const {
        return std::hash<int>()(s.x) ^ (std::hash<int>()(s.y) << 1);
    }
};

struct PrivateEntry {
    std::string name;
    int number;
    std::unordered_set<std::string> tags;
};

class UNode {
public:
    int value = 0;
    UNode* parent = this;
    int rank = -1;
    int represent = -1; // This variable is used for updating our tree index

    // Constructor
    UNode() {}
    UNode(int val) : value(val) {}
    //UNode(int val) : value(val), parent(this), rank(0), represent(val) {}
};

class TNode {
//private:
    int core = -1;
    std::set<int> nodeSet;
    std::set<int> coreSet;
    std::vector<TNode*> childList;

public:
    // Constructor
    TNode(int core) : core(core) {}

    // Getter and setter for core
    int getCore() const {
        return core;
    }
    void setCore(int core) {
        this->core = core;
    }

    // Getter and setter for nodeSet
    std::set<int>& getNodeSet() {
        return nodeSet;
    }
    void setNodeSet(const std::set<int>& newNodeSet) {
        nodeSet = newNodeSet;
    }

    // Getter and setter for coreSet
    void setCoreSet(const std::set<int>& nodes) { coreSet = nodes; }
    const std::set<int>& getCoreSet() const { return coreSet; }

    // Getter and setter for childList
    std::vector<TNode*>& getChildList() {
        return childList;
    }
    void setChildList(const std::vector<TNode*>& newChildList) {
        childList = newChildList;
    }

    void update(const std::set<int>& nodes) {
        coreSet.insert(nodes.begin(), nodes.end());
    }

    // Update coreSet from subtree
    void updateCoreSetFromSubtree() {
        // 初始化 coreSet 为当前节点的 nodeSet
        this->setCoreSet(this->getNodeSet());

        // 使用队列实现层级遍历
        std::queue<TNode*> queue;
        queue.push(this); // 将当前节点作为起点加入队列

        // 遍历子树中所有节点
        while (!queue.empty()) {
            TNode* current = queue.front();
            queue.pop();

            if (!current) continue; // 跳过空节点

            // 合并当前节点的 nodeSet 到根节点的 coreSet
            this->update(current->getNodeSet());

            // 将子节点加入队列
            for (TNode* child : current->getChildList()) {
                if (child) {
                    queue.push(child);
                }
            }
        }
    }

};

class UnionFind {
public:
    void makeSet(UNode* x) {
        x->parent = x;
        x->rank = 0;
        x->represent = x->value; // Initialize as itself for tree index
    }

    UNode* find(UNode* x) {
        if (x->parent != x) {
            x->parent = find(x->parent); // Path compression
        }
        return x->parent;
    }

    void unionSets(UNode* x, UNode* y) {
        UNode* xRoot = find(x);
        UNode* yRoot = find(y);

        if (xRoot == yRoot) {
            return; // x and y are already in the same set
        }

        // Merge the sets based on rank
        if (xRoot->rank < yRoot->rank) {
            xRoot->parent = yRoot;
        }
        else if (xRoot->rank > yRoot->rank) {
            yRoot->parent = xRoot;
        }
        else {
            yRoot->parent = xRoot;
            xRoot->rank = xRoot->rank + 1;
        }
    }
};

// ********************************************  data io ********************************************************************

inline void load_public_graph(std::string path, std::vector<std::pair<int, int>>& data) {
    std::ifstream infile(path);
    if (!infile) {
        std::cerr << "Unable to open file data.txt";
    }

    std::string line;
    data.clear();

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int num1, num2;

        if (!(iss >> num1 >> num2)) {
            std::cerr << "Error reading line: " << line << std::endl;
            continue; // Skip lines that do not match the format
        }

        data.emplace_back(num1, num2);
    }

    infile.close();
}

inline void load_private_file(const std::string& filename, std::unordered_map<int, std::vector<PrivateEntry>>& data) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Unable to open file " << filename << std::endl;
        return;
    }

    std::string line;
    int currentKey = -1;
    int linesToRead = 0;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (line[0] != '#') {
            // This line should contain two numbers
            if (!(iss >> currentKey >> linesToRead)) {
                continue;
            }
            //std::cout << "starting process id:" << currentKey <<std::endl;
        }
        else {
            // This line should contain a data entry
            PrivateEntry entry;
            std::string numberString, tagsString;

            // Extract name
            size_t firstHash = line.find('#');
            size_t secondHash = line.find('#', firstHash + 1);
            entry.name = line.substr(firstHash + 1, secondHash - firstHash - 1);

            // Extract number
            size_t thirdHash = line.find('#', secondHash + 1);
            size_t fourthHash = line.find('#', thirdHash + 1);
            numberString = line.substr(thirdHash + 1, fourthHash - thirdHash - 1);
            entry.number = std::stoi(numberString);

            // Extract tags
            size_t fifthHash = line.find('#', fourthHash + 1);
            size_t sixthHash = line.find('#', fifthHash + 1);
            tagsString = line.substr(fifthHash + 1, sixthHash - fifthHash - 1);
            std::istringstream tagStream(tagsString);
            std::string tag;
            while (std::getline(tagStream, tag, ';')) {
                entry.tags.insert(tag);
            }

            // Add entry to the data map
            data[currentKey].push_back(entry);

            // Decrement linesToRead
            --linesToRead;
        }
    }
}


inline void count_private(const std::string& filename, int& res) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Unable to open file " << filename << std::endl;
        return;
    }

    int count = 0;

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (line[0] != '#') {
            // This line should contain two numbers
            count++;
        }
        else {
            continue;
        }
    }
    res = count;
}



// ********************************************  functions ********************************************************************

// ����ģ�庯��
template<typename T>
std::unordered_set<T> intersectHelper(const std::unordered_set<T>& set1, const std::unordered_set<T>& set2) {
    std::unordered_set<T> ret_set;
    if (set1.size() < set2.size()) {
        for (const auto& elem : set1) {
            if (set2.count(elem)) {
                ret_set.insert(elem);
            }
        }
    }
    else {
        for (const auto& elem : set2) {
            if (set1.count(elem)) {
                ret_set.insert(elem);
            }
        }
    }
    return ret_set;
}

template<typename T, typename... Sets>
std::unordered_set<T> intersectHelper(const std::unordered_set<T>& set1, const std::unordered_set<T>& set2, const Sets&... sets) {
    auto temp_result = intersectHelper(set1, set2);
    return intersectHelper(temp_result, sets...);
}


template<typename T, typename... Sets>
std::unordered_set<T> intersect(const std::unordered_set<T>& firstSet, const Sets&... sets) {
    return intersectHelper(firstSet, sets...);
}