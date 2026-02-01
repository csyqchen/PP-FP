#include "../include/Graph.h"
#include <vector>
#include <unordered_set>
#include <cassert>
#include <algorithm>

using namespace std;

Graph::Graph() {
    maxN = 1;
}

Graph::Graph(const Graph& graph) {
    maxN = graph.maxN;
    ng = graph.ng;
    nodes = graph.nodes;
    edges = graph.edges;
    //attr = graph.attr;

}

Graph::Graph(const Graph* p_graph) {
    maxN = p_graph->maxN;
    ng = p_graph->ng;
    nodes = p_graph->nodes;
    edges = p_graph->edges;
    //attr = p_graph->attr;

}

Graph::~Graph() {
}

void Graph::setNodes(std::unordered_set<int> data) {
    nodes = data;
}

void Graph::setEdges(std::unordered_set<Edge, MyEdgeHash> data) {
    edges = data;
}
void Graph::setAttrs(std::unordered_map<int, Attributes> data) {
    attr = data;
}
void Graph::setNbs(std::unordered_map<int, Neighbors> data) {
    ng = data;
}

void Graph::setData(Graph* p_graph) {
    nodes.clear();
    edges.clear();
    ng.clear();

    maxN = p_graph->getMaxN();
    nodes = p_graph->getNodes();
    edges = p_graph->getEdges();
    //attr = p_graph->getAttrs();
    ng = p_graph->getNeighbors();
}

void Graph::clear() {
    nodes.clear();
    edges.clear();
    ng.clear();
}

//const Graph::Neighbors& Graph::neighbors(int v) {
//    return ng[v];
//}

const Graph::Neighbors& Graph::neighbors(int v) const {
    static const Neighbors emp;
    if (!hasNode(v)) { return emp; }
    return ng.at(v);
}

const int Graph::neighbors_node(int v, int index) const {
    Edge e = ng.at(v)[index];
    if (e.getX() == v) { return e.getY(); }
    else return e.getX();
}


const Graph::Attributes& Graph::attributes(int v) {
    return attr[v];
}

int Graph::deg(int v) const {
    if (ng.count(v) > 0)
        return (int)ng.at(v).size();
    return 0;
}

void Graph::add_attributes(int v, Attributes atts) {
    attr[v].insert(atts.begin(), atts.end());
}

void Graph::remove_attribute(int v, Attributes atts) {
    if (attr.find(v) != attr.end()) {
        for (const auto& att : atts) {
            attr[v].erase(att);  // �ӽڵ� v �����Լ������Ƴ�ָ������
        }
    }
}

void Graph::addEdge(int x, int y) {
    if (x == y) return;
    if (x > y)
        swap(x, y);
    if (hasEdge(x, y))
        return;
    Edge e = Edge(x, y);
    edges.insert(e);
    ng[x].push_back(e);
    ng[y].push_back(e);
    if (y >= maxN)
        maxN = y + 1;
    nodes.insert(x);
    nodes.insert(y);
}

void Graph::add_edge(int x, int y) {
    if (x == y) return;
    if (x > y)
        swap(x, y);
    /*if (hasEdge(x, y))
        return;*/
    Edge e = Edge(x, y);
    edges.insert(e);
    ng[x].push_back(e);
    ng[y].push_back(e);
    if (y >= maxN)
        maxN = y + 1;
    nodes.insert(x);
    nodes.insert(y);
}

void Graph::deleteEdge(Edge e) {
    if (!hasEdge(e)) return;

    // Remove from edges
    {
        //auto it = find(edges.begin(), edges.end(), e);
        //if (it == edges.end()) return;
        edges.erase(e);
    }

    int x = e.getX(), y = e.getY();
    // Remove from ng
    {
        if (ng.find(x) != ng.end()) {
            auto it = find(ng[x].begin(), ng[x].end(), e);
            assert(it != ng[x].end());
            ng[x].erase(it);
        }

        if (ng.find(y) != ng.end()) {
            auto it = find(ng[y].begin(), ng[y].end(), e);
            assert(it != ng[y].end());
            ng[y].erase(it);
        }

    }

}

void Graph::deleteNode(int v) {
    if (!hasNode(v))
        return;

    vector<Edge> edgesToBeDeleted = neighbors(v);
    for (Edge e : edgesToBeDeleted) {
        deleteEdge(e);
    }

    nodes.erase(v);
    ng.erase(v);

}

void Graph::addEdge(Edge e) {
    addEdge(e.getX(), e.getY());
}

const unordered_set<Edge, MyEdgeHash>& Graph::getEdges() const {
    return edges;
}

bool Graph::hasEdge(int x, int y) const {
    if (x > y)
        swap(x, y);
    //return find(edges.begin(), edges.end(), Edge(x, y)) != edges.end();
    return edges.count(Edge(x, y));
}

bool Graph::hasEdge(Edge e) const {
    return hasEdge(e.getX(), e.getY());
}

bool Graph::hasNode(int v) const {
    return nodes.count(v);
}

bool Graph::addNode(int v) {
    if (v > maxN)
        maxN = v;
    // nodeIds.push_back(v);
    // idToIndex[v] = nodeIds.size() - 1; // ���ڵ� ID ӳ�䵽������
    return nodes.insert(v).second;
}

int Graph::getMaxN() const {
    return maxN;
}

void Graph::unionGraph(const Graph& g) {
    for (Edge e : g.getEdges())
        addEdge(e);
}

void Graph::generateIndices() {
    nodeIds.clear();
    idToIndex.clear();
    int index = 0;
    for (int nodeId : nodes) {
        nodeIds.push_back(nodeId); // �������ڵ� ID ��ӳ��
        idToIndex[nodeId] = index; // �ڵ� ID ��������ӳ��
        index++;
    }
}

int Graph::getN() const {
    return (int)nodes.size();
}

int Graph::getM() const {
    return (int)edges.size();
}

const std::unordered_set<int>& Graph::getNodes() const {
    return nodes;
}

std::unordered_map<int, Graph::Neighbors> Graph::getNeighbors() {
    return ng;
}

std::unordered_map<int, Graph::Attributes> Graph::getAttrs() {
    return attr;
}

// void Graph::load_graph(std::string data_path) {

//     std::ifstream infile(data_path);
//     if (!infile) {
//         std::cerr << "Unable to open file data.txt";
//     }

//     std::string line;

//     while (std::getline(infile, line)) {
//         std::istringstream iss(line);
//         int num1, num2;

//         if (!(iss >> num1 >> num2)) {
//             std::cerr << "Error reading line: " << line << std::endl;
//             continue; // Skip lines that do not match the format
//         }

//         add_edge(num1, num2);
//     }

//     infile.close();

// }

void Graph::load_graph(std::string data_path) {
    std::ifstream infile(data_path);
    if (!infile) {
        std::cerr << "Unable to open file: " << data_path << std::endl;
        return;
    }

    std::string line;

    while (std::getline(infile, line)) {

        line.erase(std::remove(line.begin(), line.end(), '#'), line.end());

        std::istringstream iss(line);
        int num1, num2;

        if (!(iss >> num1 >> num2)) {
            std::cerr << "Error reading line: " << line << std::endl;
            continue;
        }
        add_edge(num1, num2);
    }

    infile.close();
}

void Graph::load_attribute(std::string data_path) {

    std::ifstream infile(data_path);
    if (!infile) {
        std::cerr << "Unable to open file data.txt";
    }

    std::string line;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);

        int node_id;
        std::unordered_set<std::string> tags;
        std::string numberString, tagsString;

        // Extract name
        size_t firstHash = line.find('#');
        size_t secondHash = line.find('#', firstHash + 1);

        // Extract number
        size_t thirdHash = line.find('#', secondHash + 1);
        size_t fourthHash = line.find('#', thirdHash + 1);
        numberString = line.substr(thirdHash + 1, fourthHash - thirdHash - 1);
        node_id = std::stoi(numberString);

        // Extract tags
        size_t fifthHash = line.find('#', fourthHash + 1);
        size_t sixthHash = line.find('#', fifthHash + 1);
        tagsString = line.substr(fifthHash + 1, sixthHash - fifthHash - 1);
        std::istringstream tagStream(tagsString);
        std::string tag;
        while (std::getline(tagStream, tag, ';')) {
            tags.insert(tag);
        }

        // Add entry to the data map
        attr[node_id] = tags;


    }

    infile.close();

}

std::vector<int> Graph::getNeighborIDs(int v) const {
    std::vector<int> neighborIDs;
    const Neighbors& edges = neighbors(v);
    for (const Edge& edge : edges) {
        int neighborID = (edge.getX() == v) ? edge.getY() : edge.getX();
        neighborIDs.push_back(neighborID);
    }
    return neighborIDs;
}

std::vector<int> Graph::getNeighborIndices(int v) const {
    std::vector<int> neighborIndices;
    const Neighbors& edges = neighbors(v);
    for (const Edge& edge : edges) {
        int neighborID = (edge.getX() == v) ? edge.getY() : edge.getX();
        neighborIndices.push_back(idToIndex.at(neighborID));
    }
    return neighborIndices;
}

// std::vector<int> Graph::getNeighborIDs(int v) const {
//     std::vector<int> neighborIDs;

//     // ���� graphData �Ǵ洢�ڽӱ��� unordered_map
//     if (graphData.find(v) == graphData.end()) {
//         std::cerr << "Warning: Node " << v << " not found in the graph." << std::endl;
//         return neighborIDs;  // ���ؿյ��ھ��б�
//     }

//     const Neighbors& edges = neighbors(v);
//     for (const Edge& edge : edges) {
//         int neighborID = (edge.getX() == v) ? edge.getY() : edge.getX();
//         neighborIDs.push_back(neighborID);
//     }
//     return neighborIDs;
// }


