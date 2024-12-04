#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>
#include <set>
#include "../include/Graph.h"
#include "../include/Utils.h"
#include "../include/Kcore.h"
#include <chrono>
#include <thread>
#include <future>

using namespace std;

int main()
{
    ofstream total_output_file("test/Output_3_fp_5.txt");

    Graph* graph = new Graph();

    //public graph data
    printf("start loading public data...... ");
    graph->load_graph("test/public_graph.txt");
    graph->load_attribute("test/public_attribute.txt");
    printf("public graph loaded \n");

    //private graph data
    printf("start loading private index...... ");
    unordered_map<int, vector<PrivateEntry>> private_index;
    load_private_file("test/test_private.txt", private_index);
    printf("private index loaded \n");

    //public graph index
    Kcore* kcore = new Kcore(graph);
    kcore->read_attrmap("test/attribute_index.txt");
    kcore->read_treeIndex("test/tree_index.txt");

    set<Edge> tempEdges;
    unordered_map<int, unordered_set<string>> tempAttributes;
    std::vector<std::string> allResults;
    for (const auto& pair : private_index)
    {
        int query_node = pair.first;
        const vector<PrivateEntry>& query_private_entry = pair.second;

        size_t private_edge_count = query_private_entry.size();

        // std::ofstream output_file("test_icde/output_5_fp_" + std::to_string(query_node) + ".txt");  //fptree

        unordered_set<string> query_node_attributes;
        for (const auto& item : query_private_entry) {
            if (item.number == query_node) {
                query_node_attributes.insert(item.tags.begin(), item.tags.end());
            }
        }
        
        // add private attribute
        graph->add_attributes(query_node, query_node_attributes);

        const unordered_set<string>& updated_query_node_attributes = graph->attributes(query_node);
        
        for (const auto& item : query_private_entry) {

            if (item.number == query_node) {
                continue;
            }

            Edge e(query_node, item.number);
            graph->addEdge(query_node, item.number);

            unordered_set<string> filtered_attributes;
            for (const auto& tag : item.tags) {
                if (updated_query_node_attributes.find(tag) != updated_query_node_attributes.end()) {
                    filtered_attributes.insert(tag);
                }
            }

            graph->add_attributes(item.number, filtered_attributes);

            tempEdges.insert(e);
            tempAttributes[item.number].insert(filtered_attributes.begin(), filtered_attributes.end());
        }


        int k = 3;
        printf("query node is: %d, k=%d \n", query_node, k);
        kcore->set_query_node(query_node);
        kcore->set_k(k);
        
        
        auto [duration, max_attributes, community] = kcore->fptree("test/fp_tree_output_" + std::to_string(query_node) + ".txt",
            "test/Output_" + std::to_string(k) + "_fp_" + std::to_string(query_node) + ".txt", allResults);  //fptree

        std::ofstream total_output_file("test/Output_" + std::to_string(k) + "_fp_" + std::to_string(query_node) + ".txt", std::ios::app); 

        total_output_file << "Query Duration: " << duration << " seconds \n";
        total_output_file.flush(); 

        for (const auto& edge : tempEdges) {
            graph->deleteEdge(edge);
        }
        for (const auto& [node, tags] : tempAttributes) {
            graph->remove_attribute(node, tags);
        }

        tempEdges.clear();
        tempAttributes.clear();

        // output_file.close();
    }

    total_output_file.close();
    delete graph;

    return 0;
}
