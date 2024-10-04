#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>
#include <set>
#include "Graph.h"
#include "Utils.h"
#include "Kcore.h"

using namespace std;

int main()
{
    // ofstream total_output_file("dblp2017/6-core/output_6_baselines23.txt"); //baselines
    // ofstream total_output_file("dblp2015/2-core/output_2_fp_combined_baselines.txt"); //fptree
    ofstream total_output_file("case/private/case_output_35914_ppgraph_3core.txt"); //fptree
    // ofstream total_output_file("/Users/beechan/Desktop/ppkcore/scalability/160w"); //fptree //sca


    Graph* graph = new Graph(); 

    //public graph data
    printf("start loading public data...... ");
    graph->load_graph("/Users/beechan/Desktop/ppkcore/data/dblp2017/public_graph.txt");
    graph->load_attribute("/Users/beechan/Desktop/ppkcore/data/dblp2017/public_attribute.txt");
    // graph->load_graph("/Users/beechan/Desktop/ppkcore/scalability/160w/public_graph.txt");
    // graph->load_attribute("/Users/beechan/Desktop/ppkcore/scalability/160w/public_attribute.txt");
    printf("public graph loaded \n");

    //private data
    printf("start loading private index...... ");
    unordered_map<int, vector<PrivateEntry>> private_index;
    // load_private_file("dblp2015/combined_baselines_2015.txt", private_index);  //fptree
    load_private_file("/Users/beechan/Desktop/ppkcore/case/private/35914.txt", private_index);  //fptree
    // load_private_file("data/baselines_test19_1.txt", private_index); //baselines
    printf("private index loaded \n");
   
    //int query_node = 933251;
    int query_node = 35914;
    //int k = 2;
    int k = 3;
    printf("query node is: %d, k=%d \n",query_node,k);

    //construct pp graph
    if (private_index.find(query_node) == private_index.end()) {
        std::cerr << "Query node " << query_node << " not found in private_index." << std::endl;
        return -1; // 或者其他错误处理
    }

    vector<PrivateEntry> query_private_entry = private_index[query_node];
    for (auto item : query_private_entry) {
        graph->addEdge(query_node, item.number);
        graph->add_attributes(item.number, item.tags);
    }

    // set<Edge> tempEdges;
    // unordered_map<int, unordered_set<string>> tempAttributes;

    // for (const auto& pair : private_index)
    // {
    //     int query_node = pair.first;
    //     const vector<PrivateEntry>& query_private_entry = pair.second;

    //     size_t private_edge_count = query_private_entry.size();

    //     // std::ofstream output_file("dblp2013/3-core/output_3_fp_" + std::to_string(query_node) + ".txt");  //fptree
    //     // std::ofstream output_file("dblp2017/6-core/output_6_" + std::to_string(query_node) + ".txt");  //baselines

    //     for (const auto& item : query_private_entry) {
    //         Edge e(query_node, item.number);
    //         graph->addEdge(query_node, item.number);
    //         graph->add_attributes(item.number, item.tags);

    //         tempEdges.insert(e);
    //         tempAttributes[item.number].insert(item.tags.begin(), item.tags.end());
    //     }

        // int k = 4;
        printf("query node is: %d, k=%d \n", query_node, k);
        Kcore* kcore = new Kcore(graph, query_node, k);  //fptree?

        // size_t success_d1, community1, success_d2, community2;    //baselines
        // double duration1 = kcore->baseline1(output_file, success_d1, community1);  //baselines
        // double duration2 = kcore->baseline2(output_file, success_d2, community2);  //baselines
        // auto [duration, max_attributes, community] = kcore->fptree(k, "dblp2013/fp_index_new/fp_tree_output_" + std::to_string(query_node) + ".txt",
        //  "dblp2013/fp_index_new/ncl_output_" + std::to_string(query_node) + ".txt",
        //  "dblp2013/output_2_fp_" + std::to_string(query_node) + ".txt");  //fptree


        try {
            // auto [duration, max_attributes, community] = kcore->fptree(k, "dblp2015/fp_index_new/fp_tree_output_" + std::to_string(query_node) + ".txt",
            // "dblp2015/fp_index_new/ncl_output_" + std::to_string(query_node) + ".txt", "dblp2015/2-core/output_2_fp_" + std::to_string(query_node) + ".txt");  // fptree
            
            auto [duration, max_attributes, community] = kcore->fptree(k, "/Users/beechan/Desktop/ppkcore/case/private/fp_tree_output_" + std::to_string(query_node) + ".txt",
            "/Users/beechan/Desktop/ppkcore/case/private/ncl_output_" + std::to_string(query_node) + ".txt", "/Users/beechan/Desktop/ppkcore/case/private/nam_output_" + std::to_string(query_node) + ".txt");  // fptree
  
            total_output_file << "Query node: " << query_node 
                          // << ", private edges: " << private_edge_count 
                          << ", fptree cost: " << duration 
                          << " seconds, max_attributes: " << max_attributes 
                          << ", community: " << community << "\n";
            //fptree
            
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: std::out_of_range exception caught: " << e.what() << std::endl;
            // 处理异常，如日志记录或其他恢复措施
        } catch (const std::exception& e) {
            std::cerr << "Error: Exception caught: " << e.what() << std::endl;
            // 捕获所有其他可能的异常
        }

        
        // total_output_file << "Query node: " << query_node 
        //                   << ", private edges: " << private_edge_count 
        //                   << ", baseline1 cost: " << duration1 
        //                   << " seconds, max_attr1: " << success_d1 
        //                   << ", community1: " << community1 
        //                   << ", baseline2 cost: " << duration2
        //                   << " seconds, max_attr2: " << success_d2 
        //                   << ", community2: " << community2 << "\n";
        // //baselines

        
        //fptree

        delete kcore;

    //     for (const auto& edge : tempEdges){
    //         graph->deleteEdge(edge);
    //     } 
    //     for (const auto& [node, tags] : tempAttributes){
    //         graph->remove_attribute(node, tags);
    //     }

    //     tempEdges.clear();
    //     tempAttributes.clear();
    //     //public

    //     //output_file.close();

    // }
    
    total_output_file.close();
    delete graph;

    return 0;
}