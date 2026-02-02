#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include "../include/Graph.h"
#include "../include/Kcore.h"

using namespace std;

int main() {

    Graph* graph = new Graph();

    std::string year = "2013";

    cout << "Loading public graph..." << endl;
    graph->load_graph("test/public_graph_" + year + ".txt");
    graph->load_attribute("test/public_attribute_" + year + ".txt");
    cout << "Public graph loaded." << endl;

    Kcore* kcore = new Kcore(graph);
    kcore->read_treeIndex("test/tree_index_" + year + ".txt");

    vector<int> k_list = {3};

    for (int k : k_list) {
        cout << "Processing k = " << k << endl;

        ofstream fp_out(
            "test/node_output.txt");

        // write header
        fp_out << "query_node\t"
            << "time\t"
            << "community_size\t"
            << "attribute_size\n";
        
        // ofstream ab2_out(
        //     "ICDE_revision/" + year + "/k" + to_string(k) + "/improve/ab2/one_node.txt");

        // query list can be read from a file if needed
        ifstream qfile("test/node.txt");
        int q;
        
        while (qfile >> q) {

            kcore->set_query_node(q);
            kcore->set_k(k);

            // PP-FP-tree querying
            auto [fp_time, fp_attr, fp_comm, fp_iter] =
                kcore->fptree(
                    "test/fp_tree_output_" + to_string(q) + ".txt",
                    "test/output_fp_" + to_string(q) + ".txt"
                );

            fp_out << q << "\t"
                   << fp_time << "\t"
                   << fp_comm << "\t"
                   << fp_attr << "\n";

            // Ablation 2 (w/o public expansion)
            // auto [ab2_time, ab2_attr, ab2_comm, ab2_iter] =
            //     kcore->ablation2(
            //         "ICDE_revision/2017/index_250_500/fp_tree_output_" + to_string(q) + ".txt",
            //         "ICDE_revision/2017/k" + to_string(k) + "/improve/ab2/output_ab2_" + to_string(q) + ".txt"
            //     );

            // ab2_out << q << "\t"
            //         << ab2_time << "\t"
            //         << ab2_comm << "\t"
            //         << ab2_attr << "\t"
            //         << ab2_iter << "\n";
        }

        fp_out.close();
        // ab2_out.close();
    }

    delete kcore;
    delete graph;
    return 0;
}

