#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <memory>
#include <cstdlib>
#include <algorithm>

// 定义 FPNode 类
class FPNode {
public:
    std::string item;
    int count = 0;
    int total_count = 0;  // 当前节点的 total_count
    std::shared_ptr<FPNode> parent;
    std::map<std::string, std::shared_ptr<FPNode>> children;

    FPNode(const std::string& item, std::shared_ptr<FPNode> parent) : item(item), parent(parent) {}
};

// 定义 FPTree 类
class FPTree {
public:
    FPTree() : root(std::make_shared<FPNode>("query node", nullptr)) {}

    // 按照输入顺序插入交易
    void insert_transaction(const std::vector<std::string>& transaction) {
        auto current = root;
        for (const auto& item : transaction) {
            // 如果孩子中没有当前项，创建新节点
            if (current->children.find(item) == current->children.end()) {
                current->children[item] = std::make_shared<FPNode>(item, current);
                item_index[item].push_back(current->children[item]);  // 记录 item 在树中的每个节点
            }
            current = current->children[item];
            current->count += 1;  // 增加计数
        }
    }

    // 更新每个节点的 total_count，表示同一个 item 在树中所有位置的 count 总和
    void update_total_counts() {
        for (auto& pair : item_index) {
            int total = 0;
            for (auto& node : pair.second) {
                total += node->count;  // 对每个相同 item 的节点计数进行汇总
            }
            for (auto& node : pair.second) {
                node->total_count = total;  // 设置该 item 的所有节点的 total_count
            }
        }
    }

    // 打印 FP-tree，按 count 降序打印每个节点的子节点，带 total_count
    void print_tree(std::ostream& os) {
        update_total_counts();  // 更新 total_count
        print_tree_helper(root, 0, os);  // 从根节点开始打印
    }

private:
    std::shared_ptr<FPNode> root;
    std::map<std::string, std::vector<std::shared_ptr<FPNode>>> item_index;  // 索引，用于存储每个 item 对应的所有节点

    // 帮助函数：打印树结构
    void print_tree_helper(std::shared_ptr<FPNode> node, int depth, std::ostream& os) {
        for (int i = 0; i < depth; ++i) os << "  ";  // 打印缩进
        os << node->item << " (" << node->count << ", " << node->total_count << ")" << std::endl;

        // 创建一个临时的子节点 vector 用于排序
        std::vector<std::shared_ptr<FPNode>> sorted_children;
        for (const auto& child_pair : node->children) {
            sorted_children.push_back(child_pair.second);
        }

        // 按照 count 进行降序排序
        std::sort(sorted_children.begin(), sorted_children.end(),
                  [](const std::shared_ptr<FPNode>& a, const std::shared_ptr<FPNode>& b) {
                      return a->count > b->count;  // 按照 count 降序排列
                  });

        // 递归打印子节点
        for (const auto& child : sorted_children) {
            print_tree_helper(child, depth + 1, os);
        }
    }
};

// 读取输入文件
void read_input_file(const std::string& filename, std::unordered_map<std::string, std::vector<std::string>>& transactions) {
    std::ifstream infile(filename);
    std::string line;
    
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string attribute;
        std::getline(iss, attribute, ':');
        std::vector<std::string> items;
        std::string temp;

        while (iss >> temp) {
            temp.erase(std::remove(temp.begin(), temp.end(), '#'), temp.end());  // 删除 # 字符
            items.push_back(temp);  // 将交易项目添加到列表中
        }
        transactions[attribute] = items;  // 存储交易记录
    }
}

// 处理文件
void process_file(const std::string& path, const std::string& output_filename) {
    std::unordered_map<std::string, std::vector<std::string>> transactions;
    read_input_file(path, transactions);

    FPTree fptree;
    for (const auto& transaction : transactions) {
        fptree.insert_transaction(transaction.second);  // 按输入顺序插入交易
    }

    std::ofstream outfile(output_filename);
    if (!outfile) {
        std::cerr << "Unable to open output file: " << output_filename << std::endl;
        return;
    }

    outfile << "FP-Tree:" << std::endl;
    fptree.print_tree(outfile);

    std::cout << "FP-Tree has been written to " << output_filename << std::endl;
}

// 主函数
int main() {
    std::string directory = "/Users/beechan/Desktop/ppkcore/case/private/";
    std::string command = "ls " + directory + "nam_output_*.txt > temp_file_list.txt";
    system(command.c_str());  // 生成文件列表

    std::ifstream file_list("temp_file_list.txt");
    std::string filename;

    while (std::getline(file_list, filename)) {
        // 修正路径和文件名的提取
        std::string query_node = filename.substr(filename.find("nam_output_") + 11);
        query_node = query_node.substr(0, query_node.find(".txt"));

        std::string output_filename = directory + "fp_tree_output_" + query_node + ".txt";

        std::cout << "Processing file: " << filename << " for query node: " << query_node << std::endl;
        process_file(filename, output_filename);
    }

    // 删除临时文件
    std::remove("temp_file_list.txt");

    return 0;
}
