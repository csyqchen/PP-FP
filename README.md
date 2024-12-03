# PP-FP

Input:
1.	public_graph.txt
- Contains the connection relationships in the public graph.
2.	public_attribute.txt
- Contains the attributes of nodes in the public graph.
3.	test_private.txt
- Represents a private graph associated with a query node q.
4.	attribute_index_2016.txt
- Stores the attribute index of the public graph.
5.	tree_index_2016.txt
- Stores the coreness index of the public graph.
6.	fp_tree_input_q.txt
- Contains the FP-tree index generated for a query node q.
Output:
1.	fp_k_output.txt
- Provides the following details:
- Query node q;
- Query duration time;
- Size of the k-core community;
- Maximum of common attributes shared by the community.
2.	Output_k_fp_q.txt
- List the nodes in the k-core community.

The provided code can be compiled with g++. Ensure all input files are correctly placed in the working directory before execution.
