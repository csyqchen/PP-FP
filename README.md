# PP-FP search algorithm
**Introduction**

The provided code can be compiled with g++. Ensure all input files are correctly placed in the working directory before execution. The dataset provided is the example shown in Fig. 1 of the paper.

**Input:**

**1.	public_graph.txt**
- Each line is two authors which stands for a co-author relationship in the public graph. (#author_1_id# #author_2_id#)
- Sample:
```html
#1# #2#
#1# #4#
```
**2.	public_attribute.txt**
- Each line is an author with his/her author research interests in the public graph. (#author_name# #author_id# #research_interest#)
- Sample:
```html
#xxx# #1# #Dynamic;Embedding;#
#xxx# #2# #Database;Embedding;#
```
**3.	test_private.txt**
- Represents a private graph associated with a query author q.

  (1) The first line contains two elements. (#query_id# #n#)

  (2) The next n lines represent (n-1) edges in the private graph, each connecting an author id to the query id, along with the private attributes of the query node itself. (#author_name# #author_id# #author_private_interests) 
- Sample:
```html
5 5
#xxx# #5# #Core;Efficient;#
#xxx# #6# #Efficient;#
#xxx# #7# #Core;Efficient;#
#xxx# #8# #Core;Efficient;#
#xxx# #9# #Efficient;#
```
**4.	attribute_index.txt**
- Stores the attribute index of the public graph.
- Each line is a research interest with corresponding authors in the public graph. (research_interest, author_1, author_2, ...)
- Sample:
```html
Database: 2 5 6 7 8 9
Embedding: 1 2 3 4 11 13
```
**5.	tree_index.txt**
- Stores the coreness index of the public graph.
- Each line is a node set with coreness=k. (coreness, author_size, author_1, author_2, ...)
- Sample:
```html
k=0 size=1 nodes: 13
  k=1 size=4 nodes: 1 2 3 5
```
**6.	fp_tree_input_q.txt**
- Contains the PP-FP-tree index generated for a query author q.
- Sample: 
```html
query node (0, 0) {}
  7 (3, 3) {Core, Database, Efficient}
   8 (3, 3) {Core, Database, Efficient}
```
**Output:**

**1.	Output_k_fp_q.txt**
- Provides the following details:
  
   (1) List the authors in the k-core community related to the query author q;

   (2) Size of the k-core community;

   (3) Maximum of common attributes shared by the community;
   
   (4) Query duration time.
- Sample:
```html
Community Node ID: 9 8 7 6 5 
Community Size: 5
Max attribute: 2
Query Duration: 0.000289 seconds 
```
