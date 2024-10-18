# CCAS

This repository contains C++ codes for the paper:

> Practically and Theoretically Efficient Algorithm for 𝑘-Clique Densest Subgraph Discovery

## Introduction

In this paper, we study the problem of k-clique densest subgraph (CDS). k-clique densest subgraph (CDS) problem aims to detect a subgraph from a graph, such that the ratio of the number of k-cliques over its vertices is maximized. This problem has received plenty of attention in the literature, and is widely used in identifying larger ``near-cliques''. We develop an efficient approximation algorithm, by combining the counting version of the best theoretical algorithm and a non-trival graph reduction technique. We have performed extensive experimental evaluation on 12 real-world large graphs and the results demonstrate the high efficiency of our algorithms.

## Environment

The codes of CCAS are implemented and tested under the following development environment:

- Hardware : Intel(R) Xeon(R) Platinum 8358 CPU @ 2.60GHz and 1TB of Memory
- Operation System : Ubuntu 20.04.4 LTS (GNU/Linux 5.13.0-40-generic x86_64)
- C++ Version: 14

## Datasets
We use 12 real-world datasets from different domains, which are downloaded from the [Stanford Network Analysis Platform](http://snap.stanford.edu/data/), [Laboratory of Web Algorithmics](http://law.di.unimi.it/datasets.php), [Network Repository](https://networkrepository.com/network-data.php), and [Networks](http://konect.cc/networks/). 

The statistical information is as follows, and you can click on the `Name` to navigate to the download link.

| Name | \|V\| | \|E\| | Kmax |
| :----: | :----: | :----: | :----: |
| [bio-SC-GT](https://networkrepository.com/bio-SC-GT.php) | 1,716 | 31,564 | 48 |
| [econ-beacxc](https://networkrepository.com/econ-beacxc.php) | 507 | 42,176 | 87 |
| [WikiTalk](https://snap.stanford.edu/data/wiki-Talk.html) | 120,834 | 237,551 | 27 |
| [Slashdot](http://konect.cc/networks/slashdot-zoo) | 77,360 | 469,180 | 26 |
| [DBLP](http://snap.stanford.edu/data/com-DBLP.html) | 317,080 | 1,049,866 | 114 |
| [HepTh](https://networkrepository.com/cit-HepTh.php) | 22,908 | 2,444,798 | 561 |
| [Hollywood](https://networkrepository.com/ca-hollywood-2009.php) | 1,069,126 | 56,306,652 | 2,209 |
| [zhishi-baidu](http://konect.cc/networks/zhishi-all) | 7,827,193 | 62,246,014 | 268 |
| [UK-2002](https://law.di.unimi.it/webdata/uk-2002/) | 18,483,190 | 261,787,260 | 944 |
| [Arabic-2005](https://law.di.unimi.it/webdata/arabic-2005/) | 22,743,892 | 553,903,073 | 3,248 |
| [IT-2004](https://law.di.unimi.it/webdata/it-2004/) | 41,290,648 | 1,027,474,947 | 3,224 |
| [Friendster](https://snap.stanford.edu/data/com-Friendster.html) | 124,836,180 | 1,806,067,135 | 1 |


## How to Run the Codes


The implementations of both `CCAS` and `CCAS-A` algorithms are included in the `ccas.cpp` and `ccas-a.cpp`. 

You can compile `ccas.cpp` in `CCAS/`:

`g++ -Ofast -mcmodel=medium -march=native -std=c++14 -mavx ccas.cpp graph.cpp pava.cpp sct.cpp heap.cpp -o  ccas`

You can compile `ccas-a.cpp` in `CCAS/`:

`g++ -Ofast -mcmodel=medium -march=native -std=c++14 -mavx ccas-a.cpp graph.cpp pava.cpp sct.cpp heap.cpp -o  ccas-a`

and run those file by :

`./ccas <graph-path> <k> <max_iter> <eps> <kmax>`

`./ccas-a <graph-path> <k> <max_iter> <eps> <kmax>`

The application of these arguments is as follows.

`graph-path` is the path to the file that contains the graph data. The graph file should be in txt format and organized as follows:

The first line of this file should contain two integers $n$ (the number of vertices) and $m$ (the number of edges) separated by a space. The index of vertices should be ranged from $0$ to $n - 1$.  Each of the following $m$ lines contains one edge by two integers separated by a space. There should be no multiple edges and self-loops in this file.

For example, here is a correct graph file with 4 vertices and 5 edges 

```
4 5
0 1
1 2
2 3
3 0
1 3
```

`k` is the k selected for k-cliques, for `CCAS-A` it will caculate CDS from 3 to k.

`max_iter` is a number of iterations the program takes. If you want the program to stop by a given $\epsilon$, please set this as 0.

`eps` is the threshold for $\epsilon$. The program will terminate itself once it reaches the $\epsilon$. This should be set as a float number greater than 0. To use this argument, please set `max_iter` as zero.

`kmax` is a parameter used for our graph reduction techniques `HCGR`, it can be set as any non-zero integer values not greater than the maximum clique size of the given graph. To achieve the best performance, please set this as the maximum clique size. The maximum clique size for all datasets we used in our experiment have been included in the table above.

For example, if we want to run the CCAS algorithm to approximate 5-cliques densest subgraph from the graph in `dblp.txt` taking 10 iterations. The command should be:

`./ccas dblp.txt 5 10 0 114`

If we want to run the CCAS-A algorithm to approximate cliques densest subgraph for all k values from the graph in `dblp.txt` by an $\epsilon$ threshold of 0.01 . The command should be:

`./ccas-a dblp.txt 114 0 0.01 114`
