# Concurrent B-Trees with Buffered Partition Arrays

## Content organization

Relevant code sections:

- `btree_tests`: Scripts for running correctness and performance tests, including microbenchmarks and YCSB workloads.
- `tlx-plain/container/`: B+-tree data structure.
- `tlx-leafds/container/`: BP-tree data structure.
- `btree_tests/ParallelTools`: Submodule for locking mechanisms used in the B-tree and BP-tree.

## Requirements
g++-11, jemalloc

The repo was recently run and compiled on a machine running Ubuntu 20.04.

### Running benchmarks/setup

See instructions here: https://docs.google.com/document/d/1GqQBpYTzSixPAQMfuK8pHyDVL7OjqDHb-Q0u406x894/edit?usp=sharing

If you use BP-tree in an academic context or publication, please cite it as

```
@article{10.14778/3611479.3611502,
author = {Xu, Helen and Li, Amanda and Wheatman, Brian and Marneni, Manoj and Pandey, Prashant},
title = {BP-Tree: Overcoming the Point-Range Operation Tradeoff for In-Memory B-Trees},
year = {2023},
issue_date = {July 2023},
publisher = {VLDB Endowment},
volume = {16},
number = {11},
issn = {2150-8097},
url = {https://doi.org/10.14778/3611479.3611502},
doi = {10.14778/3611479.3611502},
abstract = {B-trees are the go-to data structure for in-memory indexes in databases and storage systems. B-trees support both point operations (i.e., inserts and finds) and range operations (i.e., iterators and maps). However, there is an inherent tradeoff between point and range operations since the optimal node size for point operations is much smaller than the optimal node size for range operations. Existing implementations use a relatively small node size to achieve fast point operations at the cost of range operation throughput.We present the BP-tree, a variant of the B-tree, that overcomes the decades-old point-range operation tradeoff in traditional B-trees. In the BP-tree, the leaf nodes are much larger in size than the internal nodes to support faster range scans. To avoid any slowdown in point operations due to large leaf nodes, we introduce a new insert-optimized array called the buffered partitioned array (BPA) to efficiently organize data in leaf nodes. The BPA supports fast insertions by delaying ordering the keys in the array. This results in much faster range operations and faster point operations at the same time in the BP-tree.Our experiments show that on 48 hyperthreads, on workloads generated from the Yahoo! Cloud Serving Benchmark (YCSB), the BP-tree supports similar or faster point operation throughput (between .94\texttimes{}-1.2\texttimes{} faster) compared to Masstree and OpenBw-tree, two state-of-the-art in-memory key-value (KV) stores. On a YCSB workload with short scans, the BP-tree is about 7.4\texttimes{} faster than Masstree and 1.6\texttimes{} faster than OpenBw-tree. Furthermore, we extend the YCSB to add large range workloads, commonly found in database applications, and show that the BP-tree is 30\texttimes{} faster than Masstree and 2.5\texttimes{} faster than OpenBw-tree.We also provide a reference implementation for a concurrent B+-tree and find that the BP-tree supports faster (between 1.03\texttimes{}-1.2\texttimes{} faster) point operations when compared to the best-case configuration for B+-trees for point operations while supporting similar performance (about .95\texttimes{} as fast) on short range operations and faster (about 1.3\texttimes{} faster) long range operations.},
journal = {Proc. VLDB Endow.},
month = {jul},
pages = {2976–2989},
numpages = {14}
}
```
