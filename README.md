Performance of C++ DiGraph class based vs CSR based PageRank (pull method).

This experiment was for comparing the performance between:
1. Find pagerank using C++ `DiGraph` **class** directly.
2. Find pagerank using **CSR** representation of DiGraph.

Both these approaches were tried on a number of different graphs, running
each approach 5 times per graph to get a good time measure. Using a **CSR**
(Compressed Sparse Row) representation has the potential for performance
improvement for both the methods due to information on vertices and edges
being stored contiguously. See ["pagerank-push-vs-pull"] for a discussion
on *push* vs *pull* method. The input data used for this experiment is
available at ["graphs"] (for small ones), and the
[SuiteSparse Matrix Collection].

```bash
$ g++ -O3 main.cxx
$ ./a.out ~/data/min-1DeadEnd.mtx
$ ./a.out ~/data/min-2SCC.mtx
$ ...

# Loading graph /home/subhajit/data/min-1DeadEnd.mtx ...
# order: 5 size: 6 {}
# order: 5 size: 6 {} (transposeWithDegree)
# [00000.004 ms; 016 iters.] [0.0000e+00 err.] pagerankClass
# [00000.002 ms; 016 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/min-2SCC.mtx ...
# order: 8 size: 12 {}
# order: 8 size: 12 {} (transposeWithDegree)
# [00000.011 ms; 039 iters.] [0.0000e+00 err.] pagerankClass
# [00000.006 ms; 039 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/min-4SCC.mtx ...
# order: 21 size: 35 {}
# order: 21 size: 35 {} (transposeWithDegree)
# [00000.031 ms; 044 iters.] [0.0000e+00 err.] pagerankClass
# [00000.015 ms; 044 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/min-NvgraphEx.mtx ...
# order: 6 size: 10 {}
# order: 6 size: 10 {} (transposeWithDegree)
# [00000.005 ms; 023 iters.] [0.0000e+00 err.] pagerankClass
# [00000.003 ms; 023 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 {}
# order: 281903 size: 2312497 {} (transposeWithDegree)
# [00801.153 ms; 062 iters.] [0.0000e+00 err.] pagerankClass
# [00404.061 ms; 062 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/web-BerkStan.mtx ...
# order: 685230 size: 7600595 {}
# order: 685230 size: 7600595 {} (transposeWithDegree)
# [01480.407 ms; 063 iters.] [0.0000e+00 err.] pagerankClass
# [00884.922 ms; 063 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/web-Google.mtx ...
# order: 916428 size: 5105039 {}
# order: 916428 size: 5105039 {} (transposeWithDegree)
# [03014.658 ms; 060 iters.] [0.0000e+00 err.] pagerankClass
# [01495.231 ms; 060 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/web-NotreDame.mtx ...
# order: 325729 size: 1497134 {}
# order: 325729 size: 1497134 {} (transposeWithDegree)
# [00426.022 ms; 057 iters.] [0.0000e+00 err.] pagerankClass
# [00216.445 ms; 057 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/soc-Slashdot0811.mtx ...
# order: 77360 size: 905468 {}
# order: 77360 size: 905468 {} (transposeWithDegree)
# [00134.193 ms; 054 iters.] [0.0000e+00 err.] pagerankClass
# [00090.104 ms; 054 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/soc-Slashdot0902.mtx ...
# order: 82168 size: 948464 {}
# order: 82168 size: 948464 {} (transposeWithDegree)
# [00149.581 ms; 055 iters.] [0.0000e+00 err.] pagerankClass
# [00099.437 ms; 055 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/soc-Epinions1.mtx ...
# order: 75888 size: 508837 {}
# order: 75888 size: 508837 {} (transposeWithDegree)
# [00109.280 ms; 053 iters.] [0.0000e+00 err.] pagerankClass
# [00080.475 ms; 053 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/coAuthorsDBLP.mtx ...
# order: 299067 size: 1955352 {}
# order: 299067 size: 1955352 {} (transposeWithDegree)
# [00544.151 ms; 044 iters.] [0.0000e+00 err.] pagerankClass
# [00236.277 ms; 044 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/coAuthorsCiteseer.mtx ...
# order: 227320 size: 1628268 {}
# order: 227320 size: 1628268 {} (transposeWithDegree)
# [00398.993 ms; 047 iters.] [0.0000e+00 err.] pagerankClass
# [00189.835 ms; 047 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/soc-LiveJournal1.mtx ...
# order: 4847571 size: 68993773 {}
# order: 4847571 size: 68993773 {} (transposeWithDegree)
# [23022.873 ms; 050 iters.] [0.0000e+00 err.] pagerankClass
# [11676.849 ms; 050 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/coPapersCiteseer.mtx ...
# order: 434102 size: 32073440 {}
# order: 434102 size: 32073440 {} (transposeWithDegree)
# [03421.116 ms; 050 iters.] [0.0000e+00 err.] pagerankClass
# [02182.782 ms; 050 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/coPapersDBLP.mtx ...
# order: 540486 size: 30491458 {}
# order: 540486 size: 30491458 {} (transposeWithDegree)
# [03733.150 ms; 048 iters.] [0.0000e+00 err.] pagerankClass
# [02082.860 ms; 048 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/indochina-2004.mtx ...
# order: 7414866 size: 194109311 {}
# order: 7414866 size: 194109311 {} (transposeWithDegree)
# [26345.023 ms; 060 iters.] [0.0000e+00 err.] pagerankClass
# [18797.385 ms; 060 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/italy_osm.mtx ...
# order: 6686493 size: 14027956 {}
# order: 6686493 size: 14027956 {} (transposeWithDegree)
# [09790.888 ms; 062 iters.] [0.0000e+00 err.] pagerankClass
# [04075.941 ms; 062 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/great-britain_osm.mtx ...
# order: 7733822 size: 16313034 {}
# order: 7733822 size: 16313034 {} (transposeWithDegree)
# [14802.853 ms; 066 iters.] [0.0000e+00 err.] pagerankClass
# [06222.597 ms; 066 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/germany_osm.mtx ...
# order: 11548845 size: 24738362 {}
# order: 11548845 size: 24738362 {} (transposeWithDegree)
# [23849.188 ms; 064 iters.] [0.0000e+00 err.] pagerankClass
# [10025.490 ms; 064 iters.] [0.0000e+00 err.] pagerankCsr
#
# Loading graph /home/subhajit/data/asia_osm.mtx ...
# order: 11950757 size: 25423206 {}
# order: 11950757 size: 25423206 {} (transposeWithDegree)
# [17691.252 ms; 062 iters.] [0.0000e+00 err.] pagerankClass
# [07642.405 ms; 062 iters.] [0.0000e+00 err.] pagerankCsr
```

<br>
<br>


## References

- [PageRank Algorithm, Mining massive Datasets (CS246), Stanford University][this lecture]
- [SuiteSparse Matrix Collection]

<br>
<br>

[![](https://i.imgur.com/MwsC9Av.jpg)](https://www.youtube.com/watch?v=GRvZnN0iiwo)

[this lecture]: http://snap.stanford.edu/class/cs246-videos-2019/lec9_190205-cs246-720.mp4
["pagerank-push-vs-pull"]: https://github.com/puzzlef/pagerank-push-vs-pull
["graphs"]: https://github.com/puzzlef/graphs
[SuiteSparse Matrix Collection]: https://suitesparse-collection-website.herokuapp.com
