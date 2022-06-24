Performance of **sequential** execution based vs **OpenMP** based PageRank ([pull], [CSR]).

This experiment was for comparing the performance between:
1. Find pagerank using a single thread (**sequential**).
2. Find pagerank accelerated using **OpenMP**.

Both techniques were attempted on different types of graphs, running each technique 5 times per graph to get a good time measure. Number of threads for this experiment (using `OMP_NUM_THREADS`) was varied from `2` to `48`. **OpenMP** does seem to provide a **clear benefit** for most graphs (except for the smallest ones). This speedup is definitely not directly proportional to the number of threads, as one would normally expect (Amdahl's law).

Note that there is still room for improvement with **OpenMP** by using sequential versions of certain routines instead of OpenMP versions because not all calculations benefit from multiple threads (ex. ["multiply-sequential-vs-openmp"]). Also note that neither approach makes use of *SIMD instructions* which are available on all modern hardware.

All outputs are saved in [out](out/) and a small part of the output is listed here. Some [charts] are also included below, generated from [sheets]. The input data used for this experiment is available at ["graphs"] (for small ones), and the [SuiteSparse Matrix Collection]. This experiment was done with guidance from [Prof. Dip Sankar Banerjee] and [Prof. Kishore Kothapalli].

<br>

```bash
$ g++ -O3 main.cxx
$ export OMP_NUM_THREADS=4
$ ./a.out ~/data/min-1DeadEnd.mtx
$ ./a.out ~/data/min-2SCC.mtx
$ ...

# ...
#
# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 {}
# order: 281903 size: 2312497 {} (transposeWithDegree)
# [00404.625 ms; 062 iters.] [0.0000e+00 err.] pagerankSeq
# [00128.902 ms; 062 iters.] [2.8684e-07 err.] pagerankOmp
#
# ...
#
# Loading graph /home/subhajit/data/soc-LiveJournal1.mtx ...
# order: 4847571 size: 68993773 {}
# order: 4847571 size: 68993773 {} (transposeWithDegree)
# [14021.927 ms; 050 iters.] [0.0000e+00 err.] pagerankSeq
# [04175.903 ms; 050 iters.] [2.0568e-03 err.] pagerankOmp
#
# ...
```

<br>

```bash
$ g++ -O3 main.cxx
$ export OMP_NUM_THREADS=48
$ ./a.out ~/data/min-1DeadEnd.mtx
$ ./a.out ~/data/min-2SCC.mtx
$ ...

# ...
#
# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 {}
# order: 281903 size: 2312497 {} (transposeWithDegree)
# [00404.945 ms; 062 iters.] [0.0000e+00 err.] pagerankSeq
# [00050.575 ms; 062 iters.] [2.1697e-08 err.] pagerankOmp
#
# ...
#
# Loading graph /home/subhajit/data/soc-LiveJournal1.mtx ...
# order: 4847571 size: 68993773 {}
# order: 4847571 size: 68993773 {} (transposeWithDegree)
# [11737.315 ms; 050 iters.] [0.0000e+00 err.] pagerankSeq
# [00847.098 ms; 050 iters.] [2.0586e-03 err.] pagerankOmp
#
# ...
```

[![](https://i.imgur.com/Quuaqnv.gif)][sheets]

<br>
<br>


## References

- [PageRank Algorithm, Mining massive Datasets (CS246), Stanford University](https://www.youtube.com/watch?v=ke9g8hB0MEo)
- [SuiteSparse Matrix Collection]

<br>
<br>

[![](https://i.imgur.com/5vdxPZ3.jpg)](https://www.youtube.com/watch?v=rKv_l1RnSqs)
[![DOI](https://zenodo.org/badge/366356464.svg)](https://zenodo.org/badge/latestdoi/366356464)


[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://cstar.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://suitesparse-collection-website.herokuapp.com
["graphs"]: https://github.com/puzzlef/graphs
["multiply-sequential-vs-openmp"]: https://github.com/puzzlef/multiply-sequential-vs-openmp
[pull]: https://github.com/puzzlef/pagerank-push-vs-pull
[CSR]: https://github.com/puzzlef/pagerank-class-vs-csr
[charts]: https://photos.app.goo.gl/Bd8bwdZbppkdUQTU9
[sheets]: https://docs.google.com/spreadsheets/d/1Mzmo9KYunJ9yv2ZNwFv73qPjf9VYNaP5YXJT0HVZgpo/edit?usp=sharing
