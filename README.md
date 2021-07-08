Performance of **uniform-OpenMP** based vs **hybrid-OpenMP** based PageRank
([pull], [CSR], [OpenMP]).

This experiment was for comparing the performance between:
1. Find pagerank using **uniform** OpenMP (*all* routines use OpenMP).
2. Find pagerank using **hybrid** OpenMP (*some* routines are *sequential*).

Both techniques were attempted on different types of graphs, running each
technique 5 times per graph to get a good time measure. Number of threads
for this experiment (using `OMP_NUM_THREADS`) was varied from `2` to `48`.
It appears that **hybrid** approach performs **worse** in most cases, and only
slightly better than *uniform* approach in a few cases. I am not sure why
that is the case, possibly there could be some correlation between execution
time and some other parameter. Note that neither approach makes use of
*SIMD instructions* which are available on all modern hardware.

All outputs are saved in [out](out/) and a small part of the output is listed
here. Some [charts] are also included below, generated from [sheets]. The input
data used for this experiment is available at ["graphs"] (for small ones), and
the [SuiteSparse Matrix Collection]. This experiment was done with guidance
from [Prof. Dip Sankar Banerjee] and [Prof. Kishore Kothapalli].

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
# [00125.334 ms; 062 iters.] [0.0000e+00 err.] pagerankUniform
# [00135.899 ms; 062 iters.] [0.0000e+00 err.] pagerankHybrid
#
# ...
#
# Loading graph /home/subhajit/data/soc-LiveJournal1.mtx ...
# order: 4847571 size: 68993773 {}
# order: 4847571 size: 68993773 {} (transposeWithDegree)
# [04092.295 ms; 050 iters.] [0.0000e+00 err.] pagerankUniform
# [04471.955 ms; 050 iters.] [0.0000e+00 err.] pagerankHybrid
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
# [00065.371 ms; 062 iters.] [0.0000e+00 err.] pagerankUniform
# [00051.436 ms; 062 iters.] [0.0000e+00 err.] pagerankHybrid
#
# ...
#
# Loading graph /home/subhajit/data/soc-LiveJournal1.mtx ...
# order: 4847571 size: 68993773 {}
# order: 4847571 size: 68993773 {} (transposeWithDegree)
# [00805.687 ms; 050 iters.] [0.0000e+00 err.] pagerankUniform
# [01108.808 ms; 050 iters.] [2.5879e-08 err.] pagerankHybrid
#
# ...
```

[![](https://i.imgur.com/pUgQHVw.gif)][sheets]

<br>
<br>


## References

- [PageRank Algorithm, Mining massive Datasets (CS246), Stanford University](https://www.youtube.com/watch?v=ke9g8hB0MEo)
- [SuiteSparse Matrix Collection]

<br>
<br>

[![](https://i.imgur.com/wzUtVOY.jpg)](https://www.youtube.com/watch?v=rKv_l1RnSqs)

[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://cstar.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://suitesparse-collection-website.herokuapp.com
["graphs"]: https://github.com/puzzlef/graphs
[pull]: https://github.com/puzzlef/pagerank-push-vs-pull
[CSR]: https://github.com/puzzlef/pagerank-class-vs-csr
[OpenMP]: https://github.com/puzzlef/pagerank-sequential-vs-openmp
[charts]: https://photos.app.goo.gl/rpK3yyxDNgoUXuzN9
[sheets]: https://docs.google.com/spreadsheets/d/1jE4hwrq9EjO1oPjTdrRiJSkq8Rvxsb8IQ270WsAZG60/edit?usp=sharing
