Comparison of ordered vs unordered vertex processing in [OpenMP]-based
[PageRank algorithm] for [link analysis].

**Unordered PageRank** is the *standard* approach of PageRank computation (as
described in the original paper by Larry Page et al. [(1)]), where *two*
*different rank vectors* are maintained; one representing the *current* ranks of
vertices, and the other representing the *previous* ranks. On the other hand,
**ordered PageRank** uses *a single rank vector*, representing the current ranks
of vertices [(2)]. This is similar to barrierless non-blocking implementations of
the PageRank algorithm by Hemalatha Eedi et al. [(3)]. As ranks are updated in
the same vector (with each iteration), the order in which vertices are processed
*affects* the final result (hence the adjective *ordered*). However, as PageRank
is an iteratively converging algorithm, results obtained with either approach
are *mostly the same*.

In this experiment, we compare the performance of **ordered** and **unordered**
**OpenMP-based PageRank** (and compare it alongside *ordered* and *unordered*
*sequential PageRank*). A *schedule* of `dynamic, 2048` is used for
*OpenMP-based PageRank* as obtained in [(4)]. We use the follwing PageRank
parameters: damping factor `α = 0.85`, tolerance `τ = 10^-6`, and limit the
maximum number of iterations to `L = 500.` The error between the current and the
previous iteration is obtained with *L1-norm*, and is used to detect
convergence. *Dead ends* in the graph are handled by always teleporting any
vertex in the graph at random (*teleport* approach [(5)]). Error in ranks
obtained for each approach is measured relative to the *unordered sequential*
*approach* using *L1-norm*.

From the results, we observe that the **ordered OpenMP-based approach is**
**somewhat faster** than the unordered approach **in terms of time**, and follows
a trend similar to that of sequential PageRank. However, the **ordered**
**approach** (both OpenMP-based and sequential) **converges in significantly fewer**
**iterations** than the unordered approach. This indicates that the ordered
approach could have been quite a bit faster, but is not, because of *some*
overhead (possibly *cache coherence* overhead due to parallel read-write access
to the same vector). In any case, **ordered PageRank** is indeed **faster than**
**unordered Pagerank**.

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli],
[Prof. Dip Sankar Banerjee], and [Prof. Sathya Peri].

<br>

```bash
$ g++ -std=c++17 -O3 -fopenmp main.cxx
$ ./a.out ~/data/web-Stanford.mtx
$ ./a.out ~/data/web-BerkStan.mtx
$ ...

# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 [directed] {}
# order: 281903 size: 2312497 [directed] {} (transposeWithDegree)
# OMP_NUM_THREADS=12
# [00423.006 ms; 063 iters.] [0.0000e+00 err.] pagerankSeqUnordered
# [00392.738 ms; 033 iters.] [2.6483e-06 err.] pagerankSeqOrdered
# [00048.959 ms; 063 iters.] [2.7737e-07 err.] pagerankOmpUnordered
# [00043.985 ms; 034 iters.] [2.5929e-06 err.] pagerankOmpOrdered
#
# Loading graph /home/subhajit/data/web-BerkStan.mtx ...
# order: 685230 size: 7600595 [directed] {}
# order: 685230 size: 7600595 [directed] {} (transposeWithDegree)
# OMP_NUM_THREADS=12
# [00974.498 ms; 064 iters.] [0.0000e+00 err.] pagerankSeqUnordered
# [00551.588 ms; 035 iters.] [2.3710e-06 err.] pagerankSeqOrdered
# [00107.979 ms; 064 iters.] [3.6591e-06 err.] pagerankOmpUnordered
# [00076.706 ms; 036 iters.] [4.4420e-06 err.] pagerankOmpOrdered
#
# ...
```

[![](https://i.imgur.com/Nl4ijut.png)][sheetp]
[![](https://i.imgur.com/zWBKedV.png)][sheetp]

<br>
<br>


## References

- [An Efficient Practical Non-Blocking PageRank Algorithm for Large Scale Graphs; Hemalatha Eedi et al. (2021)](https://ieeexplore.ieee.org/document/9407114)
- [PageRank Algorithm, Mining massive Datasets (CS246), Stanford University](https://www.youtube.com/watch?v=ke9g8hB0MEo)
- [The PageRank Citation Ranking: Bringing Order to the Web; Larry Page et al. (1998)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [What's the difference between "static" and "dynamic" schedule in OpenMP?](https://stackoverflow.com/a/10852852/1413259)
- [OpenMP Dynamic vs Guided Scheduling](https://stackoverflow.com/a/43047074/1413259)

<br>
<br>


[![](https://i.imgur.com/50yaKL7.jpg)](https://www.youtube.com/watch?v=g2tMcMQqSbA)<br>
[![DOI](https://zenodo.org/badge/530082702.svg)](https://zenodo.org/badge/latestdoi/530082702)


[(1)]: https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427
[(2)]: https://github.com/puzzlef/pagerank-ordered-vs-unordered
[(3)]: https://ieeexplore.ieee.org/document/9407114
[(4)]: https://github.com/puzzlef/pagerank-openmp-adjust-schedule
[(5)]: https://gist.github.com/wolfram77/94c38b9cfbf0c855e5f42fa24a8602fc
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[Prof. Sathya Peri]: https://people.iith.ac.in/sathya_p/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[OpenMP]: https://en.wikipedia.org/wiki/OpenMP
[PageRank algorithm]: https://en.wikipedia.org/wiki/PageRank
[link analysis]: https://en.wikipedia.org/wiki/Network_theory#Link_analysis
[gist]: https://gist.github.com/wolfram77/e58f99e6be8bba28ce5d6a9b45ce276f
[charts]: https://imgur.com/a/jl1Wrkc
[sheets]: https://docs.google.com/spreadsheets/d/1sFMdzsATcjGc9WcS_NxOdOgrMLZyWRiYXtxHibQ4jHw/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vSl_OY4NGqwv9r7a4yXB1-RnkTwlVFrYplaaKGtBk_2Il2dgWhV5sngVHZ3KEj5MDNwFRxRDQ_amNBx/pubhtml
