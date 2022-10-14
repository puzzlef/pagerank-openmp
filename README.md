Effect of varying amounts of random thread sleep with barrier-free iterations in
[OpenMP]-based ordered [PageRank algorithm] for [link analysis].

**Unordered PageRank** is the *standard* approach of PageRank computation (as
described in the original paper by Larry Page et al. [(1)]), where *two*
*different rank vectors* are maintained; one representing the *current* ranks of
vertices, and the other representing the *previous* ranks. On the other hand,
**ordered PageRank** uses *a single rank vector*, representing the current ranks
of vertices [(2)]. This is similar to barrierless non-blocking implementations
of the PageRank algorithm by Hemalatha Eedi et al. [(3)]. As ranks are updated
in the same vector (with each iteration), the order in which vertices are
processed *affects* the final result (hence the adjective *ordered*). However,
as PageRank is an iteratively converging algorithm, results obtained with either
approach are *mostly the same*. **Barrier-free PageRank** is an *ordered*
*PageRank* where each thread processes a subset of vertices in the graph
independently, *without* waiting (with a barrier) for other threads to complete an
iteration. This minimizes unnecessary waits and allows each thread to be on a
*different iteration number* (which may or may not be beneficial for convergence)
[(3)].

In this experiment, we seek to observe the effect of random thread sleeps on the
performance of different approaches of PageRank computation. Sleep can occur
before computing that rank of any vertex (in an iteration) with a certain
probability. This sleep probability is varied from `0.0` to `1.0` in steps of
`0.1`. We also adjust the duration of each sleep from `1 ms` to `1000 ms`. We
perform two different approaches of barrier-free iterations of *OpenMP-based*
*ordered PageRank*; one in which each thread detects convergence by measuring the
difference between the previous and the current ranks of all the vertices
(**full**), and the other in which the difference is measured between the
previous and current ranks of only the subset of vertices being processed by
each thread (**part**). We compare them with OpenMP-based unordered and ordered
PageRank. We use a damping factor of `α = 0.85`, a tolerance of `τ = 10^-10`,
and limit the maximum number of iterations to `L = 500`. Convergence of ranks is
determined based of `L∞-norm` between the ranks of the previous and the current
iteration. The error between the approaches is calculated with *L1-norm*. The
*sequential unordered* approach is considered to be the *gold standard* (wrt to
which error is measured). *Dead ends* in the graph are handled by adding
self-loops to all vertices in the graph (*loopall* approach [(4)]).

From the results, we observe the following. **OpenMP-based unordered PageRank**
is **always faster or the same performance** as barrier-free (partial/full error
measurement) PageRank. A **random thread sleep with low probability** (`0.1`)
**can significantly increase the time required for convergence**, but increasing
this sleep probability does not significantly affect the convergence time
further. With respect to the **number of iterations**, **OpenMP-based**
**unordered/ordered approaches** are **not impacted at all by random thread**
**sleeps**, but barrier-free approaches and are affected (usually peaking at sleep
probability 0.1 and dropping and the probability increases).

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
# order: 281903 size: 2594400 [directed] {} (selfLoopAllVertices)
# order: 281903 size: 2594400 [directed] {} (transposeWithDegree)
# OMP_NUM_THREADS=12
# [00139.434 ms; 092 iters.] [2.9434e-08 err.] pagerankOmpUnordered       {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00202.637 ms; 087 iters.] [3.5569e-08 err.] pagerankOmpOrdered         {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00213.090 ms; 090 iters.] [3.5820e-08 err.] pagerankBarrierfreeFullOmp {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00157.038 ms; 068 iters.] [6.3337e-08 err.] pagerankBarrierfreePartOmp {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00143.502 ms; 092 iters.] [2.9434e-08 err.] pagerankOmpUnordered       {sleep_prob: 0.1, sleep_dur: 0001 ms}
# [00203.554 ms; 087 iters.] [3.5567e-08 err.] pagerankOmpOrdered         {sleep_prob: 0.1, sleep_dur: 0001 ms}
# [00200.269 ms; 080 iters.] [3.7927e-08 err.] pagerankBarrierfreeFullOmp {sleep_prob: 0.1, sleep_dur: 0001 ms}
# [00159.348 ms; 068 iters.] [6.5084e-08 err.] pagerankBarrierfreePartOmp {sleep_prob: 0.1, sleep_dur: 0001 ms}
# [00144.523 ms; 092 iters.] [2.9434e-08 err.] pagerankOmpUnordered       {sleep_prob: 0.2, sleep_dur: 0001 ms}
# [00204.447 ms; 087 iters.] [3.5567e-08 err.] pagerankOmpOrdered         {sleep_prob: 0.2, sleep_dur: 0001 ms}
# [00207.622 ms; 079 iters.] [3.7733e-08 err.] pagerankBarrierfreeFullOmp {sleep_prob: 0.2, sleep_dur: 0001 ms}
# [00173.303 ms; 068 iters.] [6.0156e-08 err.] pagerankBarrierfreePartOmp {sleep_prob: 0.2, sleep_dur: 0001 ms}
# ...
# [00132.218 ms; 092 iters.] [2.9434e-08 err.] pagerankOmpUnordered       {sleep_prob: 0.0, sleep_dur: 0005 ms}
# [00194.745 ms; 087 iters.] [3.5571e-08 err.] pagerankOmpOrdered         {sleep_prob: 0.0, sleep_dur: 0005 ms}
# [00198.631 ms; 086 iters.] [3.6316e-08 err.] pagerankBarrierfreeFullOmp {sleep_prob: 0.0, sleep_dur: 0005 ms}
# [00151.742 ms; 068 iters.] [6.3871e-08 err.] pagerankBarrierfreePartOmp {sleep_prob: 0.0, sleep_dur: 0005 ms}
# [00181.450 ms; 092 iters.] [2.9434e-08 err.] pagerankOmpUnordered       {sleep_prob: 0.1, sleep_dur: 0005 ms}
# [00233.075 ms; 087 iters.] [3.5567e-08 err.] pagerankOmpOrdered         {sleep_prob: 0.1, sleep_dur: 0005 ms}
# [00249.260 ms; 081 iters.] [3.8378e-08 err.] pagerankBarrierfreeFullOmp {sleep_prob: 0.1, sleep_dur: 0005 ms}
# [00200.137 ms; 070 iters.] [1.1288e-07 err.] pagerankBarrierfreePartOmp {sleep_prob: 0.1, sleep_dur: 0005 ms}
# ...
# [58360.914 ms; 092 iters.] [2.9434e-08 err.] pagerankOmpUnordered       {sleep_prob: 1.0, sleep_dur: 1000 ms}
# [57011.207 ms; 087 iters.] [3.5569e-08 err.] pagerankOmpOrdered         {sleep_prob: 1.0, sleep_dur: 1000 ms}
# [96603.320 ms; 089 iters.] [5.0634e-08 err.] pagerankBarrierfreeFullOmp {sleep_prob: 1.0, sleep_dur: 1000 ms}
# [91684.117 ms; 082 iters.] [3.2389e-07 err.] pagerankBarrierfreePartOmp {sleep_prob: 1.0, sleep_dur: 1000 ms}
#
# Loading graph /home/subhajit/data/web-BerkStan.mtx ...
# order: 685230 size: 7600595 [directed] {}
# order: 685230 size: 8285825 [directed] {} (selfLoopAllVertices)
# order: 685230 size: 8285825 [directed] {} (transposeWithDegree)
# OMP_NUM_THREADS=12
# [00242.956 ms; 089 iters.] [3.4431e-08 err.] pagerankOmpUnordered       {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00260.992 ms; 093 iters.] [4.3849e-08 err.] pagerankOmpOrdered         {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00287.529 ms; 096 iters.] [4.4476e-08 err.] pagerankBarrierfreeFullOmp {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00193.964 ms; 071 iters.] [1.1348e-07 err.] pagerankBarrierfreePartOmp {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00243.095 ms; 089 iters.] [3.4431e-08 err.] pagerankOmpUnordered       {sleep_prob: 0.1, sleep_dur: 0001 ms}
# [00267.031 ms; 093 iters.] [4.3864e-08 err.] pagerankOmpOrdered         {sleep_prob: 0.1, sleep_dur: 0001 ms}
# [00300.417 ms; 095 iters.] [4.3921e-08 err.] pagerankBarrierfreeFullOmp {sleep_prob: 0.1, sleep_dur: 0001 ms}
# [00210.418 ms; 071 iters.] [1.1081e-07 err.] pagerankBarrierfreePartOmp {sleep_prob: 0.1, sleep_dur: 0001 ms}
# ...
```

[![](https://i.imgur.com/U4RORcK.png)][sheetp]
[![](https://i.imgur.com/0OGRgGS.png)][sheetp]
[![](https://i.imgur.com/OrZ8uq0.png)][sheetp]
[![](https://i.imgur.com/bCTjZu4.png)][sheetp]

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


[![](https://i.imgur.com/K4dU43e.jpg)](https://www.youtube.com/watch?v=GMT18TMNQbY)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/542447553.svg)](https://zenodo.org/badge/latestdoi/542447553)


[(1)]: https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427
[(2)]: https://github.com/puzzlef/pagerank-ordered-vs-unordered
[(3)]: https://ieeexplore.ieee.org/document/9407114
[(4)]: https://gist.github.com/wolfram77/94c38b9cfbf0c855e5f42fa24a8602fc
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[Prof. Sathya Peri]: https://people.iith.ac.in/sathya_p/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[OpenMP]: https://en.wikipedia.org/wiki/OpenMP
[PageRank algorithm]: https://en.wikipedia.org/wiki/PageRank
[link analysis]: https://en.wikipedia.org/wiki/Network_theory#Link_analysis
[gist]: https://gist.github.com/wolfram77/ad0dc8fa38ef6719a93ea9b830ee9ca8
[charts]: https://imgur.com/a/H7ICGOL
[sheets]: https://docs.google.com/spreadsheets/d/1ZkTeBDk-pbQRoj9x8-EPge_Uo8BSVg0dkOiYUSgCyd4/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vScCmQLPcsbW-Zlhf57BuJ50gu5h4zVFwWzwbmkvb2zQOciFuS-tIID3LKG-zbiKu0T3DKncha0IZof/pubhtml
