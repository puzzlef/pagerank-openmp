Effect of using different values of tolerance with [OpenMP]-based ordered
[PageRank algorithm] for [link analysis].

**Unordered PageRank** is the *usual* method of computing PageRank (as given in
the original PageRank paper by Larry Page et al. [(1)]), where *two rank*
*vectors* are maintained; one denotes the *current* ranks of vertices, and the
other denotes the *previous* ranks. On the contrary, **ordered PageRank** uses
only *one rank vector*, denoting the current ranks [(2)]. This is similar to
barrierless non-blocking PageRank implementations by Hemalatha Eedi et al.
[(3)]. As ranks are updated in the same vector (with each iteration), the order
in which ranks of vertices are computed *affects* the final consequence (hence
the adjective *ordered*). However, PageRank is an iteratively converging
algorithm, rand thus results obtained with either approach are *generally the*
*same*.

In this experiment, we perform *OpenMP-based ordered PageRank* while adjusting
the tolerance `τ` from `10^-1` to `10^-14` with three different tolerance
functions: `L1-norm`, `L2-norm`, and `L∞-norm`. We also compare it with
unordered PageRank (both OpenMP-based and sequential) for the same tolerance and
tolerance function. We use a damping factor of `α = 0.85` and limit the maximum
number of iterations to `L = 500`. The error between the approaches is
calculated with *L1-norm*. The *sequential unordered* approach is considered to
be the *gold standard* (wrt to which error is measured). *Dead ends* in the
graph are handled by always teleporting any vertex in the graph at random
(*teleport* approach [(4)]). The teleport contribution to all vertices is
calculated *once* (for all vertices) at the begining of each iteration.

From the results, we observe that **OpenMP-based ordered PageRank** only
converges **faster** than the unordered approach **below a tolerance of**
`τ = 10^-6`. This may be due to *cache coherence overhead* associated with the
ordered approach, which can exceed the benefit provided by ordered approach with
loose tolerance values. In terms of the number of iterations, we interestingly
observe that iterations of OpenMP-based unordered/ordered approaches are higher
than with sequential approaches. We currently do not have an explanation for
this.

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
# [00031.683 ms; 005 iters.] [0.0000e+00 err.] pagerankSeqUnordered {tol_norm: L1, tolerance: 1e-01}
# [00041.039 ms; 004 iters.] [6.4879e-02 err.] pagerankSeqOrdered   {tol_norm: L1, tolerance: 1e-01}
# [00005.483 ms; 005 iters.] [6.4879e-02 err.] pagerankOmpUnordered {tol_norm: L1, tolerance: 1e-01}
# [00006.194 ms; 004 iters.] [1.1537e-02 err.] pagerankOmpOrdered   {tol_norm: L1, tolerance: 1e-01}
# [00074.416 ms; 012 iters.] [0.0000e+00 err.] pagerankSeqUnordered {tol_norm: L1, tolerance: 1e-02}
# [00070.314 ms; 007 iters.] [1.2687e-02 err.] pagerankSeqOrdered   {tol_norm: L1, tolerance: 1e-02}
# [00010.354 ms; 012 iters.] [1.2687e-02 err.] pagerankOmpUnordered {tol_norm: L1, tolerance: 1e-02}
# ...
# [03071.110 ms; 500 iters.] [0.0000e+00 err.] pagerankSeqUnordered {tol_norm: Li, tolerance: 1e-13}
# [04773.107 ms; 500 iters.] [1.6800e-07 err.] pagerankSeqOrdered   {tol_norm: Li, tolerance: 1e-13}
# [00360.719 ms; 500 iters.] [3.0579e-07 err.] pagerankOmpUnordered {tol_norm: Li, tolerance: 1e-13}
# [00588.157 ms; 500 iters.] [2.7384e-07 err.] pagerankOmpOrdered   {tol_norm: Li, tolerance: 1e-13}
# [03050.820 ms; 500 iters.] [0.0000e+00 err.] pagerankSeqUnordered {tol_norm: Li, tolerance: 1e-14}
# [04781.842 ms; 500 iters.] [1.6800e-07 err.] pagerankSeqOrdered   {tol_norm: Li, tolerance: 1e-14}
# [00357.327 ms; 500 iters.] [3.0517e-07 err.] pagerankOmpUnordered {tol_norm: Li, tolerance: 1e-14}
# [00581.229 ms; 500 iters.] [2.7051e-07 err.] pagerankOmpOrdered   {tol_norm: Li, tolerance: 1e-14}
#
# Loading graph /home/subhajit/data/web-BerkStan.mtx ...
# order: 685230 size: 7600595 [directed] {}
# order: 685230 size: 7600595 [directed] {} (transposeWithDegree)
# OMP_NUM_THREADS=12
# [00070.303 ms; 005 iters.] [0.0000e+00 err.] pagerankSeqUnordered {tol_norm: L1, tolerance: 1e-01}
# [00058.766 ms; 004 iters.] [9.5628e-02 err.] pagerankSeqOrdered   {tol_norm: L1, tolerance: 1e-01}
# [00010.088 ms; 005 iters.] [9.5627e-02 err.] pagerankOmpUnordered {tol_norm: L1, tolerance: 1e-01}
# [00008.914 ms; 004 iters.] [7.0768e-02 err.] pagerankOmpOrdered   {tol_norm: L1, tolerance: 1e-01}
# [00165.733 ms; 012 iters.] [0.0000e+00 err.] pagerankSeqUnordered {tol_norm: L1, tolerance: 1e-02}
# [00115.817 ms; 008 iters.] [1.9341e-02 err.] pagerankSeqOrdered   {tol_norm: L1, tolerance: 1e-02}
# [00020.469 ms; 012 iters.] [1.9340e-02 err.] pagerankOmpUnordered {tol_norm: L1, tolerance: 1e-02}
# [00015.728 ms; 008 iters.] [1.7850e-02 err.] pagerankOmpOrdered   {tol_norm: L1, tolerance: 1e-02}
# ...
```

[![](https://i.imgur.com/DXQtjER.png)][sheetp]
[![](https://i.imgur.com/f4OtE4b.png)][sheetp]

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


[![](https://i.imgur.com/qp7YIhe.jpg)](https://www.youtube.com/watch?v=69-J2m_GyhI)<br>
[![DOI](https://zenodo.org/badge/530790127.svg)](https://zenodo.org/badge/latestdoi/530790127)


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
[gist]: https://gist.github.com/wolfram77/12aef056d47ecadf05dee0fec4918021
[charts]: https://imgur.com/a/AC6gRGQ
[sheets]: https://docs.google.com/spreadsheets/d/1cXSdBuMwdhIlN1ufXdl_tvzeU2DeVIAI9efU0zQRBPM/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vSeuJwLt-RUwAn70YXDFx5soAjY2ikgySDYoFx8vOeB49d7INSRECTJnEsrLWQXstyQwE_lCA3aPQVL/pubhtml
