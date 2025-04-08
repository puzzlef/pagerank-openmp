Design of **OpenMP-based** *PageRank algorithm* for link analysis.

<br>


### Comparing with Ordered approach

**Unordered PageRank** is the *standard* approach of PageRank computation (as
described in the original paper by Larry Page et al. [(1)][page]), where *two*
*different rank vectors* are maintained; one representing the *current* ranks of
vertices, and the other representing the *previous* ranks. On the other hand,
**ordered PageRank** uses *a single rank vector*, representing the current ranks
of vertices [(2)][pagerank]. This is similar to barrierless non-blocking
implementations of the PageRank algorithm by Hemalatha Eedi et al. [(3)][eedi].
As ranks are updated in the same vector (with each iteration), the order in
which vertices are processed *affects* the final result (hence the adjective
*ordered*). However, as PageRank is an iteratively converging algorithm, results
obtained with either approach are *mostly the same*.

In this experiment ([approach-ordered]), we compare the performance of
**ordered** and **unordered OpenMP-based PageRank** (and compare it alongside
*ordered* and *unordered sequential PageRank*). A *schedule* of `dynamic, 2048`
is used for *OpenMP-based PageRank* as obtained in [(4)][pagerank-openmp]. We
use the follwing PageRank parameters: damping factor `α = 0.85`, tolerance
`τ = 10^-6`, and limit the maximum number of iterations to `L = 500.` The error
between the current and the previous iteration is obtained with *L1-norm*, and
is used to detect convergence. *Dead ends* in the graph are handled by always
teleporting any vertex in the graph at random (*teleport* approach [(5)][teleport]).
Error in ranks obtained for each approach is measured relative to the *unordered*
*sequential approach* using *L1-norm*.

From the results, we observe that the **ordered OpenMP-based approach is**
**somewhat faster** than the unordered approach **in terms of time**, and follows
a trend similar to that of sequential PageRank. However, the **ordered**
**approach** (both OpenMP-based and sequential) **converges in significantly fewer**
**iterations** than the unordered approach. This indicates that the ordered
approach could have been quite a bit faster, but is not, because of *some*
overhead (possibly *cache coherence* overhead due to parallel read-write access
to the same vector). In any case, **ordered PageRank** is indeed **faster than**
**unordered Pagerank**.

[approach-ordered]: https://github.com/puzzlef/pagerank-openmp/tree/approach-ordered

<br>


### Adjusting Tolerance (Ordered approach)

In this experiment ([adjust-tolerance-ordered]), we perform *OpenMP-based*
*ordered PageRank* while adjusting the tolerance `τ` from `10^-1` to `10^-14`
with three different tolerance functions: `L1-norm`, `L2-norm`, and `L∞-norm`.
We also compare it with unordered PageRank (both OpenMP-based and sequential)
for the same tolerance and tolerance function. We use a damping factor of
`α = 0.85` and limit the maximum number of iterations to `L = 500`. The error between
the approaches is calculated with *L1-norm*. The *sequential unordered* approach
is considered to be the *gold standard* (wrt to which error is measured). *Dead ends*
in the graph are handled by always teleporting any vertex in the graph at
random (*teleport* approach [(4)]). The teleport contribution to all vertices is
calculated *once* (for all vertices) at the begining of each iteration.

From the results, we observe that **OpenMP-based ordered PageRank** only
converges **faster** than the unordered approach **below a tolerance of**
`τ = 10^-6`. This may be due to *cache coherence overhead* associated with the
ordered approach, which can exceed the benefit provided by ordered approach with
loose tolerance values. In terms of the number of iterations, we interestingly
observe that iterations of OpenMP-based unordered/ordered approaches are higher
than with sequential approaches. We currently do not have an explanation for
this.

[adjust-tolerance-ordered]: https://github.com/puzzlef/pagerank-openmp/tree/adjust-tolerance-ordered

<br>


### Adjusting OpenMP schedule

In this experiment ([adjust-schedule]), we compare performance obtained for
*OpenMP-based PageRank* for various *schedules*. Each thread is assigned a
certain number of *vertices* to process. The **schedule kind** is adjusted among
`static` / `dynamic` / `guided` / `auto`, and the **chunk size** is adjusted
from `1` to `65536`. We do this for the **rank computation step**. PageRank
factors, contributions, and teleport contribution computation is calculated with
suitable OpenMP schedule (`auto`). We use the follwing PageRank parameters:
damping factor `α = 0.85`, tolerance `τ = 10^-6`, and limit the maximum number
of iterations to `L = 500`. The error between the current and the previous
iteration is obtained with *L1-norm*, and is used to detect convergence.

From the results, we observe that a **dynamic schedule with a chunk size of**
**2048** appears to perform the **best**. This however may change based on the
size of graphs in the dataset, or the system used. In such cases `auto` schedule
may be used as a fallback. We also observe that the *difference in ranks*
obtained from sequential and OpenMP-based approach is *relatively high*
(`< 10^-3`) on *large directed graphs*. This may be due to the fact that parallel
reduce performed for teleport contibution calculation differs from sequential
reduce due to *inaccuracies associated with 32-bit floating point format*
*(float)*, and can be avoided by using *64-bit floating point format (double)*.

[adjust-schedule]: https://github.com/puzzlef/pagerank-openmp/tree/adjust-schedule

<br>


### Comparision with Hybrid approach

This experiment ([approach-hybrid]) was for comparing the performance between
finding pagerank using **uniform** OpenMP (*all* routines use OpenMP), or using
**hybrid** OpenMP (*some* routines are *sequential*). Both techniques were
attempted on different types of graphs, running each technique 5 times per graph
to get a good time measure. Number of threads for this experiment (using
`OMP_NUM_THREADS`) was varied from `2` to `48`.

It appears that **hybrid** approach performs **worse** in most cases, and only
slightly better than *uniform* approach in a few cases. I am not sure why
that is the case, possibly there could be some correlation between execution
time and some other parameter. Note that neither approach makes use of
*SIMD instructions* which are available on all modern hardware.

[approach-hybrid]: https://github.com/puzzlef/pagerank-openmp/tree/approach-hybrid

<br>


### Comparision with Sequential implementation

This experiment ([compare-sequential]) was for comparing the performance between
finding pagerank using a single thread (**sequential**), or finding pagerank
accelerated using **OpenMP**. Both techniques were attempted on different types
of graphs, running each technique 5 times per graph to get a good time measure.
Number of threads for this experiment (using `OMP_NUM_THREADS`) was varied from
`2` to `48`.

**OpenMP** does seem to provide a **clear benefit** for most graphs (except for
the smallest ones). This speedup is definitely not directly proportional to the
number of threads, as one would normally expect (Amdahl's law). Note that there
is still room for improvement with **OpenMP** by using sequential versions of
certain routines instead of OpenMP versions because not all calculations benefit
from multiple threads (ex. [vector-multiplication-openmp]). Also note that
neither approach makes use of *SIMD instructions* which are available on all
modern hardware.

[![](https://i.imgur.com/Quuaqnv.gif)][sheets]

[compare-sequential]: https://github.com/puzzlef/pagerank-openmp/tree/compare-sequential

<br>
<br>


## References

- [An Efficient Practical Non-Blocking PageRank Algorithm for Large Scale Graphs; Hemalatha Eedi et al. (2021)](https://ieeexplore.ieee.org/document/9407114)
- [PageRank Algorithm, Mining massive Datasets (CS246), Stanford University](https://www.youtube.com/watch?v=ke9g8hB0MEo)
- [The PageRank Citation Ranking: Bringing Order to the Web; Larry Page et al. (1998)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427)
- [Ranking nodes in growing networks: When PageRank fails; Mariani et al. (2015)](https://www.nature.com/articles/srep16181)
- [Local Approximation of PageRank and Reverse PageRank; Bar-Yossef et al. (2008)](https://static.googleusercontent.com/media/research.google.com/en//pubs/archive/34455.pdf)
- [PageRank Beyond the Web; David F. Gleich (2015)](https://www.cs.purdue.edu/homes/dgleich/publications/Gleich%202015%20-%20prbeyond.pdf)
- [Some methods of speeding up the convergence of iteration methods; Boris T. Polyak (1964)](https://www.researchgate.net/publication/243648538_Some_methods_of_speeding_up_the_convergence_of_iteration_methods)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [When to Stop Slowly Convergent Iteration?; Prof. W. Kahan](https://people.eecs.berkeley.edu/~wkahan/Math128/SlowIter.pdf)
- [Simple Trick to Dramatically Improve Speed of Convergence; Vincent Granville](https://www.datasciencecentral.com/simple-trick-to-dramatically-improve-speed-of-convergence/)
- [What's the difference between "static" and "dynamic" schedule in OpenMP?](https://stackoverflow.com/a/10852852/1413259)
- [OpenMP Dynamic vs Guided Scheduling](https://stackoverflow.com/a/43047074/1413259)
- [Block Compressed Row Format (BSR)](https://scipy-lectures.org/advanced/scipy_sparse/bsr_matrix.html)
- [Aitken's delta-squared process](https://en.wikipedia.org/wiki/Aitken%27s_delta-squared_process)
- [Fixed-point iteration](https://en.wikipedia.org/wiki/Fixed-point_iteration)
- [Steffensen's method](https://en.wikipedia.org/wiki/Steffensen%27s_method)
- [Rate of convergence](https://en.wikipedia.org/wiki/Rate_of_convergence)

<br>
<br>


[![](https://i.imgur.com/5vdxPZ3.jpg)](https://www.youtube.com/watch?v=rKv_l1RnSqs)
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/366356464.svg)](https://zenodo.org/badge/latestdoi/366356464)
![](https://ga-beacon.deno.dev/G-KD28SG54JQ:hbAybl6nQFOtmVxW4if3xw/github.com/puzzlef/pagerank-openmp)

[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://cstar.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://suitesparse-collection-website.herokuapp.com
[graphs]: https://github.com/puzzlef/graphs
[vector-multiplication-openmp]: https://github.com/puzzlef/vector-multiplication-openmp
[charts]: https://photos.app.goo.gl/Bd8bwdZbppkdUQTU9
[sheets]: https://docs.google.com/spreadsheets/d/1Mzmo9KYunJ9yv2ZNwFv73qPjf9VYNaP5YXJT0HVZgpo/edit?usp=sharing
[page]: https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427
[pagerank]: https://github.com/puzzlef/pagerank
[eedi]: https://ieeexplore.ieee.org/document/9407114
[pagerank-openmp]: https://github.com/puzzlef/pagerank-openmp/tree/adjust-schedule
[teleport]: https://gist.github.com/wolfram77/94c38b9cfbf0c855e5f42fa24a8602fc
