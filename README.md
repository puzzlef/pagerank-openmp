Effect of varying amounts of random thread sleep with helper variant of
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
approach are *mostly the same*.

**Barrier-free PageRank** is an approach where each thread processes a subset of
vertices in the graph independently, *without* waiting (with a barrier) for
other threads to complete an iteration. This minimizes unnecessary waits and
allows each thread to be on a *different iteration number* (which may or may not
be beneficial for convergence) [(3)]. We have already observed that barrier-free
PageRank with static-scheduling (each thread processing a fixed set of vertices)
is not suitable for performance. This is most likely due to workload imbalance,
which is also likely to affect convergence rate as certian regions of the graph
are iterated over faster than others. A dynamically scheduled barrier-free
PageRank however has good workload balance, and is likely to solve these
problems.

Let us now consider the processing of massive graphs on shared memory systems.
Due to high CPU temperatures, memory access delays, and presence of
heterogeneous memory types (such as [NVRAM]); random thread delays are likely
during the computation of parallel PageRank. The standard approach of parallel
PageRank computation is for each thread to process a subset of vertices of the
graph and wait until all other threads have completed the iteration, before
proceeding onto the next iteration (join). This implicit barrier at the end
approach can usually be a non-issue when we do not anticipate any random thread
sleeps (assuming the chunk size for dynamic schedling is appropriately set).
However, when a random thread sleep occurs in one of the threads, the other
threads must wait for the lazy thread.

One approach to solve this problem is to remove this implicit barrier, and allow
threads to run in different iterations at the same time. This is rarely affects
convergence rate, as OpenMP dynamic scheduling assigns more work to faster
threads (in chunks of certain size), thus keeping threads in the same iteration
in most cases. However, a thread may get stuck (sleep for a long time) and not
process the set of vertices assigned to it (chunk size). Through work-stealing,
other threads can take up work from such stuck threads and complete the
iteration. Work-stealing thus allows all threads to be in the same iteration
even without a barrier, and can be beneficial, when a thread can get stuck.
However, random long-duration thread sleeps can be highly unlikely on most
systems, and work-stealing has additional associated cost. In addition, even
without work-stealing, other (not stuck) threads can still all vertices of the
graph for the next iteration. Given below is details, pros-cons of different
scheduling approaches:

**With OpenMP static scheduling (+barrier)**:
- Each thread is assigned a specific subset of vertices to process (good for cache).
- All threads join at the end of each iteration.
- Even without thread sleep, work will be unbalanced.

**With OpenMP dynamic scheduling (+barrier)**: (used here)
- Each thread picks work from a global pool (bad for cache).
- All threads join at the end of each iteration.
- Without thread sleep, work will be balanced.
- With thread sleep, other threads will be waiting for the lazy thread.

**With Barrier-free/No-sync**:
- Each thread is assigned a specific subset of vertices to process (good for cache).
- Threads do not join at the end of each iteration.
- Each thread can be in a different iteration (covergence rate affected).
- Without thread sleep, work will be initially balanced (not at the end).
- With thread sleep, processing of a subset of vertices delays (convergence affected).

**With OpenMP dynamic scheduling with work-stealing (+barrier)**:
- Each thread picks work from a global pool (bad for cache).
- All threads join at the end of each iteration.
- Free threads steal work from busy/slow/stuck threads.
- Without thread sleep, work will be balanced.
- With thread sleep, lazy threads work will be done by others.
- However, other threads will still wait for lazy thread to wake up (join).
- Additional work stealing overhead involved.

**With dynamic scheduling Barrier-free/No-sync**: (used here)
- Each thread picks work from a global pool (bad for cache).
- Threads do not join at the end of each iteration.
- Each thread can be in a different iteration (covergence rate affected).
- Without thread sleep, work will be balanced.
- With thread sleep, work will mostly be balanced (not at extreme end).

**With dynamic scheduling Barrier-free/No-sync with work-stealing**: (used here)
- Each thread picks work from a global pool (bad for cache).
- Threads do not join at the end of each iteration.
- Each thread will be in the same iteration (due to work stealing).
- Without/with thread sleep, work will be balanced.
- Additional work stealing overhead involved.

In this experiment, we seek to observe the effect of random thread sleeps on the
performance of **plain OpenMP dynamic scheduling (with barrier)**, **dynamic**
**scheduling Barrier-free/No-sync**, and **dynamic scheduling**
**Barrier-free/No-sync with work-stealing** approaches. Sleep can occur before
computing that rank of any vertex (in an iteration) with a certain probability.
This sleep probability is varied from `0.0` to `1.0` in steps of `0.2`. We also
adjust the duration of each sleep from `1 ms` to `1000 ms`. We use a damping
factor of `α = 0.85`, a tolerance of `τ = 10^-10`, and limit the maximum number
of iterations to `L = 500`. Convergence of ranks is determined based of
`L∞-norm` between the ranks of the previous and the current iteration. The error
between the approaches is calculated with *L1-norm*. The *sequential unordered*
approach is considered to be the *gold standard* (wrt to which error is
measured). *Dead ends* in the graph are handled by adding self-loops to all
vertices in the graph (*loopall* approach [(4)]).

From the results, we observe the following. As the sleep duration and sleep
probability is increased, the time taken by **plain OpenMP dynamic scheduling**
**(with barrier)** increases significantly. The time taken by both barrier-free
approaches however increases by a smaller amount. Ordered variants (using a
single rank vector) for both the barrier-free approaches is faster to converge
than unordered variants. **Dynamic** **scheduling Barrier-free/No-sync**
approach is faster than **dynamic scheduling Barrier-free/No-sync with**
**work-stealing** approach, most likely due to the added cost of work-stealing. We
therefore recommend use of **dynamic scheduling Barrier-free/No-sync** approach
in almost all cases.

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
# [00137.761 ms; 092 iters.] [2.9434e-08 err.] pagerankOmp                     {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00156.173 ms; 083 iters.] [5.6706e-09 err.] pagerankBarrierfreeOmp          {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00179.083 ms; 081 iters.] [3.6353e-08 err.] pagerankBarrierfreeOneOmp       {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00195.832 ms; 083 iters.] [5.6975e-09 err.] pagerankHelperBarrierfreeOmp    {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00203.032 ms; 082 iters.] [3.5946e-08 err.] pagerankHelperBarrierfreeOneOmp {sleep_prob: 0.0, sleep_dur: 0001 ms}
# [00145.174 ms; 092 iters.] [2.9434e-08 err.] pagerankOmp                     {sleep_prob: 0.2, sleep_dur: 0001 ms}
# [00157.046 ms; 083 iters.] [5.6560e-09 err.] pagerankBarrierfreeOmp          {sleep_prob: 0.2, sleep_dur: 0001 ms}
# [00187.969 ms; 082 iters.] [3.5937e-08 err.] pagerankBarrierfreeOneOmp       {sleep_prob: 0.2, sleep_dur: 0001 ms}
# [00208.096 ms; 083 iters.] [5.6975e-09 err.] pagerankHelperBarrierfreeOmp    {sleep_prob: 0.2, sleep_dur: 0001 ms}
# [00228.240 ms; 082 iters.] [3.5939e-08 err.] pagerankHelperBarrierfreeOneOmp {sleep_prob: 0.2, sleep_dur: 0001 ms}
# ...
# [60560.363 ms; 092 iters.] [2.9434e-08 err.] pagerankOmp                     {sleep_prob: 1.0, sleep_dur: 1000 ms}
# [08389.505 ms; 095 iters.] [3.2062e-08 err.] pagerankBarrierfreeOmp          {sleep_prob: 1.0, sleep_dur: 1000 ms}
# [08239.562 ms; 089 iters.] [3.5533e-08 err.] pagerankBarrierfreeOneOmp       {sleep_prob: 1.0, sleep_dur: 1000 ms}
# [09050.817 ms; 092 iters.] [2.9434e-08 err.] pagerankHelperBarrierfreeOmp    {sleep_prob: 1.0, sleep_dur: 1000 ms}
# [08269.615 ms; 091 iters.] [3.5544e-08 err.] pagerankHelperBarrierfreeOneOmp {sleep_prob: 1.0, sleep_dur: 1000 ms}
#
# ...
```

[![](https://i.imgur.com/Kb2YVo3.png)][sheetp]
[![](https://i.imgur.com/ZLqieI5.png)][sheetp]
[![](https://i.imgur.com/vhMxBOU.png)][sheetp]
[![](https://i.imgur.com/7jTPhrt.png)][sheetp]

[![](https://i.imgur.com/k60LaFC.png)][sheetp]
[![](https://i.imgur.com/w27vAOf.png)][sheetp]
[![](https://i.imgur.com/VjxUx0P.png)][sheetp]
[![](https://i.imgur.com/WtU75XE.png)][sheetp]

<br>
<br>


## References

- [An Efficient Practical Non-Blocking PageRank Algorithm for Large Scale Graphs; Hemalatha Eedi et al. (2021)](https://ieeexplore.ieee.org/document/9407114)
- [PageRank Algorithm, Mining massive Datasets (CS246), Stanford University](https://www.youtube.com/watch?v=ke9g8hB0MEo)
- [The PageRank Citation Ranking: Bringing Order to the Web; Larry Page et al. (1998)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [Iterating over the array in random order](https://stackoverflow.com/a/18994414/1413259)
- [std::atomic | compare_exchange_weak vs. compare_exchange_strong](https://stackoverflow.com/a/15366555/1413259)
- [Parallel OpenMP loop with break statement](https://stackoverflow.com/a/9836422/1413259)
- [OpenMP: are local variables automatically private?](https://stackoverflow.com/a/30026383/1413259)
- [What's the difference between "static" and "dynamic" schedule in OpenMP?](https://stackoverflow.com/a/10852852/1413259)
- [OpenMP Dynamic vs Guided Scheduling](https://stackoverflow.com/a/43047074/1413259)

<br>
<br>


[![](https://i.imgur.com/K4dU43e.jpg)](https://www.youtube.com/watch?v=GMT18TMNQbY)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)


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
[NVRAM]: https://en.wikipedia.org/wiki/Non-volatile_random-access_memory
[gist]: https://gist.github.com/wolfram77/8a12536a176ba851d8cf42e5030fb01b
[charts]: https://imgur.com/a/4yssfd9
[sheets]: https://docs.google.com/spreadsheets/d/1McPZou_GVoz-vvX0EffnP-QVRWAiis2e9C2AoFmEN1E/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vTFR_tu_dHW4KFZ97GJS_IAJG2lm09aKNdoov2a216PYTCVVKrxdBDUVDwF8-zSqA0dNjoPEgiyFXLo/pubhtml
