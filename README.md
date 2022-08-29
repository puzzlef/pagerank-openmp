Comparing various schedules in [OpenMP]-based [PageRank algorithm] for
[link analysis].

In this experiment, we compare performance obtained for *OpenMP-based PageRank*
for various *schedules*. Each thread is assigned a certain number of *vertices*
to process. The **schedule kind** is adjusted among `static` / `dynamic` /
`guided` / `auto`, and the **chunk size** is adjusted from `1` to `65536`. We do
this for the **rank computation step**. PageRank factors, contributions, and
teleport contribution computation is calculated with suitable OpenMP schedule
(`auto`). We use the follwing PageRank parameters: damping factor `α = 0.85`,
tolerance `τ = 10^-6`, and limit the maximum number of iterations to `L = 500`.
The error between the current and the previous iteration is obtained with
*L1-norm*, and is used to detect convergence.

From the results, we observe that a **dynamic schedule with a chunk size of**
**2048** appears to perform the **best**. This however may change based on the
size of graphs in the dataset, or the system used. In such cases `auto` schedule
may be used as a fallback. We also observe that the *difference in ranks*
obtained from sequential and OpenMP-based approach is *relatively high*
(`< 10^-3`) on *large directed graphs*. This may be due to the fact that parallel
reduce performed for teleport contibution calculation differs from sequential
reduce due to *inaccuracies associated with 32-bit floating point format*
*(float)*, and can be avoided by using *64-bit floating point format (double)*.

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli] and
[Prof. Dip Sankar Banerjee].

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
# [00498.733 ms; 063 iters.] [0.0000e+00 err.] pagerankSeq
# [00116.266 ms; 063 iters.] [2.7737e-07 err.] pagerankOmp {sch_kind: static, chunk_size: 1}
# [00128.507 ms; 063 iters.] [2.7737e-07 err.] pagerankOmp {sch_kind: static, chunk_size: 2}
# [00083.807 ms; 063 iters.] [2.7737e-07 err.] pagerankOmp {sch_kind: static, chunk_size: 4}
# ...
# [00057.358 ms; 063 iters.] [2.7737e-07 err.] pagerankOmp {sch_kind: auto, chunk_size: 16384}
# [00055.914 ms; 063 iters.] [2.7737e-07 err.] pagerankOmp {sch_kind: auto, chunk_size: 32768}
# [00055.303 ms; 063 iters.] [2.7935e-07 err.] pagerankOmp {sch_kind: auto, chunk_size: 65536}
#
# Loading graph /home/subhajit/data/web-BerkStan.mtx ...
# order: 685230 size: 7600595 [directed] {}
# order: 685230 size: 7600595 [directed] {} (transposeWithDegree)
# OMP_NUM_THREADS=12
# [00969.713 ms; 064 iters.] [0.0000e+00 err.] pagerankSeq
# [00266.877 ms; 064 iters.] [3.6591e-06 err.] pagerankOmp {sch_kind: static, chunk_size: 1}
# [00227.865 ms; 064 iters.] [3.6591e-06 err.] pagerankOmp {sch_kind: static, chunk_size: 2}
# [00207.485 ms; 064 iters.] [3.6591e-06 err.] pagerankOmp {sch_kind: static, chunk_size: 4}
# ...
```

[![](https://i.imgur.com/ao04hOr.png)][sheetp]

[![](https://i.imgur.com/MjFYtAR.png)][sheetp]
[![](https://i.imgur.com/SWWzK2n.png)][sheetp]
[![](https://i.imgur.com/btVvbY3.png)][sheetp]
[![](https://i.imgur.com/UCD1JKV.png)][sheetp]
[![](https://i.imgur.com/nHP67ZI.png)][sheetp]
[![](https://i.imgur.com/UNbGaHP.png)][sheetp]
[![](https://i.imgur.com/gIo3FOJ.png)][sheetp]
[![](https://i.imgur.com/hM9QE88.png)][sheetp]
[![](https://i.imgur.com/hCqT5CW.png)][sheetp]
[![](https://i.imgur.com/NJ33mp5.png)][sheetp]
[![](https://i.imgur.com/WiBxzPv.png)][sheetp]
[![](https://i.imgur.com/nhOMpIh.png)][sheetp]
[![](https://i.imgur.com/EEvC3Ng.png)][sheetp]
[![](https://i.imgur.com/SHTYwOn.png)][sheetp]
[![](https://i.imgur.com/XkPS5GP.png)][sheetp]
[![](https://i.imgur.com/zc6BQag.png)][sheetp]

<br>
<br>


## References

- [PageRank Algorithm, Mining massive Datasets (CS246), Stanford University](https://www.youtube.com/watch?v=ke9g8hB0MEo)
- [The PageRank Citation Ranking: Bringing Order to the Web; Larry Page et al. (1998)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [What's the difference between "static" and "dynamic" schedule in OpenMP?](https://stackoverflow.com/a/10852852/1413259)
- [OpenMP Dynamic vs Guided Scheduling](https://stackoverflow.com/a/43047074/1413259)

<br>
<br>


[![](https://i.imgur.com/0XKZ240.jpg)](https://www.bleepingcomputer.com/review/gaming/minecraft-story-mode-is-fun-for-the-whole-family/)<br>


[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[OpenMP]: https://en.wikipedia.org/wiki/OpenMP
[PageRank algorithm]: https://en.wikipedia.org/wiki/PageRank
[link analysis]: https://en.wikipedia.org/wiki/Network_theory#Link_analysis
[gist]: https://gist.github.com/wolfram77/c448cba0c0e371e25c3effa673ea82f2
[charts]: https://imgur.com/a/SwaOufq
[sheets]: https://docs.google.com/spreadsheets/d/1gG5hGrc3o8ztt_eBB7qddy6YLLmOKaDpuQUF3kNzGTc/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vTIiozzISmsE7cm1kHOC-_woCGY6GYh3xcfjovjvfZxMna-Fs9t2vCJ_aq7JIF8-AYlSkhudW7zo3lu/pubhtml
