#include <string>
#include <vector>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include "src/main.hxx"

using namespace std;




// You can define datatype with -DTYPE=...
#ifndef TYPE
#define TYPE double
#endif
// You can define number of threads with -DMAX_THREADS=...
#ifndef MAX_THREADS
#define MAX_THREADS 12
#endif




template <class G, class H>
void runPagerank(const G& x, const H& xt, int repeat) {
  using T = TYPE;
  enum NormFunction { L0=0, L1=1, L2=2, Li=3 };
  vector<T> *init = nullptr;
  T damping   = 0.85f;
  T tolerance = 1e-10f;

  for (int sleepDurationMs=1; sleepDurationMs<=1000; sleepDurationMs*=10) {
    for (float sleepProbability=0.0f; sleepProbability<1.01f; sleepProbability+=0.2f) {
      PagerankOptions<T> p = {1,      false, sleepProbability, sleepDurationMs, damping, Li, tolerance};
      PagerankOptions<T> o = {repeat, false, sleepProbability, sleepDurationMs, damping, Li, tolerance};
      // Find pagerank using a single thread for reference (unordered, no dead ends).
      auto a0 = pagerankMonolithicSeq(x, xt, init, p);
      // Find pagerank accelerated with OpenMP (no dead ends, sleep, no help).
      auto a1 = pagerankMonolithicOmp<false, true, false>(x, xt, init, o);
      auto e1 = l1NormOmp(a1.ranks, a0.ranks);
      printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmp                     {sleep_prob: %.1f, sleep_dur: %04d ms}\n", a1.time, a1.iterations, e1, sleepProbability, sleepDurationMs);
      // Find pagerank accelerated with OpenMP (no dead ends, sleep, no help, two vectors).
      auto a2 = pagerankMonolithicBarrierfreeOmp<false, true, false, false>(x, xt, init, o);
      auto e2 = l1NormOmp(a2.ranks, a0.ranks);
      printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreeOmp          {sleep_prob: %.1f, sleep_dur: %04d ms}\n", a2.time, a2.iterations, e2, sleepProbability, sleepDurationMs);
      // Find pagerank accelerated with OpenMP (no dead ends, sleep, no help, one vector).
      auto a3 = pagerankMonolithicBarrierfreeOmp<false, true, false, true>(x, xt, init, o);
      auto e3 = l1NormOmp(a3.ranks, a0.ranks);
      printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreeOneOmp       {sleep_prob: %.1f, sleep_dur: %04d ms}\n", a3.time, a3.iterations, e3, sleepProbability, sleepDurationMs);
      // Find pagerank accelerated with OpenMP (no dead ends, sleep, with help/work-stealing, two vectors).
      auto a4 = pagerankMonolithicBarrierfreeOmp<false, true, true, false>(x, xt, init, o);
      auto e4 = l1NormOmp(a4.ranks, a0.ranks);
      printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankHelperBarrierfreeOmp    {sleep_prob: %.1f, sleep_dur: %04d ms}\n", a4.time, a4.iterations, e4, sleepProbability, sleepDurationMs);
      // Find pagerank accelerated with OpenMP (no dead ends, sleep, with help/work-stealing, one vector).
      auto a5 = pagerankMonolithicBarrierfreeOmp<false, true, true, true>(x, xt, init, o);
      auto e5 = l1NormOmp(a5.ranks, a0.ranks);
      printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankHelperBarrierfreeOneOmp {sleep_prob: %.1f, sleep_dur: %04d ms}\n", a5.time, a5.iterations, e5, sleepProbability, sleepDurationMs);
    }
  }
}


int main(int argc, char **argv) {
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  printf("Loading graph %s ...\n", file);
  OutDiGraph<int64_t> x;
  readMtxW(x, file); println(x);
  auto fl = [](auto u) { return true; };
  selfLoopU(x, None(), fl); print(x); printf(" (selfLoopAllVertices)\n");
  auto xt = transposeWithDegree(x); print(xt); printf(" (transposeWithDegree)\n");
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  runPagerank(x, xt, repeat);
  printf("\n");
  return 0;
}
