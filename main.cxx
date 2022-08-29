#include <string>
#include <vector>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include "src/main.hxx"

using namespace std;




// You can define datatype with -DTYPE=...
#ifndef TYPE
#define TYPE float
#endif
// You can define number of threads with -DMAX_THREADS=...
#ifndef MAX_THREADS
#define MAX_THREADS 12
#endif




template <class G, class H>
void runPagerank(const G& x, const H& xt, int repeat) {
  using T = TYPE;
  vector<T> *init = nullptr;

  // Find pagerank using a single thread.
  auto a1 = pagerankMonolithicSeq<false>(x, xt, init, {repeat});
  auto e1 = l1Norm(a1.ranks, a1.ranks);
  printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankSeq\n", a1.time, a1.iterations, e1);

  // Find pagerank accelerated with OpenMP (static schedule).
  for (int chunkSize=1; chunkSize<=65536; chunkSize*=2) {
    omp_set_schedule(omp_sched_static, chunkSize);
    auto a2 = pagerankMonolithicOmp<false>(x, xt, init, {repeat});
    auto e2 = l1Norm(a2.ranks, a1.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmp {sch_kind: static, chunk_size: %d}\n", a2.time, a2.iterations, e2, chunkSize);
  }

  // Find pagerank accelerated with OpenMP (dynamic schedule).
  for (int chunkSize=1; chunkSize<=65536; chunkSize*=2) {
    omp_set_schedule(omp_sched_dynamic, chunkSize);
    auto a2 = pagerankMonolithicOmp<false>(x, xt, init, {repeat});
    auto e2 = l1Norm(a2.ranks, a1.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmp {sch_kind: dynamic, chunk_size: %d}\n", a2.time, a2.iterations, e2, chunkSize);
  }

  // Find pagerank accelerated with OpenMP (guided schedule).
  for (int chunkSize=1; chunkSize<=65536; chunkSize*=2) {
    omp_set_schedule(omp_sched_guided, chunkSize);
    auto a2 = pagerankMonolithicOmp<false>(x, xt, init, {repeat});
    auto e2 = l1Norm(a2.ranks, a1.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmp {sch_kind: guided, chunk_size: %d}\n", a2.time, a2.iterations, e2, chunkSize);
  }

  // Find pagerank accelerated with OpenMP (auto schedule).
  for (int chunkSize=1; chunkSize<=65536; chunkSize*=2) {
    omp_set_schedule(omp_sched_auto, chunkSize);
    auto a2 = pagerankMonolithicOmp<false>(x, xt, init, {repeat});
    auto e2 = l1Norm(a2.ranks, a1.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmp {sch_kind: auto, chunk_size: %d}\n", a2.time, a2.iterations, e2, chunkSize);
  }
}


int main(int argc, char **argv) {
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  printf("Loading graph %s ...\n", file);
  auto x  = readMtxOutDiGraph(file); println(x);
  // auto fl = [](auto u) { return true; };
  // selfLoopU(x, None(), fl); print(x); printf(" (selfLoopAllVertices)\n");
  auto xt = transposeWithDegree(x);  print(xt); printf(" (transposeWithDegree)\n");
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  runPagerank(x, xt, repeat);
  printf("\n");
  return 0;
}
