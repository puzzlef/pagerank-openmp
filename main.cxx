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

  // Find pagerank using a single thread (unordered).
  auto a1 = pagerankMonolithicSeq<false>(x, xt, init, {repeat});
  auto e1 = l1Norm(a1.ranks, a1.ranks);
  printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankSeqUnordered\n", a1.time, a1.iterations, e1);

  // Find pagerank using a single thread (ordered).
  auto a2 = pagerankMonolithicSeq<true>(x, xt, init, {repeat});
  auto e2 = l1Norm(a2.ranks, a1.ranks);
  printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankSeqOrdered\n", a2.time, a2.iterations, e2);

  // Find pagerank accelerated with OpenMP (unordered).
  auto a3 = pagerankMonolithicOmp<false>(x, xt, init, {repeat});
  auto e3 = l1Norm(a3.ranks, a1.ranks);
  printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpUnordered\n", a3.time, a3.iterations, e3);

  // Find pagerank accelerated with OpenMP (ordered).
  auto a4 = pagerankMonolithicOmp<true>(x, xt, init, {repeat});
  auto e4 = l1Norm(a4.ranks, a1.ranks);
  printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpOrdered\n", a4.time, a4.iterations, e4);
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
