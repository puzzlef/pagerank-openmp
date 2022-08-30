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
  enum NormFunction { L0=0, L1=1, L2=2, Li=3 };
  vector<T> *init = nullptr;
  float damping = 0.85;

  // Use L1-norm for convergence check.
  for (float tolerance=1e-1; tolerance>=1e-15; tolerance/=10) {
    // Find unordered pagerank using a single thread.
    auto a0 = pagerankMonolithicSeq<false>(x, xt, init, {repeat, L1, damping, tolerance});
    auto e0 = l1Norm(a0.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankSeqUnordered {tol_norm: L1, tolerance: %.0e}\n", a0.time, a0.iterations, e0, tolerance);
    // Find ordered pagerank using a single thread.
    auto a1 = pagerankMonolithicSeq<true>(x, xt, init, {repeat, L1, damping, tolerance});
    auto e1 = l1Norm(a1.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankSeqOrdered   {tol_norm: L1, tolerance: %.0e}\n", a1.time, a1.iterations, e1, tolerance);
    // Find unordered pagerank accelerated with OpenMP.
    auto a2 = pagerankMonolithicOmp<false>(x, xt, init, {repeat, L1, damping, tolerance});
    auto e2 = l1Norm(a2.ranks, a1.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpUnordered {tol_norm: L1, tolerance: %.0e}\n", a2.time, a2.iterations, e2, tolerance);
    // Find ordered pagerank accelerated with OpenMP.
    auto a3 = pagerankMonolithicOmp<true>(x, xt, init, {repeat, L1, damping, tolerance});
    auto e3 = l1Norm(a3.ranks, a1.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpOrdered   {tol_norm: L1, tolerance: %.0e}\n", a3.time, a3.iterations, e3, tolerance);
  }

  // Use L2-norm for convergence check.
  for (float tolerance=1e-1; tolerance>=1e-15; tolerance/=10) {
    // Find unordered pagerank using a single thread.
    auto a0 = pagerankMonolithicSeq<false>(x, xt, init, {repeat, L2, damping, tolerance});
    auto e0 = l1Norm(a0.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankSeqUnordered {tol_norm: L2, tolerance: %.0e}\n", a0.time, a0.iterations, e0, tolerance);
    // Find ordered pagerank using a single thread.
    auto a1 = pagerankMonolithicSeq<true>(x, xt, init, {repeat, L2, damping, tolerance});
    auto e1 = l1Norm(a1.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankSeqOrdered   {tol_norm: L2, tolerance: %.0e}\n", a1.time, a1.iterations, e1, tolerance);
    // Find unordered pagerank accelerated with OpenMP.
    auto a2 = pagerankMonolithicOmp<false>(x, xt, init, {repeat, L2, damping, tolerance});
    auto e2 = l1Norm(a2.ranks, a1.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpUnordered {tol_norm: L2, tolerance: %.0e}\n", a2.time, a2.iterations, e2, tolerance);
    // Find ordered pagerank accelerated with OpenMP.
    auto a3 = pagerankMonolithicOmp<true>(x, xt, init, {repeat, L2, damping, tolerance});
    auto e3 = l1Norm(a3.ranks, a1.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpOrdered   {tol_norm: L2, tolerance: %.0e}\n", a3.time, a3.iterations, e3, tolerance);
  }

  // Use Li-norm for convergence check.
  for (float tolerance=1e-1; tolerance>=1e-15; tolerance/=10) {
    // Find unordered pagerank using a single thread.
    auto a0 = pagerankMonolithicSeq<false>(x, xt, init, {repeat, Li, damping, tolerance});
    auto e0 = l1Norm(a0.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankSeqUnordered {tol_norm: Li, tolerance: %.0e}\n", a0.time, a0.iterations, e0, tolerance);
    // Find ordered pagerank using a single thread.
    auto a1 = pagerankMonolithicSeq<true>(x, xt, init, {repeat, Li, damping, tolerance});
    auto e1 = l1Norm(a1.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankSeqOrdered   {tol_norm: Li, tolerance: %.0e}\n", a1.time, a1.iterations, e1, tolerance);
    // Find unordered pagerank accelerated with OpenMP.
    auto a2 = pagerankMonolithicOmp<false>(x, xt, init, {repeat, Li, damping, tolerance});
    auto e2 = l1Norm(a2.ranks, a1.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpUnordered {tol_norm: Li, tolerance: %.0e}\n", a2.time, a2.iterations, e2, tolerance);
    // Find ordered pagerank accelerated with OpenMP.
    auto a3 = pagerankMonolithicOmp<true>(x, xt, init, {repeat, Li, damping, tolerance});
    auto e3 = l1Norm(a3.ranks, a1.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpOrdered   {tol_norm: Li, tolerance: %.0e}\n", a3.time, a3.iterations, e3, tolerance);
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
