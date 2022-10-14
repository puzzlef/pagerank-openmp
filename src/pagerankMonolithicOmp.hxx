#pragma once
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "transpose.hxx"
#include "dynamic.hxx"
#include "pagerank.hxx"
#include "pagerankSeq.hxx"
#include "pagerankOmp.hxx"

using std::vector;
using std::swap;




// PAGERANK-LOOP
// -------------

template <bool DEAD=false, bool SLEEP=false, bool HELP=false, class K, class T>
int pagerankMonolithicOmpLoopU(vector<T>& a, vector<T>& r, vector<T>& c, const vector<T>& f, const vector<K>& vfrom, const vector<K>& efrom, const vector<K>& vdata, vector<PagerankThreadWork*>& works, K i, K n, K N, T p, T E, int L, int EF, float SP, int SD) {
  int l = 0;
  while (l<L) {
    T c0 = DEAD? pagerankTeleportOmp(r, vdata, N, p) : (1-p)/N;
    if (HELP) pagerankCalculateHelperOmpW<SLEEP>(a, c, vfrom, efrom, i, n, c0, SP, SD, works);
    else pagerankCalculateOmpW<SLEEP>(a, c, vfrom, efrom, i, n, c0, SP, SD, works);  // update ranks of vertices
    multiplyValuesOmpW(c, a, f, i, n);        // update partial contributions (c)
    T el = pagerankErrorOmp(a, r, i, n, EF);  // compare previous and current ranks
    swap(a, r); ++l;                          // final ranks in (r)
    if (el<E) break;                          // check tolerance
  }
  return l;
}




// PAGERANK (STATIC / INCREMENTAL)
// -------------------------------

// Find pagerank using multiple threads (pull, CSR).
// @param x  original graph
// @param xt transpose graph (with vertex-data=out-degree)
// @param q  initial ranks (optional)
// @param o  options {damping=0.85, tolerance=1e-6, maxIterations=500}
// @returns {ranks, iterations, time}
template <bool DEAD=false, bool SLEEP=false, bool HELP=false, class G, class H, class T=float>
PagerankResult<T> pagerankMonolithicOmp(const G& x, const H& xt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  using K = typename G::key_type;
  K    N  = xt.order();  if (N==0) return PagerankResult<T>::initial(xt, q);
  auto ks = pagerankVertices(x, xt, o, C);
  return pagerankOmp(xt, ks, K(0), N, pagerankMonolithicOmpLoopU<DEAD, SLEEP, HELP, K, T>, q, o);
}

template <bool DEAD=false, bool SLEEP=false, bool HELP=false, class G, class T=float>
PagerankResult<T> pagerankMonolithicOmp(const G& x, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  auto xt = transposeWithDegree(x);
  return pagerankMonolithicOmp<DEAD, SLEEP, HELP>(x, xt, q, o, C);
}




// PAGERANK (DYNAMIC)
// ------------------

template <bool DEAD=false, bool SLEEP=false, bool HELP=false, class G, class H, class T=float>
PagerankResult<T> pagerankMonolithicOmpDynamic(const G& x, const H& xt, const G& y, const H& yt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  using K = typename G::key_type;
  K     N = yt.order();                                        if (N==0) return PagerankResult<T>::initial(yt, q);
  auto [ks, n] = pagerankDynamicVertices(x, xt, y, yt, o, C);  if (n==0) return PagerankResult<T>::initial(yt, q);
  return pagerankOmp(yt, ks, K(0), n, pagerankMonolithicOmpLoopU<DEAD, SLEEP, HELP, K, T>, q, o);
}

template <bool DEAD=false, bool SLEEP=false, bool HELP=false, class G, class T=float>
PagerankResult<T> pagerankMonolithicOmpDynamic(const G& x, const G& y, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  auto xt = transposeWithDegree(x);
  auto yt = transposeWithDegree(y);
  return pagerankMonolithicOmpDynamic<DEAD, SLEEP, HELP>(x, xt, y, yt, q, o, C);
}
