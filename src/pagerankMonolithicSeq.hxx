#pragma once
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "transpose.hxx"
#include "dynamic.hxx"
#include "pagerank.hxx"
#include "pagerankSeq.hxx"

using std::vector;
using std::swap;




// PAGERANK-LOOP
// -------------

template <bool DEAD=false, bool SLEEP=false, class K, class T>
int pagerankMonolithicSeqLoopU(vector<T>& a, vector<T>& r, vector<T>& c, const vector<T>& f, const vector<K>& vfrom, const vector<K>& efrom, const vector<K>& vdata, PagerankThreadWork *work, K i, K n, K N, T p, T E, int L, int EF, float SP, int SD) {
  int l = 0;
  while (l<L) {
    T c0 = DEAD? pagerankTeleport(r, vdata, N, p) : (1-p)/N;
    pagerankCalculateW<SLEEP>(a, c, vfrom, efrom, i, n, c0, SP, SD, work);  // update ranks of vertices
    multiplyValuesW(c, a, f, i, n);        // update partial contributions (c)
    T el = pagerankError(a, r, i, n, EF);  // compare previous and current ranks
    swap(a, r); ++l;                       // final ranks in (r)
    if (el<E) break;                       // check tolerance
  }
  return l;
}




// PAGERANK (STATIC / INCREMENTAL)
// -------------------------------

// Find pagerank using a single thread (pull, CSR).
// @param x  original graph
// @param xt transpose graph (with vertex-data=out-degree)
// @param q  initial ranks (optional)
// @param o  options {damping=0.85, tolerance=1e-6, maxIterations=500}
// @returns {ranks, iterations, time}
template <bool DEAD=false, bool SLEEP=false, class G, class H, class T=float>
PagerankResult<T> pagerankMonolithicSeq(const G& x, const H& xt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  using K = typename G::key_type;
  K    N  = xt.order();  if (N==0) return PagerankResult<T>::initial(xt, q);
  auto ks = pagerankVertices(x, xt, o, C);
  return pagerankSeq(xt, ks, K(0), N, pagerankMonolithicSeqLoopU<DEAD, SLEEP, K, T>, q, o);
}

template <bool DEAD=false, bool SLEEP=false, class G, class T=float>
PagerankResult<T> pagerankMonolithicSeq(const G& x, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  auto xt = transposeWithDegree(x);
  return pagerankMonolithicSeq<DEAD, SLEEP>(x, xt, q, o, C);
}




// PAGERANK (DYNAMIC)
// ------------------

template <bool DEAD=false, bool SLEEP=false, class G, class H, class T=float>
PagerankResult<T> pagerankMonolithicSeqDynamic(const G& x, const H& xt, const G& y, const H& yt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  using K = typename G::key_type;
  K    N  = yt.order();                                        if (N==0) return PagerankResult<T>::initial(yt, q);
  auto [ks, n] = pagerankDynamicVertices(x, xt, y, yt, o, C);  if (n==0) return PagerankResult<T>::initial(yt, q);
  return pagerankSeq(yt, ks, K(0), n, pagerankMonolithicSeqLoopU<DEAD, SLEEP, K, T>, q, o);
}

template <bool DEAD=false, bool SLEEP=false, class G, class T=float>
PagerankResult<T> pagerankMonolithicSeqDynamic(const G& x, const G& y, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  auto xt = transposeWithDegree(x);
  auto yt = transposeWithDegree(y);
  return pagerankMonolithicSeqDynamic<DEAD, SLEEP>(x, xt, y, yt, q, o, C);
}
