#pragma once
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "transpose.hxx"
#include "dynamic.hxx"
#include "pagerank.hxx"
#include "pagerankSeq.hxx"
#include "pagerankOmp.hxx"
#include "pagerankMonolithicSeq.hxx"

using std::vector;
using std::swap;
using std::min;




// PAGERANK-LOOP
// -------------

template <bool O, bool D, class K, class T, bool F=false>
int pagerankMonolithicBarrierfreeOmpLoopU(vector<T>& a, vector<T>& r, vector<T>& c, const vector<T>& f, const vector<K>& vfrom, const vector<K>& efrom, const vector<K>& vdata, vector<default_random_engine*>& rnds, K i, K n, K N, T p, T E, int L, int EF, float SP, int SD, K EI=K(), K EN=K()) {
  float l = 0;
  if (!O) return 0;
  if (F) EI = 0;
  if (F) EN = N;
  // Ordered approach
  int TS = int(min(K(omp_get_max_threads()), n));
  K   DN = ceilDiv(n, K(TS));
  #pragma omp parallel for schedule(static, 1) reduction(+:l)
  for (int t=0; t<TS; t++) {
    K    i1 = i+t*DN, I1 = min(i1+DN, i+n), n1 = max(I1-i1, K(0));
    if  (n1==0) continue;
    int  l1 = pagerankMonolithicSeqLoopU<O, D>(a, r, c, f, vfrom, efrom, vdata, rnds[t], i1, n1, N, p, E, L, EF, SP, SD, EI, EN);
    l += l1 * float(n1)/n;
  }
  return int(l + 0.5f);
}




// PAGERANK (STATIC / INCREMENTAL)
// -------------------------------

// Find pagerank using multiple threads (pull, CSR).
// @param x  original graph
// @param xt transpose graph (with vertex-data=out-degree)
// @param q  initial ranks (optional)
// @param o  options {damping=0.85, tolerance=1e-6, maxIterations=500}
// @returns {ranks, iterations, time}
template <bool O, bool D, bool F, class G, class H, class T=float>
PagerankResult<T> pagerankMonolithicBarrierfreeOmp(const G& x, const H& xt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  using K = typename G::key_type;
  K    N  = xt.order();  if (N==0) return PagerankResult<T>::initial(xt, q);
  auto ks = pagerankVertices(x, xt, o, C);
  return pagerankOmp(xt, ks, K(0), N, pagerankMonolithicBarrierfreeOmpLoopU<O, D, K, T, F>, q, o);
}

template <bool O, bool D, bool F, class G, class T=float>
PagerankResult<T> pagerankMonolithicBarrierfreeOmp(const G& x, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  auto xt = transposeWithDegree(x);
  return pagerankMonolithicBarrierfreeOmp<O, D, F>(x, xt, q, o, C);
}




// PAGERANK (DYNAMIC)
// ------------------

template <bool O, bool D, bool F, class G, class H, class T=float>
PagerankResult<T> pagerankMonolithicBarrierfreeOmpDynamic(const G& x, const H& xt, const G& y, const H& yt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  using K = typename G::key_type;
  K     N = yt.order();                                        if (N==0) return PagerankResult<T>::initial(yt, q);
  auto [ks, n] = pagerankDynamicVertices(x, xt, y, yt, o, C);  if (n==0) return PagerankResult<T>::initial(yt, q);
  return pagerankOmp(yt, ks, K(0), n, pagerankMonolithicBarrierfreeOmpLoopU<O, D, K, T, F>, q, o);
}

template <bool O, bool D, bool F, class G, class T=float>
PagerankResult<T> pagerankMonolithicBarrierfreeOmpDynamic(const G& x, const G& y, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  auto xt = transposeWithDegree(x);
  auto yt = transposeWithDegree(y);
  return pagerankMonolithicBarrierfreeOmpDynamic<O, D, F>(x, xt, y, yt, q, o, C);
}
