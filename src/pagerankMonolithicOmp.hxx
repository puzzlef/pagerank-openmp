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

template <bool O, class T>
int pagerankMonolithicOmpLoopU(vector<T>& a, vector<T>& r, vector<T>& c, const vector<T>& f, const vector<int>& vfrom, const vector<int>& efrom, const vector<int>& vdata, int i, int n, int N, T p, T E, int L, int EF) {
  int l = 0;
  // Unordered approach
  while (!O && l<L) {
    T c0 = pagerankTeleportOmp(r, vdata, N, p);
    pagerankCalculateOmpW(a, c, vfrom, efrom, i, n, c0);  // update ranks of vertices
    multiplyValuesOmpW(c, a, f, i, n);                    // update partial contributions (c)
    T el = pagerankErrorOmp(a, r, i, n, EF);              // compare previous and current ranks
    swap(a, r); ++l;                                      // final ranks in (r)
    if (el<E) break;                                      // check tolerance
  }
  // Ordered approach
  while (O && l<L) {
    T c0 = pagerankTeleportOmp(r, vdata, N, p);
    pagerankCalculateOrderedOmpU(a, r, f, vfrom, efrom, i, n, c0);  // update ranks of vertices
    T el = pagerankErrorOmp(a, i, n, EF); ++l;            // compare previous and current ranks
    if (el<E) break;                                      // check tolerance
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
template <bool O, class G, class H, class T=float>
PagerankResult<T> pagerankMonolithicOmp(const G& x, const H& xt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  int  N  = xt.order();  if (N==0) return PagerankResult<T>::initial(xt, q);
  auto ks = vertexKeys(xt);
  return pagerankOmp<O>(xt, ks, 0, N, pagerankMonolithicOmpLoopU<O, T>, q, o);
}

template <bool O, class G, class T=float>
PagerankResult<T> pagerankMonolithicOmp(const G& x, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  auto xt = transposeWithDegree(x);
  return pagerankMonolithicOmp<O>(x, xt, q, o);
}




// PAGERANK (DYNAMIC)
// ------------------

template <bool O, class G, class H, class T=float>
PagerankResult<T> pagerankMonolithicOmpDynamic(const G& x, const H& xt, const G& y, const H& yt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  int  N = yt.order();                             if (N==0) return PagerankResult<T>::initial(yt, q);
  auto [ks, n] = dynamicInVertices(x, xt, y, yt);  if (n==0) return PagerankResult<T>::initial(yt, q);
  return pagerankOmp<O>(yt, ks, 0, n, pagerankMonolithicOmpLoopU<O, T>, q, o);
}

template <bool O, class G, class T=float>
PagerankResult<T> pagerankMonolithicOmpDynamic(const G& x, const G& y, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  auto xt = transposeWithDegree(x);
  auto yt = transposeWithDegree(y);
  return pagerankMonolithicOmpDynamic<O>(x, xt, y, yt, q, o);
}
