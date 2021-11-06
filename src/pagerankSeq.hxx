#pragma once
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "vertices.hxx"
#include "edges.hxx"
#include "csr.hxx"
#include "pagerank.hxx"

using std::vector;
using std::swap;




template <class T>
void pagerankFactor(vector<T>& a, const vector<int>& vdata, int i, int n, T p) {
  for (int u=i; u<i+n; u++) {
    int d = vdata[u];
    a[u] = d>0? p/d : 0;
  }
}


template <class T>
void pagerankCalculate(vector<T>& a, const vector<T>& c, const vector<int>& vfrom, const vector<int>& efrom, int i, int n, T c0) {
  for (int v=i; v<i+n; v++)
    a[v] = c0 + sumAt(c, sliceIter(efrom, vfrom[v], vfrom[v+1]));
}




// PAGERANK-ERROR
// --------------

template <class T>
T pagerankError(const vector<T>& x, const vector<T>& y, int i, int N, int EF) {
  switch (EF) {
    case 1:  return l1Norm(x, y, i, N);
    case 2:  return l2Norm(x, y, i, N);
    default: return liNorm(x, y, i, N);
  }
}




// PAGERANK-LOOP
// --------------

template <class T>
int pagerankMonolithicSeqLoop(vector<T>& a, vector<T>& r, vector<T>& c, const vector<T>& f, const vector<int>& vfrom, const vector<int>& efrom, int i, int n, int N, T p, T E, int L, int EF) {
  T  c0 = (1-p)/N;
  int l = 0;
  for (; l<L; ++l) {
    multiply(c, r, f, i, n);
    pagerankCalculate(a, c, vfrom, efrom, i, n, c0);
    T el = pagerankError(a, r, i, n, EF);
    if (el < E) { ++l; break; }
    swap(a, r);
  }
  return l;
}


// Find pagerank using a single thread (pull, CSR).
// @param xt transpose graph, with vertex-data=out-degree
// @param q initial ranks (optional)
// @param o options {damping=0.85, tolerance=1e-6, maxIterations=500}
// @returns {ranks, iterations, time}
template <class H, class T=float>
PagerankResult<T> pagerankMonolithicSeq(const H& xt, const vector<T> *q=nullptr, PagerankOptions<T> o={}) {
  T    p  = o.damping;
  T    E  = o.tolerance;
  int  L  = o.maxIterations, l = 0;
  int  EF = o.toleranceNorm;
  auto ks    = vertices(xt);
  auto vfrom = sourceOffsets(xt);
  auto efrom = destinationIndices(xt);
  auto vdata = vertexData(xt);
  int  N     = xt.order();
  vector<T> a(N), r(N), c(N), f(N), qc;
  if (q) qc = compressContainer(xt, *q);
  float t = measureDurationMarked([&](auto mark) {
    if (q) copy(r, qc);
    else fill(r, T(1)/N);
    copy(a, r);
    mark([&] { pagerankFactor(f, vdata, 0, N, p); });
    mark([&] { l = pagerankMonolithicSeqLoop(a, r, c, f, vfrom, efrom, 0, N, N, p, E, L, EF); });
  }, o.repeat);
  return {decompressContainer(xt, a), l, t};
}
