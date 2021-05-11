#pragma once
#include "_main.hxx"
#include "pagerank.hxx"




template <class C, class G, class T>
T pagerankTeleport(const C& r, const G& xt, T p) {
  int N = xt.order();
  T a = (1-p)/N;
  for (int u : xt.vertices())
    if (xt.vertexData(u) == 0) a += p*r[u]/N;
  return a;
}

template <class C, class G, class T>
void pagerankFactor(C& a, const G& xt, T p) {
  int N = xt.order();
  for (int u : xt.vertices()) {
    int d = xt.vertexData(u);
    a[u] = d>0? p/d : 0;
  }
}

template <class C, class G, class T>
void pagerankClassOnce(C& a, const C& c, const G& xt, T c0) {
  for (int v : xt.vertices())
    a[v] = c0 + sumAt(c, xt.edges(v));
}

template <class C, class G, class T>
int pagerankClassLoop(C& a, C& r, const C& f, C& c, const G& xt, T p, T E, int L) {
  int l = 0;
  T e0 = T();
  for (; l<L; l++) {
    T c0 = pagerankTeleport(r, xt, p);
    multiply(c, r, f);
    pagerankClassOnce(a, c, xt, c0);
    T e1 = absError(a, r);
    if (e1 < E || e1 == e0) break;
    swap(a, r);
    e0 = e1;
  }
  return l;
}

template <class C, class G, class T>
int pagerankClassCore(C& a, C& r, C& f, C& c, const G& xt, const C *q, T p, T E, int L) {
  int N = xt.order();
  if (q) copy(r, *q);
  else fillAt(r, T(1)/N, xt.vertices());
  pagerankFactor(f, xt, p);
  return pagerankClassLoop(a, r, f, c, xt, p, E, L);
}


// Find pagerank using C++ DiGraph class directly.
// @param xt transpose graph, with vertex-data=out-degree
// @param q initial ranks (optional)
// @param o options {damping=0.85, tolerance=1e-6, maxIterations=500}
// @returns {ranks, iterations, time}
template <class G, class T=float>
PagerankResult<T> pagerankClass(const G& xt, const vector<T> *q=nullptr, PagerankOptions<T> o={}) {
  T    p = o.damping;
  T    E = o.tolerance;
  int  L = o.maxIterations, l;
  auto a = xt.vertexContainer(T());
  auto r = xt.vertexContainer(T());
  auto f = xt.vertexContainer(T());
  auto c = xt.vertexContainer(T());
  float t = measureDuration([&]() { l = pagerankClassCore(a, r, f, c, xt, q, p, E, L); }, o.repeat);
  return {a, l, t};
}
