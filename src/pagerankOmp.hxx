#pragma once
#include <vector>
#include <algorithm>
#include <chrono>
#include <omp.h>
#include "_main.hxx"
#include "vertices.hxx"
#include "edges.hxx"
#include "csr.hxx"
#include "pagerank.hxx"
#include "pagerankSeq.hxx"

using std::vector;
using std::chrono::milliseconds;
using std::swap;




// PAGERANK-FACTOR
// ---------------
// For contribution factors of vertices (unchanging).

template <class K, class T>
void pagerankFactorOmpW(vector<T>& a, const vector<K>& vdata, K i, K n, T p) {
  if (n<SIZE_MIN_OMPM) { pagerankFactorW(a, vdata, i, n, p); return; }
  #pragma omp parallel for schedule(auto)
  for (K u=i; u<i+n; u++) {
    K  d = vdata[u];
    a[u] = d>0? p/d : 0;
  }
}




// PAGERANK-TELEPORT
// -----------------
// For teleport contribution from vertices (inc. dead ends).

template <class K, class T>
T pagerankTeleportOmp(const vector<T>& r, const vector<K>& vdata, K N, T p) {
  if (N<SIZE_MIN_OMPR) return pagerankTeleport(r, vdata, N, p);
  T a = (1-p)/N;
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (K u=0; u<N; u++)
    if (vdata[u] == 0) a += p*r[u]/N;
  return a;
}




// PAGERANK-CALCULATE
// ------------------
// For rank calculation from in-edges.

template <bool SLEEP=false, class K, class T>
void pagerankCalculateOmpW(vector<T>& a, const vector<T>& c, const vector<K>& vfrom, const vector<K>& efrom, K i, K n, T c0, vector<PagerankThreadWork*>& works) {
  double sp = double(SP)/n;
  milliseconds sd(SD);
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K v=i; v<i+n; v++) {
    int t = omp_get_thread_num();
    if (SLEEP) randomSleepFor(sd, sp, works[t]->rnd);
    a[v] = c0 + sumValuesAt(c, sliceIterable(efrom, vfrom[v], vfrom[v+1]));
  }
}




// PAGERANK-CALCULATE-HELPER
// -------------------------

template <class K>
inline bool pagerankIsChunkStart(K off, int chunkSize=2048) {
  return off && (chunkSize-1) == 0;
}

template <class R>
inline int pagerankBusyThread(vector<PagerankThreadWork*>& works, R& rnd) {
  int N = works.size();
  int p = randomPrime(N+1, 10*N, rnd);
  for (int i=0, q=0; i<N; ++i) {
    q = (q + p) % N;
    if (!works[q]->empty()) return q;
  }
  return -1;
}

template <bool ALL=false>
inline bool pagerankStealWorkRange(PagerankThreadWork& me, PagerankThreadWork& victim, size_t begin, size_t end) {
  if (begin>=end) return false;
  if (!victim.end.compare_exchange_strong(end, ALL? victim.begin : begin)) return false;
  victim.stolen = true;
  me.updateRange(begin, end);
  return true;
}

template <bool ALL=false>
inline bool pagerankStealWorkFrom(PagerankThreadWork& me, PagerankThreadWork& victim, size_t begin) {
  return pagerankStealWorkRange(me, victim, begin, victim.end);
}

template <bool ALL=false>
inline bool pagerankStealWorkSize(PagerankThreadWork& me, PagerankThreadWork& victim, size_t n) {
  size_t end = victim.end, begin = end - n;
  return pagerankStealWorkRange(me, victim, begin, end);
}

template <class F>
inline bool pagerankIsThreadStuck(PagerankThreadWork& victim, size_t oldBegin) {
  return !victim.empty() && victim.begin==oldBegin;
}


template <bool SLEEP=false, class K, class T>
void pagerankCalculateHelperOmpW(vector<T>& a, const vector<T>& c, const vector<K>& vfrom, const vector<K>& efrom, K i, K n, T c0, vector<PagerankThreadWork*>& works) {
  const int chunkSize = 2048;
  double sp = double(SP)/n;
  milliseconds sd(SD);
  #pragma omp parallel
  {
    int t = omp_get_thread_num();
    PagerankThreadWork& me = *(works[t]);
    #pragma omp for schedule(dynamic, 2048) nowait
    for (K v=i; v<i+n; v++) {
      if (pagerankIsChunkStart(v-i, chunkSize)) me.updateRange(v, v+chunkSize);
      if (SLEEP) randomSleepFor(sd, sp, me.rnd);
      a[v] = c0 + sumValuesAt(c, sliceIterable(efrom, vfrom[v], vfrom[v+1]));
      ++me.begin;  // work done on current vertex
    }
    for (;;) {
      int b = pagerankBusyThread(works, me.rnd);
      if (b<0) break;  // all threads free?
      PagerankThreadWork& victim = *(works[b]);
      if (victim.size() > chunkSize/2) {
        if (!pagerankStealWorkSize(me, victim, chunkSize/2)) continue;
        pagerankCalculateHelperW(a, c, vfrom, efrom, me.begin, me.size(), works[t]);
      }
      else {
        size_t begin = victim.begin;
        pagerankCalculateW(a, c, vfrom, efrom, begin, 1, c0, works[t]);
        if (!pagerankIsThreadStuck(victim, begin)) continue;
        if (!pagerankStealWorkFrom<true>(me, victim, begin+1)) continue;
        pagerankCalculateHelperW(a, c, vfrom, efrom, me.begin(), me.size(), works[t]);
      }
    }
  }
}




// PAGERANK-ERROR
// --------------
// For convergence check.

template <class K, class T>
T pagerankErrorOmp(const vector<T>& x, const vector<T>& y, K i, K N, int EF) {
  switch (EF) {
    case 1:  return l1NormOmp(x, y, i, N);
    case 2:  return l2NormOmp(x, y, i, N);
    default: return liNormOmp(x, y, i, N);
  }
}




// PAGERANK
// --------
// For Monolithic / Componentwise PageRank.

template <class H, class J, class M, class FL, class T=float>
PagerankResult<T> pagerankOmp(const H& xt, const J& ks, size_t i, const M& ns, FL fl, const vector<T> *q, const PagerankOptions<T>& o) {
  using K  = typename H::key_type;
  K     N  = xt.order();
  T     p  = o.damping;
  T     E  = o.tolerance;
  int   L  = o.maxIterations, l = 0;
  int   EF = o.toleranceNorm;
  float SP = o.sleepProbability;
  int   SD = o.sleepDurationMs;
  int   TH = omp_get_max_threads();
  auto vfrom = sourceOffsetsAs(xt, ks, K());
  auto efrom = destinationIndicesAs(xt, ks, K());
  auto vdata = vertexData(xt, ks);
  auto works = pagerankThreadWorks(TH);
  vector<T> a(N), r(N), c(N), f(N), qc;
  if (q) qc = compressContainer(xt, *q, ks);
  float t = measureDuration([&]() {
    if (q) copyValuesOmpW(r, qc);  // copy old ranks (qc), if given
    else fillValueOmpU(r, T(1)/N);
    pagerankFactorOmpW(f, vdata, K(0), N, p); multiplyValuesOmpW(c, r, f, 0, N);       // calculate factors (f) and contributions (c)
    l = fl(a, r, c, f, vfrom, efrom, vdata, works, K(i), ns, N, p, E, L, EF, SP, SD);  // calculate ranks of vertices
  }, o.repeat);
  return {decompressContainer(xt, r, ks), l, t};
}
