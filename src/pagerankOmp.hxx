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
using std::chrono::microseconds;
using std::this_thread::sleep_for;
using std::swap;
using std::min;
using std::max;




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


template <class K, class T>
T pagerankTeleportBarrierfreeOmp(const vector<T>& r, const vector<K>& vdata, K N, T p) {
  T a = (1-p)/N;
  #pragma omp for schedule(auto) reduction(+:a)
  for (K u=0; u<N; u++)
    if (vdata[u] == 0) a += p*r[u]/N;
  return a;
}




// PAGERANK-CALCULATE
// ------------------
// For rank calculation from in-edges.

template <bool SLEEP=false, class K, class T>
void pagerankCalculateOmpW(vector<T>& a, const vector<T>& c, const vector<K>& vfrom, const vector<K>& efrom, K i, K n, T c0, float SP, int SD, vector<PagerankThreadWork*>& works) {
  double sp = double(SP)/n;
  milliseconds sd(SD);
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K v=i; v<i+n; v++) {
    int t = omp_get_thread_num();
    pagerankCalculateRankW<SLEEP>(a, c, vfrom, efrom, v, c0, sp, sd, works[t]->rnd);
  }
}


template <bool SLEEP=false, class K, class T>
void pagerankCalculateBarrierfreeOmpW(vector<T>& a, const vector<T>& r, const vector<T>& f, const vector<K>& vfrom, const vector<K>& efrom, K i, K n, T c0, float SP, int SD, vector<PagerankThreadWork*>& works) {
  double sp = double(SP)/n;
  milliseconds sd(SD);
  double err = 0;
  int t = omp_get_thread_num();
  PagerankThreadWork& me = *works[t];
  #pragma omp for schedule(dynamic, 2048) nowait
  for (K v=i; v<i+n; v++) {
    T e = pagerankCalculateRankDeltaW<SLEEP>(a, r, f, vfrom, efrom, v, c0, sp, sd, me.rnd);
    err = max(err, e);  // li-norm
  }
  me.error = err;
}




// PAGERANK-CALCULATE-HELPER
// -------------------------

template <class K>
inline bool pagerankIsChunkStart(K off, int chunkSize=2048) {
  return off && (chunkSize-1) == 0;
}

template <bool ITER=false, class R>
inline int pagerankBusyThread(vector<PagerankThreadWork*>& works, R& rnd, int l=0) {
  int N = works.size();
  int p = randomPrime(N+1, 10*N, rnd);
  for (int i=0, q=0; i<N; ++i) {
    q = (q + p) % N;
    if (ITER && works[q]->iteration > l) continue;
    if (!works[q]->empty()) return q;
  }
  return -1;
}

template <bool ALL=false>
inline bool pagerankStealWork(PagerankThreadWork& me, PagerankThreadWork& victim, size_t begin, size_t end) {
  if (!victim.end.compare_exchange_strong(end, ALL? size_t(victim.begin) : begin)) return false;
  victim.stolen = true;
  me.updateRange(begin, end);
  return true;
}

template <bool ALL=false>
inline bool pagerankStealWorkSize(PagerankThreadWork& me, PagerankThreadWork& victim, size_t n) {
  size_t begin = victim.begin, end = victim.end;
  if (begin+n>=end) return false;
  pagerankStealWork<ALL>(me, victim, min(end-n, end), end);
  return true;
}

template <bool ALL=false, class F>
inline bool pagerankStealWorkIfSlow(PagerankThreadWork& me, PagerankThreadWork& victim, size_t n, F fn) {
  size_t begin = victim.begin, end = victim.end;
  if (begin==end) return false;
  fn (begin);
  if (victim.empty() || victim.begin!=begin) return false;
  pagerankStealWork<ALL>(me, victim, min(begin+n, end), end);
  return true;
}


template <bool SLEEP=false, class K, class T>
void pagerankCalculateHelperOmpW(vector<T>& a, const vector<T>& c, const vector<K>& vfrom, const vector<K>& efrom, K i, K n, T c0, float SP, int SD, vector<PagerankThreadWork*>& works) {
  const int chunkSize = 2048;
  const int stealSize = chunkSize/4;
  double sp = double(SP)/n;
  milliseconds sd(SD);
  // 0. Reset thread works.
  for (int i=0; i<works.size(); ++i)
    works[i]->clearRange();
  #pragma omp parallel
  {
    int t = omp_get_thread_num();
    PagerankThreadWork& me = *works[t];
    // 1. Perform work assigned to me.
    #pragma omp for schedule(dynamic, 2048) nowait
    for (K v=i; v<i+n; ++v) {
      if (me.stolen)  continue;
      if (me.empty()) me.updateRange(v, min(v+chunkSize, i+n));
      pagerankCalculateRankW<SLEEP>(a, c, vfrom, efrom, v, c0, sp, sd, me.rnd);
      ++me.begin;
    }
    while (true) {
      // 2. Perform remaining/stolen work.
      while (!me.empty())
        pagerankCalculateRankW<SLEEP>(a, c, vfrom, efrom, K(me.begin++), c0, sp, sd, me.rnd);
      // 3. Find a busy thread (victim), who has work.
      int b = pagerankBusyThread(works, me.rnd);
      if (b<0) break;
      // 4. Steal work from victim.
      PagerankThreadWork& victim = *works[b];
      if (pagerankStealWorkSize(me, victim, stealSize)) continue;
      pagerankStealWorkIfSlow<true>(me, victim, 1, [&](size_t begin) {  // or if stuck
        pagerankCalculateRankW(a, c, vfrom, efrom, K(begin), c0, sp, sd, me.rnd);
        sleep_for(microseconds(4));
      });
    }
  }
}


template <bool SLEEP=false, class K, class T>
void pagerankCalculateHelperBarrierfreeOmpW(vector<T>& a, const vector<T>& r, const vector<T>& f, const vector<K>& vfrom, const vector<K>& efrom, K i, K n, T c0, float SP, int SD, vector<PagerankThreadWork*>& works) {
  const int chunkSize = 2048;
  const int stealSize = chunkSize/4;
  double sp = double(SP)/n;
  milliseconds sd(SD);
  double err = 0;
  // 0. Reset thread works.
  int t = omp_get_thread_num();
  PagerankThreadWork& me = *works[t];
  int l = me.iteration;
  me.clearRange();
  // 1. Perform work assigned to me.
  #pragma omp for schedule(dynamic, 2048) nowait
  for (K v=i; v<i+n; ++v) {
    if (me.stolen)  continue;
    if (me.empty()) me.updateRange(v, min(v+chunkSize, i+n));
    T e = pagerankCalculateRankDeltaW<SLEEP>(a, r, f, vfrom, efrom, v, c0, sp, sd, me.rnd);
    err = max(err, e);  // li-norm
    ++me.begin;
  }
  while (true) {
    // 2. Perform remaining/stolen work.
    while (!me.empty()) {
      T e = pagerankCalculateRankDeltaW<SLEEP>(a, r, f, vfrom, efrom, K(me.begin++), c0, sp, sd, me.rnd);
      err = max(err, e);  // li-norm
    }
    // 3. Find a busy thread (victim), who has work.
    int b = pagerankBusyThread<true>(works, me.rnd, l);
    if (b<0) break;
    // 4. Steal work from victim.
    PagerankThreadWork& victim = *works[b];
    if (pagerankStealWorkSize(me, victim, stealSize)) continue;
    pagerankStealWorkIfSlow<true>(me, victim, 1, [&](size_t begin) {  // or if stuck
      T e = pagerankCalculateRankDeltaW(a, r, f, vfrom, efrom, K(begin), c0, sp, sd, me.rnd);
      err = max(err, e);  // li-norm
      sleep_for(microseconds(4));
    });
  }
  me.error = err;
  ++me.iteration;
}




// PAGERANK-ERROR
// --------------
// For convergence check.

template <class K, class T>
inline T pagerankErrorOmp(const vector<T>& x, const vector<T>& y, K i, K N, int EF) {
  switch (EF) {
    case 1:  return l1NormOmp(x, y, i, N);
    case 2:  return l2NormOmp(x, y, i, N);
    default: return liNormOmp(x, y, i, N);
  }
}

inline double pagerankErrorBarrierfreeOmp(vector<PagerankThreadWork*>& works) {
  double a = 0;
  for (int t=0; t<works.size(); ++t)
    a = max(a, double(works[t]->error));  // li-norm
  return a;
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
  forEach(works, [](PagerankThreadWork *w) { delete w; });
  return {decompressContainer(xt, r, ks), l, t};
}
