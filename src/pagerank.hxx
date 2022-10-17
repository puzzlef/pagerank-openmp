#pragma once
#include <atomic>
#include <vector>
#include <utility>
#include <random>
#include "_main.hxx"
#include "components.hxx"

using std::atomic;
using std::vector;
using std::random_device;
using std::default_random_engine;
using std::move;




// PAGERANK-OPTIONS
// ----------------

template <class T>
struct PagerankOptions {
  int   repeat;
  bool  splitComponents;
  float sleepProbability;
  int   sleepDurationMs;
  T     damping;
  int   toleranceNorm;
  T     tolerance;
  int   maxIterations;

  PagerankOptions(int repeat=1, bool splitComponents=false, float sleepProbability=0.0f, int sleepDurationMs=0, T damping=0.85, int toleranceNorm=1, T tolerance=1e-6, int maxIterations=500) :
  repeat(repeat), splitComponents(splitComponents), sleepProbability(sleepProbability), sleepDurationMs(sleepDurationMs), damping(damping), toleranceNorm(toleranceNorm), tolerance(tolerance), maxIterations(maxIterations) {}
};




// PAGERANK-RESULT
// ---------------

template <class T>
struct PagerankResult {
  vector<T> ranks;
  int   iterations;
  float time;

  PagerankResult(vector<T>&& ranks, int iterations=0, float time=0) :
  ranks(ranks), iterations(iterations), time(time) {}

  PagerankResult(vector<T>& ranks, int iterations=0, float time=0) :
  ranks(move(ranks)), iterations(iterations), time(time) {}


  // Get initial ranks (when no vertices affected for dynamic pagerank).
  template <class G>
  static PagerankResult<T> initial(const G& x, const vector<T>* q=nullptr) {
    int  N = x.order();
    auto a = q? *q : createContainer(x, T());
    if (!q) fillValueAtU(a, x.vertexKeys(), T(1)/N);
    return {a, 0, 0};
  }
};




// PAGERANK-DATA
// -------------
// Using Pagerank Data for performance!

template <class G>
struct PagerankData {
  using K = typename G::key_type;
  vector2d<K> components;
  G blockgraph;
  G blockgraphTranspose;
};

template <class G, class K>
auto blockgraphD(const G& x, const vector2d<K>& cs, const PagerankData<G> *D) {
  return D? D->blockgraph : blockgraph(x, cs);
}

template <class G>
auto blockgraphTransposeD(const G& b, const PagerankData<G> *D) {
  return D? D->blockgraphTranspose : transpose(b);
}

template <class G, class H>
auto componentsD(const G& x, const H& xt, const PagerankData<G> *D) {
  return D? D->components : components(x, xt);
}




// PAGERANK-THREAD-WORK
// --------------------

struct PagerankThreadWork {
  random_device dev;          // used for random sleeps
  default_random_engine rnd;  // used for random sleeps
  volatile int    iteration;  // current iteration
  volatile double error;      // rank error wrt previous iteration for this thread
  volatile bool   stolen;     // indicates if a thread has stolen work
  atomic<size_t>  begin;      // vertex being processed
  atomic<size_t>  end;        // 1 + last vertex to be processed

  PagerankThreadWork(size_t begin=0, size_t end=0) :
  dev(), rnd(dev()), iteration(0), error(1), stolen(false), begin(begin), end(end) {}

  inline size_t size() const { return end -  begin; }
  inline bool empty()  const { return end <= begin; }
  inline void updateRange(size_t _begin, size_t _end) { begin = _begin; end = _end; }
  inline void clearRange() { stolen = false; begin = 0; end = 0; }
  inline void clear()      { iteration = 0;  error = 1; clearRange(); }
};


inline auto pagerankThreadWorks(int n) {
  vector<PagerankThreadWork*> a;
  for (int i=0; i<n; ++i)
    a.push_back(new PagerankThreadWork());
  return a;
}
