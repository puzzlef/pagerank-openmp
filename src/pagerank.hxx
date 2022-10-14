#pragma once
#include <vector>
#include <utility>
#include "_main.hxx"
#include "components.hxx"

using std::vector;
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
