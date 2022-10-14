#pragma once
#include <cmath>
#include <random>

using std::ceil;
using std::sqrt;
using std::uniform_int_distribution;




// COALESCE
// --------
// Similar to JavaScript coalescing || operator.

template <class T>
T coalesce(T x, T d=T()) {
  return x!=T()? x : d;
}




// CEIL-DIV
// --------
// For kernel launch calculation.

template <class T>
T ceilDiv(T x, T y) { return (x + y-1) / y; }
template <>
float ceilDiv<float>(float x, float y) { return ceil(x/y); }
template <>
double ceilDiv<double>(double x, double y) { return ceil(x/y); }




// SGN
// ---
// https://stackoverflow.com/a/4609795/1413259

template <typename T>
int sgn(T x) {
  return (T() < x) - (x < T());
}




// POW-2
// -----

template <class T>
constexpr bool isPow2(T x) noexcept {
  return !(x & (x-1));
}


template <class T>
constexpr T prevPow2(T x) noexcept {
  return 1 << T(log2(x));
}


template <class T>
constexpr T nextPow2(T x) noexcept {
  return 1 << T(ceil(log2(x)));
}




// PRIME
// -----

template <class T>
bool isPrime(T x) {
  // 1. 2, 3 are prime
  if (x<=3) return x>1;
  // 2. Multiples of 2, 3 not prime
  if (x % 2==0 || x % 3==0) return false;
  // 3. Factor of 6k-1 or 6k+1 => not prime
  for (T i=6, I=sqrt(x)+1; i<=I; i+=6)
    if (x % (i-1)==0 || x % (i+1)==0) return false;
  return true;
}


template <class T, class R>
T randomPrime(T begin, T end, R& rnd) {
  uniform_int_distribution<T> dis(begin, end);
  for (int i=128; i>0; --i) {
    T a = dis(rnd);
    if (isPrime(a)) return a;
  }
  return end-1;
}
