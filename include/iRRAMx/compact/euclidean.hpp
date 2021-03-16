#pragma once
// iRRAM extension for compact subsets in euclidean
// author: Jiman Hwang(molehair a.t kaist.ac.kr)

#include <iRRAM/lib.h>
#include "plot.hpp"
#include "etc.hpp"


using namespace iRRAM;

namespace iRRAM {
typedef int prec_t;
typedef size_t hash_t;

// Euclidean space using alias declaration
template <int N>
using IR = std::array<REAL, N>;

// Hypercube as a (mental) subset
template <int N>
using HyperCube = IR<N>;

// just alias
template <int N>
using Point = IR<N>;

// Dyadic point
template <int N>
using DyadicPoint = std::array<DYADIC, N>;

// integer coordinate
template <int N>
class IntegerPoint {
  std::array<INTEGER, N> coord;
public:
  INTEGER& operator[](int i) { return coord[i]; }
  hash_t hash() {
    std::string st = "";
    for(INTEGER &x : coord) st += swrite(x) + "|";
    return std::hash<std::string>{}(st);
  }
};

// return 2^n
DYADIC scale(int n) {
  if(n > 0) return (INTEGER(1) << n);
  else if(n == 0) return 1;
  else return DYADIC(1) / DYADIC(INTEGER(1) << -n);
}


//
inline DYADIC toDYADIC(sizetype &err) {
  sizetype_normalize(err);
  return scale(err.exponent) * int(err.mantissa);
}

//
inline DYADIC errSize(REAL &r) {
  sizetype err;
  r.geterror(err);
  return toDYADIC(err);
}

//
void exactify(REAL &r) {
  sizetype err;
  sizetype_exact(err);
  r.seterror(err);
}

//
template <int N>
void exactify(HyperCube<N> &hc) {
  for(int i=0;i<N;i++) exactify(hc[i]);
}


template <int N>
class Cell {
public:
  HyperCube<N> H;

  Cell() {
    // set this instance to [0,1]^N
    sizetype err;
    err = sizetype_power2(-1);
    for(int j=0;j<N;j++) {
      H[j] = RATIONAL(1,2);
      H[j].seterror(err);
    }
  }

  // 
  inline Cell halfCell(Acc &acc) {
    sizetype err;
    Cell<N> HH;
    for(int i=0;i<N;i++) {
      H[i].geterror(err);
      sizetype_half(err,err);
      HH[i] = H[i];
      exactify(HH[i]);
      HH[i] = HH[i] + toDYADIC(err) * ((acc[i]<<1)-1);
      HH[i].seterror(err);
    }
    return HH;
  }
  
  REAL& operator[](int i) { return H[i]; }
};
}