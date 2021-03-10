#pragma once
// iRRAM extension for compact subsets in euclidean
// author: Jiman Hwang(molehair a.t kaist.ac.kr)

#include <iRRAM/lib.h>
#include "plot.hpp"


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

template <int N>
class Cell {
public:
  prec_t prec;
  IntegerPoint<N> coord;

  // // the squared distance between two cells
  // static DYADIC distance2(Cell<N> &x, Cell<N> &y) {
  //   DYADIC dist2 = 0;
  //   DYADIC t1, t2;
    
  //   for(int j=0;j<N;j++) {
  //     t1 = scale(-x.prec-1) * ((x.coord[j]<<1)+1);
  //     t2 = scale(-y.prec-1) * ((y.coord[j]<<1)+1);
  //     dist2 = dist2 + (t1-t2)*(t1-t2);
  //   }
  //   return dist2;
  // }

  // set this instance to [0,1]^N
  void makeUnitCell() {
    prec = 0;
    for(int j=0;j<N;j++) coord[j] = 0;
  }

  INTEGER& operator[](int i) { return coord[i]; }
};


//
DYADIC errSize(REAL &r) {
  sizetype err;
  r.geterror(err);
  sizetype_normalize(err);
  return scale(err.exponent) * int(err.mantissa);
}

//
template <int N>
void exactify(HyperCube<N> &hc) {
  sizetype err;
  sizetype_exact(err);
  for(int i=0;i<N;i++) hc[i].seterror(err);
}

}