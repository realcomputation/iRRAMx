#pragma once
// iRRAM extension for compact subsets in euclidean
// author: Jiman Hwang(molehair a.t kaist.ac.kr)

#include <vector>

// accumulator
class Acc {
  int len;
  int base;
  std::vector<int> acc;

public:
  Acc(int len, int base=2) {
    this->len = len;
    this->base = base;
    acc.resize(len, 0);
  }

  // Increment the accumulator
  // ex base 2: 0000 -> 0001 -> 0010 -> ... -> 1111 -> 0000
  // return: true   if 1111 -> 0000
  //         false  otherwise
  bool inc() {
    int i=len-1;
    do {
      acc[i] = (acc[i]+1) % base;
      i--;
    } while(i>=0 && !acc[i+1]);
    return ((i<0) & !acc[0]);
  }

  int &operator[](int i) { return acc[i]; }
};