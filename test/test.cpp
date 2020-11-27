#include <iRRAM.h>
#include "iRRAMx/utility.hpp"
#include "iRRAMx/wiener.hpp"
// #include<pair.h>
using namespace iRRAM;


std::function <REAL(REAL)> join (REAL x, std::function<REAL(REAL)> f, std::function<REAL(REAL)> g)
{
    return ([=](REAL m)->REAL { return iRRAM::glue (m<x, f(m), g(m)); });
}


void compute () {
  WIENER w1;
  cout<<w1.compute(0.1)<<"\n";
  cout<<w1.compute(0.2)<<"\n";
  cout<<w1.compute(0.3)<<"\n";
}
