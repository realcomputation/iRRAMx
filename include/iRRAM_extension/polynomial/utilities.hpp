/*
the header file define some basic utilities (funcitons) which is no included in a native iRRAM dist.
*/
#pragma once

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
// #include "intervals.h"
#include "iRRAM_extension/polynomial.hpp"

#include <utility>
using namespace iRRAM;
namespace iRRAM{

INTEGER factorial(int );
COMPLEX power(COMPLEX, int);
void print(COMPLEX );
DYADIC minimum(DYADIC , DYADIC );
DYADIC maximum(DYADIC , DYADIC );
RATIONAL minimum(RATIONAL , RATIONAL );
RATIONAL maximum(RATIONAL , RATIONAL );
// RATIONAL abs(RATIONAL);
// void print(OPENINTERVAL);
// void print(RATIONALINTERVAL);
// void print(INTERVALCOMPONENT);
// void print(RATIONALCOMPONENT);
void print(POLYNOMIAL );
}
