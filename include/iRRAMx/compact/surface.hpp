#pragma once

#include <array>

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"
#include "homotopy.hpp"

using namespace iRRAM;

/** \addtogroup surface
*  @{
*/

namespace iRRAM {

// [0,1]^2 -> R^N
template <int N>
class Surface : public Homotopy<2,N> {
public:
  /** Init with a homotopy
   * 
   * @param a homotopy \p f :[0,2] → ℝⁿ
   */
  Surface(std::function<Point<N>(Point<2>)> f) : Homotopy<2,N>(f) {}
  Surface(std::function<Point<N>(REAL,REAL)> f) : Homotopy<2,N>(wrapper(f)) {}

private:
  std::function<Point<N>(Point<2>)> wrapper(std::function<Point<N>(REAL,REAL)> &f) {
    return [=](Point<2> x) -> Point<N> { return f(x[0],x[1]); };
  }
};

}

/** @} */