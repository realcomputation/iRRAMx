#pragma once

#include <array>

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"
#include "homotopy.hpp"

using namespace iRRAM;

/** \addtogroup path
*  @{
*/

namespace iRRAM {

// [0,1] -> R^N
template <int N>
class Path : public Homotopy<1,N> {
public:
  /** Init with a homotopy
   * 
   * @param a homotopy \p f :[0,1] → ℝⁿ
   */
  Path(std::function<Point<N>(Point<1>)> f) : Homotopy<1,N>(f) {}
  Path(std::function<Point<N>(REAL)> f) : Homotopy<1,N>(wrapper(f)) {}

private:
  std::function<Point<N>(Point<1>)> wrapper(std::function<Point<N>(REAL)> &f) {
    return [=](Point<1> x) -> Point<N> { return f(x[0]); };
  }
};

}

/** @} */