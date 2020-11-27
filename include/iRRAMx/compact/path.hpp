#pragma once

#include <array>

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"
#include "compact.hpp"
#include "euclidean.hpp"
#include "plot.hpp"

using namespace iRRAM;

/** \addtogroup path
*  @{
*/

// [0,1] -> R^N
template <int N>
class Path : public Compact<N> {
public:
  // original function
  std::function<Point<N>(REAL)> f;

  /** Init with a homotopy
   * 
   * @param a homotopy \p f :[0,1] → ℝⁿ
   */
  Path(std::function<Point<N>(REAL)> f) { this->f = f; }

  // whenever |x-z| < 2^-pArg, |f(x)-f(z)| < 2^-p for all x,z
  // will be increased when higher precision is requested
  int p=INT_MIN, pArg=INT_MIN;

  /** Increase the current precision and find the corresponding \p pArg
   * 
   * @param p target precision(higher \p p → higher precision)
   */
  void increasePrecision(int p) {
    // ignore lower or equal precision
    if(this->p >= p) return;

    // find the pArg
    // must |f(u)-f(z)| < (2^-p)/sqrt(2) to include the box with a ball
    // hence find find pArg such that whenever |x-z| < 2^-pArg, |f(x)-f(z)| < 2^(-p-1) for all x,z
    // This makes the distance between two consecutive centers of balls 2^(-p-1)*sqrt(2) at most.
    std::function<Point<N>(Point<1>)> ff = [&](Point<1> x) -> Point<N> {
      return this->f(x[0]);
    };
    this->pArg = module2<1,N>(ff, p+1);

    // update the current precision
    this->p = p;

    // update the current characteristic function
    this->cfun = [&](Point<N> pt, int p) -> bool{
      single_valued code;

      // check the membership with previously found p and pArg
      // centers of balls: f(step/2 + i*step, step/2 + j*step)
      // radius of a ball: 2^-p
      // maximum distance between two consecutive balls(centers): 2^(-p-1)*sqrt(2)   (check increasePrecision())
      // When checking inclusiveness, we run (d < radius) and (d > radius/2) in parallel for decidability.
      // For any point on the path, there exists a ball that contains the point.
      RATIONAL step = Exp(-pArg);   // 2^-pArg
      RATIONAL radius = Exp(-p);     // 2^-p, the radius of a ball
      RATIONAL radiusHalf = radius/INTEGER(2);     // 2^(-p-1)
      REAL d;
      RATIONAL u,v;
      for(u=step/2 ; u<=ONE ; u+=step) {
        d = IR_d<N>(pt, this->f(u));
        if(choose(d < radius, d > radiusHalf) == 1) return true;
      }
      return false;
    };
  }
  
  /** membership test for point with precision 2⁻ᴾ
   * 
   * If the current precision is not enough, it will increase, and then do membership test.
   * 
   * @param point a point in ℝⁿ
   * @param p target precision(higher \p p → higher precision)
   */
  bool member (Point<N> point, int p) {
    // cout << "member in path p = " << p << "\n";

    // check if previously found pArg is viable
    // if not, increase the precision
    if(this->p < p) this->increasePrecision(p);
    
    return this->cfun(point, p);
  }
};

/** @} */