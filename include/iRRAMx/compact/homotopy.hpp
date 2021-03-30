#pragma once

#include <array>
#include <stack>

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"

using namespace iRRAM;

#include "compact.hpp"
#include "euclidean.hpp"
#include "plot.hpp"
#include "etc.hpp"

// 
#define MAX_ERROR_SIZE sizetype_power2(-5)
// #define MAX_ERROR_SIZE scale(-5)

namespace iRRAM {

/** \addtogroup homotopy
*  @{
*/

// f: [0,1]^M -> R^N
template <int M, int N>
class Homotopy : public Compact<N> {
  // original function
  std::function<Point<N>(Point<M>)> orig_f;

  // temporary storage for plot
  Palette *pal = NULL;
  DYADIC x1,y1,x2,y2;

public:
  /** Init with a homotopy
   * 
   * @param a homotopy \p f :[0,1]ᵐ → ℝⁿ
   */
  Homotopy(std::function<Point<N>(Point<M>)> f) {
    // set the homotopy
    orig_f = f;
  }

  std::function<Point<N>(Point<M>)> get_f() { return orig_f; }
  
  /** membership test for point with precision 2⁻ᴾ
   * 
   * If the current precision is not enough, it will increase, and then do membership test.
   * 
   * @param point a point in ℝⁿ
   * @param p target precision(higher \p p → higher precision)
   */
  bool member(Point<N> x, prec_t p) {
    // increase the computable precision
    increaseAvailablePrecision(p);

    // start with unit cell
    Cell<M> unitCell;
    return member_(x, p, unitCell);
  }
  
  /** Save the 2D graph to an .png file
   * 
  *  Area to draw: [x1, y1] X [x2, y2]
  * 
  *  Height will be determined automatically.
  *
  *  @param width the image width in pixel
  *  @param p precision(higher \p p → higher precision)
  * 
  *  @warning REQUIRE: \p x1 < \p x2, \p y1 < \p y2
  *  @warning Previously computed data will be reset for the exact drawing.
  */
  void plot2D(Palette &pal, int width, REAL x1, REAL y1, REAL x2, REAL y2) {
    // plane only
    assert(N==2);
    
    // need an uninitialized palette
    assert(pal.width==0);

    // image height
    int height = ((y2-y1)/(x2-x1)*width+RATIONAL(1,2)).as_INTEGER();

    // eval
    if(pal.width == 0) pal.init(width, height);
    this->pal = &pal;
    this->x1 = x1.as_DYADIC();
    this->x2 = x2.as_DYADIC();
    this->y1 = y1.as_DYADIC();
    this->y2 = y2.as_DYADIC();
    Cell<M> unitCell;
    _plot2D(unitCell);
  }

private:
  // check if x is in f(H) with respect to the precision p
  bool member_(Point<N> &x, prec_t p, Cell<M> &H) {
    if(sizetype_less(H[0].geterror(), MAX_ERROR_SIZE)) {
      // get the image of the cell
      HyperCube<N> fH = f(H);

      // get the distances
      // d: the distance between x and the center of H
      // r: the half length of diagonal of H
      REAL d=INTEGER(0), r=INTEGER(0);
      sizetype exactErr;
      sizetype_exact(exactErr);
      for(int i=0;i<N;i++) {
        // r
        DYADIC tmp = errSize(fH[i]);
        r += tmp*tmp;

        // d
        fH[i].seterror(exactErr);
        d += (x[i]-fH[i])*(x[i]-fH[i]);
      }
      r = root(r, N);
      d = root(d, N);
      
      // check the membership
      DYADIC b1 = scale(-p);
      DYADIC b2 = 2*b1;
      DYADIC offset = b1/4;
      bool isOut = (choose(d-r>b1, d-r<b1+offset) == 1);
      bool isIn = (choose(d+r<b2, d+r>b2-offset) == 1);

      // return!
      if(isOut) return false;
      if(isIn) return true;
    }

    // delve into sub-hypercubes of H
    Acc b(M, 2);
    Cell<M> HH;
    do {
      HH = H.halfCell(b);
      if(member_(x, p, HH)) return true;
    } while(!b.inc());
    return false;
  }

  // plot f(H)
  void _plot2D(Cell<M> &H) {
    if(sizetype_less(H[0].geterror(), MAX_ERROR_SIZE)) {
      HyperCube<N> fH = f(H);

      // compute the boundary of fH
      // fH is a subset of [u1,u2]*[v1,v2]
      DYADIC xHalfSize = errSize(fH[0]), yHalfSize = errSize(fH[1]);
      exactify<2>(fH);
      REAL u1 = fH[0]-xHalfSize, u2 = fH[0]+xHalfSize;
      REAL v1 = fH[1]-yHalfSize, v2 = fH[1]+yHalfSize;

      // need to draw?
      // Check if [x1,x2]*[y1,y2] overlaps [u1,u2]*[v1,v2].
      // If not, don't draw.
      bool isOut = (choose(x2 < u1, u1-xHalfSize/2 < x2) == 1);
      isOut |= (choose(u2 < x1, x1 < u2+xHalfSize/2) == 1);
      isOut |= (choose(y2 < v1, v1-yHalfSize/2 < y2) == 1);
      isOut |= (choose(v2 < y1, y1 < v2+yHalfSize/2) == 1);
      if(isOut) return;

      // fH is small enough?
      if(xHalfSize*2*pal->width <= x2-x1 && yHalfSize*2*pal->width <= x2-x1) {
        // reference point to color
        int j = ((fH[0]-x1)/(x2-x1)*pal->width).as_INTEGER();
        int i = ((y2-fH[1])/(x2-x1)*pal->width).as_INTEGER();

        // color (j,i), (j-1,i), (j,i-1), (j-1,i-1)
        if(j < pal->width && i < pal->height) pal->setColor(j, i, PLOT_COLOR_R, PLOT_COLOR_G, PLOT_COLOR_B);
        if(j && i < pal->height) pal->setColor(j-1, i, PLOT_COLOR_R, PLOT_COLOR_G, PLOT_COLOR_B);
        if(j < pal->width && i) pal->setColor(j, i-1, PLOT_COLOR_R, PLOT_COLOR_G, PLOT_COLOR_B);
        if(i && j) pal->setColor(j-1, i-1, PLOT_COLOR_R, PLOT_COLOR_G, PLOT_COLOR_B);
        return;
      }
    }

    // delve into sub-hypercubes of H
    Acc b(M, 2);
    Cell<M> HH;
    do {
      HH = H.halfCell(b);
      _plot2D(HH);
    } while(!b.inc());
  }

  Point<N> f(Point<M> &x) { return orig_f(x); }
  HyperCube<N> f(Cell<M> &H) { return orig_f(H.H); }
};




/** @} */
}   // namespace iRRAM