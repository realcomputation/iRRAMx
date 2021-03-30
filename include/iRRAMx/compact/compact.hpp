#pragma once
// iRRAM extension for compact subsets in euclidean
// author: Jiman Hwang(molehair a.t kaist.ac.kr)

#include <iRRAM/lib.h>
#include <vector>
#include "euclidean.hpp"
#include "plot.hpp"


using namespace iRRAM;


#define PLOT_COLOR_R    0x00
#define PLOT_COLOR_G    0x00
#define PLOT_COLOR_B    0x00

namespace iRRAM {

/** \addtogroup compact
*  @{
*/
// R^N
template <int N>
class Compact {
private:
  // disjunction list
  std::vector<Compact<N>> disj;
  
  // characteristic function
  std::function< bool (Point<N> , int) > cfun;

public:
  /*! @brief init with the emptyset
  */
  Compact() { this->cfun = [=](Point<N>, int) -> bool { return false; }; }


  /*! @brief init with characteristic func
  *
  *  @param cfun the characteristic function that defines a compact set
  */
  Compact(std::function<bool(Point<N>, int)> cfun) { this->cfun = cfun; }


  /*! @brief membership test for point with precision 2⁻ᴾ
  *
  *  @param point a point in ℝⁿ
  *  @param p precision(higher p → higher precision)
  *
  *  @return the result of calling the characteristic function with \p point and \p p
  */
  virtual bool member(Point<N> point, int p) {
    // belong to the current characteristic function?
    if(this->cfun(point, p)) return true;

    // and unions?
    for(Compact<N> &com : disj) {
      if(com.member(point, p)) return true;
    }

    // no membership
    return false;
  }


  /** Save the 2D graph to an .png file
   * 
  *  Area to draw: [x1, y1] X [x2, y2]
  * 
  *  Height will be determined automatically.
  *
  *  @param pal image container
  *  @param width the image width in pixel
  *  @param p precision(higher \p p → higher precision)
  * 
  *  @warning REQUIRE: \p pal must be uninitialized.
  *  @warning REQUIRE: \p x1 < \p x2, \p y1 < \p y2
  */
  virtual void plot2D(Palette &pal, int width, DYADIC x1, DYADIC y1, DYADIC x2, DYADIC y2) {
    // only plane
    assert(N==2);
    
    // need an uninitialized palette
    assert(pal.width==0);

    // some values
    REAL pixelSize = (x2-x1)/REAL(width);      // single pixel size as a rect
    int height = ceil(((y2-y1)/pixelSize).as_double());               // image height

    // plot the current chracteristic function
    pal.init(width, height);
    _plot2D_this(pal, x1, y1, x2, y2);

    // plot unions for each 'disj' and merge(?)
    for(Compact<N> &com : disj) {
      Palette palDisj;
      com.plot2D(palDisj, width, x1, y1, x2, y2);
      pal.union_with(palDisj);
    }
  }

  /** Update this compact set in disjunction with another one
   * 
  *  @param com the compact set with which union
  */
  void union_with(Compact<N> &com) {
    disj.push_back(com);
  }

private:
  // plot this->cfunc
  // manually check membership for each pixel
  virtual void _plot2D_this(Palette &pal, DYADIC &x1, DYADIC &y1, DYADIC &x2, DYADIC &y2) {
    // some values..
    REAL pixelSize = (x2-x1)/REAL(pal.width);      // single pixel size as a rect

    // precision
    int p = floor((RATIONAL(1,2) - log(pixelSize)/ln2()).as_double());

    // plot each pixel to palette
    // start at the bottom row, from left to right.. row += 1 .. repeat
    // variable point stores the coordinate of the center of the current pixel
    Point<N> point = {REAL(0), y1 - pixelSize/2};      // init with coordinate outside image
    for(int i=pal.height-1 ; i>=0 ; i--) {
      point[0] = x1 + pixelSize/2;    // x; the left most pixel
      point[1] += pixelSize;          // y; +1 row

      for(int j=0 ; j<pal.width ; j++) {
        if(member(point, p)) {
          pal.setColor(j, i, PLOT_COLOR_R, PLOT_COLOR_G, PLOT_COLOR_B);
        }
        point[0] += pixelSize;
      }
    }
  }
};
/** @} */



  
/** \addtogroup compact
*  @{
*/


/** pointwise disjunction for two compact subsets
 */
template <int N>
Compact<N> disjunction(Compact<N> &com1, Compact<N> &com2) {
  auto cfun = [&] (Point<N> pt, int p) -> bool {
    return com1.member(pt, p) || com2.member(pt, p);
  };
  return Compact<N>(cfun);
}

/** @} */
}
