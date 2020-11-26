#pragma once
// iRRAM extension for compact subsets in euclidean
// author: Jiman Hwang(molehair a.t kaist.ac.kr)

#include <iRRAM/lib.h>
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
public:
  // characteristic function
  std::function< bool (Point<N> , int) > cfun;


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
  virtual bool member(Point<N> point, int p) { return this->cfun(point, p); }


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
  */
  void plot2D(const char *filename, int width, REAL x1, REAL y1, REAL x2, REAL y2) {
    // only plane
    if(N != 2) return;

    // some values..
    REAL pixelSize = (x2-x1)/REAL(width);      // single pixel size as a rect
    int height = ceil(((y2-y1)/pixelSize).as_double());               // image height

    // precision
    // Define p such that a ball centered at the center of pixel cover the pixel
    // That is, pixelSize/2*sqrt(2)  <  2^-p (radius of ball)
    // The drawn path will not be broken.
    int p = floor((REAL(0.5) - log(pixelSize)/ln2()).as_double());       // precision

    // plot each pixel to palette
    // start at the bottom row, from left to right.. row += 1 .. repeat
    // variable point stores the coordinate of the center of the current pixel
    Point<N> point = {REAL(0), y1 - pixelSize/2};      // init with coordinate outside image
    Palette pal(width, height);
    for(int i=height-1 ; i>=0 ; i--) {
      point[0] = x1 + pixelSize/2;    // x; the left most pixel
      point[1] += pixelSize;          // y; +1 row

      for(int j=0 ; j<width ; j++) {
        if(member(point, p)) {
          pal.setColor(j, i, PLOT_COLOR_R, PLOT_COLOR_G, PLOT_COLOR_B);
        }
        point[0] += pixelSize;
      }

      cout << (height-i) << " / " << height << " row done\n";
    }

    // to image
    writeImage(filename, pal);
  }

};
/** @} */



  
/** \addtogroup compact
*  @{
*/
/** pointwise conjunction for two compact subsets
 */
template <int N>
Compact<N> conjunction(Compact<N> &com1, Compact<N> &com2) {
  auto cfun = [&] (Point<N> pt, int p) -> bool {
    return com1.member(pt, p) && com2.member(pt, p);
  };
  return Compact<N>(cfun);
}

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
