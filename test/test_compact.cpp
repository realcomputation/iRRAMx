#include <iRRAM.h>
#include "iRRAM_extension/compact.hpp"

using namespace iRRAM;


void compute () {
  int width = 50;

  // compact: [0,2*pi]*[0,3]
  auto f = [](Point<2> point, int p) -> bool {
    REAL x = point[0], y = point[1];
    REAL x0=0, x1=2*pi(), y0=0, y1=3;   // boundary
    RATIONAL eps = Exp(-p);
    bool xTest = (choose(x > x0 && x < x1, x < x0+eps || x > x1-eps)==1);
    if(!xTest) return false;
    bool yTest = (choose(y > y0 && y < y1, y < y0+eps || y > y1-eps)==1);
    return yTest;
  };
  Compact<2> com(f);
  com.plot2D("sampleCompact.png", width, 0, -3, 4*pi(), 3);
  cout << "Compact drawn\n\n";

  // path: [0,1] -> R^2 (prolate cycloid)
  REAL r = 2, r0 = 3*r/2, v = 4*pi()*r;
  Path<2> path([=](REAL t) -> Point<2> {
    REAL theta = v*t/r;
    return Point<2>({v*t+r0*sin(theta), r0*cos(theta)});
  });
  path.plot2D("samplePath.png", width, 0, -3, 4*pi(), 3);
  cout << "path drawn\n\n";
  
  // surface: [0,1]^2 -> R^2 (a trapezoid, (u,v) -> (2u, 2v-2u))
  Surface<2> surface([=](REAL u, REAL v) -> Point<2> {
    return Point<2>({2*u, 2*v-2*u});
  });
  // surface.plot2D("sampleSurface.png", 20, 0, -2, 2, 2);
  surface.plot2D("sampleSurface.png", width, 0, -3, 4*pi(), 3);
  cout << "surface drawn\n\n";

  // disjunction
  Compact<2> disj = disjunction(com, path);
  // disj.plot2D("disjunctionSample.png", 100, 0, -r0, v*1, r0);
  disj.plot2D("sampleDisjunction.png", width, 0, -3, 4*pi(), 3);
  cout << "disjunction drawn\n\n";

  // conjucntion
  Compact<2> conj = conjunction(com, surface);
  // conj.plot2D("conjunctionSample.png", 20, 0, 0, v*1, r0);
  conj.plot2D("sampleConjunction.png", width, 0, -3, 4*pi(), 3);
  cout << "conjunction drawn\n\n";
}
