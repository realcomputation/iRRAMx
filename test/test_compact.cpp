#include <iRRAM.h>
#include "iRRAMx/compact.hpp"

using namespace iRRAM;



void compute () {
  int width = 600;

  // TODO: piano curve
  // // compact: [0,2*pi]*[0,3]
  // auto f = [](Point<2> point, int p) -> bool {
  //   REAL x = point[0], y = point[1];
  //   REAL x0=0, x1=2*pi(), y0=0, y1=3;   // boundary
  //   RATIONAL eps = Exp(-p);
  //   bool xTest = (choose(x > x0 && x < x1, x < x0+eps || x > x1-eps)==1);
  //   if(!xTest) return false;
  //   bool yTest = (choose(y > y0 && y < y1, y < y0+eps || y > y1-eps)==1);
  //   return yTest;
  // };
  // Compact<2> com(f);
  // com.plot2D("sampleCompact.png", width, 0, -3, 4*pi(), 3);
  // cout << "Compact drawn\n";



  // // ttt
  // Path<2> path1([=](REAL t) -> Point<2> {
  //   // return Point<2>({t, log(1+t)});
  //   return Point<2>({t, sin(t)});
  //   // return Point<2>({t, t});
  // });
  // // cout << path1.member({0,0}, 100) << "\n";
  // // cout << path1.member({RATIONAL(1,256),log(1+RATIONAL(1,256))}, 100) << "\n";
  // path1.plot2D("samplePath1.png", width, 0, 0, 1, 1);
  // cout << "samplePath1.png is generated\n";

  // path1: [0,1] -> R^2 (prolate cycloid)
  REAL r = 2, r0 = 3*r/2, v = 8*pi()*r;
  Path<2> path1([=](REAL t) -> Point<2> {
    REAL theta = v*t/r;
    return Point<2>({v*t+r0*sin(theta), r0*cos(theta)});
  });
  path1.plot2D("samplePath1.png", width, 0, -3.25, 8*pi(), 3.25);
  cout << "samplePath1.png is generated\n";

  // path2: [0,1] -> R^2 (whirl)
  Path<2> path2([=](Point<1> x) -> Point<2> {
    REAL v = 4*pi()*x[0];
    return Point<2>({v*cos(v), v*sin(v)});
  });
  path2.plot2D("samplePath2.png", width, -14, -14, 14, 14);
  cout << "samplePath2.png is generated\n";
  
  // surface: [0,1]^2 -> R^2 (filled sin function)
  Surface<2> surface1([=](REAL u, REAL v) -> Point<2> {
    REAL t = 6*pi()*u;
    return Point<2>({t, v*sin(t)});
  });
  surface1.plot2D("sampleSurface1.png", width, 0, -1, 6*pi(), 1);
  cout << "sampleSurface1.png is generated\n";
  
  // surface: [0,1]^2 -> R^2
  Surface<2> surface2([=](Point<2> x) -> Point<2> {
    REAL v = 4*pi()*x[0];
    REAL th = (x[1]*0.4+0.6)*v;
    return Point<2>({v*cos(th), v*sin(th)});
  });
  // cout << surface2.member({1,1}, 200) << "\n";
  surface2.plot2D("sampleSurface2.png", width, -14, -14, 14, 14);
  cout << "sampleSurface2.png is generated\n";



//   // disjunction
//   Compact<2> disj = disjunction(com, path);
//   // disj.plot2D("disjunctionSample.png", 100, 0, -r0, v*1, r0);
//   disj.plot2D("sampleDisjunction.png", width, 0, -3, 4*pi(), 3);
//   cout << "disjunction drawn\n\n";
// 이거 내부에서 member안 쓰고 하도록 해야함

//   // conjucntion
//   Compact<2> conj = conjunction(com, surface);
//   // conj.plot2D("conjunctionSample.png", 20, 0, 0, v*1, r0);
//   conj.plot2D("sampleConjunction.png", width, 0, -3, 4*pi(), 3);
//   cout << "conjunction drawn\n\n";
}

