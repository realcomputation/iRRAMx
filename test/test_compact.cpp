#include <iRRAM.h>
#include "iRRAMx/compact.hpp"
#include "iRRAMx/plot.hpp"

using namespace iRRAM;

void test_compact() {
  int width = 600;        // image width
  int p = 100;            // precision for membership

  // definition: unit rectangle
  Compact<2> com([](Point<2> point, int p) -> bool {
    REAL &x=point[0], &y=point[1];
    DYADIC offset = scale(-p);

    // check if it's in the unitcell or close to it
    bool xBounded = choose(       // -2^-(p+1) < x < 1+2^-(p+1) ?
      -REAL(offset)*2 < x && x < 1+REAL(offset)*2,        // in
      x < -REAL(offset) || REAL(INTEGER(1))+offset < x             // out
      ) == 1;
    bool yBounded = choose(       // -2^-(p+1) < y < 1+2^-(p+1) ?
      -REAL(offset)*2 < y && y < 1+REAL(offset)*2,        // in
      y < -REAL(offset) || REAL(INTEGER(1))+offset < y             // out
      ) == 1;
    if(!xBounded || !yBounded) return false;

    // check if close to the border
    bool x0 = choose(       // is x close to 0?
      -REAL(offset)*2 < x && x < REAL(offset)*2,        // in
      x < -REAL(offset) || REAL(offset) < x             // out
      ) == 1;
    bool x1 = choose(       // is x close to 1?
      REAL(INTEGER(1))-offset*2 < x && x < REAL(INTEGER(1))+offset*2,        // in
      x < REAL(INTEGER(1))-offset || REAL(INTEGER(1))+offset < x             // out
      ) == 1;
    bool y0 = choose(       // is y close to 0?
      -REAL(offset)*2 < y && y < REAL(offset)*2,        // in
      y < -REAL(offset) || REAL(offset) < y             // out
      ) == 1;
    bool y1 = choose(       // is y close to 1?
      REAL(INTEGER(1))-offset*2 < y && y < REAL(INTEGER(1))+offset*2,        // in
      y < REAL(INTEGER(1))-offset || REAL(INTEGER(1))+offset < y             // out
      ) == 1;
    return (x0 | x1 | y0 | y1);
  });

  // membership test
  cout << "membership: " << com.member({1,1.1}, p) << "\n";

  // image
  Palette pal;
  com.plot2D(pal, width, -0.2, -0.2, 1.2, 1.2);
  writeImage("sampleCompact.png", pal);
  cout << "sampleCompact.png is generated.\n";
}

void test_path1() {
  int width = 600;        // image width
  int p = 100;            // precision for membership

  // path1: [0,1] -> R^2 (prolate cycloid)
  REAL r = 2, r0 = 3*r/2, v = 8*pi()*r;
  Path<2> path([=](REAL t) -> Point<2> {
    REAL theta = v*t/r;
    return Point<2>({v*t+r0*sin(theta), r0*cos(theta)});
  });

  // membership test
  cout << "membership: " << path.member({0,0}, p) << "\n";

  // image
  Palette pal;
  path.plot2D(pal, width, 0, -3.25, 8*pi(), 3.25);
  writeImage("samplePath1.png", pal);
  cout << "samplePath1.png is generated.\n";
}

void test_path2() {
  int width = 600;        // image width
  int p = 100;            // precision for membership

  // path2: [0,1] -> R^2 (whirl)
  Path<2> path([=](Point<1> x) -> Point<2> {
    REAL v = 4*pi()*x[0];
    return Point<2>({v*cos(v), v*sin(v)});
  });

  // membership test
  cout << "membership: " << path.member({0,0}, p) << "\n";

  // image
  Palette pal;
  path.plot2D(pal, width, -14, -14, 14, 14);
  writeImage("samplePath2.png", pal);
  cout << "samplePath2.png is generated.\n";
}

void test_surface1() {
  int width = 600;        // image width
  int p = 10;            // precision for membership

  // surface: [0,1]^2 -> R^2 (filled sin function)
  Surface<2> surface([=](REAL u, REAL v) -> Point<2> {
    REAL t = 6*pi()*u;
    return Point<2>({t, v*sin(t)});
  });

  // membership test
  // 이거 p=100일때 개오래걸림
  // disjunction test 하세
  cout << "membership: " << surface.member({0,0}, p) << "\n";

  // image
  Palette pal;
  surface.plot2D(pal, width, 0, -1, 6*pi(), 1);
  writeImage("sampleSurface1.png", pal);
  cout << "sampleSurface1.png is generated.\n";
}

void test_surface2() {
  int width = 600;        // image width
  int p = 100;            // precision for membership

  // surface: [0,1]^2 -> R^2 (whirl2)
  Surface<2> surface([=](Point<2> x) -> Point<2> {
    REAL v = 4*pi()*x[0];
    REAL th = (x[1]*0.4+0.6)*v;
    return Point<2>({v*cos(th), v*sin(th)});
  });

  // membership test
  cout << "membership: " << surface.member({1,1}, p) << "\n";

  // image
  Palette pal;
  surface.plot2D(pal, width, -14, -14, 14, 14);
  writeImage("sampleSurface2.png", pal);
  cout << "sampleSurface2.png is generated.\n";
}

void compute () {
  cout << "testing compact\n";
  test_compact();
  cout << "\n";
  return;

  cout << "testing path 1\n";
  test_path1();
  cout << "\n";

  cout << "testing path 2\n";
  test_path2();
  cout << "\n";

  cout << "testing surface 1\n";
  test_surface1();
  cout << "\n";

  cout << "testing surface 2\n";
  test_surface2();
  cout << "\n";

  // // // ttt
  // // Path<2> path1([=](REAL t) -> Point<2> {
  // //   // return Point<2>({t, log(1+t)});
  // //   return Point<2>({t, sin(t)});
  // //   // return Point<2>({t, t});
  // // });
  // // // cout << path1.member({0,0}, 100) << "\n";
  // // // cout << path1.member({RATIONAL(1,256),log(1+RATIONAL(1,256))}, 100) << "\n";
  // // path1.plot2D("samplePath1.png", width, 0, 0, 1, 1);
  // // cout << "samplePath1.png is generated\n";

  
  


  // // disjunction
  // Compact<2> disj = disjunction(com, path);
  // // disj.plot2D("disjunctionSample.png", 100, 0, -r0, v*1, r0);
  // disj.plot2D("sampleDisjunction.png", width, 0, -3, 4*pi(), 3);
  // cout << "disjunction drawn\n\n";
}

