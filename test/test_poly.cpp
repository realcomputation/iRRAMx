#include <iRRAM.h>
#include "iRRAMx/poly.hpp"

using namespace iRRAM;
void compute () {
    // 1D real polynomial
    POLYNOMIAL<REAL, 1> p1;
    p1.set_coefficient({2}, 2.5);
    p1.set_coefficient({1}, 1);
    cout << to_string(p1, {'x'}, 10) << std::endl;
    cout << p1({-2}) << std::endl;
    // 2D real polynomial
    POLYNOMIAL<REAL, 2> p2;
    p2.set_coefficient({1,1}, 2);
    p2.set_coefficient({3,0}, 1);
    p2.set_coefficient({0,0}, -1);
    cout << to_string(p2, {'x', 'y'}, 10) << std::endl;
    cout << p2({1,2}) << std::endl;
    // another 2D real polynomial
    POLYNOMIAL<REAL, 2> p3;
    p3.set_coefficient({5,2}, 1);
    cout << to_string(p3, {'x', 'y'}, 10) << std::endl;
    cout << p3({1,1}) << std::endl;
    // addition
    auto p4 = p2+p3;
    cout << to_string(p4, {'x', 'y'}, 10) << std::endl;
    cout << p4({2,1}) << std::endl;
    // multiplication
    auto p5 = p2*p3;
     cout << to_string(p5, {'x', 'y'}, 10) << std::endl;
    cout << p5({2,4}) << std::endl;
    cout << p2({2,4})*p3({2,4}) << std::endl;
    // derivative
    auto p6 = derive(p5,{1,3});
    cout << to_string(p6, {'x', 'y'}, 10) << std::endl;
    cout << p6({2,55}) << std::endl;
}