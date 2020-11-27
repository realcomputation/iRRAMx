#pragma once
#include "iRRAMx/fix.hpp"

#include <iRRAM/lib.h>

using namespace iRRAM;
namespace iRRAM{

REAL prec(int p); // p => 2^p
std::string to_string(REAL, unsigned int);
std::string to_string(COMPLEX, unsigned int);



}
