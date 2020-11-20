#pragma once

#include <iRRAM/lib.h>
#include "iRRAM_extension/random/random-core.hpp"
#include "iRRAM_extension/random/random-real.hpp"

using namespace iRRAM;
namespace iRRAM{


COMPLEX uniform_complex();
COMPLEX uniform_complex(COMPLEX, REAL);
COMPLEX gaussian_complex();
}
