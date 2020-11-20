#pragma once

#include <iRRAM/lib.h>
#include "iRRAM_extension/random/random-core.hpp"

using namespace iRRAM;
namespace iRRAM{

REAL uniform_real();
REAL uniform_real(REAL, REAL);

REAL gaussian_real();
REAL gaussian_real(REAL, REAL);

REAL linear_real();
}
