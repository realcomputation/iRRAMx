#pragma once

#include <iRRAM/lib.h>
#include "iRRAMx/random/random-core.hpp"
#include "iRRAMx/random/random-real.hpp"

using namespace iRRAM;
namespace iRRAM{




class REALRAND
{
  private:
    int bitlength;
    std::string bits;
    REAL randreal;

  public:
    REALRAND();
    ~REALRAND();
    REAL asREAL();
};
}
