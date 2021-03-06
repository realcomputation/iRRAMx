#include "iRRAMx/random.hpp"

#include <utility>
#include <stdlib.h>
#include <random>
#include <vector>
#include <string>

using namespace iRRAM;

namespace iRRAM{
INTEGER pow(int a, int b)
{
  INTEGER result = 1;
  for(int i=0; i< b; i++)
    result = result * a;
  return result;
}


REALRAND::REALRAND()
  {
    int prec = actual_stack().actual_prec;

    sizetype err;
    std::string r;
    int s;
    INTEGER result_integer;
    int needblock, neededblock;
    std::string result_bits;
    int result_size ;



    // random int generator (uniformly distributed [0,2^16-1])
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, 65535);


    if (actual_stack().inlimit == 0)
    {
      if(get_cached(result_bits)
        && get_cached(result_size))
      {
        r = result_bits;
        s = result_size;
      }
      else
      {
        result_bits = std::to_string(dis(gen));
        result_size = 16;
        put_cached(result_bits);
        put_cached(result_size);
        r = result_bits;
        s = result_size;
      }

    INTEGER result_integer = r;
    int needblock = (-prec)/16 +1;
    int neededblock = needblock - s/16;


    for(int i=0; i<neededblock; i++)
      result_integer = result_integer*pow(2,16) + dis(gen);


    bits = swrite(result_integer);
    bitlength = needblock * 16;

    modify_cached(bits);
    modify_cached(bitlength);


    sizetype_set(err,1,-bitlength);
    randreal = RATIONAL(INTEGER(bits), pow(2,bitlength));
    randreal.adderror(err);
  }


//  in continuous section, return (inconsistent) random number with no using chache
  else
  {
    result_bits = std::to_string(dis(gen));
    result_size = 16;
    r = result_bits;
    s = result_size;

    INTEGER result_integer = r;
    needblock = (-prec)/16 +1;
    neededblock = needblock - s/16;


    for(int i=0; i<neededblock; i++)
      result_integer = result_integer * pow(2,16) + dis(gen);


    bits = swrite(result_integer);
    bitlength = needblock * 16;

    sizetype_set(err,1,-bitlength);
    randreal = RATIONAL(INTEGER(bits), pow(2,bitlength));
    randreal.adderror(err);
  }
}


REALRAND::~REALRAND()
  {
  }

REAL REALRAND::asREAL()
  {
    return randreal;
  }


}
