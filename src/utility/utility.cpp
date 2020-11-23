

#include "iRRAM_extension/utility.hpp"

using namespace iRRAM;
namespace iRRAM{

REAL prec(int p){
  return scale(REAL(1), p);
}

std::string to_string(REAL r, unsigned int p){
  return (swrite(r, p));

}

std::string to_string(COMPLEX c, unsigned int p){
  return "("+to_string(real(c), p) + " + i"+to_string(imag(c), p)+")";
}

}
