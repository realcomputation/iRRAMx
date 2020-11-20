#include <iRRAM.h>
#include "iRRAM_extension/linear.hpp"

using namespace iRRAM;
void compute () {

  REALMATRIX A(4,4);
  for (int i = 0; i < 4; i ++)
    for (int j= 0; j < 4; j ++)
      A(i, j) = i + j;

  cout <<"The 2^-100 approximations of the eigenvalues of the matrix:"<<"\n";
  cout <<to_string_double(A) <<"\n====";
  std::vector<REAL> eig = symm_eig(A, 100);
  for (REAL& r: eig)
    cout << r << "\n";

}
