// #include <iRRAM.h>
#include <cassert>

#include "iRRAMx/linear.hpp"

using namespace iRRAM;
namespace iRRAM{



double rtod (REAL r){
  return std::stod(swrite(r, 50));
}
std::vector<std::vector<double> > real_mat_to_double_mat(REALMATRIX T){
    int d = T.maxrow;
    int r = T.maxcolumn;

    std::vector<std::vector<double> > approx;

    for (int i = 0; i < d; i ++){
      approx.push_back({});
      for (int j = 0; j < r; j ++){
        approx[i].push_back(0);
      }
    }

    for (int i = 0; i < d; i ++){
      for (int j = 0; j < r; j ++){
        approx[i][j] = rtod(T(i, j));
      }
    }

    return approx;
  }





  //*********



// print double matrix in .mat format
std::string to_string_d_mat (std::vector<std::vector<double> > approx){
  std::string res = "";
  int d = approx.size();
  int r = approx[0].size();
  res +="[";
  for (int i = 0; i < d; i++){
    for (int j = 0; j < r; j ++){
      res += std::to_string(approx[i][j]);
      res += " ";
    }
    if (i != d - 1){
      res += ";";
    }
  }

  res += "]";
  return res;
}


std::string to_string_double(const REALMATRIX M){
  return to_string_d_mat(real_mat_to_double_mat(M));
}


}
