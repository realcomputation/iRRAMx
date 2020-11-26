#pragma once
#include "iRRAM_extension/fix.hpp"
// iRRAM extension for linear algebra.


#include <iRRAM/lib.h>
#include <vector>

#include "iRRAM_extension/utility.hpp"

using namespace iRRAM;
namespace iRRAM{

typedef REALMATRIX REALVECTOR;

// utilities
std::string to_string_double (REALMATRIX M);


// level 1 where no computation is required
REALMATRIX rswap(REALMATRIX, int, int);
REALMATRIX cswap(REALMATRIX, int, int);
REALVECTOR cvec(REALMATRIX, int);
REALVECTOR cvec(REALMATRIX, int, int);
REALMATRIX transpose(REALMATRIX M);
REALMATRIX block_matrix(REALMATRIX A, REALMATRIX B, REALMATRIX C, REALMATRIX D); // [A B ; C D]
REALMATRIX block_diag_matrix(REALMATRIX A, REALMATRIX B); // [A 0 ; 0 B]
REALMATRIX givens (unsigned int n, unsigned int i, unsigned int j, REAL  c, REAL s);
REALMATRIX rconcat(REALMATRIX A, REALMATRIX B); //[A B]
REALMATRIX cconcat(REALMATRIX A, REALMATRIX B); //[A ; B]
REALMATRIX basis_vec(int i, int n);
REALMATRIX submatrix(REALMATRIX M, int r, int c);

// level 2 where some computation is required
REAL inner(REALVECTOR, REALVECTOR);
REAL trace(REALMATRIX M);
REAL det(REALMATRIX M); // defined via limit on det_approx. it is a total function.
REAL Enorm(REALMATRIX M);
REAL Fnorm(REALMATRIX M);
REALVECTOR projection(REALVECTOR u, REALVECTOR v);
REALMATRIX block_decompose(REALMATRIX A, unsigned int i, unsigned int j);
void block_decompose(REALMATRIX A, unsigned int i, unsigned int j, REALMATRIX* L);
REALMATRIX normalize(REALVECTOR u);
REALMATRIX householder(REALVECTOR u);

// level 3
REALMATRIX strassen (REALMATRIX, REALMATRIX); // strassen fast multiplication
REALMATRIX linear_sys(REALMATRIX M, REALVECTOR b); // when M is regular
REALMATRIX inv(REALMATRIX M);
REALMATRIX kernel (REALMATRIX, int); // orthogonal kernel basis where the second arg. is the dimension

std::pair<REALMATRIX, REALMATRIX> QR(REALMATRIX); // QR decomposition
std::pair<REALMATRIX, REALMATRIX> QR_H(REALMATRIX); // QR decomposition for Hessenberg matrix
REALMATRIX hessenberg_reduction(REALMATRIX M, int p); // reduces to Hessenberg with eig perturbation 2^p


// level 4
REALVECTOR symm_eig (REALMATRIX); //  eigenvalues of symmetric matrix
REALMATRIX eigen_vec(REALMATRIX A, REAL e, int num); //


}
