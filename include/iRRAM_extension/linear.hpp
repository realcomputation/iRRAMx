#pragma once
// iRRAM extension for linear algebra.


#include <iRRAM/lib.h>
#include <vector>

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
REALMATRIX Givens (unsigned int n, unsigned int i, unsigned int j, REAL  c, REAL s);
REALMATRIX concat(REALMATRIX A, REALMATRIX B);
REALMATRIX basisVector(int i, int n);
REALMATRIX subMatrix(REALMATRIX M, int r, int c);

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
REALMATRIX Householder(REALVECTOR u);

// level 3
REALMATRIX strassen (REALMATRIX, REALMATRIX); // strassen fast multiplication
REALMATRIX linearSys(REALMATRIX M, REALVECTOR b);
REALMATRIX inv(REALMATRIX M);
REALMATRIX kernel (REALMATRIX, int); // orthogonal kernel basis

std::pair<REALMATRIX, REALMATRIX> QR(REALMATRIX); // QR decomposition
std::pair<REALMATRIX, REALMATRIX> QR_H(REALMATRIX); // QR decomposition for Hessenberg matrix
REALMATRIX hessenberg_reduction(REALMATRIX M, int p); // reduces to Hessenberg with eig perturbation 2^p


// level 4
std::vector<REAL> symm_eig (REALMATRIX, int); // approximate eigenvalues of symmetric matrix
REALMATRIX eigenVector(REALMATRIX A, REAL e, int num); //


}
