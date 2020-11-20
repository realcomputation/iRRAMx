#pragma once
// iRRAM extension for linear algebra.


#include <iRRAM/lib.h>
#include <vector>

using namespace iRRAM;
namespace iRRAM{

// utilities
std::string to_string_double (const REALMATRIX M);


// elementary ops
REALMATRIX transpose(REALMATRIX M);
REAL trace(REALMATRIX M);

// constructions of special matrices
REALMATRIX block_matrix(REALMATRIX A, REALMATRIX B, REALMATRIX C, REALMATRIX D); // [A B ; C D]
REALMATRIX block_diag_matrix(REALMATRIX A, REALMATRIX B); // [A 0 ; 0 B]

// matrix algorithms:
REALMATRIX strassen (REALMATRIX, REALMATRIX); // strassen fast multiplication


// Decompositions
REALMATRIX QRDecomposition(REALMATRIX X); // when X is regular

// Hessenberg matrices
REALMATRIX hessenberg_reduction(REALMATRIX M, int p); // reduces a matrix
std::pair<REALMATRIX, REALMATRIX> Hessenberg_QR_decomposition(REALMATRIX H);

// linear system
// - Regular matrix:
REALMATRIX gelim (REALMATRIX); // gaussian elimination for regular matrix
REAL determinant(REALMATRIX M);
REALMATRIX linearSys(REALMATRIX M, REALMATRIX b);
REALMATRIX inv(REALMATRIX M);

// - General matrix:
REALMATRIX kernel (REALMATRIX, int); // orthogonal kernel basis
REALMATRIX gelim (REALMATRIX, int); // gaussian elimination for irrelgualr matrix


// Symmetric Eigenproblem
std::vector<REAL> symm_eig (REALMATRIX, int); // approximate eigenvalues of symmetric matrix
REALMATRIX eigenVector(REALMATRIX A, REAL e, int num); //



}
