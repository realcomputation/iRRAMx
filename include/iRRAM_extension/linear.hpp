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

/** @defgroup linalg Linear Algebra
 *  This is the collection of function specifications related to
 *  linear algebra
 *  @{
 */


// level 1 where no computation is required

/*! @brief elementary op. for swaping rows
 *
 *  @return a matrix which is identical to A but i'th and j'th rows are swapped.
 */

/// @brief This method adds two integers.
/// @param a First integer to add.
/// @param b Second integer to add.
/// @return The sum of both parameters.
REALMATRIX rswap(REALMATRIX A, int i, int j);

/*! @brief elementary op. for swaping columns
 *
 *  @param i a row index of A
 *  @param j a column index of A
 *  @return a matrix which is identical to A but i'th and j'th columns are swapped.
 */
REALMATRIX cswap(REALMATRIX A, int i, int j);

/*! @brief from an input matrix, get a column
 *
 *  @param i a column index of A
 *  @return a column vector which is the i'th column of A
 */
REALVECTOR cvec(REALMATRIX A, int i);

/*! @brief from an input matrix, get a column below an index
 *
 *  @param i a row index of A
 *  @param j a column index of A
 *  @return a column vector which is the column below, including, (i, j)
 */
REALVECTOR cvec(REALMATRIX A, int i, int j);

/*! @brief transpose of a matrix
 *
 */
REALMATRIX transpose(REALMATRIX M);

/*! @brief create a matrix from given submatrices
 *
 *  If A is of dimension n * m, B should be of dimension n * k for some k.
 *  And, C should be of dimension of l * m for some l and D should be of
 *  dimension l * k
 *
 *  @return a matrix which is [A B ; C D]
 */
REALMATRIX block_matrix(REALMATRIX A, REALMATRIX B, REALMATRIX C, REALMATRIX D);

/*! @brief create a matrix from given submatrices
 *
 *  @return a matrix which is [A 0 ; 0 B]
 */
REALMATRIX block_diag_matrix(REALMATRIX A, REALMATRIX B); // [A 0 ; 0 B]

/*! @brief create Given's rotation
 *
 *  return the Given's n-dimentional rotation G where
 *  G(i,i) = G(j,j) = c and
 *  G(j,i) = s = - G(i, j)
 *
 *  @param n the dimension of Given's rotation
 *  @param i 0 <= a row index < n
 *  @param j 0 <= a column index < n
 *
 *  @warning r^2 + c^2 should be 1
 *
 *  @return a matrix which is [A B ; C D]
 */
REALMATRIX givens (unsigned int n, unsigned int i, unsigned int j, REAL  c, REAL s);

/*! @brief concatnate the rows of two matrices
 *
 *  When A and B have the same number of rows, return [A B]
 *
 *  @return a matrix which is [A B ; C D]
 */
REALMATRIX rconcat(REALMATRIX A, REALMATRIX B); //[A B]

/*! @brief concatnate the columns of two matrices
 *
 *  When A and B have the same number of columns, return [A ; B]
 *
 *  @return matrix which is [A ; B]
 */
REALMATRIX cconcat(REALMATRIX A, REALMATRIX B); //[A ; B]

/*! @brief create a basis vector
 *
 *  @return a n-dimensional column vector where e(i) = 1
 */
REALMATRIX basis_vec(int i, int n);

/*! @brief extract a submatrix below an index
 *
 *  @param r a row index
 *  @param c a column index
 *
 *  @return the submatrix of M(r..., c...)
 */
REALMATRIX submatrix(REALMATRIX M, int r, int c);

// level 2 where some computation is required

/*! @brief compute the inner product of two column vectors
 *
 *  @param w a column vector
 *  @param v a column vector with the same dimension of w
 *
 *  @return the inner product of v and w
 */
REAL inner(REALVECTOR v, REALVECTOR w);

/*! @brief compute the trace of a square matrix
 *
 *  @param M a square matrix
 *
 *  @return the trace of M
 */
REAL trace(REALMATRIX M);

/*! @brief compute the determinant of a square matrix
 *
 *  @param M a square matrix
 *
 *  @return the determinant of M
 */
REAL det(REALMATRIX M); // defined via limit on det_approx. it is a total function.

/*! @brief compute the Euclidean norm of a matrix
 *
 *  @return the Euclidean norm of M
 */
REAL Enorm(REALMATRIX M);

/*! @brief compute the Frobenius norm of a matrix
 *
 *  @return the Frobenius norm of M
 */
REAL Fnorm(REALMATRIX M);
REALVECTOR projection(REALVECTOR u, REALVECTOR v);
REALMATRIX block_decompose(REALMATRIX A, unsigned int i, unsigned int j);
void block_decompose(REALMATRIX A, unsigned int i, unsigned int j, REALMATRIX* L);

/*! @brief normalize a column vector w.r.t. Euclidean norm
 *
 *  @warning this runs forever if u is the zero vector
 *
 *  @param  u column vector which is nonzero
 *  @return u normalized w.r.t. Euclidean norm
 */
REALMATRIX normalize(REALVECTOR u);

/*! @brief create Householder reflector from u
 *
 *  @warning this runs forever if u is the zero vector
 *
 *  @param  u unit column vector which is nonzero
 *  @return Householder reflector I - 2 * u * u^T
 */
REALMATRIX householder(REALVECTOR u);

// level 3

/*! @brief multiply two square matrices using Strassen algorithm
 *
 *
 *  @param  A a square matrix whose dimension is of B
 *  @param  B a square matrix whose dimension is of A
 *
 *  @return A*B
 */
REALMATRIX strassen (REALMATRIX A, REALMATRIX B); // strassen fast multiplication

/*! @brief Solve a linear system
 *
 *  @warning this runs forever if M is singular matrix
 *
 *  @param  M a regular square matrix
 *  @param  b a column vector
 *
 *  @return a column vector x where M * x = b
 */
REALMATRIX linear_sys(REALMATRIX M, REALVECTOR b); // when M is regular

/*! @brief compute the inverse of a square regular matrix
 *
 *  @warning this runs forever if M singular matrix
 *
 *  @param  M a regular square matrix
 *
 *  @return the inverse of M
 */
REALMATRIX inv(REALMATRIX M);

/*! @brief compute the basis of the kernel
 *
 *  @param  M a regular square matrix
 *  @param  k the dimension of the kernel of M
 *
 *  @return a matrix whose columns are orthogonal and span the kernel
 */
REALMATRIX kernel (REALMATRIX M, int k); // orthogonal kernel basis where the second arg. is the dimension

/*! @brief QR Factorization on a square matrix
 *
 *  @warning this runs forever if M is singular
 *
 *  @param  M a regular square matrix
 *
 *  @return a pair (Q, R) where Q is an orthogonal matrix and R is an upper-triangular matrix s.t. M = Q*R
 */
std::pair<REALMATRIX, REALMATRIX> QR(REALMATRIX M); // QR decomposition

/*! @brief QR Factorization on a square Hassenberg matrix
 *
 *  @warning this runs forever if M is singular
 *
 *  @param  M a regular square Hassenberg matrix
 *
 *  @return a pair (Q, R) where Q is an orthogonal matrix and R is an upper-triangular matrix s.t. M = Q*R
 */
std::pair<REALMATRIX, REALMATRIX> QR_H(REALMATRIX); // QR decomposition for Hessenberg matrix

/*! @brief approximately reduce a square matrix to Hessenberg matrix
 *
 *  For a square matrix M, compute a Hessenberg matrix H.
 *  It holds that the eigenvalues of H are perturbed from the eigenvalues of M
 *  by at most 2^p
 *
 *  @param M a square matrix
 *  @param p negative integer
 *
 *  @return a Hessenberg matrix whose eigenvalues approximates the eigenvalues of M
 */
REALMATRIX hessenberg_reduction(REALMATRIX M, int p); // reduces to Hessenberg with eig perturbation 2^p


// level 4

/*! @brief compute the eigenvalues of a symmetric matrix
 *
 *  @param M a symmetric matrix
 *
 *  @return a column vector which consists of the eigenvalues of M
 */
REALVECTOR symm_eig (REALMATRIX M); //  eigenvalues of symmetric matrix

/*! @brief compute the egenvectors of a symmetric matrix associated with an eigenvalue
 *
 *  @param M a symmetric matrix
 *  @param e an eigenvalue of M
 *  @param num the multiplicity of e
 *
 *  @return a matrix whose columns are orthogonal and span the eigenspace associated with e
 */
REALMATRIX eigen_vec(REALMATRIX M, REAL e, int num); //

/** @} */ 

}
