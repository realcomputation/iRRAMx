// #include <iRRAM.h>
#include <cassert>

#include "iRRAMx/linear.hpp"

using namespace iRRAM;
namespace iRRAM{


REALMATRIX rswap(REALMATRIX A, int i, int j){
  REAL tmp;
  REALMATRIX X = A;
  for(unsigned int c = 0; c<X.maxcolumn; c++){
    tmp = X(i, c);
    X(i, c) = X(j, c);
    X(j, c) = tmp;
  }
  return X;
}

REALMATRIX cswap(REALMATRIX A, int i, int j){
  REAL tmp;
  REALMATRIX X = A;
  for(unsigned int c = 0; c<X.maxrow; c++){
    tmp = X(c, i);
    X(c, i) = X(c, j);
    X(c, j) = tmp;
  }
  return X;
}

REALVECTOR cvec(REALMATRIX u, int idx)
{
	REALMATRIX v = REALMATRIX(u.maxrow, 1);
	for (int i=0; i<(int) u.maxrow; i++)
	{
		v(i,0) = u(i,idx);
	}
	return v;
}
//
// REALVECTOR cvec(REALMATRIX u, int r, int c)
// {
// 	REALMATRIX v = REALMATRIX(u.maxrow - r, 1);
// 	for (int i=r; i<(int) u.maxrow; i++)
// 	{
// 		v(i,0) = u(i,c);
// 	}
// 	return v;
// }

// return colum matrix below some index (i,j)
REALMATRIX cvec(REALMATRIX M, int i, int j)
{
    REALMATRIX A = REALMATRIX(M.maxrow - i, 1);
    for (int k = i; k < (int) M.maxrow; k++)
	A(k-i,0) = M(k, j);
    return A;
}

// [A B]
// [C D]
REALMATRIX block_matrix(REALMATRIX A, REALMATRIX B, REALMATRIX C, REALMATRIX D)
{
    int row = (int) (A.maxrow + C.maxrow);
    int col = (int) (A.maxcolumn + B.maxcolumn);
    REALMATRIX S = REALMATRIX(row, col);
    for(int i = 0; i < row; i ++)
    {
	for(int j = 0; j < col; j++)
	{
	    if (i < (int) A.maxrow && j < (int) A.maxcolumn)
		S(i,j) = A(i,j);
	    if (i < (int) A.maxrow && j >= (int) A.maxcolumn)
		S(i,j) = B(i,j - A.maxcolumn);
	    if (i >= (int) A.maxrow && j < (int) A.maxcolumn)
		S(i,j) = C(i - A.maxrow,j);
	    if (i >= (int) A.maxrow && j >= (int) A.maxcolumn)
		S(i,j) = D(i - A.maxrow, j - A.maxcolumn);
	}
    }
    return S;
}

// [A 0]
// [0 B]
REALMATRIX block_diag_matrix(REALMATRIX A, REALMATRIX B)
{
    int n = (int) (A.maxrow + B.maxrow);
    REALMATRIX S = REALMATRIX(n, n);
    for(int i = 0; i < n; i ++)
    {
	for(int j = 0; j < n; j++)
	{
	    if (i < (int) A.maxrow && j < (int) A.maxcolumn)
		S(i,j) = A(i,j);
	    else
	    {
		if (i >= (int) A.maxrow && j >= (int) A.maxcolumn)
		    S(i,j) = B(i - A.maxrow, j - A.maxcolumn);
		else
		    S(i,j) = 0;
	    }
	}
    }
    return S;
}

// [A B C]
// [E F G]
// [H I J] => F (M(ii) ~ M(j-1,j-1))
REALMATRIX block_decompose(REALMATRIX A, unsigned int i, unsigned int j)
{
    REALMATRIX M = REALMATRIX(j-i, j-i);
    for (unsigned int n = i; n < j; n++)
	for (unsigned int m = i; m < j; m++)
	    M(n-i,m-i) = A(n,m);
    return M;
}

void block_decompose(REALMATRIX A, unsigned int i, unsigned int j, REALMATRIX* L)
{
    int r = A.maxrow;
    int c = A.maxcolumn;
    L[0] = REALMATRIX(i,j);
    L[1] = REALMATRIX(i,c-j);
    L[2] = REALMATRIX(r-i,j);
    L[3] = REALMATRIX(r-i,c-j);
    for(int n=0; n<r; n++)
	for(int m=0; m<c; m++)
	    if(n<i)
		if(m<j)
		    L[0](n,m)=A(n,m);
		else
		    L[1](n,m-j)=A(n,m);
	    else
		if(m<j)
		    L[2](n-i,m)=A(n,m);
		else
		    L[3](n-i,m-j)=A(n,m);
}

REALMATRIX transpose(REALMATRIX M)
{
    REALMATRIX N = REALMATRIX(M.maxcolumn, M.maxrow);
    for (unsigned int i = 0; i < M.maxcolumn; i++)
	for (unsigned int j = 0; j < M.maxrow; j++)
	    N(i,j) = M(j,i);
    return N;
}



/* QR decomposition and QR iterations
 * As matrices we deal with is now tridiagonal matrices,
 * we use Givens rotations to decompose the matrixes.
 */
// Constructs Givens matrix
// i < j <= n, c^2 + s^2 = 1
REALMATRIX givens(unsigned int n, unsigned int i, unsigned int j, REAL  c, REAL s)
{
    REALMATRIX G = eye(n);
    G(i,i) = c;
    G(j,j) = c;
    G(j,i) = s;
    G(i,j) = - s;
    return G;
}



REALMATRIX basis_vec(int i, int n)
{
	REALMATRIX M = REALMATRIX(n, 1);
	M(i,0) = 1;
	return M;
}

REALMATRIX submatrix(REALMATRIX M, int r, int c)
{
	REALMATRIX N = REALMATRIX(M.maxrow-r, M.maxcolumn-c);
	for(int i=r; i <(int)M.maxrow; i++)
		for(int j=c; j<(int)M.maxcolumn; j++)
			N(i-r, j-c) = M(i,j);
	return N;

}


REALMATRIX rconcat(REALMATRIX A, REALMATRIX B)
{
	if(A.maxrow == 0)
		return B;
	if(B.maxrow == 0)
		return A;

	REALMATRIX M(A.maxrow, A.maxcolumn + B.maxcolumn);
	for(int r = 0; r<(int)A.maxrow; r++)
	{
		for(int c = 0; c<(int)A.maxcolumn+(int)B.maxcolumn; c++)
		{
			if (c<(int)A.maxcolumn)
			{
				M(r,c) = A(r,c);
			}
			else
			{
				M(r,c) = B(r,c-A.maxcolumn);
			}
		}
	}
	return M;
}


REALMATRIX cconcat(REALMATRIX A, REALMATRIX B)
{
	if(A.maxcolumn == 0)
		return B;
	if(B.maxcolumn == 0)
		return A;

	REALMATRIX M(A.maxrow + B.maxrow, A.maxcolumn);
	for(int c = 0; c<(int)A.maxcolumn; c++)
	{
		for(int r = 0; r<(int)A.maxrow+(int)B.maxrow; r++)
		{
			if (r<(int)A.maxrow)
			{
				M(r,c) = A(r,c);
			}
			else
			{
				M(r,c) = B(r-A.maxrow,c);
			}
		}
	}
	return M;
}


}
