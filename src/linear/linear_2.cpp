// #include <iRRAM.h>
#include <cassert>

#include "iRRAM_extension/linear.hpp"

using namespace iRRAM;
namespace iRRAM{

REAL prec(int p){
  return scale(REAL(1), p);
}

REAL inner(REALVECTOR u, REALVECTOR v)
{
	assert(u.maxrow == v.maxrow);
	assert(u.maxcolumn == 1);
	assert(v.maxcolumn == 1);

	REAL sum = 0;
	for(int i=0; i<(int)u.maxrow; i++)
	{
		sum += u(i,0) * v(i,0);
	}
	return sum;
}



REALVECTOR projection(REALVECTOR u, REALVECTOR v)
{
	return u * (inner(v,u)/inner(u,u));
}


REAL Enorm(REALMATRIX M){
  REAL sum = 0;
  for (unsigned int i = 0; i<M.maxrow; i++)
    for (unsigned int j = 0; j<M.maxcolumn; i++)
      sum = sum + M(i, j) * M(i, j);
  sum = sqrt(sum);
  return sum;
}


REAL Fnorm(REALMATRIX M)
{
    REAL norm = 0;
    for (unsigned int i = 0; i < M.maxrow; i++)
	for (unsigned int j = 0; j < M.maxcolumn; j++)
	    norm += M(i,j)*M(i,j);
    return sqrt(norm);
}


// householder reflector
REALMATRIX Householder(REALVECTOR u){
//u is a column vector
  REALMATRIX X(u.maxrow, u.maxrow);
  for(unsigned int i = 0; i < u.maxrow; i++)
    for(unsigned int j = 0; j < u.maxrow; j++)
      if (i == j)
        X(i, j) = 1 - 2 * u(i, 0) * u(j, 0);
      else
        X(i, j) = - 2 * u(i, 0) * u(j, 0);
  return X;
}

REALMATRIX normalize(REALVECTOR u){
  REALMATRIX X(u.maxrow, 1);
  REAL sum = 0;
  for (unsigned int i = 0; i<u.maxrow; i++)
    sum = sum + abs(u(i, 0) * u(i, 0));
  sum = sqrt(sum);
  for (unsigned int i = 0; i<u.maxrow; i++)
    X(i, 0) = u(i, 0) / sum;
  return X;
}



REALMATRIX neg(REALMATRIX M)
{
    return zeroes(M.maxrow, M.maxcolumn) - M;
}


REAL trace(REALMATRIX mat)
{
	REAL tr = 0;
	for (int i=0; i<(int)mat.maxrow; i++)
		tr += mat(i,i);
	return tr;
}


REALMATRIX VSimilar(REALMATRIX M)
{
	REALMATRIX V = REALMATRIX(M.maxrow, M.maxcolumn);
	for(int r = 0; r<(int)M.maxrow; r ++)
	{
		for(int c = 0; c<(int)M.maxcolumn; c++)
		{
			V(r,c) = power(REAL(r+1), c);
		}
	}
	return (inv(V)*M)*V;
}



}
