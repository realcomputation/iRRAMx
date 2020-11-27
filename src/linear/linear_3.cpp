// #include <iRRAM.h>
#include <cassert>

#include "iRRAMx/linear.hpp"

using namespace iRRAM;
namespace iRRAM{


// using hadamard bound. p is negative
REAL det_approx(int p, REALMATRIX M)
{
	REALMATRIX W = M;
	REAL det = 1;
	REAL tmp = 0;
	REAL tmpV;
  int tr, tc;
	for(int i=0; i<(int)W.maxrow; i++)
	{

    tmp = 0;
    for(unsigned int l = i; l < W.maxrow; l ++)
      for(unsigned int m = i; m < W.maxrow; m ++)
        tmp = maximum(abs(W(l, m)), tmp);

    // test if the det can be neglected..
    int b = upperbound(det)+2;
    REAL eps = prec(p + b);
    if (choose(tmp < eps, tmp > eps / 2) == 1)
      return 0;

    // there is nonzero index!
    for(unsigned int l = i; l < W.maxrow; l ++)
      for(unsigned int m = i; m < W.maxrow; m ++)
        if(choose(abs(W(l, m)) > tmp / 2, abs(W(l, m)) < tmp) == 1){
          tr = l;
          tc = m;
        }

    if (tr != i)
      det = - det;
    if (tc != i)
      det = - det;

    W = rswap(W, i, tr);
    W = cswap(W, i, tc);

		for(int j=i+1; j<(int)W.maxrow; j++)
		{
			tmpV = W(j,i)/W(i,i);
			for(int l=0; l<(int)W.maxcolumn; l++)
				W(j,l) = W(j,l) - W(i,l)*tmpV;
		}
    det = det * W(i, i);
  }
  return det;
}

REAL det(REALMATRIX M)
{
  return limit(det_approx, M);
}


// Strassen



void square_block_decompose(const REALMATRIX& x, REALMATRIX* blocks)
{
  int n = x.maxrow;
  REALMATRIX A = REALMATRIX(n/2,n/2);
  REALMATRIX B = REALMATRIX(n/2,n/2);
  REALMATRIX C = REALMATRIX(n/2,n/2);
  REALMATRIX D = REALMATRIX(n/2,n/2);
  for(int i=0; i<n; i++)
  {
    for(int j=0; j<n; j++)
    {
      if(i < n/2)
	if(j < n/2)
	  A(i,j) = x(i,j);
        else
	  B(i, j - n/2) = x(i,j);
      else
	if(j < n/2)
	  C(i - n/2, j) = x(i,j);
	else
	  D(i-n/2,j-n/2) = x(i,j);
    }
  }
  blocks[0] = A;
  blocks[1] = B;
  blocks[2] = C;
  blocks[3] = D;
}

REALMATRIX square_block_compose(REALMATRIX A, REALMATRIX B, REALMATRIX C, REALMATRIX D)
{
  int n = A.maxrow + C.maxrow;
  REALMATRIX composed = REALMATRIX(n,n);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      if(i < n/2)
	if(j < n/2)
	  composed(i,j) = A(i,j);
	else
	  composed(i,j) = B(i,j-n/2);
      else
	if(j<n/2)
	  composed(i,j) = C(i-n/2,j);
	else
	  composed(i,j) = D(i-n/2,j-n/2);
  return composed;
}

REALMATRIX add_padding(const REALMATRIX& x)
{
  REALMATRIX padded = zeroes(x.maxrow+1, x.maxcolumn+1);
  for(int i=0; i<x.maxrow; i++)
    for(int j=0; j<x.maxcolumn; j++)
      padded(i,j) = x(i,j);
  return padded;
}

REALMATRIX remove_padding(const REALMATRIX& x)
{
  REALMATRIX unpadded = REALMATRIX(x.maxrow-1, x.maxcolumn-1);
  for(int i=0; i<x.maxrow-1; i++)
    for(int j=0; j<x.maxcolumn-1; j++)
      unpadded(i,j) = x(i,j);
  return unpadded;
}

// Strassen algorithm for where x and y are square matrices
// Currently, phase transform happens for dim = 101
// Later, it should be reduced... (n = 20 sounds good)
REALMATRIX mult (const REALMATRIX& x, const REALMATRIX& y)
{
  int n = x.maxcolumn;
  if(n < 101)
    return x*y;
  bool is_odd = false;
  REALMATRIX A, B;
  if (n % 2 == 1)
  {
    is_odd = true;
    A = add_padding(x);
    B = add_padding(y);
  }
  else
  {
    A = x;
    B = y;
  }

  REALMATRIX As[4];
  REALMATRIX Bs[4];
  square_block_decompose(A, As);
  square_block_decompose(B, Bs);

  REALMATRIX M1, M2, M3, M4, M5, M6, M7;
  M1 = mult(As[0] + As[3], Bs[0] + Bs[3]);
  M2 = mult(As[2] + As[3], Bs[0]);
  M3 = mult(As[0], Bs[1] - Bs[3]);
  M4 = mult(As[3], Bs[2] - Bs[0]);
  M5 = mult(As[0] + As[1], Bs[3]);
  M6 = mult(As[2] - As[0], Bs[0] + Bs[1]);
  M7 = mult(As[1] - As[3], Bs[2] + Bs[3]);

  REALMATRIX C = square_block_compose(M1+M4-M5+M7, M3+M5, M2+M4, M1-M2+M3+M6);

  if(is_odd)
  {
    C = remove_padding(C);
  }
  return C;
}

REALMATRIX strassen (const REALMATRIX& x, const REALMATRIX& y){
  return mult(x, y);
}


//
REALMATRIX linear_sys(REALMATRIX M, REALMATRIX b)
{
	REALMATRIX W = M;
	REALMATRIX w = b;
	REAL tmpfactor;
	REAL tmp = 0;
	REAL tmpV;
	for(int i=0; i<(int)W.maxrow; i++)
	{
		// finding half maximum index
		tmp = 0;
		for(int t = i; t<(int)W.maxrow; t++)
		{
			tmp = maximum(abs(W(t,i)), tmp);
		}

		for(int t = i; t<(int)W.maxrow; t++)
		{
			if(choose(abs(W(t,i)) > tmp / 2, abs(W(t,i)) < tmp) == 1)
			{
				for(int tc=0; tc<(int)W.maxcolumn; tc++)
				{
					tmpV = W(i,tc);
					W(i,tc) = W(t,tc);
					W(t,tc) = tmpV;
				}
				for(int tc=0; tc<(int)w.maxcolumn; tc++)
				{
					tmpV = w(i,tc);
					w(i,tc) = w(t,tc);
					w(t,tc) = tmpV;
				}

				tmpfactor = 1/W(i,i);

				for(int tc=0; tc<(int)W.maxcolumn; tc++)
				{
					W(i,tc) = W(i,tc)*tmpfactor;
				}
				for(int tc=0; tc<(int)w.maxcolumn; tc++)
				{
					w(i,tc) = w(i,tc)*tmpfactor;
				}

				for(int j=0; j<(int)W.maxrow; j++)
				{
					if(i!= j)
					{
						tmpfactor = 0-W(j,i)/W(i,i);

						for(int tc=0; tc<(int)W.maxcolumn; tc++)
						{
							W(j,tc) = W(j,tc)+ W(i,tc)*tmpfactor;
						}
						for(int tc=0; tc<(int)w.maxcolumn; tc++)
						{
							w(j,tc) = w(j,tc)+ w(i,tc)*tmpfactor;
						}

					}
				}
			break;
			}
		}
	}
	return w;
}



REALMATRIX inv(REALMATRIX M)
{

	REALMATRIX INV = REALMATRIX(M.maxcolumn,M.maxcolumn);
	for (int i = 0 ; i < (int)M.maxcolumn; i ++)
		INV(i,i) = 1;
	REALMATRIX ins = linear_sys(M,INV);
	return ins;
}




// Orthonormalize columnvectors of M
REALMATRIX GramSchmidt(REALMATRIX M)
{
	REALMATRIX Q;
	REALMATRIX v;
	for(int i=0; i<(int) M.maxcolumn; i++)
	{
		v = cvec(M,i);
		for(int j=0; j < i; j++)
		{
			v = v - projection(cvec(Q, j), cvec(M, i));
		}
		v = v /sqrt(inner(v,v));
		Q = rconcat(Q, v);
	}
	return Q;
}

REALMATRIX eigen_vec(REALMATRIX A, REAL eigenValue, int nulldim)
{

	REAL max, pivotvalue, scale, temp;
	int i,j,k, temp_index, flg, ROW, COL;
	int pivoti = 0;
	int pivotj = 0;
	ROW =  (int)A.maxrow;
	COL = (int)A.maxcolumn;
	int rank = ROW - nulldim;

	REALMATRIX solutionMatrix = REALMATRIX(rank, COL - rank);
	REAL tmpV;

	for (int i=0; i<COL; i++)
	{
		A(i,i) = A(i,i)-eigenValue;
	}

	int record[COL];
    for (int i=0; i<COL; i++) record[i] = i;


	for (i=0; i < rank; i++)
	{


		max = 0;
		for(int ti = i; ti < ROW; ti ++)
		{
			for(int tj = i; tj <COL; tj ++)
			{
				max = maximum(abs(A(ti,tj)), max);
			}
		}

		flg = 0;
    	for (k=i; k<COL; k++)
    	{
    		for (j=i; j<ROW; j++)
    		{

    			if(1==choose(abs(A(j,k)) > max / 2, abs(A(j,k)) < max))
    			{
					pivotvalue = A(j,k);
    				pivoti = j;
    				pivotj = k;
    				flg = 1;
    				break;
    			}
    		}
    		if (flg==1) break;
    	}

    	if (pivoti != i)
    	{
			for(int tc=0; tc<(int)A.maxcolumn; tc++)
			{
				tmpV = A(i,tc);
				A(i,tc) = A(pivoti,tc);
				A(pivoti,tc) = tmpV;
			}
    	}

		if (pivotj != i)
		{

			for(int tc=0; tc<(int)A.maxrow; tc++)
			{
				tmpV = A(tc,i);
				A(tc,i) = A(tc,pivotj);
				A(tc,pivotj) = tmpV;
			}
            temp_index = record[i];
            record[i] = record[pivotj];
            record[pivotj] = temp_index;
		}

		for(int tc=0; tc<(int)A.maxcolumn; tc++)
		{
			A(i,tc) = A(i,tc) / pivotvalue;
		}

		for (j=0; j<ROW; j++)
		{
			if (j == i)
				continue;
			scale = A(j,i);


			for(int tc=0; tc<(int)A.maxcolumn; tc++)
			{
				A(j,tc) = A(j,tc) -scale*A(i,tc);
			}
		}
	}


	// copy ** to solution matrix.
	for(i = 0; i<rank; i++)
	{
		for(j=0; j<COL - rank; j++)
		{
			solutionMatrix(i,j) = A(i,j + rank);
		}
	}


	REALMATRIX B = REALMATRIX(ROW,nulldim);
	for(int i=0; i<COL - rank; i++)
	{
		for(int j=0; j<rank; j++)
		{
			B(record[j],i) =  solutionMatrix(j,i);
		}
		for (int k=0; k<COL - rank; k++)
		{
			if(i == k)
				B(record[rank + k], i) =  REAL(-1);
			else
				B(record[rank + k], i) =  REAL(0);
		}
	}
	B = GramSchmidt(B);

	return B;
}


REALMATRIX kernel(REALMATRIX A, int nulldim){
  return eigen_vec(A, 0, nulldim);
}




}
