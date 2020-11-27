// #include <iRRAM.h>
#include <cassert>

#include "iRRAMx/linear.hpp"

using namespace iRRAM;
namespace iRRAM{


/* Reducing a symmetric matrix into a Tridiagonal matrix.
 * First one requires the matrix to be special.
 * Second one works for every symmetric matrices however,
 * it returns a tridiagonal matrix which approximates
 * eigenvalues of the original one. it is not a similarity transformation.
 * The natural question arises: if a hessenberg reduction computable...
 */
REALMATRIX hessenberg_reduction(REALMATRIX M)
{
  REALMATRIX c, u, P;
  REAL s;

  int n = (int) M.maxrow;
  for(int k = 0; k < n-2; k++)
  {

  	c = cvec(M, k+1, k);
  	s = Fnorm(c);
  	c(0, 0) = c(0, 0) - s;
  	s = Fnorm(c);
  	u = c / s;
  	P = eye(n - k - 1) - 2 * u * transpose(u);
  	P = block_diag_matrix(eye(k+1), P);
  	M = strassen(P, strassen(M, P));
  }
    return M;
}

// output: Hessenberg matrix H whose eigenvalues are perturbed
// at most by 2^p (p is considered as a negative number)
REALMATRIX hessenberg_reduction(REALMATRIX M, int p)
{
  REALMATRIX c, u, P;
  REALMATRIX L[4];
  REAL s, q;
  REAL approx = prec(p - 1);
  int n = (int) M.maxrow;
  REAL err;
  for(int k = 0; k < n-2; k++)
  {
  	c = cvec(M, k+1, k);
  	s = Fnorm(c);
  	//    q = sqrt(s*s - c(0,0)*c(0,0));
  	err = 0;
  	err = abs(c(1,0));
  	for (int i=1; i< c.maxrow; i++)
  	    err = maximum(err, abs(c(i,0)));

  	if (choose(err < approx / n / n, err > approx / 2 / n / n ) == 2)
  	{
      c(0, 0) = c(0, 0) - s;
      s = Fnorm(c);
      u = c / s;
      P = eye(n - k - 1) - 2 * u * transpose(u);
      //    P = block_diag_matrix(eye(k+1), P);
      block_decompose (M, k+1, k+1, L);
      L[1] = L[1]*P;
      L[2] = P*L[2];
      L[3] = strassen(P, strassen(L[3], P));
      for(int i=0; i<n; i++)
    		for(int j=0; j<n; j++)
    		    if(i < k+1)
    			if(j > k)
    			    M(i,j) = L[1](i,j-k-1);
    			else
    			    continue;
    		    else
    			if(j < k+1)
    			    M(i,j) = L[2](i-k-1, j);
    			else
    			    M(i,j) = L[3](i-k-1, j-k-1);
  	}
  	else
  	{
  	    for(int i = k+2; i < n; i++)
  	    {
  		M(k, i) -= M(i, k);
  		M(i, k) = 0;
  	    }

  	}
  }
    return M;
}

// 2^{-p} reduce the hessenberg matrix M
std::vector<REALMATRIX> split(REALMATRIX M, int p)
{
  REAL abs_value;
  int n = (int) M.maxrow;
  unsigned int prev = 0;
  REAL approx = prec(p -1) / n;
  std::vector<REALMATRIX> matrices;
  for(int i=1; i<(int) M.maxrow; i++)
  {
  	abs_value = abs(M(i, i-1));
  	if(choose(abs_value >  approx / 2, abs_value < approx) == 2)
  	{
	    matrices.push_back(block_decompose(M, prev, i));
	    prev = i;
  	}
  }
  matrices.push_back(block_decompose(M, prev, n));
  return matrices;
}


// QR Decomposition on Hessenberg matrix via Givens rotation
std::pair<REALMATRIX, REALMATRIX> Hessenberg_QR_decomposition(REALMATRIX H)
{
  int n = H.maxrow;
  REALMATRIX Q = eye(n);
  REAL r, c, s;
  REAL tmp1, tmp2;
  for (int i = 0; i < n-1; i++)
  {
  	r = sqrt(H(i,i)*H(i,i) + H(i+1, i)*H(i+1, i));
  	c = H(i,i) / r;
  	s = H(i+1,i) / r;
  	for(int l=0; l<n; l++)
  	{
	    tmp1 = Q(l,i); tmp2 = Q(l,i+1);
	    Q(l,i) = tmp1*c + tmp2*s;
	    Q(l,i+1) = - tmp1 * s + tmp2 * c;

	    tmp1 = H(i,l); tmp2 = H(i+1,l);
	    H(i,l) = tmp1*c +tmp2*s;
	    H(i+1,l) = - tmp1*s + tmp2*c;
  	}
  }
  return std::pair<REALMATRIX, REALMATRIX>(Q, H);
}

int sign(REAL x, REAL p)
{
    if(choose(x < p, x > - p) == 1)
	return - 1;
    return 1;
}

REALMATRIX Hessenberg_QR_step(REALMATRIX M)
{
    int n = M.maxrow-1;

    // Soft Wilkinson shift with k = 1/8
    REAL delta = (M(n-1,n-1) - M(n,n))/2;
    REAL s = M(n,n) - sign(delta, M(n,n-1) / 16) * M(n,n-1) *  M(n,n-1) / (abs(delta) + sqrt(delta*delta + M(n,n-1)*M(n,n-1)));

    for(int i=0; i<n+1; i++)
    {
	M(i,i) = M(i,i) - s;
    }

    std::pair<REALMATRIX, REALMATRIX> QR = Hessenberg_QR_decomposition(M);

    REALMATRIX H = strassen(QR.second, QR.first);

    for(int i=0; i<n+1; i++)
	H(i,i) = H(i,i) + s;

    return H;
}

REALVECTOR tridiagonal_QR_algorithm (REALMATRIX M, int p)
{
    int n = M.maxrow;
    cout << n <<"\n";

    if(n == 1){
      REALVECTOR E(1, 1);
      E(0, 0) = M(0,0);
      return E;
    }
    if(n == 2){
      REAL tr = M(0,0) + M(1, 1);
      REAL det = M(0,0) * M(1, 1) - M(1, 0) * M(0, 1);

      REAL eig1, eig2;

      eig1 = ( tr - sqrt(tr * tr - 4*det)) / 2;
      eig2 = ( tr + sqrt(tr * tr - 4*det)) / 2;
      REALVECTOR E(2, 1);
      E(0, 0) = eig1;
      E(1, 0) = eig2;
      return E;
    }



    std::vector<REAL> eigens = std::vector<REAL>();
    eigens.reserve(n);
    for (int j = 0; j < n; j++)
	eigens.push_back(0);

    REAL error = power(2, p - 1) / n;
    // RATIONAL error = err (-p+1) / n;
    REALMATRIX H = M;
    int m = n - 1;
    for (int i = 0; i < n-2; i++)
    {
	while(choose(abs(H(m,m-1)) > error / 2, abs(H(m,m-1)) < error) == 1)
	{
	    H = Hessenberg_QR_step(H);
	}
	eigens[i] = H(m,m);
	H = block_decompose(H, 0, m);
	m -= 1;
    }


    REAL a1 = H(0,0);
    REAL a2 = H(1,1);
    REAL b = H(1,0);
    REAL tmp = sqrt((a1-a2)*(a1-a2) + 4*b*b);
    REAL l1 = (a1 + a2 + tmp)/2;
    REAL l2 = (a1 + a2 - tmp)/2;
    eigens[n-2] = l1;
    eigens[n-1] = l2;


    REALVECTOR E(n, 1);
    for (int i = 0; i < n; i++)
      E(i, 0) = eigens[i];
    return E;
    //
    // return eigens;
}



REALVECTOR mysort(int p, REALVECTOR X){
  REALVECTOR A = X;
  int n = A.maxrow;
  REAL tmp;
  REAL eps = prec(p);
  for (int i = 1; i < n; i ++){
    for (int j = 0; j < i; j ++){
      if (choose(A(j, 0) > A(i, 0) - eps, A(i, 0) > A(j, 0) - eps) == 1 ){
      	tmp = A(i, 0);
      	A(i, 0) = A(j, 0);
      	A(j, 0) = tmp;
      }
    }
  }
  return A;
}

std::vector<REAL> mysort( std::vector<REAL> A, int p){

  int n = A.size();
  REAL tmp;
  REAL eps = prec(p);
  for (int i = 1; i < n; i ++){
    for (int j = 0; j < i; j ++){
      if (choose(A[j] > A[i] - eps, A[i] > A[j] - eps) == 1 ){
	tmp = A[i];
	A[i] = A[j];
	A[j] = tmp;
      }
    }
  }
  return A;
}


REALMATRIX symm_eig_approx(int p, REALMATRIX M)
{

  REALMATRIX eigen;
  int n = M.maxrow;


  M = hessenberg_reduction(M, p - 3);

  std::vector<REALMATRIX> T = split(M, p - 3);

  eigen = tridiagonal_QR_algorithm(T[0], p - 5);

  for(unsigned int k = 1; k < T.size(); k++)
    eigen = cconcat(eigen, tridiagonal_QR_algorithm(T[k], p - 5));

  return mysort(p-5, eigen);
}


REALVECTOR symm_eig(REALMATRIX M)
{
  return limit(symm_eig_approx, M);
}





REALMATRIX power(REALMATRIX mat, int k)
{
	REALMATRIX tmp = mat;
	for(int i=0; i<k; i++)
		tmp = tmp * mat;
	return tmp;

}









// QR decomposition for a regular matrix X
// return a normalized Q, orthogonal matrix
REALMATRIX QRDecomposition(REALMATRIX X)
{

	int n = (int) X.maxrow;
	REALMATRIX Y = X;
	REAL alpha;
	REALMATRIX u;
	REALMATRIX v;
	REALMATRIX Q[n];
	REALMATRIX x;
	REAL rr[n];
	for (int i=0; i < n-1; i ++)
	{
		x = cvec(Y, 0);
		alpha = sqrt(inner(x, x));
		u = x - alpha * basis_vec(0, n - i);
		v = u / sqrt(inner(u,u));
		Q[i] = eye(n-i) - 2 * v * transpose(v);
		Y = Q[i]*Y;
		rr[i] = Y(0,0);
		Y = submatrix(Y, 1,1);
	}
	rr[n-1] = Y(0,0);

	REALMATRIX T;
	for (int i=0; i< n-1; i++)
	{
		T = eye(n);
		for(int j=i; j<n;j++)
		{
			for(int k=i;k<n;k++)
			{
				T(j,k) = Q[i](k-i,j-i);
			}
		}
		Q[i] = T;
	}
	T = Q[0];
	for (int i=1;i<n-1;i++)
		T = T * Q[i];


	REALMATRIX lambda = REALMATRIX(n,n);
	for(int i=0;i<n;i++)
		if(rr[i] > 0)
			lambda(i,i) = 1;
		else
			lambda(i,i) = -1;

	return T*lambda;
}



//







std::pair<REALMATRIX, REALMATRIX> QR (REALMATRIX M){

  std::vector<REALMATRIX> reflectors;
  REALMATRIX A = M;
  for (unsigned int i = 0; i < A.maxcolumn; i ++ ){
    REALMATRIX y(A.maxrow, 1);
    REALMATRIX e(A.maxrow, 1);
    for(unsigned int j = 0; j < A.maxcolumn; j ++){
      if(j<i)
        y(j,0) = 0;
      else
        y(j,0) = A(j, i);
      if(i==j)
        e(j,0) = 1;
      else
        e(j,0) = 0;
    }
    REALMATRIX H = householder(normalize(y - Enorm(y) * e));
    reflectors.push_back(H);
    A = H * A;
  }
  REALMATRIX H = reflectors[0];
  for(int i =1; i<reflectors.size(); i++)
    H = H * reflectors[i];

  return std::make_pair(H, A);
}

std::pair<REALMATRIX, REALMATRIX> QR_H(REALMATRIX H){
  return Hessenberg_QR_decomposition(H);
}





}
