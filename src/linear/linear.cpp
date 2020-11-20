// #include <iRRAM.h>
#include <cassert>

#include "iRRAM_extension/linear.hpp"

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


REAL prec(int p){
  return scale(REAL(1), p);
}

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


REAL Fnorm(REALMATRIX M)
{
    REAL norm = 0;
    for (unsigned int i = 0; i < M.maxrow; i++)
	for (unsigned int j = 0; j < M.maxcolumn; j++)
	    norm += M(i,j)*M(i,j);
    return sqrt(norm);
}

REAL Cnorm(REALMATRIX M)
{
    return 0;
}

REAL Rnorm(REALMATRIX M)
{
    return 0;
}



/*
 * Some matrix operation utilities
 * Some of these should be in
 * REALMATRIXs definiton
 */
REALMATRIX neg(REALMATRIX M)

{
    return zeroes(M.maxrow, M.maxcolumn) - M;
}

REALMATRIX transpose(REALMATRIX M)
{
    REALMATRIX N = REALMATRIX(M.maxcolumn, M.maxrow);
    for (unsigned int i = 0; i < M.maxcolumn; i++)
	for (unsigned int j = 0; j < M.maxrow; j++)
	    N(i,j) = M(j,i);
    return N;
}

// return colum matrix below some index (i,j)
REALMATRIX colum(REALMATRIX M, int i, int j)
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


REAL trace(REALMATRIX mat)
{
	REAL tr = 0;
	for (int i=0; i<(int)mat.maxrow; i++)
		tr += mat(i,i);
	return tr;
}


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

  	c = colum(M, k+1, k);
  	s = Fnorm(c);
  	c(0, 0) = c(0, 0) - s;
  	s = Fnorm(c);
  	u = c / s;
  	P = eye(n - k - 1) - 2 * u * transpose(u);
  	P = block_diag_matrix(eye(k+1), P);
  	M = mult(P, mult(M, P));
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
  	c = colum(M, k+1, k);
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
      L[3] = mult(P, mult(L[3], P));
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



/* QR decomposition and QR iterations
 * As matrices we deal with is now tridiagonal matrices,
 * we use Givens rotations to decompose the matrixes.
 */
// Constructs Givens matrix
// i < j <= n, c^2 + s^2 = 1
REALMATRIX Givens (unsigned int n, unsigned int i, unsigned int j, REAL  c, REAL s)
{
    REALMATRIX G = eye(n);
    G(i,i) = c;
    G(j,j) = c;
    G(j,i) = s;
    G(i,j) = - s;
    return G;
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

// QR step
REALMATRIX QR(REALMATRIX H)
{
    int n = H.maxrow;
    REALMATRIX Q = eye(n);
    REAL r, c, s;
    REAL tmp1, tmp2;
    REAL Givens[2*n-2];
    int tmp;
    for (int i = 0; i < n-1; i++)
    {
	r = sqrt(H(i,i)*H(i,i) + H(i+1, i)*H(i+1, i));
	c = H(i,i) / r;
	s = H(i+1,i) / r;
	Givens[i*2] = c;
	Givens[i*2+1] = s;
	tmp = i + 3;
	if(i == n-2)
	    tmp = i + 2;
	for(int l=i; l<tmp; l++)
	{
	    tmp1 = H(i,l); tmp2 = H(i+1,l);
	    H(i,l) = tmp1*c +tmp2*s;
	    H(i+1,l) = - tmp1*s + tmp2*c;
	}
    }

    for (int i=0; i<n-1; i++)
    {
	for (int j=i-1; j<i+2; j++)
	{
	    if(j == -1 || j == n+1)
		continue;
	    tmp1 = H(j, i) * Givens[2*i] + H(j,i+1) * Givens[2*i + 1];
	    H(j, i+1) = - H(j, i) * Givens[2*i+1] + H(j,i+1) *Givens[2*i];
	    H(j,i) = tmp1;
	}
    }

    return H;
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

    REALMATRIX H = mult(QR.second, QR.first);

    for(int i=0; i<n+1; i++)
	H(i,i) = H(i,i) + s;

    return H;
}

std::vector<REAL> tridiagonal_QR_algorithm (REALMATRIX M, int p)
{
    int n = M.maxrow;
    cout << n <<"\n";

    if(n == 1)
      return {M(0,0)};
    if(n == 2){
      REAL tr = M(0,0) + M(1, 1);
      REAL det = M(0,0) * M(1, 1) - M(1, 0) * M(0, 1);

      REAL eig1, eig2;

      eig1 = ( tr - sqrt(tr * tr - 4*det)) / 2;
      eig2 = ( tr + sqrt(tr * tr - 4*det)) / 2;

      return {eig1, eig2};


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
    return eigens;
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


std::vector<REAL> symm_eig(REALMATRIX M, int p)
{

  std::vector<REAL> eigen, eigens;
  int n = M.maxrow;
  int ind = 0;
  eigens.reserve(n);
  for(int i=0; i<n; i++)
  	eigens.push_back(0);

  M = hessenberg_reduction(M, p - 3);

  std::vector<REALMATRIX> T = split(M, p - 3);

  for(unsigned int k = 0; k < T.size(); k++)
  {

  	eigen = tridiagonal_QR_algorithm(T[k], p - 3);
  	for(unsigned int i = 0; i < T[k].maxrow; i++)
  	{
  	    eigens[ind] = eigen[i];
  	    ind++;
  	}
  }
  return mysort(eigens, p);
}















//




//



REALMATRIX power(REALMATRIX mat, int k)
{
	REALMATRIX tmp = mat;
	for(int i=0; i<k; i++)
		tmp = tmp * mat;
	return tmp;

}




REALMATRIX concat(REALMATRIX A, REALMATRIX B)
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


REAL inner(REALMATRIX u, REALMATRIX v)
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

REALMATRIX colVector(REALMATRIX u, int idx)
{
	REALMATRIX v = REALMATRIX(u.maxrow, 1);
	for (int i=0; i<(int) u.maxrow; i++)
	{
		v(i,0) = u(i,idx);
	}
	return v;

}

REALMATRIX colVector(REALMATRIX u, int r, int c)
{
	REALMATRIX v = REALMATRIX(u.maxrow - r, 1);
	for (int i=r; i<(int) u.maxrow; i++)
	{
		v(i,0) = u(i,c);
	}
	return v;

}

REALMATRIX projection(REALMATRIX u, REALMATRIX v)
{
	return u * (inner(v,u)/inner(u,u));
}

// Orthonormalize columnvectors of M
REALMATRIX GramSchmidt(REALMATRIX M)
{
	REALMATRIX Q;
	REALMATRIX v;
	for(int i=0; i<(int) M.maxcolumn; i++)
	{
		v = colVector(M,i);
		for(int j=0; j < i; j++)
		{
			v = v - projection(colVector(Q, j), colVector(M, i));
		}
		v = v /sqrt(inner(v,v));
		Q = concat(Q, v);
	}
	return Q;
}


REALMATRIX linearSys(REALMATRIX M, REALMATRIX b)
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


REAL determinant(REALMATRIX M)
{
	REALMATRIX W = M;
	REAL det = 1;
	REAL tmp = 0;
	REAL tmpV;
	for(int i=0; i<(int)W.maxrow; i++)
	{
		tmp = 0;
		// finding half maximum index
		for(int t = i; t<(int)W.maxrow; t++)
			tmp = maximum(abs(W(t,i)), tmp);
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

				det *= -1;

				for(int j=i+1; j<(int)W.maxrow; j++)
				{
					tmpV = W(j,i)/W(i,i);
					for(int tc=0; tc<(int)W.maxcolumn; tc++)
					{
						W(j,tc) = W(j,tc) - W(i,tc)*tmpV;
					}
				}
			break;
			}
		}
	}
	for(int i=0; i<(int)M.maxrow; i++)
		det *= W(i,i);
	return det;
}


REALMATRIX inv(REALMATRIX M)
{

	REALMATRIX INV = REALMATRIX(M.maxcolumn,M.maxcolumn);
	for (int i = 0 ; i < (int)M.maxcolumn; i ++)
		INV(i,i) = 1;
	REALMATRIX ins = linearSys(M,INV);
	return ins;
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


REALMATRIX eigenVector(REALMATRIX A, REAL eigenValue, int nulldim)
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
  return eigenVector(A, 0, nulldim);
}


REALMATRIX basisVector(int i, int n)
{
	REALMATRIX M = REALMATRIX(n, 1);
	M(i,0) = 1;
	return M;
}

REALMATRIX subMatrix(REALMATRIX M, int r, int c)
{
	REALMATRIX N = REALMATRIX(M.maxrow-r, M.maxcolumn-c);
	for(int i=r; i <(int)M.maxrow; i++)
		for(int j=c; j<(int)M.maxcolumn; j++)
			N(i-r, j-c) = M(i,j);
	return N;

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
		x = colVector(Y, 0);
		alpha = sqrt(inner(x, x));
		u = x - alpha * basisVector(0, n - i);
		v = u / sqrt(inner(u,u));
		Q[i] = eye(n-i) - 2 * v * transpose(v);
		Y = Q[i]*Y;
		rr[i] = Y(0,0);
		Y = subMatrix(Y, 1,1);
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
//



REALMATRIX gelim(REALMATRIX A, int k)
{

  int n = A.maxrow;
  int pi, pj;
  REAL m;
  REAL temp;
  REAL div[n];
  for(int i=0; i<k-1; i++)
    {
      m = 0;
      for(int j=i; j<n; j++)
	for(int k=i; k<n; k++)
	  m = maximum(m, abs(A(j,k)));
      for(int j=i;j<n;j++)
	for(int k=i; k<n; k++)
	  if (choose(abs(A(j,k)) > m / 2, abs(A(j,k)) < m) == 1)
	    { pi = j; pj = k;}

      for(int j=0; j<n; j++)
	{
	  temp = A(j,i);
	  A(j,i) = A(j, pj);
	  A(j,pj) = temp;
	}
      for(int j=0; j<n; j++)
	{
	  temp = A(i,j);
	  A(i,j) = A(pi, j);
	  A(pi, j) = temp;
	}
      for(int j=i+1; j<n; j++)
	div[j] = A(i,j)/A(i,i);
      div[i] = 1;

      for(int j=i+1; j<n; j++)
	{
	  for(int k=i+1; k<n; k++)
	    A(j,k) = A(j,k) - div[k] * A(j,i);
	  A(j,i) = 0;
	}

    }
  return A;
}





REALMATRIX gelim(REALMATRIX A)
{
  return gelim(A, A.maxrow);
}
}
