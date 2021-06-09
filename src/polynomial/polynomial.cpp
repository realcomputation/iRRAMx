#include "iRRAMx/polynomial.hpp"
#include "iRRAMx/polynomial/utilities.hpp"

#include "iRRAM/lib.h"
#include "iRRAM/core.h"

#include <stack>
#include <queue>
#include <utility>
#include <vector>
#include <cmath>
using namespace iRRAM;
namespace iRRAM{

// void REALMATRIX::geterror (sizetype& error) const
// {
//   int i;
//   sizetype lerror;
//   ELEMENT((*this),0,0).geterror(error);
//   for (i=0;i<(*this).maxcolumn;i++)
//   for (j=0;j<(*this).maxrow;j++) {
//       ELEMENT((*this),i,j).geterror(lerror);
//       sizetype_max(error,error,lerror);
//   }
// }

// should be sup norm induced metric
int geterror_exp(const POLYNOMIAL & P)
{
	sizetype error;
  	P.coef[0].geterror(error);
  	int ex = error.exponent;
  	int tmp;

	for(int i=1; i< P.degree; i++)
	{
		P.coef[0].geterror(error);
		tmp = error.exponent;
	  	if (tmp >ex)
	  		ex = tmp;
	}

	return tmp;
}



POLYNOMIAL::POLYNOMIAL(int d, COMPLEX* c)
{
	std::vector<COMPLEX> co;
	co.reserve(d+1);
	for (int i = 0; i<d+1; i++)
	{
		co.emplace_back(c[i]);
	}
	degree = d;
	coef = co;
}


POLYNOMIAL::POLYNOMIAL(int d, const std::vector<COMPLEX> &v)
{
	std::vector<COMPLEX> co;
	co.reserve(d+1);
	for (int i = 0; i<d+1; i++)
	{
		co.emplace_back(v[i]);
	}
	degree = d;
	coef = co;
}

POLYNOMIAL::POLYNOMIAL(COMPLEX x)
{
	std::vector<COMPLEX> co;
	co.reserve(1);
	co.emplace_back(std::move(x));
	degree = 0;
	coef = co;
}



POLYNOMIAL::POLYNOMIAL()
{
}

POLYNOMIAL::~POLYNOMIAL()
{
}


COMPLEX POLYNOMIAL::operator () (const COMPLEX& z) const
{
	COMPLEX fx = COMPLEX(0,0);
	for(int i=0; i < degree + 1; i++)
		fx = fx + power(z,i) * coef[i];
	return fx;
}

orstream &operator<<(orstream &ors, const POLYNOMIAL &p) {
    for(int i = 0; i <= p.degree; i++) {
        ors << real(p.coef[i]) << "\t+\t" << imag(p.coef[i]) << "\ti\tx^" << i;
        if(i != p.degree) ors << "\n";
    }
    return ors;
}

POLYNOMIAL operator + (const POLYNOMIAL& P , const POLYNOMIAL& Q)
{
	int d;
	if ( P.degree > Q.degree)
		d =  P.degree;
	else
		d = Q.degree;

	COMPLEX c[d+1];
	for(int i=0; i< d+1; i++)
	{
		if(i < P.degree +1 && i < Q.degree+1)
			c[i] = P.coef[i] + Q.coef[i];
		else if (i > P.degree)
			c[i] = Q.coef[i];
		else if (i > Q.degree)
			c[i] = P.coef[i];
	}

	POLYNOMIAL R = POLYNOMIAL(d, c);
	return R;
}

POLYNOMIAL operator - (const POLYNOMIAL& P)
{

	COMPLEX c[P.degree+1];
	for(int i=0; i< P.degree+1; i++)
		c[i] =  - P.coef[i];

	POLYNOMIAL R = POLYNOMIAL(P.degree, c);
	return R;
}
POLYNOMIAL operator - (const POLYNOMIAL& P, const POLYNOMIAL& Q)
{
	int d;
	int p = P.degree;
	int q = Q.degree;

	if (p > q)
		d = p;
	else
		d = q;

	COMPLEX c[d+1];

	for(int i=0; i< d+1; i++)
	{
		if(i < p +1 && i < q+1)
			c[i] = P.coef[i] - Q.coef[i];
		else if (i > p)
			c[i] = COMPLEX(0,0) - Q.coef[i];
		else if (i > q)
			c[i] = P.coef[i];
	}

	POLYNOMIAL R = POLYNOMIAL(d, c);
	return R;
}


POLYNOMIAL operator * (const POLYNOMIAL& P, const POLYNOMIAL& Q)
{
	int d = P.degree + Q.degree;
	int p = P.degree;
	int q = Q.degree;
	COMPLEX c[d+1];

	for (int i=0; i < p + 1; i++)
	{
		for(int j=0; j < q +1; j++)
		{
			c[i+j] = c[i+j] +  P.coef[i] * Q.coef[j];
		}
	}
	POLYNOMIAL R = POLYNOMIAL(d, c);

	return R;
}


POLYNOMIAL operator * (const COMPLEX &c, const POLYNOMIAL &p) {
    POLYNOMIAL res(p.degree, p.coef);
    for(auto &now : res.coef) {
        now = now * c;
    }
    return res;
}


POLYNOMIAL operator * (const POLYNOMIAL &p, const COMPLEX &c) {
    POLYNOMIAL res(p.degree, p.coef);
    for(auto &now : res.coef) {
        now = now * c;
    }
    return res;
}


COMPLEX evaluate(const POLYNOMIAL& P, const COMPLEX& z)
{
	COMPLEX fx = COMPLEX(0,0);
	for(int i=0;i<P.degree + 1;i++)
		fx = fx + power(z,i) * P.coef[i];
	return fx;

}

POLYNOMIAL deriv(POLYNOMIAL P, int k)
{
	int d;
	if (k == 0)
		return POLYNOMIAL(P.degree, P.coef);

	if (k > P.degree)
	{
		d = 0;
		COMPLEX C[1];
		C[0] = COMPLEX(0,0);
		return POLYNOMIAL(d, C);

	}
	else
	{
		d = P.degree - k;
		INTEGER tmp;
		COMPLEX C[d+1];
		for(int i=0; i< d+1; i++)
		{
			tmp = 1;
			for (int j=0; j< k; j++)
			{
				tmp = tmp * (j+i+1);
			}
			C[i] = P.coef[k+i] * REAL(tmp);
		}

		return POLYNOMIAL(d, C);
	}
}

// kth taylor coefficient of P at z
COMPLEX CoefAt(POLYNOMIAL P, int k, const COMPLEX& z)
{
	COMPLEX tmp[P.degree+1];
	for(int i=0; i<P.degree + 1; i++)
		tmp[i] = P.coef[i];

	for(int j=0; j<k; j++)
	{
		for(int i=0;i<P.degree-j;i++)
		{
			tmp[i] = tmp[i+1] * REAL(i+1);
		}
	}

	COMPLEX fx = COMPLEX(0,0);
	for(int i=0;i < P.degree + 1 - k; i++)
	{
		fx = fx + power(z,i) * tmp[i];
	}
	return  fx / REAL(factorial(k));
}

// translation and dilation by g(x) = f(a + bx)
POLYNOMIAL translation(const POLYNOMIAL& P, const REAL& a, const COMPLEX& b)
{
	COMPLEX C[P.degree +1];
	for(int i=0; i<P.degree + 1; i++)
	{
		C[i] = power(a, i) * deriv(P,i)(b) / REAL(factorial(i));
	}

	return POLYNOMIAL(P.degree, C);
}


// power of the polynomial
POLYNOMIAL power(const POLYNOMIAL &p, int n) {
    assert(n >= 0);
    if(n == 0) {
        COMPLEX c(REAL(1), REAL(0));
        return POLYNOMIAL(0, &c);
    }
    if(n == 1) return p;
    auto r = power(p, n / 2);
    return r * r * power(p, n % 2);
}


// composition of the polynomial
POLYNOMIAL composite(const POLYNOMIAL &p, const POLYNOMIAL &q) {
    COMPLEX c(REAL(0));
    POLYNOMIAL res(0, &c);
    c = COMPLEX(REAL(1));
    POLYNOMIAL Q(0, &c);
    for(int i = 0; i <= p.degree; i++) {
        res = res + p.coef[i] * Q;
        if(i != p.degree) Q = Q * q;
    }
    return res;
}

/*

POLYNOMIAL operator + (const POLYNOMIAL& P , const REAL& a)
{
	REAL c[P.degree+1];
	for(int i=0; i< P.degree+1; i++)
	{
		c[i] = P.coef[i];
	}
	c[0] += a;
	return POLYNOMIAL(P.degree, c);
}
POLYNOMIAL operator + (const REAL& a , const POLYNOMIAL& P)
{
	return P + a;
}


POLYNOMIAL operator - (const POLYNOMIAL& P, const REAL& c)
{
	return P + REAL(-c);
}
POLYNOMIAL operator - (const REAL& c, const POLYNOMIAL& P)
{
	return (-P) + c;
}

POLYNOMIAL operator * (const POLYNOMIAL& P , const REAL& a)
{
	REAL c[P.degree+1];
	for(int i=0; i< P.degree+1; i++)
	{
		c[i] = P.coef[i] * a;
	}
	return POLYNOMIAL(P.degree, c);
}
POLYNOMIAL operator * (const REAL& c, const POLYNOMIAL& P)
{
	return P*c;
}
POLYNOMIAL operator / (const POLYNOMIAL& P, const REAL& c)
{
	return P*REAL(1/c);
}

*/

}
