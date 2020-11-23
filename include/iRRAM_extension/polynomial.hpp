#pragma once

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM_extension/polynomial/rcomplex.hpp"
#include "iRRAM_extension/polynomial/polynomialroot.hpp"
#include <utility>
using namespace iRRAM;
namespace iRRAM{

class POLYNOMIAL
{
	public:
		std::vector<COMPLEX> coef;
		int degree;

		POLYNOMIAL(int, COMPLEX *);
		POLYNOMIAL(int, std::vector<COMPLEX>);
		POLYNOMIAL(COMPLEX);
		POLYNOMIAL();
		~POLYNOMIAL();
		COMPLEX operator () (const COMPLEX&);
};

POLYNOMIAL operator + (const POLYNOMIAL&, const POLYNOMIAL&);
POLYNOMIAL operator - (const POLYNOMIAL&, const POLYNOMIAL&);
POLYNOMIAL operator * (const POLYNOMIAL&, const POLYNOMIAL&);

namespace internal{template <> struct is_continuous<POLYNOMIAL> : public std::true_type{};}
// void geterror(const POLYNOMIAL &);
int geterror_exp(const POLYNOMIAL & );

// Q = deriv(P,k) := Q^{(k)}
POLYNOMIAL deriv(POLYNOMIAL, int);
// y = evaluate(P, x) := P(x)
COMPLEX evaluate(POLYNOMIAL, COMPLEX);
// Taylor coefficients: a_k = CoefAt(P, k, z) := P^{(k)}(z)/k!
COMPLEX CoefAt(POLYNOMIAL, int, COMPLEX);
// Translation and dilation by g(x) = f(ax + b)
// a : REAL, b : COMPLEX
POLYNOMIAL translation(POLYNOMIAL , REAL , COMPLEX );



class cvec_wrap{
public:
  std::vector<COMPLEX> data;
  cvec_wrap(std::vector<COMPLEX> v){ data = v; }
	cvec_wrap(){ }


	void adderror (sizetype error);
	void seterror (sizetype error);
	void geterror (sizetype& error) const;




};

cvec_wrap root_approximation_newton(int p, std::string& choice, const POLYNOMIAL& Q);

// Find root

std::vector<COMPLEX > roots(POLYNOMIAL P);


}
