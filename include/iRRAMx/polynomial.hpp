#pragma once
#include "fix.hpp"

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <utility>
using namespace iRRAM;
namespace iRRAM{

/** @defgroup poly Polynomial Roots
 *  @brief include "iRRAM_extension/polynomial.h"
 *  This is the collection of function specifications related to
 *  finding roots of polynomials
 *  @{
 */


//!  A class for complex polynomials.
/*!
  A class of polynomials with complex coefficients.
	Ring operations are overloaded.
	Need to be carefull about the degree of the result.
*/
class POLYNOMIAL
{
	public:
		//! coefficients
		/*!
			coefficients are stored as a std::vector<COMPLEX> coef.
			coef[0] is the constant. coef[n] is the coefficient of x^n
		*/
		std::vector<COMPLEX> coef;
		int degree;

		//! construction with an array of coefficients
		/*!
			@warning c[d] has to be nonzero; i.e., d really has to be the degree of c
			@param d a degree of the created polynomial
			@param c an array of coefficient where c[i] is the coef. of z^i
		*/
		POLYNOMIAL(int d, COMPLEX * c);

		//! construction with a vector of coefficients
		/*!
			@warning c[d] has to be nonzero; i.e., d really has to be the degree of c
			@param d a degree of the created polynomial
			@param c an array of coefficient where c[i] is the coef. of z^i
		*/
		POLYNOMIAL(int d, const std::vector<COMPLEX> &c);

		//! construction of a constant polynomial
		/*!
			@param c the value of the constant polynomial
		*/
		POLYNOMIAL(COMPLEX c);
		POLYNOMIAL();
		~POLYNOMIAL();
		COMPLEX operator () (const COMPLEX&);

		friend orstream &operator<<(orstream &ors, const POLYNOMIAL &p);
};

POLYNOMIAL operator + (const POLYNOMIAL&, const POLYNOMIAL&);
POLYNOMIAL operator - (const POLYNOMIAL&, const POLYNOMIAL&);
POLYNOMIAL operator * (const POLYNOMIAL&, const POLYNOMIAL&);
POLYNOMIAL operator * (const COMPLEX&, const POLYNOMIAL&);
POLYNOMIAL operator * (const POLYNOMIAL&, const COMPLEX&);

namespace internal{template <> struct is_continuous<POLYNOMIAL> : public std::true_type{};}
// void geterror(const POLYNOMIAL &);
int geterror_exp(const POLYNOMIAL & );

// Q = deriv(P,k) := Q^{(k)}

/*! @brief differentiate a polynomial
 *
 *
 *  @param k the number of times the differentiate op. is applied
 *
 *  @return p^{(k)}
 */
POLYNOMIAL deriv(POLYNOMIAL p, int k);
// y = evaluate(P, x) := P(x)

/*! @brief evaluate a polynomial at a point
 *
 *
 *
 *  @return p(x)
 */
COMPLEX evaluate(POLYNOMIAL p, COMPLEX x);
// Taylor coefficients: a_k = CoefAt(P, k, z) := P^{(k)}(z)/k!

/*! @brief compute Taylor coefficient
 *
 *
 *
 *  @return p^{(k)}(z)/k!
 */
COMPLEX CoefAt(POLYNOMIAL p, int k, const COMPLEX& z);
// Translation and dilation by g(x) = f(ax + b)
// a : REAL, b : COMPLEX
POLYNOMIAL translation(const POLYNOMIAL& , const REAL& , const COMPLEX& );


/*! @brief compute a power of the polynomial
 *
 *
 *
 *  @return p^n
 */
POLYNOMIAL power(const POLYNOMIAL &p, int n);


/*! @brief compute a composition of two polynomials
 *
 *
 *
 * @return p(q(x))
 */
POLYNOMIAL composite(const POLYNOMIAL &p, const POLYNOMIAL &q);

// Find root
/*! @brief find every roots of a polynomial
 *
 *
 *
 *  @return a vector of all roots of p
 */
std::vector<COMPLEX > roots(const POLYNOMIAL& p);

/** @} */

}
