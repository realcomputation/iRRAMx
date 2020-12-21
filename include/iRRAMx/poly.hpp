#pragma once
#include "fix.hpp"

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <array>
#include <utility>
#include "utility.hpp"
using namespace iRRAM;
namespace iRRAM{

/** @defgroup polynomial Polynomial
 *  @brief include "iRRAM_extension/poly.hpp"
 *  Class template for multivariate polynomials
 *  @{
 */


//!  A class template for multivariate polynomials 
/*!
  A class template for multivariate polynomials over some arbitrary type.
	Ring operations are overloaded.
*/
template<class T, std::size_t D>
class POLYNOMIAL{
    private:
        //! coefficients
		/*!
			coefficients are polynomials with one variable less 
		*/
        using coeff_type = POLYNOMIAL<T, D-1>;
        std::vector<coeff_type> coeffs;
    public: 
        //! construction of the constant zero polynomial
        POLYNOMIAL() : POLYNOMIAL(T(0)) {};
        //! construction of a constant polynomial
		/*!
			@param c the value of the constant polynomial
		*/
        POLYNOMIAL(const T& c); 
        //! evaluation of the polynomial
		/*!
			@param x point where the polynomial is evaluated
            @return p(x)
		*/
        T operator()(const std::array<T,D>& x) const;
        //! get coefficient size
		/*!
            @return size of coefficient vector
		*/
        size_t coeff_size() const {return coeffs.size(); };
        //! return the coefficient at a specific index
		/*!
			@param index the index of the coefficient
            @return coefficient of p at given index
		*/
        T get_coefficient(const std::array<int,D>& index) const;
        //! return the sub-polynomial at index i 
        // if i is larger than the degree, an empty polynomial is returned
		/*!
			@param i the index of the coefficient
            @return D-1 dimensional polynomial corresponding to the i-th power of the first variable
		*/
        POLYNOMIAL<T,D-1> get_coefficient(const int i) const;
        //! set coefficient at index to some D-1 variate polynomial
		/*!
			@param i the index of the coefficient
            @param p the new value of the coefficient
		*/
        void set_coefficient(const int i, const POLYNOMIAL<T,D-1>& p);
        //! set coefficient at specific index to some value
		/*!
			@param index the index of the coefficient
            @param value the new value of the coefficient
		*/
        void set_coefficient(const std::array<int, D>& index, const T& value);

        std::vector<coeff_type> get_coeffs() const {return coeffs;}

};
//!  Template specialization for zero-dimensional polynomial
/*!
    A polynomial of dimension 0 is a constant 
*/
template<class T>
class POLYNOMIAL<T,1>{
    private:
        std::vector<T> coeffs;
    public:
        POLYNOMIAL(const T& c) : coeffs({c}) {};
        POLYNOMIAL() : POLYNOMIAL(T(0)) {};
        T get_coefficient(const int i) const {
            if(i < coeffs.size())
                return coeffs[i];
            return T(0);
        }
        T get_coefficient(const std::array<int,1>& index) const {return get_coefficient(index[0]);}
        void set_coefficient(const int i, const T& c) {
            if(i >= coeffs.size()) coeffs.resize(i+1);
            coeffs[i] = c;
        }
        void set_coefficient(const std::array<int, 1>& index, const T& c) {set_coefficient(index[0],c);}
        T operator()(const std::array<T,1>& x) const {
                T ans(0);
                for(int i=coeffs.size()-1; i>=0; i--){
                    ans = ans*x[0]+coeffs[i];
                }
                return ans;
        }
        size_t coeff_size() const {return coeffs.size(); };
};


/*! @brief difference of two polynomials
 *
 *
 *  @param lhs the first polynomial
 *  @param rhs the second polynomial
 *
 *  @return p1-p2
 */
template<class T, std::size_t D>
POLYNOMIAL<T,D> operator-(const POLYNOMIAL<T,D>& lhs, const POLYNOMIAL<T,D>& rhs);

/*! @brief product of two polynomials
 *
 *
 *  @param lhs the first polynomial
 *  @param rhs the second polynomial
 *
 *  @return lhs*rhs
 */
template<class T, std::size_t D>
POLYNOMIAL<T,D> operator*(const POLYNOMIAL<T,D>& lhs, const POLYNOMIAL<T,D>& rhs);

/*! @brief scalar multiplication
 *
 *
 *  @param lhs the scalar
 *  @param rhs the  polynomial
 *
 *  @return lhs*rhs
 */
template<class T, std::size_t D>
POLYNOMIAL<T,D> operator*(const T& lhs, const POLYNOMIAL<T,D>& rhs);

template<class T, std::size_t D>
POLYNOMIAL<T,D> operator*(const POLYNOMIAL<T,D>& lhs, const T& rhs){
    return rhs*lhs;
};

/*! @brief sum of two polynomials
 *
 *
 *  @param lhs the first polynomial
 *  @param rhs the second polynomial
 *
 *  @return p1+p2
 */
template<class T, std::size_t D>
POLYNOMIAL<T,D> operator+(const POLYNOMIAL<T,D>& lhs, const POLYNOMIAL<T,D>& rhs);

/*! @brief partial derivative of a polynomial
 *
 *
 *  @param p the polynomial
 *  @param alpha the order of the partial derivative
 *
 *  @return p^{(alpha_1, ..., alpha_d)}
 */
template<class T, std::size_t D>
POLYNOMIAL<T,D> derive(const POLYNOMIAL<T,D>&, const std::array<int, D>& alpha);

/*! @brief generate a polynomial from a string 
 *
 *
 *  @param alpha the order of the partial derivative
 *
 *  @return p^{(alpha_1, ..., alpha_d)}
 */

template<class T, std::size_t D>
POLYNOMIAL<T,D> from_string(const std::string& s, const std::array<char,D>& variable_symbols);


/*! @brief output a string for a polynomial with given coefficient precision
 *
 *
 *  @param p the polynomial
 *  @param variable_symbols the symbols used for the variables in order
 *  @param prec output precision for the coefficients
 *
 *  @return the polynomial as string
 */

template<class T, std::size_t D>
std::string to_string(const POLYNOMIAL<T,D>& p, const std::array<char,D>& variable_symbols, unsigned int prec);

// Template implementations
template<class T, std::size_t D>
POLYNOMIAL<T,D>::POLYNOMIAL(const T& c){
    coeffs = {POLYNOMIAL<T,D-1>(c)};
}

template<class T, std::size_t D>
POLYNOMIAL<T,D-1> POLYNOMIAL<T,D>::get_coefficient(const int i) const{
    if(i >= coeffs.size())
        return POLYNOMIAL<T,D-1>();
    return coeffs[i];
}
template<class T, std::size_t D>
T POLYNOMIAL<T,D>::get_coefficient(const std::array<int,D>& index) const{
    std::array<int, D-1> next_index;
    std::copy(index.begin()+1, index.end(), next_index.begin());
    return coeffs[index[0]].get_coefficient(next_index);
}

template<class T, std::size_t D>
void POLYNOMIAL<T,D>::set_coefficient(const int i, const POLYNOMIAL<T,D-1>& p) {
    if(i >= coeffs.size()){
        coeffs.resize(i+1);
    }
    coeffs[i] = p;    
}

template<class T, std::size_t D>
void POLYNOMIAL<T,D>::set_coefficient(const std::array<int, D>& index, const T& value) {
    if(index[0] >= coeffs.size()){
        coeffs.resize(index[0]+1);
    }
    std::array<int, D-1> next_index;
    std::copy(index.begin()+1, index.end(), next_index.begin());
    coeffs[index[0]].set_coefficient(next_index, value);
}

// evaluation using Horner's method
template<class T, std::size_t D>
T POLYNOMIAL<T,D>::operator()(const std::array<T,D>& x) const{
    T ans(0);
    std::array<T, D-1> next_x;
    std::copy(x.begin()+1, x.end(), next_x.begin());
    for(int i=coeffs.size()-1; i>=0; i--){
        ans = ans*x[0]+coeffs[i](next_x);
    }
    return ans;
}

template<class T, std::size_t D>
POLYNOMIAL<T,D> operator+(const POLYNOMIAL<T,D>& lhs, const POLYNOMIAL<T,D>& rhs){
    POLYNOMIAL<T,D> ans;
    int d=max(lhs.coeff_size(), rhs.coeff_size());
    for(int i=0; i<d; i++){
        ans.set_coefficient(i, lhs.get_coefficient(i)+rhs.get_coefficient(i));
    }
    return ans;
}

template<class T, std::size_t D>
POLYNOMIAL<T,D> operator-(const POLYNOMIAL<T,D>& lhs, const POLYNOMIAL<T,D>& rhs){
    POLYNOMIAL<T,D> ans;
    int d=max(lhs.coeff_size(), rhs.coeff_size());
    for(int i=0; i<d; i++){
        ans.set_coefficient(i, lhs.get_coefficient(i)-rhs.get_coefficient(i));
    }
    return ans;
}


template<class T, std::size_t D>
POLYNOMIAL<T,D> operator*(const POLYNOMIAL<T,D>& lhs, const POLYNOMIAL<T,D>& rhs){
    POLYNOMIAL<T,D> ans;
    for(int i=0; i<lhs.coeff_size(); i++){
        for(int j=0; j<rhs.coeff_size(); j++){
            auto p = ans.get_coefficient(i+j)+lhs.get_coefficient(i)*rhs.get_coefficient(j);
            ans.set_coefficient(i+j, p);
        }
    }
    return ans;
}
template<class T, std::size_t D>
POLYNOMIAL<T,D> operator*(const T& lhs, const POLYNOMIAL<T,D>& rhs){
    POLYNOMIAL<T,D> ans;
    for(int i=0; i<rhs.coeff_size(); i++){
        ans.set_coefficient(i, lhs*rhs.get_coefficient(i));
    }
    return ans;
}

template<class T>
T derive(const T& p, std::array<int,0>& alpha){
    return p;
}

template<class T, std::size_t D>
POLYNOMIAL<T,D> derive(const POLYNOMIAL<T,D>& p, const std::array<int,D>& alpha){
    POLYNOMIAL<T,D> ans;
    std::array<int, D-1> next_index;
    std::copy(alpha.begin()+1, alpha.end(), next_index.begin());
    int k = alpha[0];
    REAL fact = 1;
    for(int i=2; i<=k; i++)
        fact *= i;
    for(int i=k; i<p.coeff_size(); i++){
        if(i > k)
            fact = i*fact/(i-k);
        auto r = fact*derive(p.get_coefficient(i),next_index);
        ans.set_coefficient(i-k, r);
    }
    return ans;
}

template<class T>
std::vector<std::string> to_string_rec(const T& c, const std::array<char, 0>& vs, unsigned int prec){
    return {to_string(c, prec)};
}

template<class T, std::size_t D>
std::vector<std::string> to_string_rec(const POLYNOMIAL<T,D>& p, const std::array<char,D>& variable_symbols, unsigned int prec){
    std::array<char, D-1> next_symbols;
    std::copy(variable_symbols.begin()+1, variable_symbols.end(), next_symbols.begin());
    std::vector<std::string> ans;
    for(int i=0; i<p.coeff_size(); i++){
        auto v = to_string_rec(p.get_coefficient(i), next_symbols, prec);
         std::string xs ="";
         if(i > 0) xs += variable_symbols[0];
         if(i > 1) xs += "^"+std::to_string(i);
        for(auto& s : v){
            ans.push_back(s+xs);
        }
    }
    return ans;
}

template<class T, std::size_t D>
std::string to_string(const POLYNOMIAL<T,D>& p, const std::array<char,D>& variable_symbols, unsigned int prec){
    std::string ans="";
    auto v = to_string_rec(p, variable_symbols, prec);
    for(int i=0; i<v.size(); i++){
        if(i > 0) ans += " + ";
        ans += v[i];
    } 
    return ans;
}

}