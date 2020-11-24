/*
the header file defines programs for approximating roots of a complex polynomial
*/
#pragma once

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM_extension/polynomial.hpp"
#include "iRRAM_extension/polynomial/rcomplex.hpp"

#include <utility>
using namespace iRRAM;
namespace iRRAM{

namespace internal{template <> struct is_continuous<std::vector< std::pair<REAL, int> > > : public std::true_type{};}


namespace internal{template <> struct is_continuous<std::vector<  COMPLEX > > : public std::true_type{};}

inline sizetype geterror( const std::vector<COMPLEX > & l){
	sizetype error, lerror;
	l[0].geterror(error);
	for (unsigned int i=0; i< l.size(); i++)
	{
		lerror = geterror(l[i]);
		sizetype_max(error,error,lerror);
	}
	return error;
}


inline void seterror(std::vector<COMPLEX > & l, sizetype & error)
{
	for (unsigned int i=0; i< l.size(); i++)
		l[i].seterror(error);
}

inline void adderror(std::vector<COMPLEX > & l, sizetype & error)
{
	for (unsigned int i=0; i< l.size(); i++)
		l[i].adderror(error);
}


inline sizetype geterror( const std::vector<std::pair<REAL, int> > & l){	sizetype error, lerror;
	l[0].first.geterror(error);
	for (unsigned int i=0; i< l.size(); i++)
	{
		lerror = geterror(l[i].first);
		sizetype_max(error,error,lerror);
	}
	return error;
}


inline void seterror(std::vector<std::pair<REAL, int> > & l, sizetype & error)
{
	for (unsigned int i=0; i< l.size(); i++)
		l[i].first.seterror(error);
}


class cvec_wrap{
public:
  std::vector<COMPLEX> data;
  cvec_wrap(std::vector<COMPLEX> v){ data = v; }
	cvec_wrap(){ }


	void adderror (sizetype error);
	void seterror (sizetype error);
	void geterror (sizetype& error) const;




};


//
// std::vector<R_COMPLEX >
// root_approximation_naive(int , POLYNOMIAL );
//
// std::vector<R_COMPLEX >
// root_approximation_newton(int , POLYNOMIAL );
}
