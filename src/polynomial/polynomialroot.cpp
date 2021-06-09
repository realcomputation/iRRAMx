#include "iRRAMx/polynomial/polynomialroot.hpp"
#include "iRRAMx/polynomial/complexplane.hpp"
#include "iRRAMx/polynomial.hpp"
#include "iRRAMx/polynomial/utilities.hpp"
#include "iRRAMx/utility.hpp"

#include "iRRAM/lib.h"
#include "iRRAM/core.h"

#include <stack>
#include <queue>
#include <utility>
#include <vector>
#include <cmath>
#include <time.h>

#include <string>
#include <sstream>
#include <iterator>

using namespace iRRAM;
namespace iRRAM{

// < relation to use pre-impelented sorting function and prioirty queue.
/*
struct less_than_key
{
    inline bool operator() (const std::pair<R_OPENDISC, int>& C1,
    	const std::pair<R_OPENDISC, int>& C2)
	{
		return (C1.first.center < C2.first.center);
	}
};

struct IC_LESS
{
    inline bool operator() (const INTERVALCOMPONENT& C1,
    	const INTERVALCOMPONENT& C2)
	{
		return (C1.lower + C1.upper < C2.lower + C2.upper);
	}
};
*/

struct disc_wider
{
    inline bool operator() (const std::pair<R_OPENDISC, int>& C1,
    	const std::pair<R_OPENDISC, int>& C2)
	{
		return (C1.first.radius < C2.first.radius);
	}
};

struct box_wider
{
    inline bool operator() (const std::pair<R_CLOSEDBOX, int>& C1,
    	const std::pair<R_CLOSEDBOX, int>& C2)
	{
		return (C1.first.width < C2.first.width);
	}
};

struct component_wider
{
    inline bool operator() (const COMPONENT& C1,
    	const COMPONENT& C2)
	{
		return (C1.Wc() < C2.Wc());
	}
};


/*

std::vector<std::pair<RATIONALINTERVAL, int> >
casting(std::vector<INTERVALCOMPONENT> P)
{

	std::vector<std::pair<RATIONALINTERVAL, int> > Q;

	for (int i =0; i<(int) P.size(); i++)
	{
		Q.push_back(std::pair<RATIONALINTERVAL, int> (RATIONALINTERVAL(P[i].Mc(), P[i].Wc()/2), P[i].kc));
	}
	return Q;
}
*/


/******************
GRAEFFE'S ITERATION
*******************/

POLYNOMIAL polyE(const POLYNOMIAL& P)
{
	double td = P.degree / 2;
	int d = std::floor(td);

	COMPLEX C[d+1];

	for(int i = 0; i<d+1; i++)
	{
		C[i] = P.coef[i*2];
	}
	return POLYNOMIAL(td, C);
}

POLYNOMIAL polyO(const POLYNOMIAL& P)
{
	double td = (P.degree-1) / 2;
	int d = std::floor(td);

	COMPLEX C[d+1];

	for(int i = 0; i<d+1; i++)
	{
		C[i] = P.coef[1+i*2];
	}
	return POLYNOMIAL(td, C);
}

POLYNOMIAL Giterate(const POLYNOMIAL& P, int N)
{

	POLYNOMIAL Q = POLYNOMIAL(P.degree, P.coef);
	COMPLEX c[2];
	c[0] = COMPLEX(0,0); c[1] = COMPLEX(1,0);
	int k = pow(-1, P.degree);
	POLYNOMIAL T = POLYNOMIAL(1, c);
	if(k > 0)
		for(int i=1; i<N+1; i++)
		{
			Q = (polyE(Q) * polyE(Q)) - (T * (polyO(Q) * polyO(Q)));
		}
	else
		for(int i=0; i<N; i++)
		{
			Q = (T*(polyO(Q) * polyO(Q))) - (polyE(Q) * polyE(Q));
		}
	return Q;
}




/*
COMPUTABLE PREDICATES FOR ROOT CONTAINMENT OF AN OPEN DISC
*/

bool unit_disc_test(const POLYNOMIAL& P, int k)
{
	REAL LHS;
	REAL RHS = 0;
	for (int i=0; i<P.degree+1; i++)
		if(i!=k)
			RHS += abs(P.coef[i]);
	LHS = abs(P.coef[k]);
	int ans = choose(LHS>RHS, LHS<RHS, LHS*2 < RHS*3 && 2*RHS < 3*LHS) == 1;
	return ans;
}

// Soft pallet test: \tilde{\mathcal{T}}_k(P,D)
bool softTTest(const POLYNOMIAL& P, int k, const R_OPENDISC& D)
{
	COMPLEX m = COMPLEX(D.center.real(), D.center.imag());
	REAL r = D.radius;
	REAL LHS = abs(CoefAt(P,k,m)) * power(r,k);
	REAL RHS = 0;
	bool ans;
	for(int i=0;i<P.degree+1;i++)
	{
		if(i!=k)
		{
			RHS += abs(CoefAt(P,i,m)) * power(r,i);
		}
	}
	ans = choose(LHS>RHS, LHS<RHS, LHS*2 < RHS*3 && 2*RHS < 3*LHS) == 1;
	return ans;
}

// Soft pallet test on Graeffe iteration: \tilde{\mathcal{T}}_k(G_{\log(1+\log d) + 5}(P_D),D(0,1))
// on a interval with rational endpoints
bool softGTest(const POLYNOMIAL& P, int k, const R_OPENDISC& D)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, D.radius, COMPLEX(D.center.real(), D.center.imag())), N);
	// R_OPENDISC U(R_COMPLEX(0, 0), 1);
	// return softTTest(Q, k, U);
	return unit_disc_test(Q,k);
}


// Soft pallet test on Graeffe iteration: \tilde{\mathcal{T}}_k(G_{\log(1+\log d) + 5}(P_D),D(0,1))
// on a interval with real endpoints
bool softGTest(const POLYNOMIAL& P, int k, const OPENDISC& D)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, D.radius, D.center), N);
	// R_OPENDISC U(R_COMPLEX(0, 0), 1);
	// return softTTest(Q, k, U);
	return unit_disc_test(Q,k);

}


// First i such that T_i(G_N(f_{m,r}), 0, 1) /\ T_i(G_N(f_{m,3r}), 0, 1)holds
//    0 if T_0(G_N(f_{m,r}), 0, 1)
//   -1 if - T_i(G_N(f_{m,r}), 0, 1) for all i
// -i-1 if T_i(G_N(f_{m,r}), 0, 1) /\  -T_i(G_N(f_{m,3r}), 0, 1)
int softGStarThree(const POLYNOMIAL& P, const R_OPENDISC& B, int k)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, B.radius, COMPLEX(B.center.real(), B.center.imag())), N);
	POLYNOMIAL R = Giterate(translation(P, 3*B.radius, COMPLEX(B.center.real(), B.center.imag())), N);

	// R_OPENDISC U(R_COMPLEX(0, 0), 1);

	// for(int i =0;i<k+1;i++)
	// {
	// 	if(softTTest(Q, i, U))
	// 	{
	// 		if (i == 0)
	// 			return 0;
	// 		if(softTTest(R, i, U))
	// 			return i;
	// 		else
	// 			return -i-1;
	// 	}
	// }

	for(int i =0;i<k+1;i++)
	{
		if(unit_disc_test(Q, i))
		{
			if (i == 0)
				return 0;
			if(unit_disc_test(R, i))
				return i;
			else
				return -i-1;
		}
	}



	return - 1;
}

// First i such that \tilde{T}_i(G_N(f_{m,r}), 0, 1) holds
//  -1 if all fail
int softGStar(const POLYNOMIAL& P, const R_OPENDISC& B, int k)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, B.radius, COMPLEX(B.center.real(), B.center.imag())), N);
	R_OPENDISC U(R_COMPLEX(0, 0), 1);

	for(int i =0;i<k+1;i++)
		// if(softTTest(Q, i, U))
		if(unit_disc_test(Q, i))

			return i;
	return - 1;
}






// Return a list of (D_i, d_i) such that each D_i is disjoint (1,3)--isolating for positive number of roots.
// Implementation based on simple naive subdivision algorithm with Graeffe's iteration
std::vector<R_COMPLEX >
root_approximation_naive(int prec, const POLYNOMIAL& P)
{
	int n = P.degree;
	RATIONAL min = 1;


/*
COMPUTING ROOT BOUND D
*/

	R_OPENDISC D = R_OPENDISC();
	while(!softGTest(P, n, D))
	{
		D = D.multiply(2);
	}
	// cout <<D.radius<<"\n";

	REAL epsilon = power(2, prec);
    std::priority_queue<std::pair<R_CLOSEDBOX, int> , std::vector<std::pair<R_CLOSEDBOX, int> >, box_wider> Q0;
    Q0.push(std::pair<R_CLOSEDBOX, int>(R_CLOSEDBOX(R_COMPLEX(0,0), 2*D.radius), P.degree));

	int k;


	int nroot=0;
	int f = 0;
	int ff = 0;

	std::vector<std::pair<R_OPENDISC, int> > solutionDiscs;
	solutionDiscs.reserve(P.degree);

	cout << D.radius <<"\n";

/*
MAIN LOOP
*/
	while(!Q0.empty())
	{
		std::pair<R_CLOSEDBOX, int> fp = Q0.top();
		Q0.pop();

		R_CLOSEDBOX a = fp.first;
		R_OPENDISC b = R_OPENDISC(a);


		int max_multiplicity = fp.second;

		ff = 0;
		for (int i = 0; i < solutionDiscs.size(); i++)
		{
			if (intersect(b, solutionDiscs[i].first))
			{
				ff = 1; break;
			}
		}

		cout<< a.id<<"\n";
		if (ff == 0)
		{
			k = softGStarThree(P, b, max_multiplicity);

			if (k==-1)   // failed
			{
				for(int si=1;si<5;si++)
					Q0.push(std::pair<R_CLOSEDBOX, int>(a.subdivide(si), max_multiplicity));

			}
			else if (k < 0)	// T(3D) failed but T(D) succeeded with i = -k - 1
			{
				for(int si=1;si<5;si++)
					Q0.push(std::pair<R_CLOSEDBOX, int>(a.subdivide(si), -k - 1));

			}
			else if (k != 0) // D, 3D has k > 0 roots
			{
				if(choose(b.radius < epsilon, b.radius > epsilon / 2) == 2)
				{
					for(int si=1;si<5;si++)
						Q0.push(std::pair<R_CLOSEDBOX, int>(a.subdivide(si), k));
				}
				else
				{
					solutionDiscs.emplace_back(b, k);
					nroot+=k;
				}
			}
		}
	}

	std::vector<R_COMPLEX > ansvec;
	for (int i=0; i<(int) solutionDiscs.size(); i++)
	{
		for (int j=0; j<solutionDiscs[i].second; j++)
			ansvec.push_back(solutionDiscs[i].first.center);
	}
	return ansvec;
}


/*
CYAP2016 NEWTON METHOD
*/



COMPONENT
Newton(const POLYNOMIAL& P, const COMPONENT& C, const REAL& epsilon)
{

	REAL L, R;
	COMPLEX Xprime; // pick a point which is C.wc/2 away from the boundary of the commoment
	COMPONENT D;
	R_CLOSEDBOX ans;

	RATIONAL new_left, new_below, new_width;

	R_COMPLEX R_Xc = C.right_most_box.center + R_COMPLEX(C.wc()/2, 0);
	COMPLEX Xc = COMPLEX(R_Xc.real(), R_Xc.imag());
	RATIONAL r = C.Wc() / 2;

	L = 4 * r * abs(deriv(P, 1)(Xc));
	R = abs(P(Xc));

	if( choose(L>R, L<R, 2*L< 3*R && 2*R < 3*L) != 2) // when f'(Xc) != 0
	{

		Xprime = Xc - REAL(C.kc) * P(Xc)/ deriv(P,1)(Xc);
		if (is_in(C, Xprime))
		{

			OPENDISC II = OPENDISC(Xprime, REAL( C.wc()) / 8 / C.Nc);

			if( softGTest(P, C.kc, II))
			{

				new_width = C.wc() / (2*C.Nc);

				INTEGER index_x = ( (real(Xprime) - C.left_most ) / new_width ).as_INTEGER();
				INTEGER index_y = ( (imag(Xprime) - C.below_most) / new_width ).as_INTEGER();

				if (index_x == 0)
					index_x = index_x + 1;
				if (index_x == (C.right_most - C.left_most / C.wc()) )
					index_x = index_x - 1;

				if (index_y == 0)
					index_y = index_y + 1;
				if (index_y == (C.upper_most - C.below_most / C.wc()) )
					index_y = index_y - 1;

				new_left  = C.left_most + new_width * (index_x - 1);
				new_below = C.below_most + new_width * (index_y - 1);
				ans = R_CLOSEDBOX(R_COMPLEX(new_left+new_width, new_below+new_width),2*new_width);

				D = COMPONENT();
				D.add(ans);
				D = D.split();
				D.Nc = C.Nc * C.Nc;
				D.depth = C.depth+1;
			}

		}

	}
	return D;
}


std::vector< COMPONENT >
Bisect(const POLYNOMIAL& P, const COMPONENT& Comp)
{
  COMPONENT C;
  std::vector< COMPONENT > U, J;
  R_CLOSEDBOX B;

  int flg;
  bool specialflg = true;

  C = Comp.split();  //C.depth ++

  for (int i =0; i< C.size(); i++)
    {
      B = C[i];
      if (!softGTest(P, 0, R_OPENDISC(B)))
	{
	  COMPONENT T = COMPONENT(B);
	  T.depth = C.depth;
	  print(T);

	  U.push_back(T);
	}
      // else
      // 	{
      // 	  COMPONENT T = COMPONENT(B);
      // 	  T.depth = C.depth;
      // 	  cout <<"no root in \n";
      // 	  print(T);	  
      // 	}
      
    }


  std::vector<int> rm_idx;
  for (int i = 0; i<U.size(); i++)
    for (int j = i+1; j<U.size(); j++)
      {
	if(adj(U[i], U[j]))
	  {
	    for (int k = 0; k < U[i].size(); k++)
	      U[j].add(U[i][k]);
	    rm_idx.push_back(i);
	    break;
	  }
      }



  int tind = 0;
  if (!rm_idx.empty()){
    for (int i = 0; i<U.size(); i++)
      {
      if(tind < rm_idx.size())
	{
        if(i == rm_idx[tind])
	  {
          tind+=1;
        }
        else{
          J.push_back(U[i]);
        }
      }
      else{
        J.push_back(U[i]);
      }
    }
  }
  else
    {
      for (int i = 0; i<U.size(); i++)
	{
	  J.push_back(U[i]);
	}
      }

  if ((int)J.size() ==1)
    specialflg = false;

  for (int i=0; i<(int)J.size(); i++)
    {
      if (specialflg)
	J[i].Nc = 4;
      else
	if (sqrt(C.Nc) > 4)
	  J[i].Nc = sqrt(C.Nc);
	else
	  J[i].Nc = 4;
    }

  return J;
}

template <typename Out>
void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

REAL dist(const R_COMPLEX& a, const R_COMPLEX& b){
  REAL w = REAL(a.real()) - REAL(b.real());
  REAL h = REAL(a.imag()) - REAL(b.imag());
  return sqrt(w * w + h * h);
}

// std::vector<COMPLEX >
cvec_wrap
root_approximation_newton(int p, std::string& choice, const POLYNOMIAL& Q){

  POLYNOMIAL P = Q;
  int n = P.degree;

  /*
    COMPUTING ROOT BOUND
  */
  R_OPENDISC D = R_OPENDISC(R_COMPLEX(0,0),1);
  while(!softGTest(P, P.degree, D))
    {
      D = D.multiply(2);
    }


  std::vector< COMPONENT > Q_out, Q_main;
  std::vector< COMPONENT > bisect_result;

  COMPONENT C = COMPONENT(R_CLOSEDBOX(D.center, D.radius*2));
  C.Nc = 4;
  Q_main.push_back(C);
  COMPONENT Cprime;

  REAL epsilon = prec(p);

  int pi;

  RATIONAL tmp, qtmp;
  COMPONENT tc;


  int fll;
  while(!Q_main.empty())
    {

      // Swap Q1[end -1] <-> Q1[max wc]
      pi = 0;
      qtmp = 0;
      tmp = 0;
      for(int i=0; i< (int) Q_main.size(); i++)
	{
	  qtmp = Q_main[i].Wc();
	  if (qtmp > tmp)
	    {
	      pi = i;
	      tmp = qtmp;
	    }
	}

      tc = Q_main[pi];
      Q_main[pi] = Q_main[(int) Q_main.size() - 1];
      Q_main[(int) Q_main.size() -1] = tc;
      // swap done


      C = Q_main.back();
      Q_main.pop_back();
      fll = 0;
      R_OPENDISC II = R_OPENDISC(C.Mc(),C.Rc());

      for(int i=0; i< (int) Q_main.size(); i++)
	{
	  if(intersect(R_OPENDISC(Q_main[i].Mc(), Q_main[i].Rc()), II.multiply(4)))
	    {
	      fll = 1;
	      break;
	    }
	}

      if (fll == 0)
	{
	  C.kc = softGStar(P, II, P.degree);

	  if(C.kc >0)
	    {
	      if(choose(C.Wc() > epsilon / 2, C.Wc() < epsilon ) == 1)
		{
		  Cprime= Newton(P,C,epsilon);

		  if(!Cprime.is_empty())
		    {
		      // std::cout <<" ==== \n";
		      // printr(C);
		      // std::cout <<" newton worked\n";
		      Q_main.push_back(Cprime);
		      continue;
		    }
		}
	      else if (C.Wc() <= 3 * C.wc())
		{
		  // std::cout <<" ==== \n";
		  // printr(C);
		  // std::cout <<" adding to q out\n";
				    
		  Q_out.push_back(C);
		  continue;
		}
	    }
	}
      // std::cout<<"F\n";

      // std::cout <<"bisection start!\n";
      // print(C);
      // std::cout <<"reduces to!\n";
      bisect_result = Bisect(P,C);
      for(int i=0; i<(int)bisect_result.size(); i++){
	Q_main.push_back(bisect_result[i]);
	// print(bisect_result[i]);
	// std::cout<<"---\n";
      }
      // std::cout<<"G\n";
		
		
    }

  // std::cout<<"\n\n doneeeeeeeee\n";
  std::vector< R_COMPLEX > roots;
  for (int i =0; i<(int) Q_out.size(); i++)
    for (int j=0; j< Q_out[i].kc; j++)
      {
	// print(Q_out[i]);
	roots.push_back(Q_out[i].Mc());
	// std::cout <<"haha!\n";
      }


	
  // CHECK M-cache

  std::vector< R_COMPLEX> ordered_roots;
  if (choice.empty()){
    std::string r = std::to_string(p);
    r += ",";
    for (int i = 0; i < roots.size(); i++){
      r += (swrite(numerator(roots[i].real())) +","+swrite(denominator(roots[i].real()))+",");
      r += (swrite(numerator(roots[i].imag())) +","+swrite(denominator(roots[i].imag())));
      if(i != roots.size() - 1) r += ",";
    }
    choice = r;

    ordered_roots = roots;
  }
  else{

    int cached_prec;
    std::vector< R_COMPLEX> cached_roots;
    std::vector< std::string> cached_vec = split(choice, ',');
    cached_prec = std::stoi(cached_vec[0]);
    cached_roots.reserve((cached_vec.size() - 1) / 4);
    for(int i = 0; i < (cached_vec.size() - 1) / 4; i++){
      cached_roots.emplace_back(RATIONAL(INTEGER(cached_vec[i*4+1]), INTEGER(cached_vec[i*4+2])), RATIONAL(INTEGER(cached_vec[i*4+3]), INTEGER(cached_vec[i*4+4])));
    }

    for(int i = 0; i<cached_roots.size(); i ++){
      for(int j = 0; j<roots.size(); j ++){
        if(dist(cached_roots[i], roots[j]) < prec(cached_prec) + prec(p))
	  {
	    ordered_roots.push_back(roots[j]);
	    roots.erase(roots.begin() + j);
	    break;
	  }
      }
    }
    std::string r = std::to_string(p);
    r += ",";
    for (int i = 0; i < ordered_roots.size(); i++){
      r += (swrite(numerator(ordered_roots[i].real())) +","+swrite(denominator(ordered_roots[i].real()))+",");
      r += (swrite(numerator(ordered_roots[i].imag())) +","+swrite(denominator(ordered_roots[i].imag())));
        if(i != ordered_roots.size() - 1) r += ",";
    }
    choice = r;
  }

  std::vector<COMPLEX> croots;
  croots.reserve(ordered_roots.size());
    for (int i = 0; i < ordered_roots.size(); i++)
    {
      croots.emplace_back(REAL(ordered_roots[i].real()), REAL(ordered_roots[i].imag()));
    }


  // return croots;
  return cvec_wrap(croots);

}


void cvec_wrap::adderror (sizetype error)
{
	for (unsigned int i=0; i < this->data.size(); i++)
  	this->data[i].adderror(error);
}

void cvec_wrap::seterror (sizetype error)
{
	for (unsigned int i=0; i < this->data.size(); i++)
  	this->data[i].seterror(error);
}

void cvec_wrap::geterror (sizetype& error) const
{

  sizetype lerror;
	this->data[0].geterror(error);
	for (unsigned int i=0; i< this->data.size(); i++)
	{
		this->data[i].geterror(lerror);
		error = max(error,lerror);
	}
}


namespace internal{template <> struct is_continuous<cvec_wrap > : public std::true_type{};}




std::vector<COMPLEX> roots(const POLYNOMIAL& P){
  std::string r;

  // return limit_mv<std::vector<COMPLEX> >(root_approximation_newton, r, P);
  cvec_wrap t = limit_mv_<cvec_wrap >(root_approximation_newton, r, P);
  return t.data;

}

// template <typename Out>
// void split(const std::string &s, char delim, Out result) {
//     std::istringstream iss(s);
//     std::string item;
//     while (std::getline(iss, item, delim)) {
//         *result++ = item;
//     }
// }




}
