#include "iRRAM_extension/wiener.hpp"

using namespace iRRAM;


REAL chauder(REAL x, int n, INTEGER k)
{
    if(n==0) return x;
    // INTEGER k=2*j-1;
    REAL max=  iRRAM::power(REAL(2),-1*REAL(n-1)/REAL(2));

    REAL first = REAL( (k-1)/prec(n));
    REAL second = REAL(k/prec(n));
    REAL third = REAL((k+1)/prec(n));
    // auto chauderL=[=](REAL t)  ->REAL{return  max*(t-first)*prec(n);};
    // auto chauderR=[=](REAL t)  ->REAL{return  max*(third-t)*prec(n);};
    return glue( x<first , 0, glue(x<second , max*(x-first)*prec(n) , glue( x<third, max*(third-x)*prec(n) , 0  ) ) );

    /*
    Chauder Function
           /\
          /  \
         /    \
        /      \
   ____________________
        1   2   3
  */
} 

namespace iRRAM{

    REAL WIENER::getX(int n,int k)
    {
        if(X.count(std::make_pair(n,k)) ==0) 
            X[std::make_pair(n,k)]=gaussian_real();
        return X[std::make_pair(n,k)];
    }

    WIENER::WIENER()
    {
        X_0=gaussian_real();
        curPrec= actual_stack().actual_prec;
        probPrec=50;
    }

    WIENER::WIENER(int pBound)
    {
        X_0=gaussian_real();
        curPrec= actual_stack().actual_prec;
        probPrec=pBound;
    }

    int WIENER::modulus(int prec1,int prec2)
    {
        // to be calculated.
        return 50;
    }

    void WIENER::setPrec(int p1,int p2)
    {
        curPrec=p1;
        probPrec=p2;
    }

    REAL WIENER::compute(REAL t)
    {
        REAL Y;
        int N=modulus(curPrec,probPrec);
        REAL result=X_0*t;
        //W_N construction. 
        for(int n=1;n<=N;n++)
        {
            // Only need to consider (k-1)/2^n <= t <= (k+1)/2^n
            // In other word, (2^n)*t - 1 <= k <= (2^n)*t + 1 
            INTEGER k= (prec(n)*t).as_INTEGER();
            result += getX(n,int(k-1))*chauder(t,n,k-1)  ;
            result += getX(n,int(k))*chauder(t,n,k)  ;
            result += getX(n,int(k+1))*chauder(t,n,k+1);
        }
        return result;
    }

    REAL WIENER::operator()(REAL t)
    {
        return this->compute(t);
    }
}