#pragma once
#include "iRRAM_extension/random.hpp"
#include "iRRAM_extension/utility.hpp"
#include<functional>
#include<map>
using namespace iRRAM;

namespace iRRAM{

class WIENER
{
    private:
        std::map< std::pair<INTEGER,INTEGER>, REAL > X;
        REAL X_0;
        int curPrec,probPrec;
        int modulus(int ,int );
        REAL getX(int ,int );
    public:
        WIENER();
        WIENER(int);
        void setPrec(int,int);
        REAL compute(REAL);
        REAL operator() (REAL);
};

}