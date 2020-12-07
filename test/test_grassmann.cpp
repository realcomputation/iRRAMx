#include <iRRAM.h>
#include "iRRAMx/grassmann.hpp"
using namespace iRRAM;





void compute() {
    REALMATRIX A(2, 10);
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 0;

    A(1, 0) = 0;
    A(1, 1) = 1;
    A(1, 2) = 1;



    REALMATRIX B(2, 10);
    B(0, 0) = 1;
    B(0, 1) = 0;
    B(0, 2) = 0;

    B(1, 0) = 0;
    B(1, 1) = 1;
    B(1, 2) = 0;
    
    GRASSMANN R = GRASSMANN(A).meet(GRASSMANN(B), 1);
    R.show();
    

    
}
