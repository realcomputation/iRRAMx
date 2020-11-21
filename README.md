# iRRAM_extensionn
iRRAM_extension is a static C++ library that provides functionalities for rigorous scientific computation. 
It is designed in the way that it extends iRRAM. 
iRRAM is another C++ library that provides exact real-number computation. 
iRRAM provides data type REAL for real numbers and primitive operators for the field arithmetic of real numbers. 
They are computed exactly in iRRAM such that they do not introduce any rounding errors! 
Going beyond and based on those primitive operations, this extension provides more advanced data types and operations for continuous mathematics.

## Prerequisite
Of course, being built on iRRAM, iRRAM is a prerequisite of this library. 
The users of this library first have to install iRRAM in their system. 
Our library is, at this moment, compatible with the version of iRRAM that is accessible through their [Github](https://github.com/norbert-mueller/iRRAM) in the most recent commit (_`a4d2409` on May 1, 2015_). When iRRAM is being installed, it asks whether it will be installed directly in the system or in a local directory. For now, iRRAM_extension only detects iRRAM when it is installed in system.

## Installation
iRRAM_extension is a static C++ library. Its compilation is deferred to the users.
When a compatible iRRAM is installed in the system, run `make` on the directory of the source directory creates a static library file in `lib/libiRRAM_extension.a`:

```
% mkdir ~/git_dir
% cd ~/git_dir
% git clone https://github.com/realcomputation/iRRAM_extension.git
% cd git_dir
% make
```

When building the library is complete, the library can be either included in the system. Or, the library and the header files can be linked each time when a C++ code using iRRAM_extension library is used:

```
% g++ -std=c++14 -O2 -Wall -I~/git_dir/iRRAM_extension/include -L~/git_dir/iRRAM_extension/lib -liRRAM_extension -liRRAM -lmpfr -gmp my_cpp.cpp -o my_out
```

## Using iRRAM_extension
iRRAM_extension consists of orthogonal components:
- compact: 
- linear: linear algebra including fast matrix multiplication, matrix eigenproblem, ...
- polynomial:
- random:
- plot:

When a component is going to be used, the corresponding header file should be used.
For example, when the user want to have infinite precision random number,
```
#include <iRRAM.h>
#include "iRRAM_extension/random.hpp"

void compute(){
  REAL x = Gaussian(); // A real number sampled from the normal distribution
  ...
}
```
Here, `compute` is a function that replaces `main` in iRRAM. The readers who are not familiar with iRRAM are referred to a [short tutorial](..).

## Reference Manual

### iRRAM_extension/compact.hpp

### iRRAM_extension/linear.hpp

### iRRAM_extension/polynomial.hpp

### iRRAM_extension/random.hpp

### iRRAM_extension/plot.hpp

## Contributing to iRRAM_extension
