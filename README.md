# iRRAMxn
iRRAMx is a static C++ library that provides functionalities for rigorous scientific computation.
It is designed in the way that it extends iRRAM.
(See )
iRRAM is another C++ library that provides exact real-number computation.
iRRAM provides data type REAL for real numbers and primitive operators for the field arithmetic of real numbers.
They are computed exactly in iRRAM such that they do not introduce any rounding errors!
Going beyond and based on those primitive operations, this extension provides more advanced data types and operations for continuous mathematics.

## Prerequisite
Of course, being built on iRRAM, iRRAM is a prerequisite of this library.
The users of this library first have to install iRRAM in their system.
Our library is, at this moment, compatible with the version of iRRAM that is accessible through a forked [Github repository](https://github.com/fbrausse/iRRAM) in the most recent commit (_`6545f78 on Mar 19, 2019`_). When iRRAM is being installed, it asks whether it will be installed directly in the system or in a local directory. For now, iRRAMx only detects iRRAM when it is installed in system.

Also,

```bash
$ sudo apt-get install libpng-dev
```



## Installation
iRRAMx is a static C++ library. Its compilation is deferred to the users.
When a compatible iRRAM is installed in the system, run `make` on the directory of the source directory creates a static library file in `lib/libiRRAMx.a`:

```
% mkdir ~/git_dir
% cd ~/git_dir
% git clone https://github.com/realcomputation/iRRAMx.git
% cd git_dir
% make
```

When building the library is complete, the library can be either included in the system. Or, the library and the header files can be linked each time when a C++ code using iRRAMx library is used:

```
% g++ -std=c++14 -O2 -Wall -I/path/to/git_dir/iRRAMx/include -L/path/to/git_dir/iRRAMx/lib -liRRAMx -liRRAM -lmpfr -lgmp my_cpp.cpp -o my_out
```

## Using iRRAMx
iRRAMx consists of orthogonal components:
- compact:
- linear: linear algebra including fast matrix multiplication, matrix eigenproblem, ...
- polynomial:
- random:
- plot:

When a component is going to be used, the corresponding header file should be used.
For example, when the user want to have infinite precision random number,
```
#include <iRRAM.h>
#include "iRRAMx/random.hpp"

void compute(){
  REAL x = Gaussian(); // A real number sampled from the normal distribution
  ...
}
```
Here, `compute` is a function that replaces `main` in iRRAM. The readers who are not familiar with iRRAM are referred to a [short tutorial](..).

## Reference Manual

- [iRRAMx/polynomial.h](https://github.com/realcomputation/iRRAMx/wiki/iRRAMx-linear.h)

### iRRAMx/compact.hpp

#### Compact

- `Compact()`
- `Compact(std::function<bool(Point<N>, int)> cfun)`

> *cfun*: the characteristic function that defines a compact set
>
> **returns** a compact set(defined by *cfun*)

- `Compact::plot2D(filename, width, x1, y1, x2, y2)`

> Create an image of the compact set on ℝ² in PNG format.

- `disjunction(Compact c1, Compact c2)`

> **returns** a compact set made by disjunction of *c1* and *c2*

- `conjunction(Compact c1, Compact c2)`

> **returns** a compact set made by conjunction of *c1* and *c2*

#### Path

**Path** inherits **Compact**.

- `Path()`
- `Path(std::function<Point<N>(REAL)> f)`

> *f*: a function [0,1]→ℝⁿ
>
> **returns** a path(defined by *f*)

#### Surface

**Surface** inherits **Compact**.

- `Surface()`
- `Surface(std::function<Point<N>(REAL,REAL)> f)`

> *f*: a function [0,1]²→ℝⁿ
>
> **returns** a surface(defined by *f*)


### iRRAMx/polynomial.hpp

### iRRAMx/random.hpp
#### random real numbers
- `REAL uniform_real()`

> **returns** z : REAL sampled from the uniform distribution on the unit interval _[0,1]_.

- `REAL uniform_real(REAL x, REAL y)`

> **returns** z : REAL sampled from the uniform distribution on the unit interval _[x,y]_.

- `REAL gaussian_real()`

> **returns** z : REAL sampled from the normal distribution.

- `REAL gaussian_real(REAL e, REAL s)`

> **returns** z : REAL sampled from the Gaussian distribution of an expectation e and a standard deviation s.

#### random complex numbers
- `COMPLEX uniform_complex()`

> **returns** z : COMPLEX sampled from the uniform distribution on the unit disc.

- `COMPLEX uniform_complex(COMPLEX c, REAL r)`

> **returns** z : COMPLEX sampled from the uniform distribution on the disc centered at c with radius r.

- `COMPLEX gaussian_complex()`

> **returns** z : COMPLEX sampled from the normal distribution on the complex plane.

#### random matrices
- `REALMATRIX gaussian_symmetric_matrix(unsigned int n)`

> **returns** a $n \times n$ random symmetric matrix A : REALMATRIX where each entry is sampled from the normal distribution.
- `REALMATRIX gaussian_asymmetric_matrix(unsigned int n)`

> **returns** a $n \times n$ random asymmetric matrix A : REALMATRIX where each entry is sampled from the normal distribution.
- `REALMATRIX gaussian_matrix(unsigned int n)`

> **returns** a $n \times n$ random matrix A : REALMATRIX where each entry is sampled from the normal distribution.
- `REALMATRIX haar_orthogonal_matrix(unsigned int n)`

> **returns** a $n \times n$ random orthogonal matrix which follows Haar distribution in $O(n)$. See _[Stewart, Gilbert W. "The efficient generation of random orthogonal matrices with an application to condition estimators." SIAM Journal on Numerical Analysis 17.3 (1980): 403-409.]_ for more detail.

### iRRAMx/plot.hpp

## Contributing to iRRAMx
This section is for those who want to contribute to this library.
For the purpose, we introduce the structure of this library.

This library consists of multiple sublibraries. (Here, sublibrary is not a technical term.. It means header files that the users will separately include in their applications.)
For each sublibrary, its public header file must stay in `./include/iRRAMx`, its header files that are used internally
must stay in its own directory in `./include/iRRAMx`, and its source codes must say in its own directory in `./src`.

For example, the public header file of _polynomial_ can be found in `./include/iRRAMx/polynomial.h`.
Some header files of _polynomial_ that are used only internally can be found in `./include/iRRAMx/polynomial/`.
The C++ source codes of _polynomial_ can be found in `./src/polynomial`.

New sublibraries can be freely created following the above organization. If the organizational requirement is satisfied, `make` will do every work for the compilation.


# iRRAM
iRRAM is a C++ library for exact real-number computation.
Due to the lack of any up-to-date documentation for it, we put some simple note here for the beginners.

## Installation
Though the official release can be found in the [official website](http://irram.uni-trier.de),
we (recommend to) use the recent version of it that can be found in their [github](https://github.com/fbrausse/iRRAM). Once the repository is cloned, run `./QUICKINSTALL_run_me`. It requires some autotools to be installed. (_please write which of those are required precisely_.)

Besides of those arbitrary installation related dependencies, it requires MPFR and GMP libraries installed.

If iRRAM is installed properly in the system, a cpp file using iRRAM can be compiled using the link option:
```
g++ -std=c++14 my_prog.cpp -liRRAM -lmpfr -lgmp
```

## How To Use
iRRAM provides one wrapper header file `iRRAM.h` that provides every functionalities of iRRAM.
Assuming, iRRAM is installed in the system, a default user program that uses iRRAM would look like this:
```
#include <iRRAM.h>

void compute(){
  ...
}
```

### compute()
`void compute()` is the function that replaces `main` when `iRRAM.h` is included.
There is a reason for this. If this technical reason is not so interesting,
this section can be skipped. Just regard `compute()` as the new `main()` and forget
everything about `main()`.

iRRAM achieves errorless computation via reiterations.
A real number is internally represented as an interval that is implemented using MPFR + GMP.
The idea is that when we perform interval computation repeatedly with higher and higher precision, it essentially is
equivalent to a computation without any errors.

A user will write an instruction for a computation. And, it needs to be computed repeatedly with
increasing precision. However, if the instruction is written in `main()` function, due to the
design of C++ language, it cannot be reiterated.
Therefore, what iRRAM does is as follows. `iRRAM.h` has already defined
`main` function which is designed to reiterate `compute()` function.

Now, when the user define their own `compute()` function regarding it as the `main` function,
the `main` function in `iRRAM.h` runs the `compute()` function with the additional cares: reiterate!

### REAL arithmetic
iRRAM provides a data type called REAL. It overloads the primitive operators.
For example, the following code
```
iRRAM::REAL x = 42;
iRRAM::REAL y = 20;
iRRAM::cout << x + y
```
assigns the real number 42 to the variable x and 20 to the variable y.
And, it prints x + y, which ought to be 62.

The idea is that iRRAM brings exact real numbers (REAL) as a first class citizen to C++.

### Partial comparison and Choose

### Limit
