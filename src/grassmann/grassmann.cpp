
#include "iRRAMx/grassmann.hpp"

#include <assert.h>
#include <iRRAM/REALMATRIX.h>

using namespace iRRAM;

namespace iRRAM {

GRASSMANN::GRASSMANN(const REALMATRIX &M) : basis(M) {}
int GRASSMANN::dimension() const { return rows(basis); }
int GRASSMANN::ambient_dimension() const { return columns(basis); }

GRASSMANN GRASSMANN::join(const GRASSMANN &rhs, int dim) const {
    assert(ambient_dimension() == rhs.ambient_dimension());
    //int d0 = dimension(); Unused variable
    int d1 = rhs.dimension();
    int a = ambient_dimension();

    REALMATRIX M = concat(basis, basis, rhs.basis, zeroes(d1, a));
    zassenhaus(M, dim);
    
    REALMATRIX S = submatrix(M, 0, dim - 1, 0, a - 1);
    return GRASSMANN(S);
}

GRASSMANN GRASSMANN::meet(const GRASSMANN &rhs, int dim) const {
    assert(ambient_dimension() == rhs.ambient_dimension());
    int d0 = dimension();
    int d1 = rhs.dimension();
    int a = ambient_dimension();

    REALMATRIX M = concat(basis, basis, rhs.basis, zeroes(d1, a));
    zassenhaus(M, d0 + d1 - dim);
    
    REALMATRIX S = submatrix(M, d0 + d1 - dim, rows(M)-1, a, a*2 - 1);
    return GRASSMANN(S);
}

GRASSMANN GRASSMANN::complement() const { 
    int a = ambient_dimension();
    int d = dimension();
    REALMATRIX M = concat(transpose(basis), eye(a), REALMATRIX(0,d), REALMATRIX(0,a));

    // Gaussian Elimination
    for (int i=0; i<d; i++) {
        int pi, pj;
        find_nonzero(M, pi, pj, i, rows(M)-1, i, i);

        for (int j=i; j<(int)columns(M); j++) swap(M(i, j), M(pi, j));
        for (int ii=i+1; ii<(int)rows(M); ii++) {
            REAL c = M(ii, i) / M (i, i);
            for (int j=i; j<(int)columns(M); j++)
                M(ii, j) -= M(i, j) * c;
        }
    }

    REALMATRIX S = submatrix(M, d, rows(M)-1, d, columns(M) - 1);
    return GRASSMANN(S); 
}

GRASSMANN GRASSMANN::project_onto(const GRASSMANN &space, int dim) const {
    assert(ambient_dimension() == space.ambient_dimension());
    REALMATRIX A = transpose(space.basis);
    REALMATRIX M = basis * (A / (space.basis * A)) * space.basis;

    // Gaussian Elimination
    for (int i=0; i<dim; i++) {
        int pi, pj;
        find_nonzero(M, pi, pj, i, rows(M)-1, 0, columns(M)-1);

        for (int j=i; j<(int)columns(M); j++) swap(M(i, j), M(pi, j));
        for (int ii=i+1; ii<(int)rows(M); ii++) {
            REAL c = M(ii, i) / M (i, i);
            for (int j=0; j<(int)columns(M); j++)
                M(ii, j) -= M(i, j) * c;
        }
    }

    REALMATRIX S = submatrix(M, 0, dim-1, 0, columns(M)-1);
    return GRASSMANN(S);
}

void GRASSMANN::show() const {
    cout << "----------------------------\n";
    for (int i = 0; i < dimension(); i++) {
        for (int j = 0; j < ambient_dimension(); j++) {
            cout << basis(i, j) << " ";
        }
        cout << "\n";
    }
    cout << "----------------------------\n";
}

void GRASSMANN::find_nonzero(const REALMATRIX &M, int &pi, int &pj, int Li,
                             int Hi, int Lj, int Hj) {
    REAL maxval(0);
    for (int i = Li; i <= Hi; i++)
        for (int j = Lj; j <= Hj; j++) maxval = maximum(maxval, abs(M(i, j)));

    for (int i = Li; i <= Hi; i++)
        for (int j = Lj; j <= Hj; j++)
            if (choose(0 < abs(M(i, j)), abs(M(i, j)) < maxval) == 1) {
                pi = i;
                pj = j;
                return;
            }

    // All entries within the range are zero
    assert(false);
}

void GRASSMANN::zassenhaus(REALMATRIX &M, int d) {
    assert(columns(M) % 2 == 0);

    for (int i = 0; i < (int)rows(M); i++) {
        int pi, pj;
        if (i < d)
            GRASSMANN::find_nonzero(M, pi, pj, i, rows(M) - 1, 0, columns(M) / 2 - 1);
        else
            GRASSMANN::find_nonzero(M, pi, pj, i, rows(M) - 1, columns(M) / 2,
                                    columns(M) - 1);

        // Swap two rows.
        for (int j = 0; j < (int)columns(M); j++) swap(M(i, j), M(pi, j));

        // Entries below M(i, pj) are made zeros.
        for (int ii = i + 1; ii < (int)rows(M); ii++) {
            REAL c = M(ii, pj) / M(i, pj);
            for (int j = 0; j < (int)columns(M); j++) {
                M(ii, j) -= M(i, j) * c;
            }
        }
    }
}

REALMATRIX GRASSMANN::concat(const REALMATRIX &M00,
                                    const REALMATRIX &M01,
                                    const REALMATRIX &M10,
                                    const REALMATRIX &M11) {
    assert(rows(M00) == rows(M01));
    assert(rows(M10) == rows(M11));
    assert(columns(M00) == columns(M10));
    assert(columns(M01) == columns(M11));

    int r0 = rows(M00);
    int c0 = columns(M00);
    int r1 = rows(M11);
    int c1 = columns(M11);
    
    REALMATRIX A(r0 + r1, c0 + c1);
    for (int i = 0; i < r0; i++) for (int j = 0; j < c0; j++) A(i,    j   ) = M00(i, j);
    for (int i = 0; i < r0; i++) for (int j = 0; j < c1; j++) A(i,    c0+j) = M01(i, j);
    for (int i = 0; i < r1; i++) for (int j = 0; j < c0; j++) A(r0+i, j   ) = M10(i, j);
    for (int i = 0; i < r1; i++) for (int j = 0; j < c1; j++) A(r0+i, c0+j) = M11(i, j);
    return A;
}

REALMATRIX GRASSMANN::submatrix(const REALMATRIX &M, int Li, int Hi, int Lj,
                                int Hj) {
    REALMATRIX A(Hi - Li + 1, Hj - Lj + 1);
    for (int i = 0; i < (int)rows(A); i++)
        for (int j = 0; j < (int)columns(A); j++) A(i, j) = M(Li + i, Lj + j);
    return A;
}

REALMATRIX GRASSMANN::transpose(const REALMATRIX &M) {
    REALMATRIX R(columns(M), rows(M));
    for (int i=0; i<(int)rows(R); i++)
        for (int j=0; j<(int)columns(R); j++)
            R(i, j) = M(j, i);
    return R;
}

}  // namespace iRRAM