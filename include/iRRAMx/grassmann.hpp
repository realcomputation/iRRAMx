#pragma once
#include <iRRAM/REALMATRIX.h>
namespace iRRAM {

/** @defgroup grassmann Grassmannian
 *  @brief include "iRRAMx/grassmann.h"
 *  A Grassmannian is a linear subspace of a Euclidean space.
 *  A Grassmannian is represented by a basis,
 *  and the basis is represented by linearly independent rows of a matrix.
 * @{
 */

//! The class for Grassmanians.
class GRASSMANN {
   public:
    
    //! The basis of the Grassmannian.
    /*!
        Linearly independent rows of the matrix 
        are considered the basis of the Grassmannian.
        The behavior of this class is undefined
        if the rows are not linearly independent.
        (Your program will likely to be stuck in that case, though.)
    */
    REALMATRIX basis;

    //! Constructor with a basis matrix.
    /*!
        @param basis a matrix whose linearly independent rows represent the basis of Grassmannian
    */
    GRASSMANN(const REALMATRIX &basis);

    //! The dimension of the Grassmannian.
    /*!
        This is equal to the number of rows of the basis matrix.
        @return the dimension of the Grassmannian, not the ambient dimension
    */
    int dimension() const;

    //! The ambient dimension of the Grassmannain.
    /*!
        This is equal to the number of columns of the basis matrix.
        @return the ambient dimension of the Grassmannian, not the vector space dimension
    */
    int ambient_dimension() const;

    //! Compute the Minkowski sum of two Grassmannians.
    /*!
        @param dim the dimension of the resulting Grassmannian
        It needs the additional parameter dim to be computable.
    */
    GRASSMANN join(const GRASSMANN &rhs, int dim) const;

    //! Compute the intersection of two Grassmannians.
    /*!
        @param dim the dimension of the resulting Grassmannian
        It needs the additional parameter dim to be computable.
    */
    GRASSMANN meet(const GRASSMANN &rhs, int dim) const;

    //! Compute the orthogonal complement.
    /*!

    */
    GRASSMANN complement() const;

    //! Compute the projection of one Grassmannian into the other Grassmannian.
    /*!
        @param dim the dimension of the resulting Grassmannian
        @param space the Grassmannian into which this Grassmannian is projected
        It needs the additional parameter dim to be computable.
        Call A.project_onto(B, d) when you want project A onto B where the resulting dimension is d.
    */
    GRASSMANN project_onto(const GRASSMANN &space, int dim) const;

    //! Print to iRRAM::std the numeric values of the basis matrix.
    /*!
        
    */
    void show() const;

   private:
    static void find_nonzero(const REALMATRIX &M, int &pi, int &pj,
                             int Li, int Hi, int Lj, int Hj);
    static void zassenhaus(REALMATRIX &M, int d);
    static REALMATRIX concat(const REALMATRIX &M00,
                             const REALMATRIX &M01,
                             const REALMATRIX &M10,
                             const REALMATRIX &M11);
    static REALMATRIX submatrix(const REALMATRIX &M,
                                int Li, int Lj, int Hi, int Hj);
    static REALMATRIX transpose(const REALMATRIX &M);
};

/** @} */

}  // namespace iRRAM