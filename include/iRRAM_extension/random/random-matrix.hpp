#pragma once

#include <iRRAM/lib.h>
#include "iRRAM_extension/random/random-core.hpp"
#include "iRRAM_extension/random/random-real.hpp"

using namespace iRRAM;
namespace iRRAM{

/** @defgroup randmat Matrix Sampler
 *  @brief Sampling real-valued matrix
 *  @ingroup rand
 *  @{
 */

//! Gaussian symmetric matrix sampler
/*!
    returns nXn symmetric matrix whose entry follows normal distribution i.i.d.
*/ 
REALMATRIX gaussian_symmetric_matrix(unsigned int n);

//! Gaussian assymmetric matrix sampler
/*!
    returns nXn assymetric matrix whose entry follows normal distribution i.i.d.
*/
REALMATRIX gaussian_asymmetric_matrix(unsigned int n);

//! Gaussian regular Matrix sampler
/*!
    returns nXn square matrix which is regular almost surely whose entry follows normal distribution i.i.d.
*/
REALMATRIX gaussian_matrix(unsigned int n);

//! Gaussian regular Matrix sampler
/*!
    returns uniformly distributed orthogonal matrix (that follows Haar measure in O(n)).
    See [Stewart, Gilbert W. "The efficient generation of random orthogonal matrices with an application to condition estimators." SIAM Journal on Numerical Analysis 17.3 (1980): 403-409.] for more detail.
*/
REALMATRIX haar_orthogonal_matrix(unsigned int n);

/** @} */
}