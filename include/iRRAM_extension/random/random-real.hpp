#pragma once

#include <iRRAM/lib.h>
#include "iRRAM_extension/random/random-core.hpp"

using namespace iRRAM;
namespace iRRAM{

/** @defgroup randreal Real Sampler
 *  @brief Sampling real number from unifrom or gaussian distribution
 *  @ingroup rand
 *  @{
 */

//! Uniform real sampler
/*!
    returns REAL sampled from the uniform distribution on the unit interval [0,1]
*/
REAL uniform_real();
//! Uniform real sampler
/*!
    returns REAL sampled from the uniform distribution on the closed interval [a,b].
*/
REAL uniform_real(REAL a, REAL b);
//! Gaussian real sampler
/*!
    returns REAL sampled from the gaussian distribution with mean 0 and  standard deviation 1.
*/
REAL gaussian_real();

//! Gaussian real sampler
/*!
    returns REAL sampled from the gaussian distribution with mean exp and  standard deviation std.
*/
REAL gaussian_real(REAL exp, REAL std);


REAL linear_real();

/** @} */
}

