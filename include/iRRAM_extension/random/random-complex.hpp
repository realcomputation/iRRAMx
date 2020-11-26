#pragma once

#include <iRRAM/lib.h>
#include "iRRAM_extension/random/random-core.hpp"
#include "iRRAM_extension/random/random-real.hpp"

using namespace iRRAM;
namespace iRRAM{

/** @defgroup randcomp Complex Sampler
 *  @brief Sampling Complex number
 *  @ingroup rand
 *  @{
 */


//! Uniform complex number sampler
/*!
    returns COMPLEX sampled from the uniform distribution on the unit disc.
*/
COMPLEX uniform_complex();

//! Uniform complex number sampler
/*!
    returns COMPLEX sampled from the uniform distribution on the disc centered at c with radius r.
*/
COMPLEX uniform_complex(COMPLEX, REAL);

//! Gaussian complex number sampler
/*!
    returns COMPLEX sampled from the normal distribution on the complex plane.
*/
COMPLEX gaussian_complex();

/** @} */
}
