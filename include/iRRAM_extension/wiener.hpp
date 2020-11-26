#pragma once
#include "iRRAM_extension/random.hpp"
#include "iRRAM_extension/utility.hpp"
#include<functional>
#include<map>
using namespace iRRAM;

namespace iRRAM{

/** @defgroup wiener Function Sampler
 *  @brief Sampling random continuous function on [0,1] from Wiener distribution
 *  @ingroup rand
 *  @{
*/

//!  A class for random function
/*!
    The instance of this class indicates the single sample path of wiener process which is sample from
    Wiener distribution
*/
class WIENER
{
    private:
        std::map< std::pair<INTEGER,INTEGER>, REAL > X;
        REAL X_0;
        int curPrec,probPrec;
        int modulus(int ,int );
        REAL getX(int ,int );
    public:
        //! Generate the sample path in probability
        /*!
            Generate the sample path computed with probability 1-2^(-50)
            using Levy-ciesielski construction. See [
                Mörters, P., & Peres, Y. (2010). Brownian motion (Vol. 30)
            ]
        */
        WIENER();

        //! Generate the sample path in probability
        /*!
            Generate the sample path computed with probability 1-2^(-p)
            using Levy-ciesielski construction. See [
                Mörters, P., & Peres, Y. (2010). Brownian motion (Vol. 30)
            ]
        */
        WIENER(int p);

        //! Manually set the probability bound and desired precision
        /*!
            @param p1 the desired precision of sampled path
            @param p2 the probability bound.
        */
        void setPrec(int p1,int p2);
        //! Compute the generated sample path on t
        /*!
            @param t must be in [0,1]
        */
        REAL compute(REAL t);

        //! operator overloading for compute()
        /*!
            @param t must be in [0,1]
        */
        REAL operator() (REAL t);
};

/** @} */

}