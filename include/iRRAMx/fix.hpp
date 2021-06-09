#pragma once
#include "iRRAM/lib.h"
#include "iRRAM/core.h"
using namespace iRRAM;


template <typename R, typename C, typename... Args>
[[deprecated("Use limti_mv_() instead.")]]
std::enable_if_t<internal::is_cacheable<C>::value, R>
limit_mv(R (*seq)(int prec, C &choice, const Args &... args),
         C choice, const Args &... x);

template <typename R, typename C, typename... Args>
std::enable_if_t<internal::is_cacheable<C>::value, R>
limit_mv_(R (*seq)(int prec, C &choice, const Args &... args),
        C choice, const Args &... x) {
    bool inlimit = actual_stack().inlimit != 0;

    limit_computation env;

    R lim, limnew;
    sizetype limnew_error, element_error;
    sizetype lim_error;

    int element = env.saved_prec();
    int element_step = env.saved_step();
    int firsttime = 2;

    if (!inlimit &&
        !get_cached(choice)) /* TODO: why only in case of !inlimit? */
        put_cached(choice);

    int x_error_exp = geterror_exp(x...);

    limit_debug("starting gen limit_mv");

    while (1) {
        try {
            iRRAM_DEBUG2(2, "trying to compute gen limit_mv with "
                            "precicion 2^(%d)...\n",
                         element);


            limnew = seq(element, choice, x...);


            modify_cached(choice);
            element_error = sizetype_power2(element);
            limnew.geterror(limnew_error);
            limnew_error += element_error;
            if (firsttime == 2)
                if (limnew_error.exponent >
                    env.saved_prec(-1) &&
                    limnew_error.exponent >
                    x_error_exp - env.saved_prec(-1)) {
                    iRRAM_DEBUG2(2,
                                 "computation not precise "
                                 "enough (%d*2^%d), trying "
                                 "normal p-sequence\n",
                                 limnew_error.mantissa,
                                 limnew_error.exponent);
                    element_step = 1;
                    element =
                            4 +
                            iRRAM_prec_array[element_step];
                    firsttime = 1;
                }
            if (firsttime != 0 || sizetype_less(limnew_error, lim_error)) {
                lim = limnew;
                lim_error = limnew_error;
                iRRAM_DEBUG2(2, "getting result with error %d*2^(%d)\n",
                             lim_error.mantissa, lim_error.exponent);
            } else {
                iRRAM_DEBUG1(
                        2,
                        "computation successful, but no improvement\n");
            }
            firsttime = 0;
            if (element <= env.saved_prec())
                break;
            element_step += 4;
            element = iRRAM_prec_array[element_step];

        } catch (const Iteration &it) {
            if (firsttime == 0) {
                iRRAM_DEBUG1(2, "computation failed, using "
                                "best success\n");
                break;
            } else if (firsttime == 2) {
                iRRAM_DEBUG1(2, "computation failed, trying "
                                "normal p-sequence\n");
                element_step = 1;
                element = 4 + iRRAM_prec_array[element_step];
                firsttime = 1;
            } else {
                iRRAM_DEBUG1(1, "computation of limit_gen1 "
                                "failed totally\n");
                iRRAM_REITERATE(0);
            }
        }
    }
    lim.seterror(lim_error);
    iRRAM_DEBUG2(2, "end of gen limit_mv with error %d*2^(%d)\n",
                 lim_error.mantissa, lim_error.exponent);
    return lim;
}


// [[deprecated]]
// template <typename R, typename C, typename... Args>
// std::enable_if_t<internal::is_cacheable<C>::value, R>
// limit_mv(R (*seq)(int prec, C &choice, const Args &... args),
//          C choice, const Args &... x);
