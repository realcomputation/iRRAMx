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
        C choice, const Args &... x);


// [[deprecated]]
// template <typename R, typename C, typename... Args>
// std::enable_if_t<internal::is_cacheable<C>::value, R>
// limit_mv(R (*seq)(int prec, C &choice, const Args &... args),
//          C choice, const Args &... x);
