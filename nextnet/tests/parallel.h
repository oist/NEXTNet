#pragma once

#include "nextnet/tests/stdafx.h"

#include "nextnet/random.h"

template <typename R>
std::vector<R> parallel(std::size_t n, rng_t &engine, std::function<R(rng_t &)> body)
{

    /* Pre-allocate result vector */
    std::vector<R> r(n, R());

    /* Create sub-RNGS so that each of the n work units gets their own RNG */
    sub_rngs rngs(n, engine);

    /* Execute work units, in parallel if possible */
    typedef boost::counting_iterator<std::size_t> count_it;
#if PARALLELIZE && HAVE_STD_EXECUTION
    /* Parallel execution policy should be available */
    std::transform(std::execution::par_unseq,
                   count_it(0), count_it(n), r.begin(),
                   [&body, &rngs](std::size_t i) -> R { return body(rngs[i]); });
#else
    /* Parallel execution policy likely not available */
    std::transform(count_it(0), count_it(n), r.begin(),
                   [&body, &rngs](std::size_t i) -> R { return body(rngs[i]); });
#endif

    return r;
}
