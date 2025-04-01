#pragma once

#include "stdafx.h"
#include "types.h"

template <typename T>
struct permutation
{
    typedef T index_t;

    permutation(){};

    permutation(index_t n, rng_t &engine)
    {
        /* Fill vector with 0, 1, ...  then shuffle */
        p.resize(n);
        std::iota(p.begin(), p.end(), 0);
        std::shuffle(p.begin(), p.end(), engine);
    }

    index_t operator[](index_t i) const
    {
        if ((i < 0) || (i >= (index_t)p.size()))
            return i;
        else
            return p[i];
    }

private:
    std::vector<index_t> p;
};
