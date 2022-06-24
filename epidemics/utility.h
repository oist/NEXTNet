////  utility.h
//  epidemics_update
//
//  Created by Cure Samuel Cyrus on 2022/02/19.
//
#pragma once

#include "stdafx.h"
#include "types.h"

/**
 * @brief Hashing support for std::pair, required for use in std::unorderedet_map
 * 
 * The C++ standard library weirdly enough does not provide this. Simply
 * XOR-ing the two hashes is pretty silly, for it seems sufficient for the
 * std::unordered_map declared below.
 */
struct pair_hash {
    template<typename T, typename U>
    std::size_t operator () (const std::pair<T,U> &p) const {
        return std::hash<T>{}(p.first) ^ std::hash<U>{}(p.second);
    }
};  

/**
 * @brief Approximates m with p/q under the constraint that q <= L
 * @param x value to approximate
 * @param L largest allowed denominator
 * @return the pair (p, q)
 */
std::pair<unsigned int, unsigned int> fraction(double x, unsigned int L);

template<typename T, typename ...Args>
double inverse_survival_function(double u, double precision, T f, Args... args) {
    if ((u > 1) || (u < 0))
        return NAN;
    /* Use bisection to invert f. We assume that f(0) = 1 >= u,
     * so we canstart with the left bound l = 0, but we have to find
     * a right bound r such that  phi(r) <= u */
    double l = 0;
    double f_l = 1;
    double r = 1;
    double f_r = f(r, std::forward(args)...);
    while (std::isfinite(r) && (f_r > u)) {
        r *= 2;
        f_r = f(r, std::forward(args)...);
    }
    /* Now we split the interval and pick the left or right subinterval
     * until we reach the desired precision
     */
    while ((f_l != f_r) && ((r - l) > precision)) {
        const double m = std::isfinite(r) ? (l + r) / 2 : l*2 ;
        const double f_m = f(m, std::forward(args)...);
        if (f_m >= u) {
            l = m;
            f_l = f_m;
        }
        else {
            r = m;
            f_r = f_m;
        }
    }
    /* Return mid-point of the final interval */
    return (l + r) / 2;
}
