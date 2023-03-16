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

template<typename T>
struct interger_set {
    std::size_t size() const {
        return size;
    }

    bool insert(T e) {
        /* The possible cases of existing ranges are:
         *    I: [a,     e-1]  [c > e+1, d] ==> extend [a, e-1] to [a, e]
         *   II: [a, b < e-1]  [e+1    , d] ==> extend [e+1, d] to [e, d]
         *  III: [a,     e-1]  [e+1    , d] ==> merge [a, e-1] and [e+1, d] to [a, d]
         *  IVa: [a, b < e-1]               ==> insert [e, e]
         *  IVb: [a, b < e-1]  [c > e+1, d] ==> insert [e, e]
         *   Va: [a <= e-1, b >= e]         ==> do nothing
         *   Vb: [e       , b >= e]         ==> do nothing
         */

        if ((e == std::numeric_limits<T>::min()) || (e == std::numeric_limits<T>::max()))
            throw std::range_error("element in a integer_set cannot be maximum or minimum of type's domain");

        /* Find the first range [x, y] that does not wholly lie left of e - 1, i.e. where y >= e - 1*/
        const auto i = ranges.lower_bound(e - 1);

        /* IVa or IVb. Just insert [e, e] */
        if ((i == ranges.end()) || (i->first > e + 1)) {
            const range r(e, e);
            ranges.insert(r);
            size++;
            return true;
        }

        /* I, III or Va. If the range contains e - 1, extend or merge it */
        if (i->first <= e - 1) {
            /* Va. If it already extends past e, we're done */
            if (i->last >= e)
                return false;
            assert(i->last == e - 1);

            /* III. Before we extend the range, check whether we would create two adjacent ranges */
            const s = std::next(i);
            if ((s != ranges.end()) && (s->first == e + 1)) {
                /* Next range is adjacent, i.e. we have [a, e - 1], [ e, b ]. Merge ranges. */
                const range r(i->first, s->last);
                ranges.erase(i);
                ranges.erase(s);
                ranges.insert(r);
                size++;
                return true;
            }

            /* I. Finally, extend the range */
            const range r(i->first, e);
            ranges.erase(i);
            ranges.insert(r);
            size++;
            return true;
        }

        /* Vb. The range contains e, we're done */
        if (i->first <= e) {
            assert(i->second >= e);
            return false;
        }

        /* II. Extend the range. */
        assert(i->first == e + 1);
        const range r(e, i->second);
        ranges.erase(i);
        ranges.insert(r);
        size++;
        return true;
    }

    bool erase(T e) {
        /* The possible cases are:
         *    Ia: [a < e, e]      ==> restrict to [a, e-1]
         *    Ib: [e, b > e]      ==> restrict to [e+1, b]
         *    Ic: [e, e]          ==> remove
         *    II: [a < e, b > e]  ==> split into [a, e-1], [e+1, b]
         *  IIIa: [a > e, b]      ==> do nothing
         *  IIIb: [a, b < e]      ==> do nothing
         */

        if ((e == std::numeric_limits<T>::min()) || (e == std::numeric_limits<T>::max()))
            throw std::range_error("element in a integer_set cannot be maximum or minimum of type's domain");

        /* Find the first range [x, y] that does not wholly lie left of e, i.e. where y >= e*/
        const auto i = ranges.lower_bound(e);

        /* IIIa or IIIb. Do nothing */
        if ((i == ranges.end()) || (i->first > e))
            return false;
        assert(i->second >= e);

        /* Ia or Ic. Restrict range */
        if (i->second == e) {
            ranges.erase(i);
            /* Ia. Insert restricted range if non-empty */
            if (i->first < e) {
                const range r(i->first, e-1);
                ranges.insert(r);
            }
            size--;
            return true;
        }

        /* Ib. Restrict range */
        if (i->first == e) {
            const range r(e+1, i->second);
            ranges.erase(i);
            ranges.insert(r);
            size--;
            return true;
        }

        /* II. Split range */
        assert(i->first <= e-1);
        assert(i->second >= e+1);
        const range r1(i->first, e-1);
        const range r2(e+1, i->second);
        ranges.erase(i);
        ranges.insert(r1);
        ranges.insert(r2);
        size--;
        return true;
    }

    T draw_present(rng_t& engine) {
        /* Draw random range, with probabilities prop. to the ranges' lengths */
        const double p = std::uniform_real_distribution<>(0, 1)(engine);
        const std::size_t l = 0;
        for(const auto& r: ranges) {
            /* Compute fraction q covered by ranges up to the current one */
            l += (r.last - r.first + 1);
            const double q = (double)l / size;
            /* If that fraction exceeds p, we found our range */
            if (q > p) {
                /* Pick uniformly random from the range */
                return std::uniform_int_distribution<T>(r.first, r.last)(engine);
            }
        }

        throw std::logic_error("draw_present() failed due to internal inconsistency");
    }

    T draw_absent(T lb, T ub, rng_t& engine) {
        if (!ranges.empty() && (ranges.first().first < lb))
            throw std::range_error("lower bound is larger than smallest element of integer_set");
        if (!ranges.empty() && (ranges.last().second > ub))
            throw std::range_error("upper bound is smaller than largest element of integer_set");

        /* Draw random gap between ranges, with probabilities prop. to the gaps' lengths */
        const double p = std::uniform_real_distribution<>(0, 1)(engine);
        const std::size_t size_gaps = ((st::size_t)(ub - lb) + 1 - size);
        const std::size_t l = 0;
        const T a = lb;
        auto i = ranges.begin();
        const auto e = ranges.end();
        do {
            /* Gap extends to beginning of next range or upper bound for last range */
            range gap;
            if (i != e) {
                gap = range(a, i->first - 1);
                a = i->second + 1;
                i++;
            } else {
                gap = range(a, ub);
                done = true;
            }
            /* Compute fraction q covered by gaps up to the current one */
            l += (gap.last - gap.first + 1);
            const double q = (double)l / size_gaps;
            /* If that fraction exceeds p, we found our range */
            if (q > p) {
                /* Pick uniformly random from the gap */
                return std::uniform_int_distribution<T>(gap.first, gap.second)(engine);
            }
        } while (!done);

        throw std::logic_error("draw_present() failed due to internal inconsistency");
    }

private:
    struct range {
        range(T f, T l) :first(f), last(l)
        {}

        T first;
        T last;
    };

    std::size_t size;
    std::set<range> ranges;
};
