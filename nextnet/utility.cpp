//
//  utility.cpp
//  epidemics_update
//
//  Created by Cure Samuel Cyrus on 2022/02/19.
//

#include "utility.h"

std::pair<unsigned int, unsigned int> fraction(double x, unsigned int L)
{
    // Algorithm due to Farey, see e.g. https://www.johndcook.com/blog/2010/10/20/best-rational-approximation/
    // We track two approximation, a/b <= x <= c/d
    unsigned int a = 0, b = 1;
    unsigned int c = 1, d = 1;

    while (true) {
        // If one approximation's denominator exceeds L, the other is the optional one
        if (b > L)
            return { c, d };
        else if (d > L)
            return { a, b };
        assert((b <= L) && (d <= L));

        // Next approximation is x' = p/q = (a+c) / (b+d)
        const unsigned int p = a + c, q = b + d;
        const double xp = (double)p / q;
        if (x == xp) {
            // We hit x exactly, optimal approximation is one of
            // a/b, c/d and p/q.
            if (p <= L) return { p, q };
            return (b > d) ? std::make_pair(a, b) : std::make_pair(c, d);
        } else if (x > xp) {
            a = p;
            b = q;
        } else {
            c = p;
            d = q;
        }
    }
}
