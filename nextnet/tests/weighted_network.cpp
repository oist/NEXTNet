#include "nextnet/tests/stdafx.h"

#include "nextnet/weighted_network.h"
#include "nextnet/NextReaction.h"

/**
 * @brief Test case to verify `erdos-reyni`
 */
TEST_CASE("Weighted Erdös-Reyni networks", "[weighted_network]")
{
    std::mt19937 engine;
    const int M = 100;
    for (int n : { 2, 3, 5, 10, 100 }) {
        for (double d : { 0.0, 0.1, 1.0, 2.0, 10.0, n - 2.0, n - 1.0 }) {
            // the complete graph with d = n - 1 has the largest degree possible
            if (d > n - 1)
                continue;

            // compute the mean degree
            double dsum = 0.0;
            for (int i = 0; i < M; ++i) {
                weighted_erdos_reyni nw(n, d, { 1.0, 5.0 }, { 3.0 / 4.0, 1.0 / 4.0 }, engine);
                for (int j = 0; j < n; ++j)
                    dsum += nw.outdegree(j);
            }
            const double dmean = dsum / ((double)n * (double)M);

            // the sum of degree is twice the sum of n * (n-1) / 2 indicator
            // variables which represent the presence or absence of an
            // edge and which are 1 with probability p = d / (n-1). The
            // std. dev. of a single indicator is sqrt(p * (1 - p)), the std.
            // dev of the sum of degree is thus
            //     2 * sqrt( n * (n-1) * p * (1 - p) / 2 )
            //   = sqrt(2 * n * d * (n - d - 1) / (n - 1)).
            // the std. dev. of the mean degree is then
            //      sqrt(2 * d * (n - d - 1) / n * (n - 1))
            const double sd = std::sqrt(2.0 * d * (n - d - 1) / (M * n * (n - 1)));
            REQUIRE(std::abs(d - dmean) <= 4.0 * sd);
        }
    }

    // TODO: Test the weights
}
