#include "tests/stdafx.h"

#include "graph.h"

TEST_CASE("fully connected networks", "[graph]") {
    std::mt19937 engine;
    const int M = 3;
    fully_connected nw(M, engine);
    for(int n=0; n < M; ++n) {
        REQUIRE(nw.outdegree(n) == M-1);
        int sum=0;
        for(int i=0; i < M-1; ++i) {
            const int np = nw.neighbour(n, i);
            REQUIRE(np >= 0);
            REQUIRE(np < M);
            REQUIRE(np != n);
            sum += np;
        }
        // every node except n must be a neighbour
        REQUIRE(sum == (M*(M-1)/2 - n));
    }
}

TEST_CASE("modified R0 for acyclic graphs","[graph]") {
    std::mt19937 engine; 

    for(double mean = 2; mean < 4; ++mean){

        double modified_r0 = acyclic::lambda(mean, 6);

        // std::vector<int> data({});
        std::poisson_distribution<> poisson( modified_r0 );
        
        double s = 0.0;
        int n = 0;
        for(int i=0; i < 100000; ++i){
            const int r = poisson(engine);
            if (r == 0)
                continue;

            s += r;
            n += 1;
        }
        double empirical_mean = s / n;
        REQUIRE(std::abs(mean - empirical_mean) < 0.01);
    }
}
