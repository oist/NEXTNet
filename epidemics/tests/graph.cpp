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

TEST_CASE("acyclic networks", "[graph]") {
    std::mt19937 engine;
    const int K = 3;

    const int M = 50;
    const int N = 100;
    double total_neighbours = 0.0;
    double total_nodes = 0.0;
    for(int m=0; m < M; ++m) {
        std::unordered_set<node_t> visited;
        acyclic nw(K, engine);

        // traverse first N nodes of the network in
        // breadth-first order
        int n = 0;
        std::deque<node_t> queue = { 0 };
        while (!queue.empty()) {
            // stop after visiting N nodes of the current network
            if (n >= N)
                break;
            n += 1;
            // get next node to visit
            const node_t node = queue.front();
            queue.pop_front();
            // check if we've seen it before
            REQUIRE(visited.find(node) == visited.end());
            visited.insert(node);
            // get node's degree
            const index_t deg = nw.outdegree(node);
            total_neighbours += deg;
            total_nodes += 1;
            // enqueue non-visited neighbours at the end of the queue,
            // which for all nodes except the first are neighbours 1,...,deg-1 
            if (node > 0)
                REQUIRE(visited.find(nw.neighbour(node, 0)) != visited.end());
            for(index_t i=1; i < deg; ++i) {
                const node_t q = nw.neighbour(node, i);
                REQUIRE(q > 0);
                queue.push_back(q);
            }
        }
    }

    const double sd = std::sqrt(K / total_nodes);
    REQUIRE(std::abs(K - total_neighbours / total_nodes) < 3.0 * sd);
}
