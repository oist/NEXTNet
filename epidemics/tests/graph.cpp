#include "tests/stdafx.h"

#include "graph.h"

/**
 * @brief Test case to verify `fully_connected`
 */
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

/**
 * @brief Test case to verify `erdos-reyni`
 */
TEST_CASE("Erdös-Reyni networks", "[graph]") {
    std::mt19937 engine;
    const int M = 100;
    for(int n: {2, 3, 5, 10, 100}) {
        for(double d: {0.0, 0.1, 1.0, 2.0, 10.0, n-2.0, n-1.0}) {
            // the complete graph with d = n - 1 has the largest degree possible
            if (d > n - 1)
                continue;

            // compute the mean degree
            double dsum = 0.0;
            for(int i=0; i < M; ++i) {
                erdos_reyni nw(n, d, engine);
                for(int j=0; j < n; ++j)
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
            REQUIRE(std::abs(d - dmean) <= 3.0 * sd);
          }
    }
}

/**
 * @brief Test case to verify `acyclic::lambda`
 *
 * The static member `lambda` finds the Poisson parameter lambda
 * such that the distribution conditioned on k >= 1 has the
 * specified mean.
 */
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

/**
 * @brief Test case to verify `acyclic` in the uniform-degree case
 *
 * Here, the tree-like network has the same degree distribution
 * for every node. If we count only edges pointing *away* from the
 * root, the degree of the root is thus on average one higher than
 * that of any other node (because it has no edge pointing towards
 * the root, while all other nodes have exactly one)
 */
TEST_CASE("acyclic networks (reduce_root_degree=false)", "[graph]") {
    std::mt19937 engine;
    const int K = 3;

    // number of networks to generate and test
    const int M = 1000;
    // number of nodes of each network to scan 
    const int N = 1 + K + K*(K-1);
    double total_neighbours = 0.0;
    double total_nodes = 0.0;
    for(int m=0; m < M; ++m) {
        std::unordered_set<node_t> visited;
        acyclic nw(K, false, engine);

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
            // enqueue non-visited neighbours at the end of the queue
            int v = 0;
            for(index_t i=0; i < deg; ++i) {
                const node_t q = nw.neighbour(node, i);
                REQUIRE(q >= 0);
                if (visited.find(q) == visited.end())
                    queue.push_back(q);
                else
                    v += 1;
            }
            // for non-roots, exactly one neighbour must have been already visited
            if (node == 0)
                REQUIRE(v == 0);
            else
                REQUIRE(v == 1);
            
        }
    }

    const double sd = std::sqrt(K / total_nodes);
    REQUIRE(std::abs(K - total_neighbours / total_nodes) < 3.0 * sd);
}

/**
 * @brief Test case to verify `acyclic` in the uniform-degree case
 *
 * Here, the tree-like network has the same distribution of the
 * number of edges pointing *away* from the root for every node.
 * If we count all edges, the average degree of the root is thus one
 * lower than the average degree of all other nodes (because it has no
 * edge pointing towards the root, while all other nodes have exactly one)
 */
TEST_CASE("acyclic networks (reduce_root_degree=true)", "[graph]") {
    std::mt19937 engine;
    const int K = 3;

    // number of networks to generate and test
    const int M = 1000;
    // number of nodes of each network to scan 
    const int N = 1 + (K-1) + (K-1)*(K-1);
    double tota_outgoing_neighbours = 0.0;
    double total_nodes = 0.0;
    for(int m=0; m < M; ++m) {
        std::unordered_set<node_t> visited;
        acyclic nw(K, true, engine);

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
            total_nodes += 1;
            // enqueue non-visited neighbours at the end of the queue
            int v = 0;
            for(index_t i=0; i < deg; ++i) {
                const node_t q = nw.neighbour(node, i);
                REQUIRE(q >= 0);
                if (visited.find(q) == visited.end()) {
                    tota_outgoing_neighbours += 1;
                    queue.push_back(q);
                }
                else
                    v += 1;
            }
            // for non-roots, exactly one neighbour must have been already visited
            if (node == 0)
                REQUIRE(v == 0);
            else
                REQUIRE(v == 1);
        }
    }

    const double sd = std::sqrt(K / total_nodes);
    REQUIRE(std::abs((K-1) - tota_outgoing_neighbours / total_nodes) < 3.0 * sd);
}


/**
 * @brief Test case to verify `cm::cm`
 *
 * The variance and mean degree of the resulting network should match
 * with the initial degree list that was the input.
 * 
 */
TEST_CASE("Configuration model networks","[graph]") {
    std::mt19937 engine;
    const int size = 10000;

    // Generate a Poisson graph with the configuration model
    std::poisson_distribution<> poisson(3);    
    std::vector<int> degreeList(size,0);
    std::size_t total_degree = 0;
    for (int i = 0; i < size; i++)
    {
        const int k = poisson(engine); 
        degreeList[i] = k;
        total_degree += k;
    }

    // make sure the total degree is even, otherwise no graph can exist
    while (total_degree % 2 == 1) {
        // re-generate a random degree
        const std::size_t i = std::uniform_int_distribution<>(0, size-1)(engine);
        const int d = degreeList[i];
        const int dp = poisson(engine);
        degreeList[i] = dp;
        total_degree += dp - d;
    }

    config_model nw(degreeList, engine); 
    std::size_t nw_degree = 0;
    for (node_t i = 0; i < size; i++)
        nw_degree += nw.adjacencylist[i].size();
    REQUIRE(std::abs((long)total_degree - (long)nw_degree) < 0.00001);
}