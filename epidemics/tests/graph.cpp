#include "tests/stdafx.h"

#include "graph.h"
#include "NextReaction.h"

/**
 * @brief Test case to verify `barabasi_albert`
 * 1. The size of the giant component should be the size of the network.
 * 2. The label of the node should be uncorrelated to its degree 
 * (usually the first nodes are the ones with the largest degree).
 * 3. The number of edges is +1 at every step, thus 2 * (SIZE -1) in total.
 * 4. The assortativity should be negative.
 * 
 * We repeat the above tests for various values of m.
 */
TEST_CASE("barabasi_albert", "[graph]") {
    std::mt19937 engine;
    const int SIZE = 1e4;
    const bool SHUFFLE = false;
    const double MEAN = 1;
    const double VARIANCE = 1;

    for (int m = 1; m<4 ;m++){
        barabasi_albert nw(SIZE,engine,m);
        transmission_time_gamma psi(MEAN,VARIANCE);
		simulate_next_reaction::params p;
		p.shuffle_neighbours = SHUFFLE;
		simulate_next_reaction simulation(nw, psi, nullptr, p);

        simulation.add_infections({ std::make_pair(0, 0)});
        int gc = 0;

        while(true){
            auto point = simulation.step(engine);
            if (!point )
                break;
            gc++;
        }
        REQUIRE(gc == SIZE);


        double mean_label = (SIZE - 1)/2;
        double mean_degree = 2*m;
        double emp_mean = 0;

        double numerator = 0.0;
        double den_x = 0.0;
        double den_y = 0.0;

        int number_edges = 0;

        for (node_t node = 0; node < SIZE; ++node) {
            const int k = nw.outdegree(node);
            number_edges += k;
            emp_mean += (double) k / SIZE;

            numerator += (node - mean_label) * (k - mean_degree);
            den_x += (node - mean_label)*(node - mean_label);
            den_y += (k - mean_degree) * (k - mean_degree);
        }

        double correlation_coefficient = numerator / (std::sqrt(den_x) * std::sqrt(den_y));
        
        REQUIRE(abs(correlation_coefficient) < 0.02);

        REQUIRE(number_edges ==2*m + 2* m * (SIZE-m-1));           
        REQUIRE(abs(emp_mean - mean_degree)/mean_degree < 0.01);
        REQUIRE(assortativity(nw) < -0.01);
    }
}


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
TEST_CASE("Acyclic networks with modified R0","[graph]") {
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
TEST_CASE("Acyclic networks without root degree reduction", "[graph]") {
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
TEST_CASE("Acyclic networks with root degree reduction", "[graph]") {
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

/**
 * @brief Test case to verify `cm::cm`
 *
 * The variance and mean degree of the resulting network should match
 * with the initial degree list that was the input.
 *
 */
TEST_CASE("Clustered configuration model networks (Serrano)","[graph]") {
	std::mt19937 engine;
	const int size = 100;

	// Generate a Poisson graph with the configuration model
	std::poisson_distribution<> poisson(3);
	std::vector<int> degreeList(size,0);
	std::size_t total_degree = 0;
	std::size_t max_degree = 0;
	for (int i = 0; i < size; i++)
	{
		const int k = poisson(engine);
		degreeList[i] = k;
		total_degree += k;
		max_degree = std::max(max_degree, (std::size_t)k);
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
	
	config_model_clustered_serrano nw(degreeList, 1.0, 1.0, engine);
	std::size_t nw_degree = 0;
	for (node_t i = 0; i < size; i++)
		nw_degree += nw.adjacencylist[i].size();
	REQUIRE(std::abs((long)total_degree - (long)nw_degree) < 0.00001);
}


/**
 * @brief Test case to verify `watts_strogatz::watts_strogatz`
 *
 * The number of edges should be preserved and should be equal to n
 *
 */
TEST_CASE("Watts-Strogatz model","[graph]") {
    std::mt19937 engine;
    const int size = 1000;
    const double p = 0.2;
    
    watts_strogatz network(size, p,engine);
    
    int number_of_edges = 0;
    for (std::vector<int> nei_list : network.adjacencylist){
        number_of_edges += (int) nei_list.size();
    }
    REQUIRE(std::abs(number_of_edges - 2* size ) < 0.01);
}

/**
 * @brief Test case to verify `cubic_lattice`
 *
 */
TEST_CASE("Cubic lattice","[graph]") {
    std::mt19937 engine;

    {
        cubic_lattice_2d nw(2);
        REQUIRE(nw.nodes() == 4);
        const node_t n11 = nw.node({0, 0});
        const node_t n12 = nw.node({0, 1});
        const node_t n21 = nw.node({1, 0});
        const node_t n22 = nw.node({1, 1});
        REQUIRE(nw.outdegree(n11)==2);
        REQUIRE(nw.outdegree(n12)==2);
        REQUIRE(nw.outdegree(n21)==2);
        REQUIRE(nw.outdegree(n22)==2);
        REQUIRE(nw.neighbour(n11, 0) == n21);
        REQUIRE(nw.neighbour(n11, 1) == n12);
		REQUIRE(nw.neighbour(n12, 0) == n22);
		REQUIRE(nw.neighbour(n12, 1) == n11);
		REQUIRE(nw.neighbour(n21, 0) == n11);
		REQUIRE(nw.neighbour(n21, 1) == n22);
		REQUIRE(nw.neighbour(n22, 0) == n12);
		REQUIRE(nw.neighbour(n22, 1) == n21);
    }

    {
        cubic_lattice_2d nw(3);
        REQUIRE(nw.nodes() == 9);
        const node_t n11 = nw.node({-1, -1});
        const node_t n12 = nw.node({-1,  0});
        const node_t n13 = nw.node({-1,  1});
        const node_t n21 = nw.node({ 0, -1});
        const node_t n22 = nw.node({ 0,  0});
        const node_t n23 = nw.node({ 0,  1});
        const node_t n31 = nw.node({ 1, -1});
        const node_t n32 = nw.node({ 1,  0});
        const node_t n33 = nw.node({ 1,  1});
        REQUIRE(nw.outdegree(n11)==2);
        REQUIRE(nw.outdegree(n13)==2);
        REQUIRE(nw.outdegree(n31)==2);
        REQUIRE(nw.outdegree(n33)==2);
        REQUIRE(nw.outdegree(n12)==3);
        REQUIRE(nw.outdegree(n21)==3);
        REQUIRE(nw.outdegree(n23)==3);
        REQUIRE(nw.outdegree(n32)==3);
        REQUIRE(nw.outdegree(n22)==4);
		for(int i = -1; i < 1; ++i) {
			for(int j = -1; j < 1; ++j) {
				const node_t n = nw.node({i, j});
				std::vector<node_t> nn;
				if (i == -1)
					nn.push_back(nw.node({i+1, j}));
				else if ((i > -1) && (i < 1)) {
					nn.push_back(nw.node({i+1, j}));
					nn.push_back(nw.node({i-1, j}));
				}
				else if (i == +1)
					nn.push_back(nw.node({i-1, j}));
				if (j == -1)
					nn.push_back(nw.node({i, j+1}));
				else if ((j > -1) && (j < 1)) {
					nn.push_back(nw.node({i, j+1}));
					nn.push_back(nw.node({i, j-1}));
				}
				else if (j == +1)
					nn.push_back(nw.node({i, j+1}));
				for(index_t k=0; k < nn.size(); ++k)
					REQUIRE(nw.neighbour(n, k) == nn[k]);
			}
		}
    }

    cubic_lattice_2d nw;
    const node_t n11 = nw.node({nw.coordinate_min, nw.coordinate_min});
    REQUIRE(nw.outdegree(n11)==2);
    const node_t n1k = nw.node({nw.coordinate_min, nw.coordinate_max});
    REQUIRE(nw.outdegree(n1k)==2);
    const node_t nk1 = nw.node({nw.coordinate_max, nw.coordinate_min});
    REQUIRE(nw.outdegree(nk1)==2);
    const node_t nkk = nw.node({nw.coordinate_max, nw.coordinate_max});
    REQUIRE(nw.outdegree(nkk)==2);
    const node_t ncc = nw.node({0, 0});
    REQUIRE(nw.outdegree(ncc)==4);
    REQUIRE(nw.neighbour(ncc, 0)==nw.node({  1,  0}));
    REQUIRE(nw.neighbour(ncc, 1)==nw.node({ -1,  0}));
    REQUIRE(nw.neighbour(ncc, 2)==nw.node({  0, +1}));
    REQUIRE(nw.neighbour(ncc, 3)==nw.node({  0, -1}));
}


/**
 * @brief Test case to verify `knn`
 *
 * When the network is uncorrelated knn(k) should be
 * independent of k. In the case of a ER network,
 * knn(k) should be approx constant and equal to <k^2>/<k>
 * as N-> infinity.
 */
TEST_CASE("Measuring average neighbour degree","[graph]") {
    std::mt19937 engine;
    int size = 1000000;
    double avg_degree = 3;

    erdos_reyni network(size, avg_degree,engine);
    

    double k1 = 0;
    double k2 = 0;
    for (node_t node = 0; node < size; node++){
        k1 += network.outdegree(node);
        k2 += std::pow(network.outdegree(node),2);
    }
    std::vector<double> vec = knn(network);
    double mu = k2/k1;
    
    int mu_ER = 4;

    REQUIRE(std::abs( mu-mu_ER )/mu_ER < 0.05);
    REQUIRE(std::abs( knn(network)[3]-mu )/mu < 0.05);
    REQUIRE(std::abs( knn(network)[5]-mu )/mu < 0.05);
    REQUIRE(std::abs( knn(network)[10]-mu)/mu < 0.05);

    
}


/**
 * @brief Test case to verify `assortativity`
 *
 * When the network is uncorrelated the assortativity should be
 * close to 0 as N-> infty. In the case of a ER network,
 * the assortativity should tend to 0;
 */
TEST_CASE("Measuring degree correlation","[graph]") {
    std::mt19937 engine;
    int size = 100000;
    erdos_reyni nk(size,3.0,engine);

    double r = assortativity(nk);
    REQUIRE(std::abs(r) < 0.05);


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

    r = assortativity(nw);
    REQUIRE(std::abs(r) < 0.05);
}

