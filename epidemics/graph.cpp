//
//  graph.cpp
//  epidemics
//
//  Created by Cure Samuel Cyrus on 2022/03/07.
//

#include "stdafx.h"
#include "types.h"
#include "graph.h"
#include "utility.h"


using boost::math::erfc;


/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-------------- NETWORK: BASE CLASS -----------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

node_t graph::nodes() {
    return -1;
}

graph::~graph()
{}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-------------- NETWORK: ADJACENCY LIST -------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

node_t graph_adjacencylist::nodes() {
    return (node_t)adjacencylist.size();
}    

node_t graph_adjacencylist::neighbour(node_t node, int neighbour_index) {
    const auto& n = adjacencylist.at(node);
    if ((neighbour_index < 0) || (n.size() <= (unsigned int)neighbour_index))
        return -1;
    return n[neighbour_index];
}

index_t graph_adjacencylist::outdegree(node_t node) {
    return (index_t) adjacencylist.at(node).size();
}


/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-------------- NETWORK: WATTS-STROGATZ -------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/



watts_strogatz::watts_strogatz(node_t size, int k, double p, rng_t& engine) {
    if (k <= 0)
        throw std::range_error("k must be positive for Watts-Strogatz networks");
    if (k % 2 != 0)
        throw std::range_error("k must be even for Watts-Strogatz networks");
    if (k > size)
        throw std::range_error("k cannot exceed size for Watts-Strogatz networks");

    /* First, create circular 1D lattice. For nodes labelled 0,...,n-1, each node is connected
     * to k/2 neighbours on each side, i.e. i to i-k/2,...,i-1,i+1,...,i+k/2/. We
     * actually insert a self-loop into the neighbour set here as well, that will avoid
     * checking for self-loops later, and we'll ignore them we constructing the actual
     * adjacency list.
     */
    std::vector<integer_set<node_t>> nodes_neighbours;
    nodes_neighbours.reserve(size);
    for(node_t i=0; i < size; ++i) {
        nodes_neighbours.emplace_back(0, size-1);
        for (node_t j=i-k/2; j <= i+k/2; ++j)
            nodes_neighbours[i].insert((size + j) % size);
    }

    /* Replace each edge u-v with probability p by a new edge u-w where w is sampled uniformly. */
	std::bernoulli_distribution rewire(p);
    for (node_t u = 0; u < size; u++) {
        /* Copy neighbours into a vector */
        auto& u_neighbours = nodes_neighbours[u];
        std::vector<node_t> vs;
        vs.reserve(u_neighbours.size());
		std::copy(u_neighbours.begin(), u_neighbours.end(), std::back_inserter(vs));

        /* Iterate over neighbours and rewire */
        for(node_t v: vs) {
            /* Rewire with probability p (and skip self-loops inserted above)
			 * Also skip if the node is already connected to every node
			 */
			if ((v == u) || ((node_t)u_neighbours.size() == size) || !rewire(engine))
                continue;

            /* Draw replacement w, delete u-v, add u-w.
             * The self-loops ensure that we don't draw w == u
             */
            const node_t w = u_neighbours.draw_complement(engine);
            const bool s1 = u_neighbours.erase(v); assert(s1); _unused(s1);
            const bool s2 = nodes_neighbours[v].erase(u); assert(s2); _unused(s2);
            const bool s3 = u_neighbours.insert(w).second; assert(s3); _unused(s3);
            const bool s4 = nodes_neighbours[w].insert(u).second; assert(s4); _unused(s4);
        }

        /* Remove self-loop */
        u_neighbours.erase(u);
    }
    
    /* Finally, initialise adjacency list */
    adjacencylist.resize(size);
	for (node_t i = 0; i < size; i++) {
		adjacencylist[i].reserve(nodes_neighbours[i].size());
		std::copy(nodes_neighbours[i].begin(), nodes_neighbours[i].end(),
				  std::back_inserter(adjacencylist[i]));
	}
}


/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-------------- NETWORK: ERDÖS-REYNI ----------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

erdos_reyni::erdos_reyni(int size, double avg_degree, rng_t& engine){
    /*--------------Initialisation--------------

     Construct erdos-Reyni graph and for each link we add an infection time:
     => A[i][j] is the time node i (or j) takes to infect node j (or i) once it has itself been infected by someone else.
     
     */
    const double p = avg_degree / (size - 1); // probability of an edge: if size ->infty and degree-> fixed then we get Poisson Graph.
    
    /* Using geometric distribution */
    std::geometric_distribution<> skip_edge(p); // comment: equals 0 with prob. p

    adjacencylist.resize(size);
    for (int i=0; i < size; i++) {
        int s = skip_edge(engine);
        // for very low probabilities p, skip_edge can return very large numbers
        // and so we have to be carefull not to overflow. Hence the slightly weird
        // coding of this loop; we test that we don't skip past the end of the nodes
        // before we attempt to update j, that way j cannot overflow.
        for(int j = -1; (s >= 0) && (s < (i - j - 1)); s = skip_edge(engine)) {
            j += 1 + s;
            assert(j < i);
            adjacencylist[i].push_back(j);
            adjacencylist[j].push_back(i);
        }
    }
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-------------- NETWORK: FULLY CONNECTED ------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

fully_connected::fully_connected(int size, rng_t& engine_)
    :engine(engine_)
    ,neighbours(size, std::vector<node_t>({}))
{}

node_t fully_connected::nodes() {
    return (node_t)neighbours.size();
}    

node_t fully_connected::neighbour(node_t node, int neighbour_index) {
    // index has to be in the range [0, size - 1]
    if ((neighbour_index < 0) || ((std::size_t)neighbour_index >= neighbours.size() - 1))
        return -1;
    // get neighbours of node <node>
    const auto& n = neighbours.at(node);
    if (!n.empty())
        return n.at(neighbour_index);
    // neighbours not yet generated, generate it now
    std::vector<node_t> np;
    np.reserve(neighbours.size() - 1);
    // add nodes [0,...,node-1,node+1,...,size] as neighbours
    for(int i=0; i < node; ++i)
        np.push_back(i);
    for(int i=node+1, m=(int)neighbours.size(); i < m; ++i)
        np.push_back(i);
    // shuffle
    std::shuffle(np.begin(), np.end(), engine);
    // store and return
    neighbours[node] = np;        
    return np.at(neighbour_index);
}

index_t fully_connected::outdegree(node_t node) {
    // network is fully connected, every node is every node's neighbour
    return (int)neighbours.size() - 1;
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-------------- NETWORK: ACYCLIC --------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

const node_t acyclic::incomplete_neighbours;

double acyclic::lambda(double mean, int digits) {
    /*
     * A Poisson distribution with parameter lambda conditioned
     * on strictly positive values has mean
     *   lambda / (1 - exp(-lambda)).
     * We invert this expression numerically to find lambda
     * for a given mean.
     */
    using namespace boost::math::tools;
    return newton_raphson_iterate([mean](double l) {
        return std::make_tuple(mean - l - mean * exp(-l),
                               -1 + mean * exp(-l));
    }, mean, 0.0, mean, digits);
}

acyclic::acyclic(double avg_degree, bool reduced_root_degree_, rng_t& engine_)
    :engine(engine_)
    ,degree_distribution(lambda(avg_degree, 10))
    ,reduced_root_degree(reduced_root_degree_)
    ,adjacencylist({ std::vector<node_t>{ incomplete_neighbours } })
{}

void acyclic::generate_neighbours(node_t node) {
    // check if the node is known to us
    if ((node < 0) || ((std::size_t)node >= adjacencylist.size()))
        throw std::range_error(std::string("invalid node ") + std::to_string(node));
    // get its adjacencies, possibly incomplete
    std::vector<node_t>& neighbours = adjacencylist[node];
    if (neighbours.empty() || (neighbours.back() != incomplete_neighbours))
        return;

    // neighbours not yet determined, so generate them

    // draw number of neighbours
    // we draw from a Poisson distribution conditioned on k > 0
    unsigned int k = 0;
    while (k == 0)
        k = degree_distribution(engine);
    // if reduced_root_degree is set, we treat the root node as
    // if it had an additional edge connecting it to the outside world,
    // and thus reduce it's degree here by one
    if ((node == 0) && reduced_root_degree)
        k -= 1;

    // remove the incomplete marker from the node's adjacency list
    neighbours.pop_back();

    // add neighbours
    // note that all nodes except 0 already have one neighbour when we get here
    while (neighbours.size() < k) {
        // get the first unused node index
        if (adjacencylist.size() > (std::size_t)std::numeric_limits<node_t>::max())
			throw std::range_error("maximum number of nodes exceeded");
        const node_t n = (node_t)adjacencylist.size();
        // add an entry marked incomplete for the new node to the adjacencylist
        std::vector<node_t> n_entry = { node, incomplete_neighbours };
        adjacencylist.push_back(n_entry);
        // add neighbour to node's adjacency list
        neighbours.push_back(n);
    }

    // we dont have to shuffle all neighbours, but we *do* have to ensure
    // that the neighbour that connects us to the root is at a random
    // position. Otherwise, traversals from the root will always find that
    // the neighbour at position 0 has already been traversed, while all
    // others haven't, which e.g. leads to bias towards later infection times
    // running an epidemic process on the networks that starts from node 0.
    if ((node > 0) && !neighbours.empty()) {
        auto idx_dist = std::uniform_int_distribution<std::size_t>(0, neighbours.size()-1);
        std::swap(neighbours[0], neighbours[idx_dist(engine)]);
    }
}

node_t acyclic::neighbour(node_t node, int neighbour_index) {
    generate_neighbours(node);
    const std::vector<node_t>& neighbours = adjacencylist[node];
    // index has to be in the range [0, size - 1]
    if ((neighbour_index < 0) || ((std::size_t)neighbour_index >= neighbours.size()))
        return -1;
    return neighbours[neighbour_index];
}

index_t acyclic::outdegree(node_t node) {
    generate_neighbours(node);
    const std::vector<node_t>& neighbours = adjacencylist[node];
    return (index_t)neighbours.size();
}


//----------------------------------------------------
//----------------------------------------------------
//-------- NETWORK: CONFIGURATION MODEL --------------
//----------------------------------------------------
//----------------------------------------------------

config_model::config_model(std::vector<int> degreelist, rng_t& engine){

    // Verify that all half edges can be paired
    const std::size_t size = degreelist.size();
    const std::size_t total_degree = std::accumulate(degreelist.begin(), degreelist.end(), 0);
    if (total_degree % 2 != 0)
        throw std::range_error("Total degree is odd, all edges cannot be paired."); 
    
    // Generate the list of all stubs (i.e. half edge without end point e.g. i--?)
    std::vector<node_t> stubs;
    stubs.reserve(total_degree);
    for (size_t i=0; i < size ; i++){
        for (int j=0; j < degreelist[i]; j++)
            stubs.push_back((int) i);
    }

    // Shuffle the list of stubs in order to sample without replacement
    std::shuffle(stubs.begin(), stubs.end(), engine);


    // set to keep track of multi-edges
    std::unordered_set<edge_t, pair_hash> edges = {};

    // Add the selected edges to the adjacency list.
    adjacencylist.resize(size);
    for (std::size_t i = 0; i < total_degree - 1; i+=2)
    {
        const node_t a = stubs[i];
        const node_t b = stubs[i+1];
        
        adjacencylist.at(a).push_back(b);
        adjacencylist.at(b).push_back(a);

        // keep track of number of self loops
        if (a == b)
            ++selfloops;

        // keep track of number of multi-edges
        // make undirected edges unambiguous by storing them
        // as a pair with a <= b.
        const edge_t e = (a < b) ? edge_t(a, b) : edge_t(b, a);
        const bool is_in = edges.find(e) != edges.end();
        if (!is_in)
            edges.insert(e);
        else
            ++multiedges;
    }
}

//-------- LOGNORMAL CONFIGURATION MODEL --------------
std::vector<int> lognormal_degree_list(double mean, double variance, int size, rng_t& engine){

    const double MU =2 * log(mean) - 0.5 * log( pow(mean,2.0)+ variance );
    const double SIGMA = sqrt( log( 1 + variance/pow(mean,2.0)));

    std::lognormal_distribution<double> distribution(MU,SIGMA);

    std::vector<int> degreelist({});
    int total = 0;
    while((int) degreelist.size() < size){
        
        const int k = std::ceil( distribution(engine) ) ;
        degreelist.push_back(k);

        total+= k;
    }
    // Make sure that the degree sequence is even so that all edges get connected to 2 different nodes.
    std::uniform_int_distribution<> d(0,size-1);
    while (total % 2 == 1){
        const int index = d(engine);
        total -= degreelist[index];;
        const int k = distribution(engine);
        degreelist[index] = k;
        total += k;
    }

    return degreelist;

    // double g = [](int k) {return (0.5 * erfc( (MU - log[k+1]) / (sqrt(2) * SIGMA) ) )/ ( 1- 0.5*erfc(MU/(sqrt(2)*SIGMA))) };

    // std::uniform_real_distribution<> d(0.0,1.0);
    // double r = d(engine);
    // int k = 0;
    // while (r >= g(k))
    //     k++;
    
    
}


//-------- SCALE FREE CONFIGURATION MODEL --------------
std::vector<int> powerlaw_degree_list(double g, int size, rng_t& engine){

    const double x0 = 1;
    const double x1 = size;

    std::uniform_real_distribution<double> distribution(0.0,1.0);

    std::vector<int> degreelist(size,0);
    int total = 0;
    for (int i = 0; i < size; i++){
        const double u =  distribution(engine);
        const double y = pow(   (  pow(x1,(g+1)) - pow(x0,(g+1))  ) * u + pow(x0,(g+1)),1/(g+1));
        const int k = std::floor(y) ;
        degreelist[i] = k;
        total += k;
    }
    // Make sure that the degree sequence is even so that all edges get connected to 2 different nodes.
    std::uniform_int_distribution<> d(0,size-1);
    while (total % 2 == 1){
        const int index = d(engine);
        total -= degreelist[index];;
        const int k = distribution(engine);
        degreelist[index] = k;
        total += k;
    }

    return degreelist;
   
    
}

//----------------------------------------------------
//----------------------------------------------------
//-------- NETWORK: CLUSTERED CONFIGURATION MODEL ----
//----------------------------------------------------
//----------------------------------------------------

/**
 * @brief Generates clustered random networks from a degree distribution and triangle counts
 *
 * Implements the algorithm proposed by Serrano & Boguna, 2005.
 */
struct config_model_clustered_serrano_builder {
	/**
	 * @brief Stub index
	 *
	 * Nodes contains as many stubs as their degree. These stubs are indexed by
	 * stub objects, which contain the node index and the index of the stub within
	 * the node (ranging from 0 to k-1 where k is the node's degree). Stub indices
	 * can also be empty, in which case their fields are set to -1.
	 */
	struct stub {
		node_t node = -1;
		index_t index = -1;
		
		/**
		 * @brief True if the stub index refers to an actual stub
		 * @return true or false
		 */
		bool is_filled() const { return (node != -1) && (index != -1); }

		/**
		 * @brief True if the stub index does not refer to any stub
		 * @return true or false
		 */
		bool is_empty() const { return !is_filled(); }
		
		bool operator==(const stub& other) const {
			return (node == other.node) && (index == other.index);
		}
		bool operator!=(const stub& other) const {
			return (node != other.node) || (index != other.index);
		}
	};
	
	struct stub_hasher {
		std::size_t operator()(const stub& e) const {
			return hash_combine(0, e.node, e.index);
		}
	};

	/**
	 * @brief Largest degree in the network
	 */
	const int kmax;

	/**
	 * @brief Exponent used in the definition of P(k) when
	 * degree classes are drawn randomly.
	 */
	const double beta;
	
	/**
	 * Maximum number of retries when forming triangles before we give up
	 */
	std::size_t max_triangle_attempts = 10;

	/**
	 * Random number generator
	 */
	rng_t& rng;

	/**
	 * @brief Number of remaining triangles to form for each degree
	 */
	std::unordered_map<int, std::size_t> triangles_remaining;

	/**
	 * @brief Number of triangles per degree class
	 */
	std::unordered_map<int, std::size_t> triangles_formed;

	/**
	 * @brief Factor used to convert fractional weights to integers
	 * Required because dyndist::vector_distribution only supports integral weights
	 */
	const double weightfactor = 1000000000;
	dyndist::vector_distribution<unsigned long> pk;
	
	/**
	 * @brief Network representation as a list of (empty or filled) stubs per node.
	 *
	 * Filled stubs represents edges and contain the index of the second stub forming the edge.
	 */
	std::unordered_map<node_t, std::vector<stub>> node_stubs;

	/**
	 * @brief List of degree classes, i.e. for every degree k a list of nodes of that degree
	 */
	std::unordered_map<int, std::set<node_t>> idx_k_nodes;

	/**
	 * Number of attempts since we last succesfully formed a triangle
	 */
	std::size_t triangle_attempts = 0;
	
	/**
	 * Set to true during network constructions if we fail to construct enough trianglges
	 */
	bool triangles_unsatisfied = false;
	
	/**
	 * @brief The "eligible components" of Serrano & Boguna
	 *
	 * Whereas Serrano & Boguna define "components" as being either edges or stubs, we
	 * represent both cases as a stub index. "Edges" are stub indices that *refer* to a filled
	 * stub, i.e. whose corresponding entry in node_stubs contains the index of another stub
	 * to which the first stub is connected. "Stubs" in the sense of Serrano & Boguna are stub
	 * indices of empty stubs, i.e. where the corresponding node_stubs entry is empty and which
	 * thus not connected to another stub.
	 *
	 * For every degree class, a set of stubs of nodes in that class is kept which may be
	 * still selected during triangle building
	 */
	std::unordered_map<int, drawable_set<stub, stub_hasher>> idx_k_stubs_eligible;
	
	/**
	 * @brief Build a network using the clustered configuration model of Serrano & Boguna 2005
	 * @param degrees list of node degrees
	 * @param degree_triangles number of triangles for each degree class (0, 1, 2, ...)
	 * @param _beta Parameter that determines P(k) which defines how likely nodes of different degrees are chosen
	 * @param _rng Random number generator
	 */
	config_model_clustered_serrano_builder(std::vector<node_t> degrees, const std::vector<int>& degree_triangles, double _beta, rng_t& _rng)
		:kmax(degree_triangles.size() - 1)
		,beta(_beta)
		,rng(_rng)
		,pk(kmax + 1, 0)
	{
		// I. Create stubs
		node_t next_node = 0;
		for(const int k: degrees) {
			const node_t node = next_node++;
			// Insert into set of nodes of degree k
			idx_k_nodes[k].insert(node);
			// Create stubs for current node
			node_stubs.insert({node, std::vector<stub>(k, stub())});
			// Insert stubs into indexes
			auto& ec_k = idx_k_stubs_eligible[k];
			for(index_t i=0; i < (index_t)k; ++i) {
				stub s {.node=node, .index = i};
				ec_k.insert(s);
			}
		}
		
		// II. Setup distribution P(k)
		for(std::size_t i=0; i < degree_triangles.size(); ++i) {
			if ((i < 2) && (degree_triangles[i] > 0))
				throw std::runtime_error("nodes with degrees 0 and 1 cannot have overlapping triangles");
			if (degree_triangles[i] > 0)
				triangles_remaining.insert({i, degree_triangles[i]});
		}
		for(const auto& e: triangles_remaining) {
			const int k = e.first;
			const std::size_t t = e.second;
			if ((t > 0) && (idx_k_nodes[k].size() == 0))
				throw std::runtime_error("triangles requested for empty degree class k=" + std::to_string(k));
			update_pk(k);
		}
	}
	
	void build() {
		// I. Create triangles
		create_triangles_serrano();
		
		// II. Closure of the network
		// Collect all remaining unconnected stubs
		std::vector<stub> remaining_stubs;
		for(node_t n=0; n < (node_t)node_stubs.size(); ++n) {
			const auto& ns = node_stubs[n];
			for(index_t i=0; i < (index_t)ns.size(); ++i) {
				const stub s {.node=n, .index=i};
				if (is_unconnected(s))
					remaining_stubs.push_back(s);
			}
		}
		// Do as the vanilla configuration model does
		if (remaining_stubs.size() % 2 != 0)
			throw std::logic_error("uneven number of stubs remaining");
		std::shuffle(remaining_stubs.begin(), remaining_stubs.end(), rng);
		while (!remaining_stubs.empty()) {
			const stub s1 = remaining_stubs.back();
			remaining_stubs.pop_back();
			const stub s2 = remaining_stubs.back();
			remaining_stubs.pop_back();
			// TODO: We might add self-edges here.
			// Those may be counted by add_edge() as triangles. Should we avoid that? How?
			add_edge(s1, s2, true, true);
		}
	}
	
	void create_triangles_serrano() {
		// c.f. Section B. Triangle formation in Serrano & Boguna, 2005.
		while (pk.weight() > 0) {
			++triangle_attempts;
			if (triangle_attempts > max_triangle_attempts) {
				triangles_unsatisfied = true;
				return;
			}
			
			// (i). Chose degree class, a stub s1 from a node A with that degree,
			// and a second stub s2 from the same node A.
			stub s1 = draw_stub(0);
			if (s1.is_empty())
				continue;
			stub s2 = draw_sibling_stub(s1);
			const node_t na = s1.node;

			// (ii). The two stubs are edges
			if (is_edge(s1) && is_edge(s2)) {
				const stub& s3 = find_unconnected_stub(connected_stub(s1).node);
				const stub& s4 = find_unconnected_stub(connected_stub(s2).node);
				// If the edges lead to nodes with no free stubs, can't form a triangle
				if (s3.is_empty() || s4.is_empty())
					continue;
				add_edge(s3, s4);
			}

			// (iii). One stub is an edge, one is unconnected
			else if ((is_edge(s1) && is_unconnected(s2)) || (is_edge(s2) && is_unconnected(s1))) {
				// wlog make s1 the edge, s2 the unconnected stub
				if (is_unconnected(s1)) {
					using std::swap;
					swap(s1, s2);
				}
				// Draw stub of node B connected to A via s1
				const stub& s3 = draw_sibling_stub(connected_stub(s1));
				const node_t nb = s3.node;
				assert(nb != na);
				if (is_edge(s3)) {
					// Stub is an edge, connects node B to node C. Find free stub of node C
					const stub& s4 = find_unconnected_stub(connected_stub(s3).node);
					// If C has no free stubs, can't form a triangle
					if (s4.is_empty())
						continue;
					add_edge(s2, s4);
				} else {
					// Stub is unconnected, must draw a node C with at least two free stubs
					const node_t nc = draw_stub(2).node;
					if ((nc == -1) || (nc == na) || (nc == nb))
						continue;
					const stub& s4 = find_unconnected_stub(nc);
					add_edge(s2, s4);
					const stub& s5 = find_unconnected_stub(nc);
					add_edge(s3, s5);
				}
			}

			// (iv). The two stubs are unconnected
			else if (is_unconnected(s1) && is_unconnected(s2)) {
				// Pick a node B with at least one free stub
				stub s3 = draw_stub(1);
				if (s3.is_empty())
					continue;
				const node_t nb = s3.node;
				if (nb == na)
					continue;
				// Pick second stub of node B, and swaps stubs
				// such that either both are unconnected, or
				// the first (s3) is connected and the second
				// unconnected.
				stub s4;
				assert(node_stubs[s3.node].size() >= 2);
				if (is_edge(s3))
					s4 = find_unconnected_stub(nb);
				else {
					s4 = s3;
					s3 = draw_sibling_stub(s4);
				}
				assert(is_unconnected(s4));
				if (is_unconnected(s3)) {
					// The remaining stub of node B is unconnected
					// Pick a node C with a t least two free stubs
					const node_t nc = draw_stub(2).node;
					if ((nc == -1) || (nc == na) || (nc == nb))
						continue;
					add_edge(s1, s4); // Connect A and B
					const stub& s5 = find_unconnected_stub(nc);
					add_edge(s2, s5); // Connect A and C
					const stub& s6 = find_unconnected_stub(nc);
					add_edge(s3, s6); // Connect B and C
				} else {
					// One stub of node B is connected to a node C
					const node_t nc = connected_stub(s3).node;
					if ((nc == na) || (nc == nb))
						continue;
					const stub& s5 = find_unconnected_stub(nc);
					// If node C has no free stubs, can't form a triangle
					if (s5.is_empty())
						continue;
					add_edge(s1, s4); // Connect A and B
					add_edge(s2, s5) ; // Connect A and C
				}
			}

			// Above list of cases should be exhaustive
			else throw std::logic_error("unreachable code reached");
		}
	}

	/**
	 * @brief Adds an edge connecting stubs a and b
	 *
	 * Checks for any triangles that may have been formed and updates the count
	 *
	 * @param a first stub
	 * @param b second stub
	 */
	bool add_edge(stub a, stub b, bool allow_selfedge=false, bool allow_multiedge = false) {
		if (a == b)
			throw std::logic_error("add_edge() called with two identical stubs");
		const node_t& na = a.node;
		const std::size_t ka = node_stubs[na].size();
		const node_t& nb = b.node;
		const std::size_t kb = node_stubs[nb].size();
		stub& ca = connected_stub(a);
		stub& cb = connected_stub(b);
		if (ca.is_filled() || cb.is_filled())
			throw std::logic_error("add_edge() called for connected stubs");
		
		// I. Check whether the two nodes are already connected, if so
		// do nothing unless we were asked to create multi-edges
		auto& na_stubs = node_stubs[na];
		auto& nb_stubs = node_stubs[nb];
		if (!allow_multiedge) {
			if (na_stubs.size() < nb_stubs.size()) {
				for(const stub& s: na_stubs)
					if (s.node == nb)
						return false;
			} else {
				for(const stub& s: nb_stubs)
					if (s.node == na)
						return false;
			}
		}
		
		// II. If we weren't specifically told to allow self-edges, complain about them
		if (!allow_selfedge && (na == nb))
			throw std::logic_error("add_edge() called for a self-edge but allow_selfedge is false");
		
		// III. Enumerate triangles that will be completed
		// First, build a set of neighbours of a
		std::set<node_t> a_neighbours;
		for(const stub& s: na_stubs) {
			if (!s.is_filled())
				continue;
			
			a_neighbours.insert(s.node);
		}
		// Now scan neighbours of b for overlaps
		for(const stub& s: nb_stubs) {
			if (!s.is_filled())
				continue;
			const node_t nc = s.node;
			if (a_neighbours.find(nc) == a_neighbours.end())
				continue;
			// Found a common neighbour of the two nodes to be connected,
			// i.e. a triangle that will be completed.
			
			// a. Update triangle counts, and remove stubs of nodes in satisified degree classes
			const std::size_t kc = node_stubs[nc].size();
			on_overlapping_triangle(ka, na);
			on_overlapping_triangle(kb, nb);
			on_overlapping_triangle(kc, nc);
		}
		
		// IV. Connect stubs
		ca.node = b.node;
		ca.index = b.index;
		cb.node = a.node;
		cb.index = a.index;
		
		// TODO: We should probably remove edges/nodes from the EC that cannot form triangles anymore.
		// One reasonable conditions seems to be that either
		//   a. Both nodes of an edge have free stubs, or
		//   b. one node has a free stub and an *neighbour* of the other node has a free stub.
		// However, some amount of retries during section II (Create triangles) are probably
		// needed even with this in place. Remove these edges is therefore more of a performance
		// optimization than an integral part of the algorithm.
		
		return true;
	}
	
	/**
	 * @brief Callled whenever a triangle is created for every node overlapping the triangle
	 *
	 * @param k degree class of the node
	 * @param n node index
	 */
	void on_overlapping_triangle(int k, node_t n) {
		assert(node_stubs[n].size() == (std::size_t)k);
		
		// Count triangle
		triangle_attempts = 0;
		triangles_formed[k]++;
		
		// Reduce number of triangles still to be formed for class
		bool will_satisfy_class = (triangles_remaining[k] == 1);
		if (triangles_remaining[k] > 0) {
			--triangles_remaining[k];
			update_pk(k);
		}
		
		// Remove stubs of nodes in newly satisfied class from EC
		if (will_satisfy_class) {
			auto& idx = idx_k_stubs_eligible[k];
			// Scan EC for degree class k and remove opposing stubs of connected stubs
			for(const stub& s: idx) {
				const stub& sp = connected_stub(s);
				if (sp.is_empty())
					continue;
				// Find EC tier which contains the connected node and remove
				// the connected stub from it, unless it lives in the same
				// degree class. In that case, it's removed anyway when the
				// whole class is removed below.
				const int kp = node_stubs[sp.node].size();
				if (k != kp)
					idx_k_stubs_eligible[kp].erase(sp);
			}
			// Remove EC for degree class k
			idx_k_stubs_eligible.erase(k);
		}
	}
	
	/**
	 * @brief Returns the index of an unconnected stub of the given node, or an empty index if none exists
	 *
	 * @param node node to find an unconnected stub for
	 * @return index of the stub as a stub object
	 */
	stub find_unconnected_stub(node_t node) {
		auto& ns = node_stubs[node];
		for(index_t i=0; i < (index_t)ns.size(); ++i) {
			if (ns[i].is_empty())
				return stub {.node=node, .index=i};
		}
		return stub();
	}
	
	/**
	 * @brief Returns the index of the stub that the given stub is connected to.
	 *
	 * If the given stub is not an edge, return an empty index.
	 *
	 * @param s stub index
	 * @return index of the connected stub
	 */
	stub& connected_stub(stub s) {
		assert(!s.is_empty());
		return node_stubs[s.node][s.index];
	}
	
	/**
	 * @brief Returns true if the stub indexed by s is part of an edge
	 * @param s stub index
	 * @return true or false
	 */
	bool is_edge(stub s) {
		assert(!s.is_empty());
		return connected_stub(s).is_filled();
	}

	/**
	 * @brief Returns true if the stub indexed by s is unconnected
	 * @param s stub index
	 * @return true or false
	 */
	bool is_unconnected(stub s) {
		assert(!s.is_empty());
		return connected_stub(s).is_empty();
	}

	/**
	 * @brief Draw a degree class k according to P(k)
	 * @return node degree k
	 */
	int draw_k() {
		return pk(rng);
	}
	
	/**
	 * @brief Draw a random stub whose node has at least rmin free stubs
	 * @param rmin minimal number of free stubs of node
	 * @return stub index of selected stub
	 */
	stub draw_stub(int rmin) {
		const int k = draw_k();
		return draw_stub(k, rmin);
	}
	
	/**
	 * @brief Draw a random stub whose node has degree class k and at least rmin free stubs
	 * @param k degree class of node
	 * @param rmin minimal number of free stubs of node
	 * @return stub index of selected stub
	 */
	stub draw_stub(int k, int rmin) {
		if (rmin > k)
			throw std::logic_error("impossible number of free stubs requested");
		
		// I. Fast path
		// Draw random stub
		const stub s = *idx_k_stubs_eligible[k](rng);
		// Check whether the node has rmin free stubs
		auto& ns = node_stubs[s.node];
		int r = 0;
		for(index_t i=0; i < (index_t)ns.size(); ++i) {
			if (ns[i].is_filled())
				continue;
			// Count number of free stubs, return stub if we found enough
			++r;
			if (r >= rmin)
				return s;
		}
		
		// II. Slow path, used after failing in the fast path
		std::vector<stub> stubs;
		stubs.reserve(idx_k_stubs_eligible[k].size());
		for(const auto& s: idx_k_stubs_eligible[k]) {
			// Scan sibling stubs to see whether at least rmin are free
			auto& ns = node_stubs[s.node];
			int r = 0;
			for(index_t i=0; i < (index_t)ns.size(); ++i) {
				if (ns[i].is_filled())
					continue;
				// Count number of free stubs, return stub if we found enough
				++r;
				if (r < rmin)
					continue;
				// Found a matching stub
				stubs.push_back(s);
				break;
			}
		}
		if (stubs.empty()) {
			// No matching stubs exist, return an empty stub
			return stub();
		}
			
		// Draw random valid stub
		return stubs[std::uniform_int_distribution<std::size_t>(0, stubs.size() - 1)(rng)];
	}
	
	/**
	 * @brief Draw a random sibling of a stub, i.e. a stub of the same node
	 * @param s1 index of a stub
	 * @return index of the randomly chosen sibling stub
	 */
	stub draw_sibling_stub(stub s1) {
		// Get degree of node
		const int k = node_stubs[s1.node].size();
		if (k <= 1)
			throw std::logic_error("cannot draw siblings of node with degree one");
		// Draw an stub different from s1
		int i = s1.index;
		while (i == s1.index)
			i = std::uniform_int_distribution<int>(0, k-1)(rng);
		return stub {.node=s1.node, .index=i};
	}
	
	/**
	 * @brief Update degree class distribution P(k) after modifying the number of triangkes
	 * @param k degree class for which the number of triangles was modified
	 */
	void update_pk(int k) {
		assert(k >= 2);
		pk[k] = (unsigned long)(std::pow((double)(triangles_remaining[k]), beta) * weightfactor);
	}
};

std::vector<int> config_model_clustered_serrano::triangles_binomial
 (std::function<double(int)> ck, std::vector<int> degrees, rng_t& engine)
{
	std::vector<int> triangles(*std::max_element(degrees.begin(), degrees.end()) + 1, 0);
	for(std::size_t i=0; i < degrees.size(); ++i) {
		const int k = degrees[i];
		if (k < 2)
			continue;
		const double p = ck(k);
		const double M = k * (k-1) / 2;
		const int t = std::binomial_distribution<int>(M, p)(engine);
		triangles.at(k) += t;
	}
	return triangles;
}

config_model_clustered_serrano::config_model_clustered_serrano
 (std::vector<int> degrees, double alpha, double beta, rng_t& engine)
	:config_model_clustered_serrano(degrees, [alpha](int k) -> double { return 0.5*std::pow((double)(k-1), -alpha); }, beta, engine)
{};

config_model_clustered_serrano::config_model_clustered_serrano
 (std::vector<int> degrees, std::function<double(int)> ck, double beta, rng_t& engine)
	:config_model_clustered_serrano(degrees, triangles_binomial(ck, degrees, engine), beta, engine)
{}

config_model_clustered_serrano::config_model_clustered_serrano
 (std::vector<int> degreelist, std::vector<int> degreetriangles, double beta, rng_t& engine)
{
	config_model_clustered_serrano_builder b(degreelist, degreetriangles, beta, engine);
	b.build();
	triangles_unsatisfied = b.triangles_unsatisfied;
	/* Copy network structure from builder */
	for(std::size_t i = 0; i < b.node_stubs.size(); ++i) {
		const auto& ns = b.node_stubs[i];
		// Get adjacencylist for i-th node, make sure it exists
		adjacencylist.resize(i+1);
		auto& i_al = adjacencylist[i];
		for(const auto& sc: ns) {
			if (sc.is_empty())
				throw std::logic_error("network construction failed");
			// Add edge. Reverse edge will be added later, we rely on node_stubs being symmetric
			i_al.push_back(sc.node);
		}
	}
}


//--------------------------------------
//--------SCALE FREE NETWORK------------
//--------------------------------------

scale_free::scale_free(int size, rng_t& engine, int m){

    // To generate efficiently a BA network, we used the approached used in the python library networkx.
    // -> instead of re-initialising the distribution at every step by updating the weights, 
    // we sample at uniform random from a list, but each node is duplicated k times where k is their degree.
    // This will mimic a random sample over a weighted distribution.
    
    std::vector<node_t> repeated_nodes;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    // In the BA algorithm, there is a negative correlation between the label of the node and their degree.
    // In the NR scheme, we fetch the next neighbour in order of their label. 
    // Therefore we wish to make the degree of a node independent of their label.
    // To do so, we introduce 'mixed_nodes'.

    std::vector<node_t> mixed_nodes;
    for( int i = 0; i < size; i++ )
       mixed_nodes.push_back(i);
    std::shuffle(mixed_nodes.begin(), mixed_nodes.end(), engine);


    // Initialise the adjacency list by adding the first m+1 nodes. (no need to sample a random number for the first step.)
    // We start with m0 nodes, the links between which are chosen arbitrarily,
    // as long as each node has at least one link.
    // We form a chain to ensure that there the network does not have disconnected components.
    adjacencylist.resize(size);
    for (int i = 0; i < m;i ++){
        adjacencylist[mixed_nodes[i]].push_back(mixed_nodes[i+1]);
        adjacencylist[mixed_nodes[i+1]].push_back(mixed_nodes[i]);
    }
    // adjacencylist[0].push_back(1);
    // adjacencylist[1].push_back(0);

    // There are only m+1 nodes in the network, all with degree 2 except the first and last:
    repeated_nodes.push_back(mixed_nodes[0]);
    repeated_nodes.push_back(mixed_nodes[m]);

    for (int i = 1; i < m;i ++){
        repeated_nodes.push_back(mixed_nodes[i]);
        repeated_nodes.push_back(mixed_nodes[i]);
    }
    
    int len = (int) repeated_nodes.size();

    // i is the current number of nodes in the graph
    for (int i=m+1; i<size; i++) {

        // Select next node to be added to the network;
        const node_t new_node = mixed_nodes[i];
        // std::cout << "node : " << new_node << "\n";
        // const node_t new_node = i;

        std::unordered_set<index_t> chosen_nodes = {};
        for (int j = 0; j < m; j++){
            // Determine who it will be attached to.
            const double u = distribution(engine);
            const index_t index = floor(u * len);

            // Forbid multiple-edges
            if (chosen_nodes.find(index) != chosen_nodes.end()) {
                j--;
                continue;
            }

            chosen_nodes.insert(index);
            const node_t selected_node = repeated_nodes[index];

            adjacencylist[new_node].push_back(selected_node);
            adjacencylist[selected_node].push_back(new_node);

            // Update the weights of the preferential attachement
            repeated_nodes.push_back(new_node);
            repeated_nodes.push_back(selected_node);
            len += 2;
        }

    }   
}


//--------------------------------------
//--------IMPORTED NETWORK----------
//--------------------------------------

int imported_network::file_size(std::string path_to_file){
    
    std::ifstream file(path_to_file);
    std::string line;

    int size = 0;
    while (std::getline(file, line))
        size++;
    return size;
}

imported_network::imported_network(std::string path_to_file){
    std::ifstream file(path_to_file); // Open the CSV file
    std::vector<std::vector<int>> data; // Create a vector to store the data

    std::string line;
    while (std::getline(file, line)) { // Read each line of the CSV file
        std::vector<int> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) { // Parse each cell of the line
            row.push_back(std::stoi(cell)); // Convert the cell to an integer and add it to the row vector
        }
        adjacencylist.push_back(row); // Add the row vector to the data vector
    }

    // Print the data to verify that it was read correctly
    // for (const auto& row : adjacencylist) {
    //     for (const auto& cell : row) {
    //         std::cout << cell << " ";
    //     }
    //     std::cout << std::endl;
    // }
//     std::ifstream file(path_to_file);
//     int size = file_size(path_to_file);
//     std::cout << "file size :" << size << "\n";
    
//     std::string line, value;
//     int i = 0;
//     while (std::getline(file,line)) {
//         std::vector<node_t> neighbours({});
// //        adjacencylist.push_back( );
//         size_t start;
//         size_t end = 0;
//         while ((start = line.find_first_not_of(",", end)) != std::string::npos) {
//             end = line.find(",", start);
//             node_t value = stoi(line.substr(start, end - start));

//             neighbours.push_back(value);
//             // std::cout << value << ",";
//         }
//         adjacencylist.push_back(neighbours);
//         // std::cout << "\n";
//         i++;
//     }
}
        // std::string j;
        // while (std::getline(line, j[0],","))
        // {
        //     int num = std::stoi(j);
        //     adjacencylist[i].push_back(num );
        //     adjacencylist[num].push_back( i );

        //     std::cout << num << " ";
        // }
        // std::cout << "\n";
        
        // // read an entire row and
        // // store it in a string variable 'line'
        // std::getline(fin,line);
  
        // // used for breaking words
        // std::stringstream s(line);

        
        // // read every column data of a row and
        // // store it in a string variable, 'word'
        // while (std::getline(s, word, ', ')) {
  
        //     // add all the column data
        //     // of a row to a vector
        //     adjacencylist[i].push_back(std::stoi(word));

//---------------------------------------------------
//-----Measure edge multiplicity in a network--------
//---------------------------------------------------

std::vector<std::vector<double>> edge_multiplicity(graph_adjacencylist& nw){

    const int SIZE = (int) nw.adjacencylist.size();

    int kmax = 0;
    for (node_t node = 0; node < SIZE; node++){
        const int k0 = nw.outdegree(node);
        kmax = std::max(k0,kmax);
    }

    std::vector<double> T2(kmax+1,0);
    std::vector<std::vector<double>> T1(kmax+1,std::vector<double>(kmax+1,0));
    std::vector<std::vector<double>> e_kk(kmax+1,std::vector<double>(kmax+1,0));
    std::vector<std::vector<std::vector<double>>> T(kmax+1, std::vector<std::vector<double>>(kmax+1, std::vector<double>(kmax+1,0)));
    
    for (node_t node = 0; node < SIZE; node++){
        
        const int k0 = nw.outdegree(node);
        for (node_t neigh_1 : nw.adjacencylist[node])
        {
            const int k1 = nw.outdegree(neigh_1);
            e_kk[k0][k1] ++;

            for (node_t neigh_2 : nw.adjacencylist[neigh_1]){
                if (neigh_2 == node)
                    continue;
                const int k2 = nw.outdegree(neigh_2);
                node_t small_node = (k0 <= k2) ? node : neigh_2;
                node_t large_node = (k0 <= k2) ? neigh_2 : node;

                // verify if edge exists between node and neighbour n°2
                auto it = std::find(nw.adjacencylist[small_node].begin(), nw.adjacencylist[small_node].end(), large_node);
                if (it != nw.adjacencylist[small_node].end()){
                    T[k0][k1][k2] ++;
                    T1[k0][k1]++;
                    T2[k0]++;
                }
            }
        }
    }
    std::vector<std::vector<double>> mkk(kmax+1,std::vector<double>(kmax+1,0));

    for (int k=0; k <= kmax; k++){
        for (int q=0; q <= kmax; q++){
            if (e_kk[k][q]==0)
                continue;
            mkk[k][q] = T1[k][q] / ( e_kk[k][q] );
        }
    }
    return mkk;
}

//---------------------------------------------------
//-----ADD DEGREE CORRELATION TO THE NETWORK---------
//---------------------------------------------------

// helper function to verify if edge exists.
bool edge_exists(node_t a, node_t b, const graph_adjacencylist& nw){
    auto it = find(nw.adjacencylist[b].begin(), nw.adjacencylist[b].end(), a);
    if (it == nw.adjacencylist[b].end())
        return false;
    // else
    return true;
}

void add_correlation(double r,graph_adjacencylist& nw,rng_t& engine){
    
    
    // To sample an edge a uniform random without generating the entire list of edges
    // we can sample a node from a discrete distr with their degree as the prob. weights.
    std::vector<double> degree_list({});
    for (node_t node = 0; node < (node_t)nw.adjacencylist.size(); node ++)
        degree_list.push_back(nw.adjacencylist[node].size());
    
    std::discrete_distribution<int> sample_source_node(degree_list.begin(),degree_list.end());
    

//    int swap_count = 0;
    
//    if (abs(assortativity(nw)-r)/abs(r) < 0.1)
//        return;
    for (int iteration=0; iteration < 9000; iteration++) {
        
 

//        u--v            u  v        u v
//               becomes  |  |   or    X
//        x--y            x  y        x y
//
//
        
        step_1: // Label for goto. Bad practice but here it allows to leave the nested conditions very smoothly.
                        
        //STEP 1 : Link Selection
        // Choose at random two links
        // label the 4 nodes a,b,c,d s.t:
        // deg(a)>= deg(b)>= deg(c)>= deg(d)
        
    
        // Pick one node i randomly according to its degree, then pick randomly one of its neighbours j,
        // => (i,j) is a sampled edge uniformly at random among all edges.
        
        node_t u = sample_source_node(engine);
        node_t x = sample_source_node(engine);
        
        if (u == x)
            continue; // same source, skip
        
        // Choose at uniform among their neigbours
        
        std::uniform_int_distribution<int> neighbour_u(0,(int) nw.adjacencylist[u].size()-1);
        std::uniform_int_distribution<int> neighbour_x(0,(int) nw.adjacencylist[x].size()-1);
        
        node_t v = nw.adjacencylist[u][neighbour_u(engine)];
        node_t y = nw.adjacencylist[x][neighbour_x(engine)];
        
        if (v==y)
            continue; //same target, skip
        
        if (x==v && y==u)
            continue;
        
        std::pair<int,node_t> p_u = {nw.outdegree(u) , u };
        std::pair<int,node_t> p_x = {nw.outdegree(x) , x };
        std::pair<int,node_t> p_v = {nw.outdegree(v) , v };
        std::pair<int,node_t> p_y = {nw.outdegree(y) , y };
        
        std::vector<std::pair<int,node_t>> sampled_nodes({p_u , p_x , p_v , p_y});
        std::sort(sampled_nodes.begin(), sampled_nodes.end());

        node_t a = sampled_nodes[0].second; // Highest degree
        node_t b = sampled_nodes[1].second;
        node_t c = sampled_nodes[2].second;
        node_t d = sampled_nodes[3].second; // Lowest degree
 
        
        // STEP 2 : REWIRING
        
        if (r>=0) {
            // STEP 2A: Assortative
            // pair low (c & d) and pair high (a & b)
    
            //   c  d        c d             c--d
            //   |  |   or    X    becomes
            //   a  b        a b             a--b

            
            // if highs are already together then skip.
            if (  (a==u && b==v) || (a==x && b==y) || (a==v && b==u) || (a==y && b==x) )
                goto step_1;
            
            // Check if rewiring leads to multilinks: a-b and/or c-d already exists
            if ( edge_exists(a,b,nw) || edge_exists(c,d,nw) )
                goto step_1;
  
            // delete u-v edge
            std::vector<int>::iterator v_index = find (nw.adjacencylist[u].begin(),nw.adjacencylist[u].end(), v);
            nw.adjacencylist[u].erase(v_index);
            
            std::vector<int>::iterator u_index = find (nw.adjacencylist[v].begin(),nw.adjacencylist[v].end(), u);
            nw.adjacencylist[v].erase(u_index);
            
            // delete x-y edge
            std::vector<int>::iterator x_index = find (nw.adjacencylist[y].begin(),nw.adjacencylist[y].end(), x);
            nw.adjacencylist[y].erase(x_index);
            
            std::vector<int>::iterator y_index = find (nw.adjacencylist[x].begin(),nw.adjacencylist[x].end(), y);
            nw.adjacencylist[x].erase(y_index);
            
            
            // add edge a-b and add edge c-d
            nw.adjacencylist[a].push_back(b);
            nw.adjacencylist[b].push_back(a);
            
            nw.adjacencylist[c].push_back(d);
            nw.adjacencylist[d].push_back(c);
            

            
        } else {
            // STEP 2B: Disassortative
            // Pair the highest with the lowest:
            // (a & d) and (b & c)

            //   c--d        c  d             c d
            //         or    |  |   becomes    X
            //   a--b        a  b             a b

            //        u--v            u  v        u v
            //               becomes  |  |   or    X
            //        x--y            x  y        x y
            
            // if correct pairs are already formed then skip.
            if (  (a==u && d==v) || (a==x && d==y) || (c==v && b==u) || (c==y && b==x) )
                goto step_1;
            
            // Check if rewiring leads to multilinks: a-d and/or b-c already exists
            if ( edge_exists(a,d,nw) || edge_exists(b,c,nw) )
                goto step_1;

            // delete u-v edge
            auto v_index = find (nw.adjacencylist[u].begin(),nw.adjacencylist[u].end(), v);
            nw.adjacencylist[u].erase(v_index);
            
            auto u_index = find (nw.adjacencylist[v].begin(),nw.adjacencylist[v].end(), u);
            nw.adjacencylist[v].erase(u_index);
            
            // delete x-y edge
            auto x_index = find (nw.adjacencylist[y].begin(),nw.adjacencylist[y].end(), x);
            nw.adjacencylist[y].erase(x_index);
            
            auto y_index = find (nw.adjacencylist[x].begin(),nw.adjacencylist[x].end(), y);
            nw.adjacencylist[x].erase(y_index);
            
            
            // add edge a-d and add edge b-c
            nw.adjacencylist[a].push_back(d);
            nw.adjacencylist[d].push_back(a);
            
            nw.adjacencylist[c].push_back(b);
            nw.adjacencylist[b].push_back(c);
            
        }
        
        if (iteration % 1000 == 0)
        {
            double cur_r = assortativity(nw);
            std::cout << "cur r : " << cur_r << "\r";
            double condition = abs(cur_r-r)/abs(r);
            if (condition < 0.25 || ( r > 0 && cur_r > r) || (r < 0 && cur_r < r)){
                std::cout << "converged.\n";
                break;
            }
        }
    
    }


    
}

//------------------------------------------------
//-----Measure degree correlation in a network----
//------------------------------------------------

// Average neighbour's degree for a node of degree k
std::vector<double> knn(graph_adjacencylist& nw) {

    int size = (int) nw.adjacencylist.size();
    std::vector<double> knn_degree(size, 0);
    std::vector<double> nb_with_degree(size,0);

    int kmax = 0;
    for (node_t node = 0; node < size; node++){
        
        int k = nw.outdegree(node);
        
        kmax = std::max(k,kmax);
        for (node_t neigh : nw.adjacencylist[node])
        {
            const double k_neigh = (double) nw.outdegree(neigh);
            knn_degree[k] += k_neigh;
            nb_with_degree[k] += 1;
        }

    }

    while ((int) knn_degree.size() > kmax + 1){
        knn_degree.pop_back();
        nb_with_degree.pop_back();
    }
    
    for (int k=0; k < kmax+1; k++)
    {
        if (nb_with_degree[k]!=0){
            knn_degree[k] =  knn_degree[k]/nb_with_degree[k];
        }
        // else {
        //     knn_degree[k]=0;
        // }
    }
    return knn_degree;

}

//------------------------------------------------
//-----Measure degree correlation in a network----
//------------------------------------------------


// From the definition of the assortativity, after simple manipulations we can rewrite:
// r = (<k^2 knn(k)>/<k> - mu^2 ) / (<k^3>/<k> - mu^2)
// where mu^2 = <k^2>/<k>
// However <k^2 knn(k)> is tricky to measure, instead, we measure the fraction of links that connect deg i to deg j, W(i,j)
double assortativity(graph_adjacencylist& nw)
{
    
    int size = (int) nw.adjacencylist.size();
    std::vector<double> Knn = knn(nw);
    
    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double kkWkk = 0.0;

    for( int i = 0; i < size; i++ ){
        int k = nw.outdegree(i);
        k1 += (double) k / size;
        k2 += (double) k*k / size;
        k3 += (double) k*k*k / size;
    }
    auto wkk = Wkk(nw);

    for (std::size_t i = 1; i < wkk.size(); i++){
        for (std::size_t j = 1; j < wkk[i].size(); j++)
        {
            kkWkk += i*j* wkk[i][j];
        }
    }

    double mu = k2/k1;

    double sigma2 = k3/k1 - mu*mu;
    
    double r = (kkWkk - mu*mu) / sigma2;

    return r;

}


//Fraction of links that connect a node of deg k to a node of deg k'
std::vector<std::vector<double>> Wkk(graph_adjacencylist& nw){
    
    //Determine k_max (Figured that it was the most sane option and the most readable,
    // at the cost of going through all nodes once more.
    int k_max = 0;
    double total_edges = 0.0;
    for (node_t node = 0; node < (int) nw.adjacencylist.size(); node++){
        const int nb_edges = (int) nw.adjacencylist[node].size();
        k_max = std::max(k_max,nb_edges);
        total_edges += nb_edges;
    }
    //total_edges /= 2; //Avoid double counting
    
    std::vector<std::vector<double>> w(k_max + 1, std::vector<double> (k_max + 1, 0.0));

    for (node_t node = 0; node < (int) nw.adjacencylist.size(); node++) {
        
        const int k = (int) nw.adjacencylist[node].size();
        
        for (node_t neigh : nw.adjacencylist[node])
        {
            const int k_prime = (int) nw.adjacencylist.at(neigh).size();
            w.at(k ).at( k_prime) += 1.0 / total_edges;
        }
        
    }
    
    return w;
}
