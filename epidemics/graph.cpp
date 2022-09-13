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


/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-------------- NETWORK: ADJACENCY LIST -------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

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
    if (neighbours.back() != incomplete_neighbours)
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
        for (size_t j=0; j<degreelist[i]; j++)
            stubs.push_back(i);
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
//----------------------------------------------------
//----------------------------------------------------
//-------- NETWORK: CM CORRELATED------ --------------
//----------------------------------------------------
//----------------------------------------------------

config_model_correlated::config_model_correlated(std::vector<int> degreelist, rng_t& engine,bool assortative){

    //Step 0: Generate uncorrelated network


    config_model nk(degreelist,engine);

    adjacencylist = nk.adjacencylist;

    //Step 1: Link Selection
    // Choose at random two links
    // label the 4 nodes a,b,c,d s.t:
    // deg(a)<= deg(b)<= deg(c)<= deg(d)

    std::vector<edge_t> vec_edges({});

    vec_edges.reserve(nk.edges.size());
    for (auto it = nk.edges.begin(); it != nk.edges.end(); ) {
        vec_edges.push_back(std::move(nk.edges.extract(it++).value()));
    }

    std::uniform_int_distribution<> dist(0,(int) vec_edges.size()-1);
    
    for (int iter=0; iter < 100; iter++){
        std::vector<std::pair<int,node_t>> selected_nodes({});
        
        for (int i = 0; i < 2; i++){
            const edge_t e = vec_edges[dist(engine)];
            const std::pair<int,node_t> pair_first ={adjacencylist[e.first].size(),e.first};
            const std::pair<int,node_t> pair_second ={adjacencylist[e.second].size(),e.second};

            // break the existing links
            // https://stackoverflow.com/questions/3385229/c-erase-vector-element-by-value-rather-than-by-position
            adjacencylist[e.first].erase(std::remove(adjacencylist[e.first].begin(), adjacencylist[e.first].end(), e.second), adjacencylist[e.first].end());
            adjacencylist[e.second].erase(std::remove(adjacencylist[e.second].begin(), adjacencylist[e.second].end(), e.first), adjacencylist[e.second].end());

            selected_nodes.push_back(pair_first);
            selected_nodes.push_back(pair_second);
        }
        
        std::sort(selected_nodes.begin(), selected_nodes.end());
        
        node_t a = selected_nodes[0].second;
        node_t b = selected_nodes[1].second;
        node_t c = selected_nodes[2].second;
        node_t d = selected_nodes[3].second;

        //Step 3: Rewiring

        if (assortative==true)
        {
            // we check whether the particular rewiring leads to multi-links.
            // If it does, we reject it, returning to Step 1.
            // const edge_t e = (a < b) ? edge_t(a, b) : edge_t(b, a);
            // const bool is_in = edges.find(e) != edges.end();
            // if (is_in)
            //     continue;

            adjacencylist.at(a).push_back(b);
            adjacencylist.at(b).push_back(a);

            const edge_t e1 = (a < b) ? edge_t(a, b) : edge_t(b, a);
            bool is_in = nk.edges.find(e1) != nk.edges.end();
            if (is_in)
                std::cout << 1 << "\n";

            adjacencylist.at(c).push_back(d);
            adjacencylist.at(d).push_back(c);

            const edge_t e2 = (c < d) ? edge_t(c, d) : edge_t(d, c);
            adjacencylist.at(c).push_back(b);
            adjacencylist.at(b).push_back(c);
        }
    }
}


// }
// }           is_in = nk.edges.find(e2) != nk.edges.end();
//             if (is_in)
//                 std::cout << 1 << "\n";
//         } else {
//             adjacencylist.at(a).push_back(d);
//             adjacencylist.at(d).push_back(a);

//  


//--------------------------------------
//--------SCALE FREE NETWORK------------
//--------------------------------------

scale_free::scale_free(int size, rng_t& engine){
    
    adjacencylist.resize(size);
    
    std::vector<node_t> M(2*size,0);
    
    for (node_t v=0; v<size; v++) {
        //node_t node = mixed_nodes[v];
        M[2*v]=v;
        std::uniform_int_distribution<> dist(0,2*v);
        int r = dist(engine);
        M[2*v+1]=M[r];
    }
    
    for (int i=0; i<size; i++) {
        adjacencylist[M[2*i]].push_back(M[2*i+1]);
        adjacencylist[M[2*i+1]].push_back(M[2*i]);
    }

        
}


////OLD VERSION, CORRECT AND TESTED BUT TOO SLOW: 
//scale_free::scale_free(int size, rng_t& engine){
//
//    //Initialisation: Create the first node with no links.
//
//    //To avoid biases in the nodes when doing the next reaction scheme, we shuffle the nodes
//    // so that there is no correlation in the labels.
//    std::vector<node_t> mixed_nodes;
//    for( int i = 0; i < size; i++ )
//       mixed_nodes.push_back(i);
//    std::shuffle (mixed_nodes.begin(), mixed_nodes.end(), engine);
//
//    adjacencylist.resize(size);
//
////    adjacencylist[mixed_nodes[0]].push_back(mixed_nodes[1]);
////    adjacencylist[mixed_nodes[1]].push_back(mixed_nodes[0]);
//
//
//
//    // i is the current number of nodes
//    for (int i=1; i<size; i++) {
//
//        const node_t new_node = mixed_nodes[i];
//
//
//        // find who is connected to the new node by generating a rand. num. between 0 and i-1
//        // where the neighbour will be mixed_node[ rand. num.]
//        // the probability weights are the current degree (Bara-Alb model)
//
//        std::vector<double> prob(i);
//
//        for (int j = 0; j<i; j++) {
//            const auto degree = adjacencylist[mixed_nodes[j]].size();
//            prob[j]=degree;
//        }
//
//        std::discrete_distribution<int> distr(prob.begin(),prob.end());
//
//        const int index = distr(engine);
//
//        const node_t neighbour = mixed_nodes[index];
//
//        adjacencylist[new_node].push_back(neighbour);
//        adjacencylist[neighbour].push_back(new_node);
//    }
//
//
//};



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

    std::ifstream file(path_to_file);

    int size = file_size(path_to_file);
    
    adjacencylist.resize(size);

    std::string line, value;
    int i = 0;
    while (std::getline(file,line)) {

        size_t start;
        size_t end = 0;
        while ((start = line.find_first_not_of(",", end)) != std::string::npos) {
            end = line.find(",", start);
            int value = stoi(line.substr(start, end - start));
            adjacencylist[i].push_back(value);
            // std::cout << value << ",";
        }
        // std::cout << "\n";
        i++;
    }
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

//------------------------------------------------
//-----Measure degree correlation in a network----
//------------------------------------------------

std::vector<double> knn(graph_adjacencylist& nk){

    std::vector<std::vector<int>> nn_degree({});
    nn_degree.resize(nk.adjacencylist.size()); // k_max <= size of network
    int k_max = 0;

    for (node_t node = 0; node < nk.adjacencylist.size(); node++)
    {
        const int k = nk.adjacencylist[node].size();
        k_max = std::max(k,k_max);

        for (node_t neigh = 0; neigh < k; neigh++)
        {
            const int neigh_k = nk.adjacencylist[neigh].size();
            nn_degree[k].push_back(neigh_k);
        }
        
    }

    std::vector<double> result({});
    // compute average degree of the vertices connected to a vertex of degree k
    for (int k = 0; k <= k_max; k++)
    {
        if (nn_degree[k].size() == 0)
        {
            result.push_back(0);
            continue;
        }
        double average = 0;
        for (int i = 0; i < nn_degree[k].size(); i++)
        {
            average += nn_degree[k][i];
        }
        result.push_back(average /= nn_degree[k].size());

    }

    return result;


}

double assortativity(graph_adjacencylist& nk)
{
    double num = 0.0;
    double den = 0.0;

    // number of edges
    int m = 0;
    for (std::vector<node_t> neighbours : nk.adjacencylist)
    {
        m += neighbours.size();
    }

  
    for (int i = 0; i < nk.adjacencylist.size(); i++)
    {
        for (int j = 0; j < i; j++)
        {
            const int Aij = std::find(nk.adjacencylist[i].begin(), nk.adjacencylist[i].end(), j) != nk.adjacencylist[i].end();
            num += Aij - nk.outdegree(i)*nk.outdegree(j)/ (2 * m);
            den += nk.outdegree(i)*Aij - nk.outdegree(i)*nk.outdegree(j)/ (2 * m);  
        }
    }
    return num/den;
}