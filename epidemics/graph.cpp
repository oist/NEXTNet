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
    std::cout << "file size :" << size << "\n";
    
    std::string line, value;
    int i = 0;
    while (std::getline(file,line)) {
        std::vector<node_t> neighbours({});
//        adjacencylist.push_back( );
        size_t start;
        size_t end = 0;
        while ((start = line.find_first_not_of(",", end)) != std::string::npos) {
            end = line.find(",", start);
            node_t value = stoi(line.substr(start, end - start));
            
//            adjacencylist[i].push_back(value);
            neighbours.push_back(value);
            // std::cout << value << ",";
        }
        adjacencylist.push_back(neighbours);
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
    for (int iteration=0; iteration < 90000; iteration++) {
        
 

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
            std::cout << "cur r : " << cur_r << "\n";
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

std::vector<double> knn(graph_adjacencylist& nw) {

    std::vector<std::vector<int>> nn_degree({});
    nn_degree.resize(nw.adjacencylist.size()); // k_max <= size of network
    int k_max = 0;

    for (node_t node = 0; node < (node_t)nw.adjacencylist.size(); node++)
    {
        const int k = (int) nw.adjacencylist[node].size();
        k_max = std::max(k,k_max);

        for (node_t neigh = 0; neigh < k; neigh++)
        {
            const int neigh_k = (int) nw.adjacencylist[neigh].size();
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
        for (int i = 0; i < (int)nn_degree[k].size(); i++)
        {
            average += nn_degree[k][i];
        }
        result.push_back(average /= nn_degree[k].size());

    }

    return result;


}

//------------------------------------------------
//-----Measure degree correlation in a network----
//------------------------------------------------




//
// assortativity a = num/den
// num = sum_kk' [ w(k,k') * k * k' ] - ( sum_k [ w(k) * k ] ) ^ 2
//den = sum_k [ w(k) * k ^ 2] - ( sum_k [ w(k) * k ] ) ^ 2

double assortativity(graph_adjacencylist& nw)
{
    std::vector<std::vector<double>> wkk = Wkk(nw);
    std::vector<double> wk = Wk(wkk);
    
    double k1 = 0.0;
    double k2 = 0.0;
    double kk = 0.0;
    for (int k =0; k<(int) wk.size(); k++) {
        k1 += wk[k] * k;
        k2 += wk[k] * k * k;
        
        for (int k_prime = 0; k_prime<(int) wk.size(); k_prime++) {
            kk += wkk[k][k_prime] * k * k_prime;
        }
    }
    
    double num = kk - k1*k1;
    double den = k2 - k1 *k1;
    double r = num/den;
    if (abs(r)>1) {
        throw std::range_error(std::string("correlation stronger than ± 1:  r=") + std::to_string(r));
    }
    return r ;
    

}
//    double num = 0.0;
//    double den = 0.0;
//
//    // number of edges
//    int m = 0;
//    for (std::vector<node_t> neighbours : nw.adjacencylist)
//    {
//        m += neighbours.size();
//    }
//
//
//    for (int i = 0; i < nw.adjacencylist.size(); i++)
//    {
//        for (int j = 0; j < i; j++)
//        {
//            const int Aij = std::find(nw.adjacencylist[i].begin(), nw.adjacencylist[i].end(), j) != nw.adjacencylist[i].end();
//            num += Aij - nw.outdegree(i)*nw.outdegree(j)/ (m);
//            den += nw.outdegree(i)*Aij - nw.outdegree(i)*nw.outdegree(j)/ (m);
//        }
//    }
//    return num/den;
//}

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

std::vector<double> Wk(std::vector<std::vector<double>>& wkk){
    std::vector<double> wk({});
    for (std::vector<double> v: wkk)
        wk.push_back(std::accumulate(v.cbegin(), v.cend(), 0.0));
    return wk;
}
