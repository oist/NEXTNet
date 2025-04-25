# Overview

The *NEXT-NET* (**N**ext-reaction-based **E**pidemics e**X**tended to **T**emporal **Net**works) C++ library contains efficient algorithms for simulating epidemics with time-varying infectiousness (so-called "non-Markovian" epidemics) on complex networks. The main algorithm provided by *NEXT-NET* is based on the next reaction method and scales to networks with millions of nodes, see our [preprint](https://arxiv.org/abs/2412.07095).

The recommended way of using *NEXT-NET* is through our wrappers [*NEXTNetR*](https://oist.github.io/NEXTNetR) and [*nextnet*](https://github.com/oist/NEXTNetPy) for R and Python. Currently, the R wrapper provides the most complete access to the features of the *NEXT-NET* library. For quick experiments the C++ library also includes a simply command-line interface for running simulations on a range of networks.

# Synopsis

To download and install the *NEXT-NET* command-line interface do
```
git clone --recurse-submodules https://github.com/oist/NEXTNet.git
cd NEXTNet
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../
make
```

To simulate an epidemic starting until time t=250 from node 1 of a Watts-Strogatz network of size N=10.000 nodes with  parameters K=6 and beta=0.1, gamma distributed infection times with mean 8 and variance 1, and log-normally distributed recovery times with mean 50 and variance 20 run

```
./nextnet -n 'watts_strogatz(10000, 4, 0.2)' -p 'gamma(8,1)' -r 'lognormal(30,20)' -i 1 -t 250 > trajectory.txt
```

To plot the resulting number of currently infected nodes over time using gnuplot, run

```
set key autotitle columnhead
set log y
plot 'trajectory.txt' using 1:9 with lines
```

The available network types and infection time distributions can be listed with 
```
./nextnet --list-networks
./nextnet --list-times
```

# Functionality

The NEXT-NET library offers various epidemic models, static and dynamic
networks, transmission time distributions, and simulation algorithms. 

## Epidemic models

* *SI* (**S**usceptible **I**nfected). Individuals start out susceptible and become infected upon transmission of the disease.
* *SIS* (**S**usceptible **I**nfected **S**usceptible). Individuals start out susceptible, become infected upon transmission of the disease, and eventually recover and become susceptible again.
* *SIR* (**S**usceptible **I**nfected **R**ecovered). Individuals start out susceptible, become infected upon transmission of the disease, and eventually recover. Recovered individuals do not become susceptible agian.

## Types of networks

### Static networks

* [*Erdős–Rényi*](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model). Random networks in two nodes are connected by an edge with a certian probabilitz. These networks exhibit no clustering or hubs.
* [*Watts-Strogatz*](https://en.wikipedia.org/wiki/Watts%E2%80%93Strogatz_model). Random small-world networks exhibiting clustering and hubs.
* [*Barabási-Albert*](https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model). Random preferential attachment networks exhibiting hubs and scale-free degree distributions.
* *Acyclic networks*. Random tree-like networks with Poissonian distribution of number of offsprings
* *Configuration mode*l. Random networks with given degree distributions
* *Clustered configuration model*. Random networks with given degree distribution and clustering coefficients.
* *Lattice*. Regular spatial networks.
* *Empirical networks*. Networks defined by an adjacency list, i.e. a list of neighbours of each node. See [*NEXTNet-EmpiricalNetworks*](https://github.com/oist/NEXTNet-EmpiricalNetworks) for a list of empirical networks that can directly be loaded into *NEXTNet*.

### Temporal networks

* *Temporal Erdös-Réyni*. Temporal networks that at any instant resembles an [*Erdős–Rényi*](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model) networks, but in which nodes randomly detach from a currently neighbour and reattach to a new, randomly chosen neighbour.
* *SIRX network*. A network version of the [SIRX model](https://www.science.org/doi/10.1126/science.abb4557) in nodes stochastically deactivate, i.e. remove their links, in response to an epidemic outbreak.
* *Activity-driven network*. A temporal network in which nodes stochastically detatch from their current neighbour and reattach elsewhere.
* *Brownian proximity network*. A temporal network with spatial structure in which nodes move randomly in two dimensions and links represent spatial proximity between nodes.
* *Epirical contact networks*. Networks defined a list *(t,i,j)* of contact times *t* between individuals *i* and *j*.

## Time distributions

Time distributions define (a) the time it takes from the infection of a node until it transmits the disease across a specific link, and (b) the time it takes for a node to recover. NEXT-NET offers the following pre-defined distributions, and allows users to easily implement additional distributions.

* *Exponential*.
* *Deterministic*. 
* *Gamma*. 
* *Lognormal*
* *Weibull*
* *Polynomial rate*.

## Simulation algorithms

* [*Next reaction method*](https://doi.org/10.1063/1.2799998). The most efficient algorithm available, in which the runtime to find the next infection depends only weakly on the number of infected nodes. Scales to networks with millions of nodes. 
* [*REGIR* (**Re**jection **Gi**llespie algorithm for non-Markovian **R**eactions *)](https://arxiv.org/abs/2212.05059). Based on a similar approximation as *nGMA*, but removes the quadratic growth of the time it takes to find the next infection. Typically slower than the *next reaction method* but scales similarly with network size.
* [*nMGA* (**n**on-**M**arkovian **G**illespie **A**lgorithm)](https://doi.org/10.1103/PhysRevE.90.042108). An approximate algorithms which works well on small networks. However, the runtime to find the next infection grows quadratically with the number of infected nodes, which limits the scalability to large networks.
