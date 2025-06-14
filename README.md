# Overview

The *NEXT-NET* (**N**ext-reaction-based **E**pidemics e**X**tended to **T**emporal **Net**works) C++ library contains efficient algorithms for simulating epidemics with time-varying infectiousness (so-called "non-Markovian" epidemics) on complex networks. The main algorithm provided by *NEXT-NET* is based on the next reaction method and scales to networks with millions of nodes, see our [preprint](https://arxiv.org/abs/2412.07095).

We recommend that users who are already familiar with Python or R use *NEXT-Net* through the packages [*NEXTNetR*](https://oist.github.io/NEXTNetR) and [*NEXTNetPy*](https://github.com/oist/NEXTNetPy) which provide convenient and flexible access to the features of the C++ library. Currently, *NEXTNetR* is the most well-tested library and provides the most complete access to *NEXTNet*'s features; however we aim to provide the same set of features in both packages eventually. For users not familiar with Python or R, the C++ library comes with a basic command-line tool called `nextnet`. This tool is less flexible than the Python and R packages when it comes to setting up and running simulations, but provides quick access to the main features of NEXT-Net. It also makes it easy to test custom modifications of the C++ code.

With *[NEXTNet-EmpiricalNetworks](https://github.com/oist/NextNet-EmpiricalNetworks)* we also provide a database of empirical networks in a format compatible with *NEXT-Net*, compiled from the SNAP ([Leskovec and Krevl, 2014](http://snap.stanford.edu/data)), ICON ([Clauset *et al.*, 2016](https://icon.colorado.edu/)), and KONECT ([Kunegis, 2013](https://doi.org/10.1145/2487788.2488173)) and other publications ([Rocha *et al.*, 2010](https://doi.org/10.1073/pnas.0914080107)).

# Installation

To download and install a prebuilt binary of the *NEXT-NET* command-line interface, go to [*Releases*](https://github.com/oist/NEXTNet/releases), select the latest release, download the binary release for your platform, and decompress; this should yield a binary called `nextnet`. Currently binary releases for Linux system (`NEXTNet-v<VERSION>-linux-x86_64.tar.gz`, supports most Linux systems) and Mac OS (`NEXTNet-v<VERSION>-mac.tar.gz`, supports Mac OS 10.13 and later) are provided. For other platforms, *NEXT-Net* must be built from source, see below for instructions.

To build *NEXT-Net* from source, and working C++ compiler with C++17 support, *[cmake](https://cmake.org/)* and *[Boost](https://www.boost.org/)* are required. If *git* is installed, the latest version of *NEXT-Net* can be downloaded with
```
git clone --recurse-submodules --branch latest-release https://github.com/oist/NEXTNet.git
```
Otherwise, please got to [*Releases*](https://github.com/oist/NEXTNet/releases), select the latest release, download `NEXTNet-v<VERSION>-source.tar.gz`, decompress into a folder `NEXTNet`, and run the following commands to produce the command-line tool `nextnet`
```
cd NEXTNet
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

# Synopsis

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
* *SIR* (**S**usceptible **I**nfected **R**ecovered). Individuals start out susceptible, become infected upon transmission of the disease, and eventually recover. Recovered individuals do not become susceptible agian (Currently only available through the *R* and *Python* packages).

## Types of networks

*NEXT-Net* supports a range of synthetic and empirical networks, both static and dynamic. In the command-line simulator, the list of available networks can be queried with `--list-networks`, and a specific network is selected with option `-n`.

### Static networks

* [*Erdős–Rényi*](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model). Random networks in two nodes are connected by an edge with a certian probability. These networks exhibit no clustering or hubs. Command-line syntax `erdos_renyi(size, avg_degree)`.
* [*Watts-Strogatz*](https://en.wikipedia.org/wiki/Watts%E2%80%93Strogatz_model). Random small-world networks exhibiting clustering and hubs. Command-line syntax `watts_strogatz(size, k, p)`.
* [*Barabási-Albert*](https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model). Random preferential attachment networks exhibiting hubs and scale-free degree distributions. Command-line syntax `barabasi_albert(size, m)`.
* *Acyclic networks*. Random tree-like networks with Poissonian distribution of number of offsprings. Command-line syntax `acyclic(avg_degree) -i 1`.
* *Configuration mode*l. Random networks with given degree distributions. Command-line syntax `config_model(degrees)`.
* *Clustered configuration model*. Random networks with given degree distribution and clustering coefficients. Command-line syntax `config_model_clustered(degrees, triangles, beta)`.
* *Lattice*. Regular spatial networks. Command-line syntax `lattice_2d(edgelength)`, `lattice_3d(edgelength)`
* *Empirical networks*. Networks defined by an adjacency list, i.e. a list of neighbours of each node. See [*NEXTNet-EmpiricalNetworks*](https://github.com/oist/NEXTNet-EmpiricalNetworks) for a list of empirical networks that can directly be loaded into *NEXTNet*. Command-line syntax `empirical(file)` for unweighted networks and `weighted_empirical(file)` for weighted networks. See the R package documentation for [empirical_network](https://oist.github.io/NEXTNetR/reference/empirical_network.html) and [empirical_weighted_network](https://oist.github.io/NEXTNetR/reference/empirical_weightednetwork.html) for a description of the file format and optional parameters. Gzip-compressed files are supported if they end in ".gz".

### Temporal networks

* *Temporal Erdös-Réyni*. Temporal networks that at any instant resembles an [*Erdős–Rényi*](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model) networks, but in which nodes randomly detach from a currently neighbour and reattach to a new, randomly chosen neighbour. Command-line syntax `temporal_erdos_renyi(size, avg_degree)`, see the R package documentation for [erdos_renyi_temporalnetwork](https://oist.github.io/NEXTNetR/reference/erdos_renyi_temporalnetwork.html) for details.
* *SIRX network*. A network version of the [SIRX model](https://www.science.org/doi/10.1126/science.abb4557) in nodes stochastically deactivate, i.e. remove their links, in response to an epidemic outbreak. Command-line syntax `temporal_sirx(nw, kappa0, kappa)`, see the R package documentation for [sirx_temporalnetwork](https://oist.github.io/NEXTNetR/reference/sirx_temporalnetwork.html) for details.
* *Activity-driven network*. A temporal network in which nodes stochastically detatch from their current neighbour and reattach elsewhere. Command-line syntax `activity_driven(activities = {1}, m, eta_sus, eta_inf, b_sus, b_inf)`, see the R package documentation for [activity_driven_temporalnetwork](https://oist.github.io/NEXTNetR/reference/activity_driven_temporalnetwork.html) for details.
* *Brownian proximity network*. A temporal network with spatial structure in which nodes move randomly in two dimensions and links represent spatial proximity between nodes. Command-line syntax `brownian_proximity(size, avg_degree, radius = 1, D0 = 1, D1 = 1)`, see the R package documentation [brownian_proximity_temporalnetwork](https://oist.github.io/NEXTNetR/reference/brownian_proximity_temporalnetwork.html) for details.
* *Epirical contact networks*. Networks defined a list *(t,i,j)* of contact times *t* between individuals *i* and *j*. Command-line syntax `empirical_contact(file, contact_kind = instantaneous, dt = 1)`, see the R package documentation for [empirical_contact_temporalnetwork](https://oist.github.io/NEXTNetR/reference/empirical_contact_temporalnetwork.html) for a description of the file format. Gzip-compressed files are supported if they end in ".gz".

## Time distributions

Time distributions define (a) the time it takes from the infection of a node until it transmits the disease across a specific link, and (b) the time it takes for a node to recover. NEXT-NET offers the following pre-defined distributions, and allows users to easily implement additional distributions. The list of avaible time distribution can be queried on the command-line with `--list-times`. Many of the common distributions listed below have an additional parameter `pinf` which specifies the probability of no infection. See the R package manual pages [](https://oist.github.io/NEXTNetR/reference/time_distributions.html) and [time_functions](https://oist.github.io/NEXTNetR/reference/time_functions.html) for additional details. The distribution can be selected on the command-line with option `-p` for the transmission time and `-r` for the recovery time.

* *Exponential*. Exponential distribution with a given rate and optionally a given probability of no infection (i.e. an infinite transmission time). Command-line syntax `exponential(lambda, pinf = 0)`. 
* *Gamma*. Gamma distribution parametrized with mean and variance and optionally the probability of no infection. Command-line syntax `gamma(mean, variance, pinf = 0)`.
* *Lognormal*. Log-normal distribution parametrized with mean and variance and optionally the probability of no infection. Command-line syntax `lognormal(mean, variance, pinf = 0)`.
* *Weibull*. Weibull distribution parametrized with shape and scale and optionally the probability of no infection. Command-line syntax `weibull(shape, scale, pinf = 0)`.
* *Polynomial rate*. Polynomial infectiousness function of arbitrary degree. Command-line syntax `polynomial_rate(coeffs)`, where coeffs is a list `{c0, c1, ...}` of coefficients for x, x squared, etc.
* *Infectiousness time-course*. Arbitrary infectionsness function defined by a list of times and corresponding infection rates. Command-line syntax `infectiousness(taus, lambdas)`, where `taus` and `lambdas` are lists `{v0, v1, ...}` of values describing time points and the corresponding lambda values.
* *Deterministic*. Fixed time of infection. Command-line syntax `deterministic(tau)`.

## Simulation algorithms

The recommonded simulation algorithm for all user-cases is currently the *NEXT-Net* algorithm. The command-line simulator allows the algorithm to be selected from the following list by specifying option `-a next|regir|nmga`.

* [*Next reaction method*](https://doi.org/10.1063/1.2799998). The most efficient algorithm available, in which the runtime to find the next infection depends only weakly on the number of infected nodes. Scales to networks with millions of nodes. Command-line syntax `next`. 
* [*REGIR* (**Re**jection **Gi**llespie algorithm for non-Markovian **R**eactions *)](https://arxiv.org/abs/2212.05059). Based on a similar approximation as *nGMA*, but removes the quadratic growth of the time it takes to find the next infection. Typically slower than the *next reaction method* but scales similarly with network size. Command-line syntax `regir`
* [*nMGA* (**n**on-**M**arkovian **G**illespie **A**lgorithm)](https://doi.org/10.1103/PhysRevE.90.042108). An approximate algorithms which works well on small networks. However, the runtime to find the next infection grows quadratically with the number of infected nodes, which limits the scalability to large networks. Command-line syntax `ngma`.
