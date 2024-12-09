# Overview

The NEXT-NET C++ library contains efficient algorithms for simulating non-Markovian epidemics on complex networks. The main algorithm provided by
NEXT-NET is based on the next reaction method, and scales to networks with millions of nodes. The 

# Functionality

## Epidemic models

* *SI*
* *SIS*
* *SIR*

## Simulation algorithms

* *Next reaction method*
* *REGIR*
* *nMGA*

## Types of networks

### Static networks

* *Erdös-Réyni*. 
* *Watts-Strogatz*.
* *Barabasi-Albert*.
* *Acyclic networks*. (Poissionan-distributed number of neighbours
* *Configuration mode*l. ()
* *Clustered configuration model*.
* *Lattice*.
* *Empirical networks*.

### Temporal networks

* *Temporal Erdös-Réyni*
* *SIRX network*
* *Activity-driven network*
* *Brownian proximity network*

## Transmission time distributions

* *Exponential*
* *Deterministic*
* *Gamma*
* *Lognormal*
* *Weibull*
* *Polynomial rate*

# Simulating epidemics in Python and R

The functionality of this library can be accessed in Python through the package
*nextnet* (https://github.com/oist/NEXTNetPy) and in R through *NEXTNetR*
(https://github.com/oist/NEXTNetR).
