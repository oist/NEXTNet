#pragma once

#include "stdafx.h"
#include "algorithm.h"
#include "network.h"

typedef std::function<simulation_algorithm *(rng_t &engine, int n)> simulation_factory_t;

// function to compute the derivative of f(x) using central finite differences method
std::vector<double> derivative(std::vector<double> &X, std::vector<double> &Y);

void knn_BA();

void exportData(std::vector<double> &trajectory, std::string filename);

void exportData(std::vector<double> &X, std::vector<double> &Y, std::string filename);

void exportData(std::vector<int> &X, std::vector<double> &Y, std::string filename);

void exportData(std::vector<int> &trajectory, std::string filename);

void export_adjacency_list(const std::vector<std::vector<node_t>> &adjacencyList, std::string filename);
void export_adjacency_matrix(const std::vector<std::vector<node_t>> &adjacencyList, std::string filename);
void export_adjacency_dot(const std::vector<std::vector<node_t>> &adjacencyList, std::string filename, bool directed = false);

void print_matrix(std::vector<std::vector<double>> &A);

double measure_running_time(adjacencylist_network &network, rng_t &engine);

void generate_data_running_time(rng_t &engine, int size, bool isNMGA);

// void average_epidemic(graph_adjacencylist& network,rng_t& engine, )
// grjijffw

void measure_runtime(rng_t &engine, simulation_factory_t factory, int SMAX, double TMAX, std::string filename);

void measure_running_time_next_reaction_ER(rng_t &engine, int SMAX, std::string filename);
void measure_running_time_next_reaction_BA(rng_t &engine, int SMAX, std::string filename);
void measure_running_time_nMGA_ER(rng_t &engine, int SMAX, std::string filename);
void measure_running_time_nMGA_BA(rng_t &engine, int SMAX, std::string filename);

template <typename Factory>
void measure_runtime(rng_t &engine, Factory factory,
                     int SMAX, double TMAX, int KMAX, std::string filename)
{
    using namespace std;

    // fill sizes
    vector<int> sizes;
    for (int k = 7; k <= KMAX; k++) {
        int i = (int)floor(pow(2, k));
        sizes.push_back(i);
    }

    // fill times
    vector<double> average_run_time;
    vector<double> std_run_time;

    for (int N : sizes) {
        cout << N << "\n";
        double mean        = 0;
        double mean_square = 0;
        for (int i = 0; i < SMAX; i++) {
            // Create simulation environment
            auto s                           = factory(engine, N);
            simulation_algorithm &simulation = *s.simulator;

            // Start measuring performance
            auto start = std::chrono::high_resolution_clock::now();

            // 10% of the network will be infected at t=0;
            const int N0 = floor(simulation.get_network().nodes() / 10);
            for (node_t node = 0; node < N0; node++) {
                simulation.add_infections({ std::make_pair(node, 0.0) });
            }

            // Run simulation, collect transmission times
            while (true) {

                auto point = simulation.step(engine);
                if (!point || (point->time > TMAX))
                    break;
            }

            auto stop = std::chrono::high_resolution_clock::now();

            auto duration     = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            const double time = duration.count();

            mean += time / SMAX;
            mean_square += time * time / SMAX;
        }

        const double std_dev = std::sqrt(mean_square - mean * mean);
        average_run_time.push_back(mean);
        std_run_time.push_back(std_dev);
    }
    // export data
    exportData(sizes, average_run_time, filename + "_AV.dat");
    exportData(sizes, std_run_time, filename + "_STD.dat");
}

template <typename Factory>
void simulate_trajectory(rng_t &engine, Factory factory, double I0, double TMAX, std::string filename)
{
    using namespace std;

    // Initialise vectors containing the trajectory points
    std::vector<double> times;
    std::vector<double> infected;
    int nb_infected = 0;

    // Create simulation environment
    auto s                           = factory(engine);
    simulation_algorithm &simulation = *s.simulator;

    // I0 % of the network will be infected at t=0;
    const int N  = simulation.get_network().nodes();
    const int N0 = floor(N * I0 / 100);
    for (node_t node = 0; node < N0; node++) {
        simulation.add_infections({ std::make_pair(node, 0.0) });
    }

    // Run simulation, collect transmission times
    while (true) {

        auto point = simulation.step(engine);
        if (!point || (point->time > TMAX))
            break;
        times.push_back(point->time);

        switch (point->kind) {
            case epidemic_event_kind::infection:
            case epidemic_event_kind::outside_infection:
                nb_infected += 1;
                break;
            case epidemic_event_kind::reset:
                nb_infected -= 1;
                break;
            default:
                throw std::logic_error("invalid event kind");
                break;
        }
        infected.push_back(nb_infected);
    }

    // export data
    exportData(times, infected, filename);
}
