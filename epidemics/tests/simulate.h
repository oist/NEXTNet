#pragma once

#include "tests/stdafx.h"
#include "tests/parallel.h"

#include "types.h"
#include "random.h"

/**
 * @brief simulate and return transmission times
 * @param engine RNG engine
 * @param psi transmission time distribution
 * @param T stopping time
 * @param Mn number of network instances to simulate on
 * @param Ms number of simulations per network instance
 * @param args parameters to pass to the network
 * @return a pair of pairs containing (1) ordered transmission times
 *         and (2) total number of infections
 */
template<typename N, typename S, typename ...Args>
std::pair<std::vector<absolutetime_t>, std::vector<absolutetime_t>>
simulate_SIR(rng_t& engine, const transmission_time& psi, absolutetime_t T, std::size_t Mn, std::size_t Ms, Args&&...args) {
    auto ts = parallel<std::vector<double>>(Mn, engine, [&psi, T, Ms, &args...](rng_t& thread_engine){
        std::vector<double> t = {};

        // Simulate Ms times on one network
        N nw(std::forward<Args>(args)..., thread_engine);
        for(std::size_t i=0; i < Ms; ++i) {
            S sim(nw, psi);
            sim.add_infections({ std::make_pair(0, 0.0)});

            // Run simulation, collect transmission times
            while (true) {
                auto point = sim.step(thread_engine);
                if (!point || (point->time > T))
                    break;
                t.push_back(point->time);
            }
        }

        return t;
    });

    // Merge results
    std::vector<double> t;
    for(const auto& tp: ts)
      std::copy(tp.begin(), tp.end(), std::back_inserter(t));

    // Sort transmission times (only necessary for Mn or Ms > 1)
    std::sort(t.begin(), t.end());

    // Create vector representing the total number of infections
    // at each time point
    std::vector<double> y = { (1.0/Mn) * (1.0/Ms) };
    for(std::size_t i=1; i < t.size(); ++i)
        y.push_back(y.back() + (1.0/Mn) * (1.0/Ms) );

    return std::make_pair(t, y);
}


/**
 * @brief simulate and return transmission times
 * @param engine RNG engine
 * @param psi transmission time distribution
 * @param T stopping time
 * @param args parameters to pass to the network
 * @return a pair of pairs containing (1) ordered transmission times
 *         and (2) total number of infections
 */
template<typename N, typename S, typename ...Args>
void
simulate_SIS(rng_t& engine, transmission_time& psi, transmission_time& rho,
             std::vector<absolutetime_t> &times, std::vector<double> &infections, std::vector<double>& infected,
             absolutetime_t T, Args&&...args)
{
    N nw(std::forward<Args>(args)..., engine);
    S sim(nw, psi, &rho);
    sim.add_infections({ std::make_pair(0, 0.0)});
    int current_infected = 0, total_infected = 0;
    // Run simulation, collect transmission times
    while (true) {

        auto point = sim.step(engine);
        if (!point || (point -> time > T))
            break;

        switch (point-> kind) {
            case event_kind::infection:
            case event_kind::outside_infection:
                ++total_infected;
                ++current_infected;
                break;
            case event_kind::reset:
                --current_infected;
                break;
            default:
                throw std::logic_error("unexpected event kind");
        }

        times.push_back(point->time);
        infections.push_back(total_infected);
        infected.push_back(current_infected);
    }
}

