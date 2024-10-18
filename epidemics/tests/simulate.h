#pragma once

#include "tests/stdafx.h"
#include "tests/parallel.h"

#include "types.h"
#include "random.h"
#include "algorithm.h"

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

    if (t.empty())
        return std::make_pair(t, std::vector<double>());

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
template<typename Factory>
void
average_trajectories(rng_t& engine, Factory factory,
					 std::vector<absolutetime_t> &times,
					 std::vector<double>& cumulative_infected,
					 std::vector<double>& current_infected,
					 absolutetime_t Tmax, int runs = 1, int output_every = -1)
{
	typedef std::pair<double, double> pair_t;
	typedef std::vector<pair_t> pairs_t;
	
	if (output_every < 0)
		output_every = runs;
	
	// Run <runs> simulations, collect a large vector of
	// (time, delta) pairs which allow easy averaging
	auto results_par = parallel<pairs_t>(runs, engine, [factory, Tmax](rng_t& thread_engine) {
		pairs_t r = {};

		// Create simulation environment
		auto s = factory(thread_engine);
		simulation_algorithm& sim = *s.simulator;
		sim.add_infections({ std::make_pair(0, 0.0)});
		
		// Run simulation, collect transmission times
		while (true) {
			auto point = sim.step(thread_engine);
			if (!point || (point->time > Tmax))
				break;
			r.push_back({point->time, delta_infected(point->kind)});
		}

		return r;
	});
	
	// Merge results
	pairs_t results;
	for(const auto& r: results_par)
	  std::copy(r.begin(), r.end(), std::back_inserter(results));

	// Sort by times
	std::sort(results.begin(), results.end(), [](const pair_t& a, pair_t& b) {
		return a.first < b.first;
	});
	
	// Output results
	int i = 0;
	double cur_inf = 0.0;
	double cum_inf = 0.0;
	for(const pair_t& p: results) {
		/* Track number of infected */
		const double t = p.first;
		const double d = p.second / (double)runs;
		cur_inf += d;
		cum_inf += std::max(d, 0.0);

		/* Only output every <output_every>-th element */
		if ((output_every > 0) && (++i != output_every))
			continue;
		i = 0;

		times.push_back(t);
		cumulative_infected.push_back(cum_inf);
		current_infected.push_back(cur_inf);
	}
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
template<typename Factory, typename DeltaInfectedFunctor>
void
average_trajectories(rng_t& engine, Factory factory, DeltaInfectedFunctor f,
					 std::vector<absolutetime_t> &times,
					 std::vector<double>& cumulative_infected,
					 std::vector<double>& current_infected,
					 absolutetime_t Tmax, int runs = 1, int output_every = -1)
{
	typedef std::pair<double, double> pair_t;
	typedef std::vector<pair_t> pairs_t;
	
	if (output_every < 0)
		output_every = runs;
	
	// Run <runs> simulations, collect a large vector of
	// (time, delta) pairs which allow easy averaging
	auto results_par = parallel<pairs_t>(runs, engine, [factory, f, Tmax](rng_t& thread_engine) {
		pairs_t r = {};

		// Create simulation environment
		auto s = factory(thread_engine);
		auto& sim = *s.simulator;
		
		// Run simulation, collect pairs
		while (true) {
			// Perform step, stop if no more events
			const auto point = sim.step(thread_engine);
			if (!point)
				break;
			// Convert event to (time, delta) pair, stop if time exceeds Tmax
			const pair_t p = f(*point);
			if (p.first > Tmax)
				break;
			// Only store points with non-zero delta
			if (p.second != 0.0)
				r.push_back(p);
		}

		return r;
	});
	
	// Merge results
	pairs_t results;
	for(const auto& r: results_par)
	  std::copy(r.begin(), r.end(), std::back_inserter(results));

	// Sort by times
	std::sort(results.begin(), results.end(), [](const pair_t& a, pair_t& b) {
		return a.first < b.first;
	});
	
	// Output results
	int i = 0;
	double cur_inf = 0.0;
	double cum_inf = 0.0;
	for(const pair_t& p: results) {
		/* Track number of infected */
		const double t = p.first;
		const double d = p.second / (double)runs;
		cur_inf += d;
		cum_inf += std::max(d, 0.0);

		/* Only output every <output_every>-th element */
		if ((output_every > 0) && (++i != output_every))
			continue;
		i = 0;

		times.push_back(t);
		cumulative_infected.push_back(cum_inf);
		current_infected.push_back(cur_inf);
	}
}

