//
//  stdafx.h
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//
// Added to follow apparent conventions of including all libraries in a stdafx file.
#pragma once

// C++ standard library includes

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <deque>
#include <optional>
#include <variant>
#include <utility>
#include <numeric>
#include <chrono>

// Boost

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wparentheses"
#endif

#include <boost/math/distributions.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/tools/roots.hpp>

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

/* Make boost::math available as simply bm */
namespace bm = boost::math;

// Priority Queue

#define STD_PRIORITY_QUEUE_DEQUE 1
#define EXT_PRIO_QUEUE 2

#ifndef NEXT_REACTION_QUEUE
#define NEXT_REACTION_QUEUE EXT_PRIO_QUEUE
#endif

#if NEXT_REACTION_QUEUE == EXT_PRIO_QUEUE
#include "../ext/prio_queue/prio_queue.hpp"
#endif

// Dynamic Distribution

/* Use C++ standard library random number generation.
 * This should work regardless of the specific RNG engine implementation
 * was selected via the RNG defines, since RNG engines in the epidemic
 * project are all assumed to comply to the C++ standard library API
 */
#include "dyndist/rng_stdcxx.h"
#include "dyndist/vector_distribution.h"
