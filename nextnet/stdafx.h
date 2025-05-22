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
#include <cstdint>
#include <cmath>
#include <memory>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <deque>
#include <optional>
#include <variant>
#include <any>
#include <utility>
#include <numeric>
#include <chrono>
#include <filesystem>

// Macros

#define STRINGIFY(v) STRINGIFY_(v)
#define STRINGIFY_(v) #v

// Boost

#if defined(__clang__)
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Weverything"
#    pragma clang diagnostic ignored "-Wdeprecated-declarations"
#    pragma clang diagnostic ignored "-Wno-pedantic"
#elif defined(__GNUC__) || defined(__GNUG__)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wall"
#    pragma GCC diagnostic ignored "-Wsign-compare"
#    pragma GCC diagnostic ignored "-Wparentheses"
#    pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#    pragma GCC diagnostic ignored "-Wno-pedantic"
#endif

#include <boost/config.hpp>

#ifdef NEXTNET_BOOST_NO_CXX98_FUNCTION_BASE
#    ifndef BOOST_NO_CXX98_FUNCTION_BASE
#        define BOOST_NO_CXX98_FUNCTION_BASE
#    endif
#endif

#ifdef NEXTNET_BOOST_NO_CXX17_HDR_EXECUTION
#    ifndef BOOST_NO_CXX17_HDR_EXECUTION
#        define BOOST_NO_CXX17_HDR_EXECUTION
#    endif
#endif

#include <boost/math/distributions.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/tag.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>

#if defined(__clang__)
#    pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
#    pragma GCC diagnostic pop
#endif

/* Make boost namespaces available under shortcuts */
namespace bm  = boost::math;
namespace bmi = boost::multi_index;

// Priority Queue

#define STD_PRIORITY_QUEUE_DEQUE 1
#define EXT_PRIO_QUEUE 2

#ifndef NEXT_REACTION_QUEUE
#    define NEXT_REACTION_QUEUE EXT_PRIO_QUEUE
#endif

#if NEXT_REACTION_QUEUE == EXT_PRIO_QUEUE
#    include "prio_queue/prio_queue.hpp"
#endif

// Dynamic Distribution

/* Use C++ standard library random number generation.
 * This should work regardless of the specific RNG engine implementation
 * was selected via the RNG defines, since RNG engines in the epidemic
 * project are all assumed to comply to the C++ standard library API
 */
#include "dyndist/rng_stdcxx.h"
#include "dyndist/vector_distribution.h"

// To fix warnings on bools being 'unused'

#define _unused(x) ((void)(x))
