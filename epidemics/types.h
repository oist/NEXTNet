//
//  types.hpp
//  epidemics
//
//  Created by Florian G. Pflug on 2022/03/08.
//

#pragma once

#include "stdafx.h"

typedef int node_t;
typedef int index_t;
typedef double interval_t;
typedef double absolutetime_t;
typedef std::pair<node_t,node_t> edge_t;

enum class event_kind {
    none = 0, infection = 1, outside_infection = 2, reset = 3
};

struct event_t {
    event_kind kind = event_kind::none;
    node_t node = -1;
    absolutetime_t time = INFINITY;
};

/******************************
 * rng_t - the RNG (random number generator) to use
 *
 * The RNG to use can either be selected from a list of supported RNGs by
 * defining RNG to be one of the RNG_* constants below. RNG_STDMT19937 for
 * example is the 32-bit Mersenne Twister from the C++11 random library.
 * 
 * Or a custom RNG can be provided by defining RNG=RNG_CUSTOM, and implementing
 * the function std::uint32_t rng_t::operator()(). 
 * Define RNG=RNG_STDMT19937 to use the Mersenne Twister from the C++ stdlib,
 * or define RNG=RNG_CUSTOM and RNG_TYPE=<the RNG engine type>.
 ******************************/

#define RNG_STDMT19937 1
#define RNG_CUSTOM 2

#if !defined(RNG) || (RNG == RNG_STDMT19937)

/* Use C++ mersenne twister (default) */
#define RNG RNG_STDMT19937
typedef std::mt19937 rng_t;

#elif RNG == RNG_CUSTOM

/* Forward-declare state type */
struct RNG_STATE_TYPE;

/* Custom RNG type. operator() must be implemented by the user */
struct rng_t {
  typedef std::uint32_t result_type;

  static constexpr result_type min() { return 0; }
  
  static constexpr result_type max() { return UINT32_MAX; }

  result_type operator()();
  
#if defined(RNG_STATE_TYPE)
  typedef RNG_STATE_TYPE state_type;
  
  typedef std::uniqie_ptr<state_type> state_ptr; 
  
  rng_t(state_ptr state_) :state(std::move(state_)) {}
  
  state_ptr state;
#endif
};
  
#else

/* Unknown RNG */
#error "unknown RNG " ## RNG

#endif
