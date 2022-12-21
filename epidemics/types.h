//
//  types.hpp
//  epidemics
//
//  Created by Florian G. Pflug on 2022/03/08.
//

#pragma once

#include "stdafx.h"

/**
 * Represents a node in a graph
 */
typedef int node_t;

/**
 * Index into the neighbours of a node
 */
typedef int index_t;

/**
 * Represents a time interval
 */
typedef double interval_t;

/**
 * Represents an absolute time
 */
typedef double absolutetime_t;

/**
 * Represents an edge in a graph connecting two nodes;
 */
typedef std::pair<node_t,node_t> edge_t;

/**
 * Type of event that occured during a simulation
 * TODO: Rename to epidemic_event_kind
 */
enum class event_kind : unsigned int {
    none = 0,
	infection = 1,
	outside_infection = 2,
	reset = 3
};
  
/**
 * Translate event kinds to their name
 */
inline const char* name(event_kind kind) {
    switch (kind) {
        case event_kind::none: return "none";
        case event_kind::infection: return "infection";
        case event_kind::outside_infection: return "outside_infection";
        case event_kind::reset: return "reset";
        default: return NULL;
    }
}

/**
 * Describes an event that occured
 * TODO: Rename to epidemic_event_t
 */
struct event_t {
    event_kind kind = event_kind::none;
    node_t node = -1;
    absolutetime_t time = INFINITY;
};

/**
 * Type of dynamic network event that occured during a simulation
 */
enum class network_event_kind : unsigned int {
	none = 0,
	neighbour_added = 1,
	neighbour_removed = 2
};
  
/**
 * Translate dynamic network event kinds to their name
 */
inline const char* name(network_event_kind kind) {
	switch (kind) {
		case network_event_kind::none: return "none";
		case network_event_kind::neighbour_added: return "neighbour_added";
		case network_event_kind::neighbour_removed: return "neighbour_removed";
		default: return NULL;
	}
}

/**
 * Describes an event that occured when evolving a dynamic network
 */
struct network_event_t {
	network_event_kind kind = network_event_kind::none;
	node_t node = -1;
	node_t neighbour = -1;
	absolutetime_t time = INFINITY;
};

/**
 * Describes an event that occured during a simulation on a dynamic network.
 * Can either be a change to the network, or an infection or recovery.
 */
typedef std::variant<event_t, network_event_t> network_or_epidemic_event_t;

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

	/* C++ standard RNG */

	/* Use C++ mersenne twister (default) */
	#define RNG RNG_STDMT19937
	typedef std::mt19937 rng_t;

#elif RNG == RNG_CUSTOM

	/* Curstom RNG */

	/* Forward-declare state type */
	struct RNG_STATE_TYPE;

	/* Custom RNG type. operator() must be implemented by the user */
	struct rng_t {
	  typedef std::uint32_t result_type;

	  static constexpr result_type min() { return 0; }
	  
	  static constexpr result_type max() { return UINT32_MAX; }

	  template<class Sseq>
	  explicit rng_t(Sseq& s)
	#if defined(RNG_STATE_TYPE)
		: state(new state_type(s)) {}
	#else
		{ throw std::logic_error("not implemented"); }
	#endif

	  result_type operator()();
	  
	#if defined(RNG_STATE_TYPE)
	  typedef RNG_STATE_TYPE state_type;
	  
	  typedef std::uniqie_ptr<state_type> state_ptr;

	  template<typename Args...>
	  rng_t(Args... &&args) :state(new state_type(std::forward<Args>(args)...)) {}
	  
	  rng_t(state_ptr state_) :state(std::move(state_)) {}
	  
	  state_ptr state;
	#else
	  rng_t() {}
	#endif
	};
  
#else

	/* Unknown RNG */
	#error "unknown RNG " ## RNG

#endif
