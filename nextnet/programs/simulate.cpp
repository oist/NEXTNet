#include "stdafx.h"
#include "random.h"
#include "analysis.h"
#include "REGIR.h"
#include "nMGA.h"
#include "NextReaction.h"

#include "popl.hpp"

using namespace std;
using namespace std::literals;
using namespace popl;

rng_t engine;

/*
 * Factories for time distribution and network objects
 */

namespace factories {

#define DECLARE_ARGUMENT(_name, _type, _defaultvalue) \
	struct _name { \
		static const char* name; \
		typedef _type value_type; \
		static const optional<value_type> defaultvalue; \
	}; \
	const char* _name::name = #_name; \
	const optional<_name::value_type> _name::defaultvalue = _defaultvalue

/**
 * Convert argument to correct type for the specified argument Arg
 */
template<typename Arg>
pair<typename Arg::value_type, int> argument(const vector<string>& vs, size_t i) {
	if (i < vs.size())
		return { boost::lexical_cast<typename Arg::value_type>(vs[i]), 1 };
	else if (Arg::defaultvalue)
		return { *Arg::defaultvalue, 0 };
	else
		throw std::runtime_error("missing value for argument "s + Arg::name);
}

/**
 * Generate description for the argument Arg
 */
template<typename Arg>
string description() {
	if (Arg::defaultvalue)
		return string(Arg::name) + " = " + boost::lexical_cast<string>(*Arg::defaultvalue);
	else
		return string(Arg::name);
}

/**
 * Convert a vector or strings into a tuple of arguments as specified by the argument pack
 */
template<typename Arg, typename ...Args>
struct argument_tuple_converter {
	typedef tuple<typename Arg::value_type, typename Args::value_type...> value_tuple_type;
	
	value_tuple_type operator()(const vector<string>& vs, std::size_t i = 0) {
		pair<typename Arg::value_type, int> a = argument<Arg>(vs, i);
		return tuple_cat(make_tuple(a.first), argument_tuple_converter<Args...>()(vs, i + a.second));
	}
};
template<typename Arg>
struct argument_tuple_converter<Arg> {
	typedef tuple<typename Arg::value_type> value_tuple_type;
	
	value_tuple_type operator()(const vector<string>& vs, std::size_t i = 0) {
		return make_tuple(argument<Arg>(vs, i).first);
	}
};

/**
 * Fill vector with descriptions of the arguments in the pack
 */
template<typename Arg, typename ...Args>
struct argument_description_converter {
	void operator()(vector<string>& v) {
		string d = description<Arg>();
		if (!d.empty())
			v.push_back(d);
		argument_description_converter<Args...>()(v);
	}
};
template<typename Arg>
struct argument_description_converter<Arg> {
	void operator()(vector<string>& v) {
		string d = description<Arg>();
		if (!d.empty())
			v.push_back(d);
	}
};


template<class T, class Tuple, std::size_t... I>
constexpr std::unique_ptr<T> make_new_from_tuple_impl(Tuple&& v, std::index_sequence<I...>)
{
	return make_unique<T>(get<I>(std::forward<Tuple>(v))...);
}

/**
 * Version of the standard libraries make_from_tuple that allocates on the heap
 */
template<class T, class Tuple>
constexpr std::unique_ptr<T> make_new_from_tuple(Tuple&& v)
{
	return make_new_from_tuple_impl<T>(std::forward<Tuple>(v),
									   make_index_sequence<tuple_size_v<remove_reference_t<Tuple>>>{});
	
}

/**
 * Factory that generated objects which derive from the common type T from strings of the form
 *   typename(arg1, arg2, ...., argN)
 */
template<typename T>
struct factory {
	/**
	 * Represents a factory function that produces type T given a list of arguments in string form
	 */
	typedef function<unique_ptr<T>(const vector<string>&)> factoryfun;
	
	/**
	 * Map of types names to factory functions
	 */
	unordered_map<std::string, factoryfun> products;

	/**
	 * List of factory descriptions
	 */
	vector<string> descriptions;

	/**
	 * Parses a string of the form "function(arg1, arg2, "quoted\ arg")" with arbitarry number of arguments
	 * into a pair comprising the function name and a vector of arguments.
	 */
	static std::pair<string, std::vector<string>> parse(std::string s)
	{
		string name;
		std::vector<string> args;
		
		enum { WS1, NAME, WS2, BRA, WS3, ARG, ARG_DQ, ARG_DQE, WS4, KET_OR_COMMA, WS5, DONE } state = WS1;
		for(std::ptrdiff_t i=0; i < s.size(); ++i) {
			const char c = s[i];
			switch (state) {
				case WS1:
				case WS2:
				case WS3:
				case WS4:
				case WS5:
					/* scan until non-whitespace, then update state */
					if (isspace(c)) continue;
					switch (state) {
						case WS1: state = NAME; --i; break;
						case WS2: state = BRA; --i; break;
						case WS3:
							switch (c) {
								case '"': state = ARG_DQ; break;
								default: state = ARG; --i; break;
							};
							args.push_back(std::string());
							break;
						case WS4: state = KET_OR_COMMA; --i; break;
						case WS5: state = DONE; --i; break;
						default: throw std::logic_error("invalid state");
					}
					break;
				case NAME:
					/* collect until non-alphanumeric character */
					if (isalnum(c)) {
						name.push_back(c);
						continue;
					}
					state = WS2;
					--i;
					break;
				case BRA:
					/* error if not opening bracket */
					if (c != '(')
						throw std::runtime_error("unable to parse '" + s + "', " +
												 "invalid symbol '" + c + "''");
					state = WS3;
					break;
				case ARG:
					/* collect until non-alphanumeric, then update state */
					if (isalnum(c) || c == '.' || c == '-' || c == '+') {
						args.back().push_back(c);
						continue;
					}
					state = WS4;
					--i;
					break;
				case ARG_DQ:
					/* scan until closing quote, then update state */
					switch (c) {
						case '"':
							state = WS4;
							break;
						case '\\':
							state = ARG_DQE;
							break;
						default:
							args.back().push_back(c);
					}
					break;
				case ARG_DQE:
					/* escaped character, collect unconditionally */
					args.back().push_back(c);
					state = ARG_DQ;
					break;
				case KET_OR_COMMA:
					switch (c) {
						case ')': state = WS5; break;
						case ',': state = WS3; break;
						default: throw std::runtime_error("unable to parse '" + s + "', " +
														  "invalid symbol '" + c + "''");
					}
					break;
				case DONE:
					throw std::runtime_error("unable to parse '" + s + "', " +
											 "invalid symbol '" + c + "''");
			}
		}
		
		return { name, args };
	}

	/**
	 * Add a certain type with the specified constructor arguments to the factory
	 */
	template<typename U, typename ...Args>
	factory add(std::string name) {
		/* Create factory function */
		std::pair<std::string, factoryfun> p = {
			name,
			[](const std::vector<std::string>& v) -> unique_ptr<T> {
				const auto vp = argument_tuple_converter<Args...>()(v);
				return make_new_from_tuple<U>(vp);
			}
		};
		products.insert(p);
		
		/* Create description */
		vector<string> ds;
		argument_description_converter<Args...>()(ds);
		stringstream d;
		d << name << "(";
		for(std::size_t i = 0; i < ds.size(); ++i)
			d << ((i > 0) ? ", " : "") << ds[i];
		d << ")";
		descriptions.push_back(d.str());
		
		return std::move(*this);
	}
	
	/**
	 * Executes an expression of the form "type(arg1, arg2, ...)" and returns an object instance
	 */
	unique_ptr<T> make(std::string s) {
		auto p = parse(s);
		auto i = products.find(p.first);
		if (i == products.end())
			throw std::runtime_error(p.first + " does not exist");
		return i->second(p.second);
	}
};

/*
 * Time distribution arguments and factors
 */

/* Arguments of time distributions */

DECLARE_ARGUMENT(pinf, double, 0.0);
DECLARE_ARGUMENT(lambda, double, nullopt);
DECLARE_ARGUMENT(mean, double, nullopt);
DECLARE_ARGUMENT(variance, double, nullopt);
DECLARE_ARGUMENT(shape, double, nullopt);
DECLARE_ARGUMENT(scale, double, nullopt);
DECLARE_ARGUMENT(tau, double, nullopt);

/* Time distribution factory */

auto time_factory = factory<transmission_time>()
	.add<transmission_time_lognormal, mean, variance, pinf>("lognormal")
	.add<transmission_time_gamma, mean, variance, pinf>("gamma")
	.add<transmission_time_exponential, lambda>("exponential")
	.add<transmission_time_weibull, shape, scale, pinf>("weibull")
	.add<transmission_time_deterministic, tau>("deterministic");

/*
 * Network arguments and factory
 */

/**
  * Implicit argument used to pass the RNG to network constructors
 */
struct rng {
	static const char* name;
	typedef reference_wrapper<rng_t> value_type;
};
const char* rng::name = nullptr;

template<>
pair<reference_wrapper<rng_t>, int> argument<rng>(const vector<string>& vs, size_t i) {
	return { engine, 0 };
}

template<>
string description<rng>() { return ""; }

/* Arguments of networks */

DECLARE_ARGUMENT(size, node_t, nullopt);
DECLARE_ARGUMENT(avg_degree, double, nullopt);
DECLARE_ARGUMENT(reduced_root_degree, bool, true);
DECLARE_ARGUMENT(m, int, nullopt);
DECLARE_ARGUMENT(k, int, nullopt);
DECLARE_ARGUMENT(p, int, nullopt);

/* Network factory */

auto network_factory = factory<network>()
	.add<fully_connected, size, rng>("fullyconnected")
	.add<acyclic, avg_degree, reduced_root_degree, rng>("acyclic")
	.add<erdos_reyni, size, avg_degree, rng>("erdos-reyni")
	.add<watts_strogatz, size, k, p, rng>("watts-strogatz")
	.add<barabasi_albert, size, rng, m>("barabasi-albert");

}

/*
 * Main method
 */

int program_simulate(int argc, const char * argv[]) {
	OptionParser op("options");
	auto help_opt   = op.add<Switch>("h", "help", "produce help message");
	auto psi_opt = op.add<Value<std::string>>("p", "transmission", "transmission time (psi)");
	auto rho_opt = op.add<Value<std::string>>("r", "recovery", "recovery time (rho)");
	auto nw_opt = op.add<Value<std::string>>("g", "network", "the network to simulat on");
	auto initial_opt = op.add<Implicit<node_t>>("i", "initial", "initial infected node", -1);
	auto tmax_opt = op.add<Value<double>>("T", "Tmax", "stop simulation at this time");
	auto list_times_opt = op.add<Switch>("", "list-times", "list distributions");
	auto list_networks_opt = op.add<Switch>("", "list-networks", "list network types");
	op.parse(argc, argv);
	
	if (help_opt->is_set()) {
		cout << op << endl;
		return 0;
	}
	
	if (list_times_opt->is_set()) {
		cout << "time distributions\n";
		cout << "------------------\n";
		for(const string& s: factories::time_factory.descriptions)
			cout << s << "\n";
	}
	
	if (list_networks_opt->is_set()) {
		cout << "networks\n";
		cout << "--------\n";
		for(const string& s: factories::network_factory.descriptions)
			cout << s << "\n";
	}
	
	unique_ptr<transmission_time> psi;
	if (psi_opt->is_set())
		psi = factories::time_factory.make(psi_opt->value());

	unique_ptr<transmission_time> rho;
	if (rho_opt->is_set())
		rho = factories::time_factory.make(rho_opt->value());

	unique_ptr<network> nw;
	if (nw_opt->is_set())
		nw = factories::network_factory.make(nw_opt->value());

	return 0;
}
