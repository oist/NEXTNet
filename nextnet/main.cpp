#include <filesystem>

#include "popl.hpp"

#include "nextnet/stdafx.h"
#include "nextnet/random.h"
#include "nextnet/REGIR.h"
#include "nextnet/nMGA.h"
#include "nextnet/NextReaction.h"
#include "nextnet/network.h"
#include "nextnet/weighted_network.h"
#include "nextnet/temporal_network.h"
#include "nextnet/brownian_proximity_network.h"

using namespace std;
using namespace std::literals;
using namespace popl;

rng_t engine;

/*
 * Errors
 */

struct program_argument_error : public runtime_error
{
    program_argument_error(const string &arg, const string &msg)
        : runtime_error(msg)
        , invalid_argument(arg)
    {
    }

    string invalid_argument;
};

/*
 * Factories for time distribution and network objects
 */

namespace factories {

struct factory_error : public runtime_error
{
    factory_error(const string &msg)
        : runtime_error(msg)
    {
    }
};

#define DECLARE_ARGUMENT_5(_name, _type, _defaultvalue, _parser, _renderer) \
    struct _name                                                            \
    {                                                                       \
        static const std::string name;                                      \
        typedef _type value_type;                                           \
        typedef std::function<value_type (const std::string&)> parser_type; \
        typedef std::function<std::string (const value_type&)> render_type; \
        static const parser_type parser;                                    \
        static const render_type renderer;                                  \
        static const optional<value_type> defaultvalue;                     \
    };                                                                      \
    const std::string _name::name                         = #_name;         \
    const optional<_name::value_type> _name::defaultvalue = _defaultvalue;  \
    const _name::parser_type _name::parser                = _parser;        \
    const _name::render_type _name::renderer              = _renderer;

#define DECLARE_ARGUMENT_3(_name, _type, _defaultvalue)             \
    DECLARE_ARGUMENT_5(_name, _type, _defaultvalue,                 \
                       (&boost::lexical_cast<_type, std::string>),  \
                       (&boost::lexical_cast<std::string, _type>))

/**
 * Convert argument to correct type for the specified argument Arg
 */
template <typename Arg>
pair<typename Arg::value_type, int> argument(const vector<string> &vs, size_t i)
{
    if (i < vs.size()) {
        try {
            return { Arg::parser(vs[i]), 1 };
        } catch (const std::exception& e) {
            throw factory_error("invalid value '"s + vs[i] + "' for argument " + Arg::name +
                                "\nreason: " + e.what());
        }
    }
    else if (Arg::defaultvalue)
        return { *Arg::defaultvalue, 0 };
    else
        throw factory_error("missing value for argument "s + Arg::name);
}

/**
 * Generate description for the argument Arg
 */
template <typename Arg>
string description()
{
    if (Arg::defaultvalue)
        return string(Arg::name) + " = " + Arg::renderer(*Arg::defaultvalue);
    else
        return string(Arg::name);
}

/**
 * Convert a vector or strings into a tuple of arguments as specified by the argument pack
 */
template <typename Arg, typename... Args>
struct argument_tuple_converter
{
    typedef tuple<typename Arg::value_type, typename Args::value_type...> value_tuple_type;

    value_tuple_type operator()(const vector<string> &vs, std::size_t i = 0)
    {
        pair<typename Arg::value_type, int> a = argument<Arg>(vs, i);
        return tuple_cat(make_tuple(a.first), argument_tuple_converter<Args...>()(vs, i + a.second));
    }
};
template <typename Arg>
struct argument_tuple_converter<Arg>
{
    typedef tuple<typename Arg::value_type> value_tuple_type;

    value_tuple_type operator()(const vector<string> &vs, std::size_t i = 0)
    {
        return make_tuple(argument<Arg>(vs, i).first);
    }
};

/**
 * Fill vector with descriptions of the arguments in the pack
 */
template <typename Arg, typename... Args>
struct argument_description_converter
{
    void operator()(vector<string> &v)
    {
        string d = description<Arg>();
        if (!d.empty())
            v.push_back(d);
        argument_description_converter<Args...>()(v);
    }
};
template <typename Arg>
struct argument_description_converter<Arg>
{
    void operator()(vector<string> &v)
    {
        string d = description<Arg>();
        if (!d.empty())
            v.push_back(d);
    }
};

template <class T, class Tuple, std::size_t... I>
constexpr std::unique_ptr<T> make_new_from_tuple_impl(Tuple &&v, std::index_sequence<I...>)
{
    return make_unique<T>(get<I>(std::forward<Tuple>(v))...);
}

/**
 * Version of the standard libraries make_from_tuple that allocates on the heap
 */
template <class T, class Tuple>
constexpr std::unique_ptr<T> make_new_from_tuple(Tuple &&v)
{
    return make_new_from_tuple_impl<T>(std::forward<Tuple>(v),
                                       make_index_sequence<tuple_size_v<remove_reference_t<Tuple>>>{});
}

/**
 * Factory that generated objects which derive from the common type T from strings of the form
 *   typename(arg1, arg2, ...., argN)
 */
template <typename T>
struct factory
{
    /**
     * Represents a factory function that produces type T given a list of arguments in string form
     */
    typedef function<unique_ptr<T>(const vector<string> &)> factoryfun;

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
        std::vector<char> bstack;
        bool qinclude = false;
        std::vector<string> args;
\
        enum { WS1,
               NAME,
               WS2,
               BRA,
               WS3,
               ARG,
               ARG_DQ,
               ARG_DQE,
               WS4,
               KET_OR_COMMA,
               WS5,
               DONE } state = WS1;
        for (std::ptrdiff_t i = 0; i < (ptrdiff_t)s.size(); ++i) {
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
                        case WS1:
                            state = NAME;
                            --i;
                            break;
                        case WS2:
                            state = BRA;
                            --i;
                            break;
                        case WS3:
                            state = ARG;
                            --i;
                            args.push_back(std::string());
                            break;
                        case WS4:
                            state = ARG;
                            --i;
                            break;
                        case WS5:
                            state = DONE;
                            --i;
                            break;
                        default:
                            throw std::logic_error("invalid state");
                    }
                    break;
                case NAME:
                    /* collect until non-alphanumeric character */
                    if (isalnum(c) || (c == '_') || (c == '-')) {
                        name.push_back(c);
                        continue;
                    }
                    state = WS2;
                    --i;
                    break;
                case BRA:
                    /* error if not opening bracket */
                    if (c != '(')
                        throw factory_error("unable to parse '" + s + "', " +
                                            "invalid symbol '" + c + "''");
                    state = WS3;
                    break;
                case ARG:
                    /* collect until space, comma or bracket */
                    if (isspace(c)) {
                        state = WS4;
                        --i;
                        break;
                    }
                    switch (c) {
                        case '(': bstack.push_back(')'); goto bopen;
                        case '{': bstack.push_back('}'); goto bopen;
                        bopen:
                            args.back().push_back(c);
                            break;
                        case ',':
                        case ')':
                        case '}':
                            if (bstack.empty()) {
                                state = KET_OR_COMMA;
                                --i;
                                break;
                            }
                            /* in a braced subexpression, braces must match, argument continues */
                            if (c != ',') {
                                if (c != bstack.back())
                                    throw factory_error("expected '"s + bstack.back() + "' but found '" + c + "'");
                                bstack.pop_back();
                            }
                            args.back().push_back(c);
                            continue;
                        case '"':
                            qinclude = !args.back().empty();
                            if (qinclude)
                                args.back().push_back(c);
                            state = ARG_DQ;
                            break;
                        default:
                            args.back().push_back(c);
                            continue;
                    }
                    break;
                case ARG_DQ:
                    /* scan until closing quote, then update state */
                    switch (c) {
                        case '"':
                            if (qinclude)
                                args.back().push_back(c);
                            state = ARG;
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
                        case ')':
                            state = WS5;
                            break;
                        case ',':
                            state = WS3;
                            break;
                        default:
                            throw factory_error("unable to parse '" + s + "', " +
                                                "invalid symbol '" + c + "''");
                    }
                    break;
                case DONE:
                    throw factory_error("unable to parse '" + s + "', " +
                                        "invalid symbol '" + c + "''");
            }
        }

        if (state != WS5)
            throw factory_error("unable to parse '" + s + "', " +
                                "incomplete expression");

        return { name, args };
    }

    /**
     * Add a certain type with the specified constructor arguments to the factory
     */
    template <typename U, typename... Args>
    factory add(std::string name)
    {
        /* Create factory function */
        std::pair<std::string, factoryfun> p = {
            name,
            [](const std::vector<std::string> &v) -> unique_ptr<T> {
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
        for (std::size_t i = 0; i < ds.size(); ++i)
            d << ((i > 0) ? ", " : "") << ds[i];
        d << ")";
        descriptions.push_back(d.str());

        return std::move(*this);
    }

    /**
     * Executes an expression of the form "type(arg1, arg2, ...)" and returns an object instance
     */
    unique_ptr<T> make(std::string s)
    {
        try {
            auto p = parse(s);
            auto i = products.find(p.first);
            if (i == products.end())
                throw factory_error(p.first + " does not exist");
            return i->second(p.second);
        } catch (const factory_error &e) {
            throw factory_error("unable to parse: "s + s + "\nreason: " + e.what());
        }
    }
};

/**
 * Base class for alogrithm factories
 */
struct algorithm
{
    typedef pair<string, string> param_t;

    virtual unique_ptr<simulation_algorithm> create(network &nw, transmission_time &psi, transmission_time *rho,
                                                    const vector<param_t> &ps) = 0;
};

/**
 * Factory for generating algorithm objects of type T
 */
template <typename T>
struct algorithm_implementation : public algorithm
{
    using algorithm::param_t;
    typedef T algorithm_type;
    typedef typename algorithm_type::params algorithm_params_type;
    typedef function<void(algorithm_params_type &p, string value)> setter_function_type;

    static algorithm_params_type default_params;

    unordered_map<string, setter_function_type> setters;

    /**
     * Add algorithm parameters
     */
    template <typename U>
    algorithm_implementation param(U algorithm_params_type::*param, string name)
    {
        setters.insert({ name,
                         [param](algorithm_params_type &p, string value) {
                             p.*param = boost::lexical_cast<U>(value);
                         } });
        return std::move(*this);
    }

    /**
     * Create algorithm instance
     */
    virtual unique_ptr<simulation_algorithm> create(network &nw, transmission_time &psi, transmission_time *rho,
                                                    const vector<param_t> &ps)
    {
        algorithm_params_type p;
        for (const param_t &pv : ps) {
            /* parameter name and value */
            const std::string &pname  = pv.first;
            const std::string &pvalue = pv.second;

            /* lookup setter */
            auto si = setters.find(pname);
            if (si == setters.end())
                throw factory_error("unknown parameter "s + pname);
            setter_function_type set = si->second;

            /* set parameter */
            set(p, pvalue);
        }

        return unique_ptr<simulation_algorithm>(new algorithm_type(nw, psi, rho, p));
    }
};

template <typename T>
typename algorithm_implementation<T>::algorithm_params_type algorithm_implementation<T>::default_params;

/*
 * Available time distributions and their arguments
 */

/* Arguments of time distributions */

DECLARE_ARGUMENT_3(pinf, double, 0.0);
DECLARE_ARGUMENT_3(lambda, double, nullopt);
DECLARE_ARGUMENT_3(mean, double, nullopt);
DECLARE_ARGUMENT_3(variance, double, nullopt);
DECLARE_ARGUMENT_3(shape, double, nullopt);
DECLARE_ARGUMENT_3(scale, double, nullopt);
DECLARE_ARGUMENT_3(tau, double, nullopt);

/* Time distribution factory */

auto time_factory = factory<transmission_time>()
                        .add<transmission_time_lognormal, mean, variance, pinf>("lognormal")
                        .add<transmission_time_gamma, mean, variance, pinf>("gamma")
                        .add<transmission_time_exponential, lambda>("exponential")
                        .add<transmission_time_weibull, shape, scale, pinf>("weibull")
                        .add<transmission_time_deterministic, tau>("deterministic");

/*
 * Available networks and their arguments
 */

/**
 * Implicit argument used to pass the RNG to network constructors
 */
struct rng
{
    static const char *name;
    typedef reference_wrapper<rng_t> value_type;
};
const char *rng::name = nullptr;

template <>
pair<reference_wrapper<rng_t>, int> argument<rng>(const vector<string> &vs, size_t i)
{
    return { engine, 0 };
}

template <>
string description<rng>() { return ""; }

/* Arguments of networks */

empirical_contact_network::edge_duration_kind parse_contact_kind(const std::string& s)
{
    if (s == "instantenous")
        return empirical_contact_network::infitesimal_duration;
    else if (s == "finite")
        return empirical_contact_network::finite_duration;
    else
        throw std::range_error("contact_kind must be 'instantenous' or 'finite'");
}

std::string render_contact_kind(const empirical_contact_network::edge_duration_kind& s)
{
    switch (s) {
        case empirical_contact_network::infitesimal_duration:
            return "instantenous";
        case empirical_contact_network::finite_duration:
            return "finite";
        default:
            return "?";
    }
}

template<typename T>
std::vector<T> parse_vector(const std::string& s)
{
    if ((s.size() < 2) || (s.front() != '{') || (s.back() != '}'))
        throw std::range_error("list must be surrounded by '{', '}'");

    std::stringstream in(std::string(++s.begin(), --s.end()));

    std::vector<T> r;
    std::string v;
    while (std::getline(in, v, ',')) {
        try{
            r.push_back(boost::lexical_cast<T>(v));
        }
        catch (const boost::bad_lexical_cast& e) {
            throw std::range_error("invalid list element '"s + v + "'");
        }
    }
    return r;
}

template<typename T>
std::string render_vector(const std::vector<T> v)
{
    std::stringstream s;
    std::size_t i = 0;
    s << '{';
    for(const T& e: v)
        (i++ ? s : s << ',') << e;
    s << '}';
    return s.str();
}

DECLARE_ARGUMENT_3(file, std::filesystem::path, nullopt);
DECLARE_ARGUMENT_5(contact_kind, empirical_contact_network::edge_duration_kind,
                   empirical_contact_network::infitesimal_duration, parse_contact_kind, render_contact_kind);
DECLARE_ARGUMENT_3(dt, double, 1.0);
DECLARE_ARGUMENT_3(size, node_t, nullopt);
DECLARE_ARGUMENT_3(avg_degree, double, nullopt);
DECLARE_ARGUMENT_5(weights, std::vector<double>, std::vector<double> { 1.0 },
                   parse_vector<double>, render_vector<double>);
DECLARE_ARGUMENT_5(probabilities, std::vector<double>, std::vector<double> { 1.0 },
                   parse_vector<double>, render_vector<double>);
DECLARE_ARGUMENT_3(timescale, double, 1.0);
DECLARE_ARGUMENT_5(activities, std::vector<double>, std::vector<double> { 1.0 },
                   parse_vector<double>, render_vector<double>);
DECLARE_ARGUMENT_3(eta_sus, double, nullopt);
DECLARE_ARGUMENT_3(eta_inf, double, nullopt);
DECLARE_ARGUMENT_3(b_sus, double, nullopt);
DECLARE_ARGUMENT_3(b_inf, double, nullopt);
DECLARE_ARGUMENT_3(k, int, nullopt);
DECLARE_ARGUMENT_3(p, int, nullopt);
DECLARE_ARGUMENT_3(m, int, nullopt);
DECLARE_ARGUMENT_5(degrees, std::vector<int>, nullopt, parse_vector<int>, render_vector<int>);
DECLARE_ARGUMENT_5(triangles, std::vector<int>, nullopt, parse_vector<int>, render_vector<int>);
DECLARE_ARGUMENT_3(beta, double, nullopt);
DECLARE_ARGUMENT_3(edgelength, int, nullopt);
DECLARE_ARGUMENT_3(radius, double, 1.0);
DECLARE_ARGUMENT_3(D0, double, 1.0);
DECLARE_ARGUMENT_3(D1, double, 1.0);
DECLARE_ARGUMENT_3(gamma, double, 1.0);
DECLARE_ARGUMENT_3(reduced_root_degree, bool, true);

/* Network factory */

auto network_factory = factory<network>()
   .add<empirical_network, file>("empirical")
   .add<empirical_contact_network, file, contact_kind, dt>("empirical-contact")
   .add<erdos_renyi, size, avg_degree, rng>("erdos-renyi")
   .add<weighted_erdos_renyi, size, avg_degree, weights, probabilities, rng>("weighted-erdos-renyi")
   .add<temporal_erdos_renyi, size, avg_degree, timescale, rng>("temporal-erdos-renyi")
   .add<activity_driven_network, activities, m, eta_sus, eta_inf, b_sus, b_inf, rng>("activity-driven")
   .add<watts_strogatz, size, k, p, rng>("watts-strogatz")
   .add<barabasi_albert, size, rng, m>("barabasi-albert")
   .add<config_model, degrees, rng>("config")
   .add<config_model_clustered_serrano, degrees, triangles, beta, rng>("config-clustered")
   .add<cubic_lattice_2d, edgelength>("lattice_2d")
   .add<cubic_lattice_3d, edgelength>("lattice_3d")
   .add<brownian_proximity_network, size, avg_degree, radius, D0, D1, gamma, rng>("brownian-proximity")
   .add<acyclic, avg_degree, reduced_root_degree, rng>("acyclic")
   .add<fully_connected, size, rng>("fullyconnected");

/*
 * Available algorithms and their parameters
 */

auto algorithm_next = algorithm_implementation<simulate_next_reaction>()
                          .param(&simulate_next_reaction::params::SIR, "SIR")
                          .param(&simulate_next_reaction::params::edges_concurrent, "edges_concurrent")
                          .param(&simulate_next_reaction::params::shuffle_neighbours, "shuffle_neighbours");

auto algorithm_nmga = algorithm_implementation<simulate_nmga>()
                          .param(&simulate_nmga::params::SIR, "SIR")
                          .param(&simulate_nmga::params::maximal_dt, "maximal_dt")
                          .param(&simulate_nmga::params::approximation_threshold, "approximation_threshold")
                          .param(&simulate_nmga::params::tau_precision, "tau_precision");

auto algorithm_regir = algorithm_implementation<simulate_regir>()
                           .param(&simulate_regir::params::SIR, "SIR")
                           .param(&simulate_regir::params::approximation_threshold, "approximation_threshold")
                           .param(&simulate_regir::params::tau_precision, "tau_precision");

unordered_map<string, algorithm &> algorithms = {
    { "next"s, algorithm_next },
    { "nmga"s, algorithm_nmga },
    { "regir"s, algorithm_regir }
};

} // namespace factories

/*
 * Main method
 */

int main(int argc, const char *argv[])
{
    unique_ptr<transmission_time> psi;
    unique_ptr<transmission_time> rho;
    unique_ptr<network> nw;
    vector<factories::algorithm::param_t> alg_params;
    unique_ptr<simulation_algorithm> alg;
    unique_ptr<simulate_on_temporal_network> sotn_alg;

    OptionParser op("options");
    auto help_opt          = op.add<Switch>("h", "help", "produce help message");
    auto psi_opt           = op.add<Value<std::string>>("p", "transmission-time", "transmission time (psi)");
    auto rho_opt           = op.add<Value<std::string>>("r", "recovery-time", "recovery time (rho)");
    auto nw_opt            = op.add<Value<std::string>>("n", "network", "network to simulat on");
    auto alg_opt           = op.add<Value<std::string>>("a", "algorithm", "simulation algorithm to use", "next");
    auto param_opt         = op.add<Value<std::string>>("s", "parameter", "set simulation parameter");
    auto initial_opt       = op.add<Value<node_t>>("i", "initial-infection", "initial infected node");
    auto ev_opt            = op.add<Value<bool>>("e", "epidemic-events", "output epidemic events", true);
    auto nv_opt            = op.add<Value<bool>>("w", "network-events", "output network events", false);
    auto tmax_opt          = op.add<Value<double>>("t", "stopping-time", "stop simulation at this time");
    auto list_times_opt    = op.add<Switch>("", "list-times", "list distributions");
    auto list_networks_opt = op.add<Switch>("", "list-networks", "list network types");

    try {
        /* Parse arguments */
        op.parse(argc, argv);

        /* Output help if requested, then stop */
        if (help_opt->is_set()) {
            cout << op << endl;
            return 0;
        }

        /* Output list of time distributions / networks if requested, then stop */

        if (list_times_opt->is_set()) {
            cout << "time distributions\n";
            cout << "------------------\n";
            for (const string &s : factories::time_factory.descriptions)
                cout << s << "\n";
            return 0;
        }

        if (list_networks_opt->is_set()) {
            cout << "networks\n";
            cout << "--------\n";
            for (const string &s : factories::network_factory.descriptions)
                cout << s << "\n";
            return 0;
        }

        /* Run simulation */

        if (psi_opt->is_set())
            psi = factories::time_factory.make(psi_opt->value());

        if (rho_opt->is_set())
            rho = factories::time_factory.make(rho_opt->value());

        if (nw_opt->is_set())
            nw = factories::network_factory.make(nw_opt->value());

        for (size_t i = 0; i < param_opt->count(); ++i) {
            const std::string &p = param_opt->value(i);
            std::size_t j        = p.find('=');
            if (j == p.npos)
                throw program_argument_error("parameter", "invalid parameter setting '"s + p + "', does not contain '='");
            const std::string pname  = p.substr(0, j);
            const std::string pvalue = p.substr(j + 1, p.npos);
            alg_params.push_back({ pname, pvalue });
        }

        {
            auto a_i = factories::algorithms.find(alg_opt->value());
            if (a_i == factories::algorithms.end())
                throw program_argument_error("algorithm", "unknown algorithm "s + alg_opt->value());
            factories::algorithm &alg_factory = a_i->second;

            if (!nw)
                throw program_argument_error("network", "no network specified");

            if (!psi)
                throw program_argument_error("transmission-time", "no transmission time distribution specified");

            alg = alg_factory.create(*nw.get(), *psi.get(), rho.get(), alg_params);
        }

        for (size_t i = 0; i < initial_opt->count(); ++i) {
            const node_t node = initial_opt->value(i) - 1;
            if ((node < 0) || (node >= nw->nodes()))
                throw program_argument_error("initial-infection", "invalid initial node "s + boost::lexical_cast<string>(initial_opt->value(i)));
            alg->add_infections({ { node, 0.0 } });
        }
        if (initial_opt->count() == 0)
            std::cerr << "WARNING: No initially infected node specified with -initial-infection, no epidemic will commence" << std::endl;

        unique_ptr<simulate_on_temporal_network> sotn_alg;
        if (dynamic_cast<temporal_network *>(nw.get()) != nullptr)
            sotn_alg.reset(new simulate_on_temporal_network(*alg.get()));
    } catch (program_argument_error &e) {
        cerr << op << endl;
        cerr << "Error in argument " << e.invalid_argument << ": " << e.what() << endl;
        return 1;
    } catch (runtime_error &e) {
        cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (exception &e) {
        cerr << "Internal error: " << e.what() << std::endl;
        return 127;
    }

    try {
        cout << "time" << '\t';
        cout << "epidemic_step" << '\t';
        cout << "network_step" << '\t';
        cout << "kind" << '\t';
        cout << "node" << '\t';
        cout << "neighbour" << '\t';
        cout << "total_infected" << '\t';
        cout << "total_reset" << '\t';
        cout << "infected" << '\n';

        const bool epidemic_events = ev_opt->value();
        const bool network_events  = nv_opt->value();
        const double tmax          = tmax_opt->is_set() ? tmax_opt->value() : INFINITY;
        size_t epidemic_step       = 0;
        size_t network_step        = 0;
        double infected            = 0;
        double total_infected      = 0;
        double total_reset         = 0;
        while (true) {
            // Execute next event
            const std::optional<network_or_epidemic_event_t> any_ev_opt =
                sotn_alg ? sotn_alg->step(engine, tmax) : alg->step(engine, tmax);

            // Stop if there are no more events
            if (!any_ev_opt)
                break;
            const network_or_epidemic_event_t any_ev = *any_ev_opt;

            // Fill columns
            double time;
            const char *kind;
            node_t node;
            node_t neighbour = -1;
            if (std::holds_alternative<epidemic_event_t>(any_ev)) {
                // Epidemic event
                const auto ev = std::get<epidemic_event_t>(any_ev);
                // Update state
                time = ev.time;
                ++epidemic_step;
                switch (ev.kind) {
                    case epidemic_event_kind::outside_infection:
                    case epidemic_event_kind::infection:
                        ++total_infected;
                        ++infected;
                        break;
                    case epidemic_event_kind::reset:
                        ++total_reset;
                        --infected;
                        break;
                    default:
                        break;
                }
                // Report event?
                if (!epidemic_events)
                    continue;
                // Fill row
                kind = name(ev.kind);
                node = ev.node + 1;
                if (ev.kind == epidemic_event_kind::infection)
                    neighbour = ev.source_node + 1;
            } else if (std::holds_alternative<network_event_t>(any_ev)) {
                // Network event
                const auto ev = std::get<network_event_t>(any_ev);
                // Update state
                time = ev.time;
                ++network_step;
                // Report event?
                if (!network_events)
                    continue;
                // Fill row
                kind      = name(ev.kind);
                node      = ev.source_node + 1;
                neighbour = ev.target_node + 1;
            } else
                throw logic_error("unknown event type");

            // Output
            cout << time << '\t';
            cout << epidemic_step << '\t';
            cout << network_step << '\t';
            cout << kind << '\t';
            cout << node << '\t';
            cout << neighbour << '\t';
            cout << total_infected << '\t';
            cout << total_reset << '\t';
            cout << infected << '\n';
        }

        return 0;
    } catch (exception &e) {
        cerr << "Internal error: " << e.what() << std::endl;
        return 127;
    }
}
