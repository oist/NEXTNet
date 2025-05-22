//
//  factories.h
//  NEXTNet
//
//  Created by Florian G. Pflug on 08.04.25.
//

#pragma once

#include "nextnet/stdafx.h"
#include "nextnet/types.h"

namespace factories {

using namespace std::literals;

struct factory_error : public std::runtime_error
{
    factory_error(const std::string &msg)
        : runtime_error(msg)
    {
    }
};

/**
 * Stores expressions of the form function(arg1,...,argN)
 */
typedef std::pair<std::string, std::vector<std::string>> parsed_expression_t;

/**
 * Parses a string of the form "function(arg1, arg2, "quoted\ arg")" with arbitarry number of arguments
 * into a pair comprising the function name and a vector of arguments.
 */
extern parsed_expression_t parse_expression(std::string s);

std::ostream &operator<<(std::ostream &o, const parsed_expression_t &expr);

#define DECLARE_ARGUMENT_5(_name, _type, _dfl, _parser, _renderer)          \
    struct _name                                                            \
    {                                                                       \
        static const std::string name;                                      \
        typedef _type value_type;                                           \
        typedef std::function<value_type(const std::string &)> parser_type; \
        typedef std::function<std::string(const value_type &)> render_type; \
        static const std::optional<value_type> defaultvalue;                \
        static const bool implicit;                                         \
        static const parser_type parser;                                    \
        static const render_type renderer;                                  \
    };                                                                      \
    const std::string _name::name                              = #_name;    \
    const std::optional<_name::value_type> _name::defaultvalue = _dfl;      \
    const bool _name::implicit                                 = false;     \
    const _name::parser_type _name::parser                     = _parser;   \
    const _name::render_type _name::renderer                   = _renderer;

#define DECLARE_ARGUMENT_3(_name, _type, _dfl)                     \
    DECLARE_ARGUMENT_5(_name, _type, _dfl,                         \
                       (&boost::lexical_cast<_type, std::string>), \
                       (&boost::lexical_cast<std::string, _type>))

/**
 * Convert argument to correct type for the specified argument Arg
 */
template <typename Arg>
std::pair<typename Arg::value_type, int> argument(const std::vector<std::string> &vs, size_t i)
{
    if (i < vs.size()) {
        try {
            return { Arg::parser(vs[i]), 1 };
        } catch (const std::exception &e) {
            throw factory_error("invalid value '"s + vs[i] + "' for argument " + Arg::name +
                                "\nreason: " + e.what());
        }
    } else if (Arg::defaultvalue)
        return { *Arg::defaultvalue, 0 };
    else
        throw factory_error("missing value for argument "s + Arg::name);
}

/**
 * Generate description for the argument Arg
 */
template <typename Arg>
std::string description()
{
    if (Arg::defaultvalue)
        return std::string(Arg::name) + " = " + Arg::renderer(*Arg::defaultvalue);
    else
        return std::string(Arg::name);
}

/**
 * Convert a vector or strings into a tuple of arguments as specified by the argument pack
 */
template <typename Arg, typename... Args>
struct argument_tuple_converter
{
    typedef std::tuple<typename Arg::value_type, typename Args::value_type...> value_tuple_type;

    value_tuple_type operator()(const std::vector<std::string> &vs, std::size_t i = 0)
    {
        std::pair<typename Arg::value_type, int> a = argument<Arg>(vs, i);
        return std::tuple_cat(std::make_tuple(a.first), argument_tuple_converter<Args...>()(vs, i + a.second));
    }
};
template <typename Arg>
struct argument_tuple_converter<Arg>
{
    typedef std::tuple<typename Arg::value_type> value_tuple_type;

    value_tuple_type operator()(const std::vector<std::string> &vs, std::size_t i = 0)
    {
        return std::make_tuple(argument<Arg>(vs, i).first);
    }
};

/**
 * Convert of a tuple of arguments into a vector of string values
 */
template <typename Arg, typename... Args>
struct argument_vector_converter
{
    template <typename Tuple, std::size_t I = 0>
    void operator()(bool named, const Tuple &vs, std::vector<std::string> &dst)
    {
        if (!Arg::implicit and named)
            dst.push_back(std::string(Arg::name) + " = " + Arg::renderer(std::get<I>(vs)));
        else if (!Arg::implicit)
            dst.push_back(Arg::renderer(std::get<I>(vs)));
        argument_vector_converter<Args...>().template operator()<Tuple, I + 1>(named, vs, dst);
    }
};
template <typename Arg>
struct argument_vector_converter<Arg>
{
    template <typename Tuple, std::size_t I = 0>
    void operator()(bool named, const Tuple &vs, std::vector<std::string> &dst)
    {
        if (!Arg::implicit and named)
            dst.push_back(std::string(Arg::name) + " = " + Arg::renderer(std::get<I>(vs)));
        else if (!Arg::implicit)
            dst.push_back(Arg::renderer(std::get<I>(vs)));
    }
};

/**
 * Fill vector with descriptions of the arguments in the pack
 */
template <typename Arg, typename... Args>
struct argument_description_converter
{
    void operator()(std::vector<std::string> &v)
    {
        if (!Arg::implicit)
            v.push_back(description<Arg>());
        argument_description_converter<Args...>()(v);
    }
};
template <typename Arg>
struct argument_description_converter<Arg>
{
    void operator()(std::vector<std::string> &v)
    {
        if (!Arg::implicit)
            v.push_back(description<Arg>());
    }
};

template <class T, class Tuple, std::size_t... I>
constexpr std::unique_ptr<T> make_new_from_tuple_impl(Tuple &&v, std::index_sequence<I...>)
{
    return std::make_unique<T>(std::get<I>(std::forward<Tuple>(v))...);
}

/**
 * Version of the standard libraries make_from_tuple that allocates on the heap
 */
template <class T, class Tuple>
constexpr std::unique_ptr<T> make_new_from_tuple(Tuple &&v)
{
    return make_new_from_tuple_impl<T>(std::forward<Tuple>(v), std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>{});
}

/**
 * Factory that generated objects which derive from the common type T from strings of the form
 *   typename(arg1, arg2, ...., argN)
 */
template <typename T>
struct factory
{
    /**
     * Unique ptr containing an object plus a unique ptr to an any containing a tuple of
     * of all the argument values. Keeping the any around until the object is destroyed
     * ensures that by-reference survive for the lieftime of the object.
     */
    typedef std::pair<std::unique_ptr<T>, std::unique_ptr<std::any>> object_holder_t;

    /**
     * Represents a factory function that produces type T given a list of arguments in string form
     */
    typedef std::function<object_holder_t(const std::vector<std::string> &)> factoryfun;

    /**
     * Represents a factory function that parses  arguments and outputs them again as strings
     */
    typedef std::function<parsed_expression_t(const std::vector<std::string> &, bool)> parsefun;

    factory(std::string name_)
        : name(name_)
    {
    }

    const std::string name;

    /**
     * Map of types names to factory functions
     */
    std::unordered_map<std::string, factoryfun> products;

    /**
     * Map of types names to parse functions
     */
    std::unordered_map<std::string, parsefun> parsers;

    /**
     * List of factory descriptions
     */
    std::vector<std::string> descriptions;

    /**
     * Add a certain type with the specified constructor arguments to the factory
     */
    template <typename U, typename... Args>
    factory add(std::string name)
    {
        /* Create factory function */
        std::pair<std::string, factoryfun> ff = {
            name,
            [](const std::vector<std::string> &v) -> object_holder_t {
                const auto vp = argument_tuple_converter<Args...>()(v);
                auto vp_ptr   = std::make_unique<std::any>(std::move(vp));
                return { make_new_from_tuple<U>(std::any_cast<decltype(vp) &>(*vp_ptr)),
                         std::move(vp_ptr) };
            }
        };
        products.insert(ff);

        /* Create factory function */
        std::pair<std::string, parsefun> pf = {
            name,
            [name](const std::vector<std::string> &v, bool named) -> parsed_expression_t {
                const auto vp = argument_tuple_converter<Args...>()(v);
                std::vector<std::string> r;
                argument_vector_converter<Args...>()(named, vp, r);
                return { name, r };
            }
        };
        parsers.insert(pf);

        /* Create description */
        std::vector<std::string> ds;
        argument_description_converter<Args...>()(ds);
        std::stringstream d;
        d << name << "(";
        for (std::size_t i = 0; i < ds.size(); ++i)
            d << ((i > 0) ? ", " : "") << ds[i];
        d << ")";
        descriptions.push_back(d.str());

        return std::move(*this);
    }

    /**
     * Executes an expression of the form "type(arg1, arg2, ...)" and returns an object instance
     *
     * The return value is a pair comprising a unique_ptr_ and a *holder*. The holder
     * is a tuple stored in an std::any that holds instances of all parameter values. For
     * most parameters this is unnecessary, but for by-reference parameters keeping
     * this any around with the unique_ptr ensures that these parameters are not
     * destructed before the object is.
     */
    object_holder_t make(std::string s)
    {
        parsed_expression_t p;
        try {
            p = parse_expression(s);
        } catch (const factory_error &e) {
            throw factory_error("unable to parse: "s + s + "\nreason: " + e.what());
        }
        try {
            auto i = products.find(p.first);
            if (i == products.end())
                throw factory_error("no "s + name + " named " + p.first);
            return i->second(p.second);
        } catch (const std::runtime_error &e) {
            throw factory_error(name + " " + p.first + ": " + e.what());
        }
    }

    parsed_expression_t parse(std::string s, bool named = true)
    {
        parsed_expression_t p;
        try {
            p = parse_expression(s);
        } catch (const factory_error &e) {
            throw factory_error("unable to parse: "s + s + "\nreason: " + e.what());
        }
        auto i = parsers.find(p.first);
        if (i == parsers.end())
            throw factory_error("no "s + name + " named " + p.first);
        return i->second(p.second, named);
    }
};

/**
 * rng. Implicit argument used to pass the RNG to network constructors
 */

extern rng_t *random_engine;

struct rng
{
    static const char *name;
    typedef std::reference_wrapper<rng_t> value_type;
    static const bool implicit;
    static const std::function<std::string(const value_type &)> renderer;
};

template <>
std::pair<std::reference_wrapper<rng_t>, int> argument<rng>(const std::vector<std::string> &vs, size_t i);

template <>
std::string description<rng>();

/**
 * @brief Argument used to pass an istream that reads from a user-specified file
 */

struct istream_ref
{
    static istream_ref parse(const std::string &path);
    static std::string render(const istream_ref &is);

    operator std::istream &() const { return *file.get(); }

    std::filesystem::path path;
    std::shared_ptr<std::istream> file;
};

/**
 * @brief Render and parse function for vector arguments
 */

template <typename T>
std::vector<T> parse_vector(const std::string &s)
{
    if ((s.size() < 2) || (s.front() != '{') || (s.back() != '}'))
        throw std::range_error("list must be surrounded by '{', '}'");

    std::stringstream in(std::string(++s.begin(), --s.end()));

    std::vector<T> r;
    std::string v;
    while (std::getline(in, v, ',')) {
        try {
            r.push_back(boost::lexical_cast<T>(v));
        } catch (const boost::bad_lexical_cast &e) {
            throw std::range_error("invalid list element '"s + v + "'");
        }
    }
    return r;
}

template <typename T>
std::string render_vector(const std::vector<T> v)
{
    std::stringstream s;
    std::size_t i = 0;
    s << '{';
    for (const T &e : v)
        (i++ ? (s << ',') : s) << e;
    s << '}';
    return s.str();
}

} // namespace factories

#include "nextnet/factories/algorithm.h"
#include "nextnet/factories/time.h"
#include "nextnet/factories/network.h"
