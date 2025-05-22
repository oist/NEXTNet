//
//  algorithm_factory.h
//  NEXTNet
//
//  Created by Florian G. Pflug on 08.04.25.
//

#pragma once

#include "nextnet/algorithm.h"

namespace factories {

/**
 * Base class for alogrithm factories
 */
struct algorithm
{
    typedef std::pair<std::string, std::string> param_t;

    virtual std::unique_ptr<simulation_algorithm> create(network &nw, transmission_time &psi, transmission_time *rho,
                                                         const std::vector<param_t> &ps) = 0;
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
    typedef std::function<void(algorithm_params_type &p, std::string value)> setter_function_type;

    static algorithm_params_type default_params;

    std::unordered_map<std::string, setter_function_type> setters;

    /**
     * Add algorithm parameters
     */
    template <typename U>
    algorithm_implementation param(U algorithm_params_type::*param, std::string name)
    {
        setters.insert({ name,
                         [param](algorithm_params_type &p, std::string value) {
                             p.*param = boost::lexical_cast<U>(value);
                         } });
        return std::move(*this);
    }

    /**
     * Create algorithm instance
     */
    virtual std::unique_ptr<simulation_algorithm> create(network &nw, transmission_time &psi,
                                                         transmission_time *rho,
                                                         const std::vector<param_t> &ps)
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

        return std::unique_ptr<simulation_algorithm>(new algorithm_type(nw, psi, rho, p));
    }
};

template <typename T>
typename algorithm_implementation<T>::algorithm_params_type algorithm_implementation<T>::default_params;

extern std::unordered_map<std::string, algorithm &> algorithms;

} /* namespace factories */
