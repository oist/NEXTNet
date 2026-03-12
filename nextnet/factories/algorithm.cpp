//
//  algorithm_factory.cpp
//  NEXTNet
//
//  Created by Florian G. Pflug on 08.04.25.
//

#include "nextnet/stdafx.h"
#include "nextnet/algorithm.h"
#include "nextnet/REGIR.h"
#include "nextnet/nMGA.h"
#include "nextnet/NextReaction.h"
#include "nextnet/factories/factories.h"

namespace factories {

algorithm_implementation<simulate_next_reaction> algorithm_next = algorithm_implementation<simulate_next_reaction>("next")
    .param(&simulate_next_reaction::params::SIR, "SIR")
    .param(&simulate_next_reaction::params::exact_reinfection, "exact_reinfection")
    .param(&simulate_next_reaction::params::edges_concurrent, "edges_concurrent")
    .param(&simulate_next_reaction::params::shuffle_neighbours, "shuffle_neighbours");

algorithm_implementation<simulate_nmga> algorithm_nmga = algorithm_implementation<simulate_nmga>("nmga")
    .param(&simulate_nmga::params::SIR, "SIR")
    .param(&simulate_nmga::params::maximal_dt, "maximal_dt")
    .param(&simulate_nmga::params::approximation_threshold, "approximation_threshold")
    .param(&simulate_nmga::params::tau_precision, "tau_precision");

algorithm_implementation<simulate_regir> algorithm_regir = algorithm_implementation<simulate_regir>("regir")
    .param(&simulate_regir::params::SIR, "SIR")
    .param(&simulate_regir::params::approximation_threshold, "approximation_threshold")
    .param(&simulate_regir::params::tau_precision, "tau_precision");

std::unordered_map<std::string, algorithm &> algorithms = {
    { algorithm_next.name, algorithm_next },
    { algorithm_nmga.name, algorithm_nmga },
    { algorithm_regir.name, algorithm_regir }
};

} /* namespace factories */
