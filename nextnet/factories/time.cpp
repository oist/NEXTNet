//
//  time.cpp
//  NEXTNet
//
//  Created by Florian G. Pflug on 08.04.25.
//

#include "nextnet/stdafx.h"
#include "nextnet/random.h"
#include "nextnet/factories/factories.h"

namespace factories {

/* Arguments of time distributions */

DECLARE_ARGUMENT_3(pinf, double, 0.0);
DECLARE_ARGUMENT_3(lambda, double, std::nullopt);
DECLARE_ARGUMENT_3(mean, double, std::nullopt);
DECLARE_ARGUMENT_3(variance, double, std::nullopt);
DECLARE_ARGUMENT_3(shape, double, std::nullopt);
DECLARE_ARGUMENT_3(scale, double, std::nullopt);
DECLARE_ARGUMENT_3(tau, double, std::nullopt);

/* Time distribution factory */

factory<transmission_time> time_factory = factory<transmission_time>()
	.add<transmission_time_lognormal, mean, variance, pinf>("lognormal")
	.add<transmission_time_gamma, mean, variance, pinf>("gamma")
	.add<transmission_time_exponential, lambda>("exponential")
	.add<transmission_time_weibull, shape, scale, pinf>("weibull")
	.add<transmission_time_deterministic, tau>("deterministic");

} /* namespace factories */
