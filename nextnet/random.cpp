//
//  randomGenerators.cpp
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#include "nextnet/stdafx.h"
#include "nextnet/random.h"
#include "nextnet/utility.h"

/*----------------------------------------------------*/
/*----------------------CLASS-------------------------*/
/*----------------TRANSMISSION TIME-------------------*/
/*----------------------------------------------------*/

transmission_time::transmission_time(double pinf)
    : pinfinity(pinf)
{
    if ((pinfinity < 0.0) || (pinfinity > 1.0) || !std::isfinite(pinfinity))
        throw std::range_error("pinfinity must lie within [0, 1]");
}

transmission_time::~transmission_time()
{
    // Nothing to do
}

/* "this ->"" is also used to ensure that if some of the functions are
 * redefined in herited classes then they are they ones being used and
 * not the by-default implementations.*/

interval_t transmission_time::sample(rng_t &rng, interval_t t, double m) const
{
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("condition t must be non-negative and finite");
    if (m < 0.0)
        throw std::range_error("m must be positive");
    if (m == 0.0)
        return INFINITY;
    // By default, sample using the quantile function
    const double u = std::uniform_real_distribution<double>(0, 1)(rng);
    return this->survivalquantile(u, t, m);
}

double transmission_time::density(interval_t tau, interval_t t, double m) const
{
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("condition t must be non-negative and finite");
    if (m < 0.0)
        throw std::range_error("m must be positive");
    // By default, compute the conditional density from the hazard rate
    // psi(tau | t, m) = m * l(t + tau) Psi(tau | t, m).
    // TODO: This is slighly convoluted since the hazard rate is by default
    // defined in terms of psi(tau) and Psi(tau), should we inline this
    // definition here?
    return m * this->hazardrate(t + tau) * this->survivalprobability(tau, t, m);
}

double transmission_time::hazardrate(interval_t tau) const
{
    // By default, compute the hazard rate from psi and Psi as
    // lambda(tau) = psi(tau) / Psi(tau).
    return this->density(tau) / this->survivalprobability(tau);
}

double transmission_time::hazardbound(interval_t tau) const
{
    return INFINITY;
}

double transmission_time::survivalprobability(interval_t tau, interval_t t, double m) const
{
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("condition t must be non-negative and finite");
    if (m < 0.0)
        throw std::range_error("m must be positive");
    // By default compute the conditional survival probability
    // using the unconditional survival function Psi(tau) based on
    // Psi(tau | t, m) = (Psi(t + tau) / Psi(t))^m
    return std::pow((this->survivalprobability(t + tau) / this->survivalprobability(t)), m);
}

interval_t transmission_time::survivalquantile(double u) const
{
    if ((u < 0) || (u > 1) || !std::isfinite(u))
        throw std::range_error("u-quantile undefined for u not in [0,1]");
    // pinfinity() is the lower bound of survivalprobability(), so the u-quantile
    // for any u < pinfinity is infinity.
    if (u < pinfinity)
        return INFINITY;
    // By default, numerically invert the survival function
    return inverse_survival_function(u, 1e-6, [&](double tau) { return this->survivalprobability(tau); });
}

interval_t transmission_time::survivalquantile(double u, interval_t t, double m) const
{
    if ((u < 0) || (u > 1) || !std::isfinite(u))
        throw std::range_error("u-quantile undefined for u not in [0,1]");
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("condition t must be non-negative and finite");
    if (m < 0.0)
        throw std::range_error("m must be positive");
    if (m == 0.0)
        return INFINITY;
    /* We consider i.i.d. tau_1, ..., tau_m and ask for the probability
     *
     *   p(tau) = P[ tau_1, ..., tau_m >= t + tau | tau_1, ..., tau_m >= t ].
     *
     * In terms of the (unconditional) CDF F of tau this yields
     *
     *                                   m
     *            (   1 - F(t + tau)   )
     *   p(tau) = ( ------------------ )
     *            (       1 - F(t)     )
     *
     * and solving p(tau) = u gives
     *
     *  t + tau = F^-1[ 1 - u^(1/m) * (1 - F(t)) ].
     *
     * In terms of the (unconditional) survival function G = 1 - F we get
     *
     *  t + tau = G^-1[ u^(1/m) * G(t) ].
     */
    const double up             = this->survivalprobability(t) * std::pow(u, 1.0 / double(m));
    const interval_t t_plus_tau = this->survivalquantile(up);
    if (t_plus_tau < t)
        throw std::logic_error("encountered invalid result when inverting the survival function");
    return (t_plus_tau - t);
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:EXPONENTIAL------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

double transmission_time_exponential::density(interval_t tau) const
{
    if (tau < 0)
        return 0.0;
    return lambda * exp(-lambda * tau);
}

double transmission_time_exponential::hazardrate(interval_t tau) const
{
    if (tau < 0)
        return 0.0;
    return lambda;
}

double transmission_time_exponential::hazardbound(interval_t) const
{
    return lambda;
}

double transmission_time_exponential::survivalprobability(interval_t tau) const
{
    if (tau < 0)
        return 1.0;
    return exp(-lambda * tau);
}

double transmission_time_exponential::survivalprobability(interval_t tau, interval_t t, double m) const
{
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("condition t must be non-negative and finite");
    if (m < 1)
        throw std::range_error("m must be positive");
    if (tau < 0)
        return 1.0;
    return exp(-lambda * m * tau);
}

double transmission_time_exponential::survivalquantile(double u) const
{
    if ((u < 0) || (u > 1) || !std::isfinite(u))
        throw std::range_error("u-quantile undefined for u not in [0,1]");
    return -mean * log(u);
}

double transmission_time_exponential::survivalquantile(double u, interval_t t, double m) const
{
    if ((u < 0) || (u > 1) || !std::isfinite(u))
        throw std::range_error("u-quantile undefined for u not in [0,1]");
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("condition t must be non-negative and finite");
    if (m < 0.0)
        throw std::range_error("m must be positive");
    if (m == 0.0)
        return INFINITY;
    return -mean / m * log(u);
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:GENERIC------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

double transmission_time_gamma::hazardbound(interval_t) const
{
    return mean / variance;
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:POLYNOMIAL RATE--------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

transmission_time_polynomial_rate::transmission_time_polynomial_rate(std::vector<double> &&_coeffs, int dummy)
    : coeffs(std::move(_coeffs))
{
    bool all_zero = true;
    for (double c : coeffs) {
        if ((c < 0.0) || (!std::isfinite(c)))
            throw std::range_error("coefficients must be non-negative");
        all_zero = all_zero && (c == 0.0);
    }
    if (all_zero)
        throw std::range_error("polynomial cannot be identically zero");
}

double transmission_time_polynomial_rate::density(interval_t tau) const
{
    if (tau < 0)
        return 0.0;
    return hazardrate(tau) * std::exp(-totalhazard(tau));
}

double transmission_time_polynomial_rate::hazardrate(interval_t tau) const
{
    if (tau < 0)
        return 0.0;
    double r = 0.0;
    for (std::size_t i = 0; i < coeffs.size(); ++i)
        r += std::pow(tau, i) * coeffs[i];
    return r;
}

double transmission_time_polynomial_rate::totalhazard(interval_t tau) const
{
    if (tau < 0)
        return 0.0;
    double s = 0.0;
    for (std::size_t i = 0; i < coeffs.size(); ++i)
        s += std::pow(tau, i + 1) * coeffs[i] / (i + 1);
    return s;
}

double transmission_time_polynomial_rate::hazardbound(interval_t) const
{
    return INFINITY;
}

double transmission_time_polynomial_rate::survivalprobability(interval_t tau) const
{
    if (tau < 0)
        return 0.0;
    return std::exp(-totalhazard(tau));
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME: INFECTIOUSNESS -------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

transmission_time_infectiousness::transmission_time_infectiousness(const std::vector<double>& taus, const std::vector<double>& lambdas)
{
    if ((taus.size() != lambdas.size()) || (taus.size() == 0))
        throw std::runtime_error("size of tau and lambda vectors must agree and be non-empty");

    /* Fill ordered map representing lambda(tau) */
    lambda_max = 0.0;
    for(std::size_t i=0,n=taus.size(); i < n; ++i) {
		lambda.emplace(taus[i], std::make_pair(lambdas[i], NAN));
        lambda_max = std::max(lambda_max, lambdas[i]);
    }

    /* Fill ordered map representing Lambda(tau) = int 0 to tau lambda(tau') dtau' */
    double lambda_cum = 0;
    double tau_last = 0;
    double lambda_last = 0;
    for(auto& i: lambda) {
        const double tau_i = i.first;
        const double lambda_i = i.second.first;
		double& Lambda_i = i.second.second;

        const double dtau = (tau_i - tau_last);
        assert(dtau >= 0.0);
        if (dtau == 0.0)
            throw std::runtime_error("tau vector must not contain duplicate values");

		/* Update Lambda(tau) */
        lambda_cum += dtau * (lambda_i + lambda_last) / 2.0;
		Lambda_i = lambda_cum;
		
		/* Update Lambda^(-1)(tau).
		 * NOTE: It is essential that we overwrite existing entries here.
		 * If lambda(tau) === 0 on an interval, the map needs to contain
		 * the endpoint, not the starting point of that interval
		 */
		lambda_inverse[lambda_cum] = std::make_pair(tau_i, lambda_i);

        tau_last = tau_i;
        lambda_last = lambda_i;
    }
}

double transmission_time_infectiousness::density(interval_t tau) const
{
    if (tau < 0)
        return 0.0;
    return hazardrate(tau) * std::exp(-totalhazard(tau));
}

double transmission_time_infectiousness::hazardrate(interval_t tau) const
{
    if (tau <= 0)
        return 0.0;

	/* Find i,j with tau_i <= tau < tau_j.
	 * First, find first tau_j > tau. If no such element exists, extrapolate.
	 * Otherwise, set i = j-1 or tau_i = 0 if tau_j is the first element, and
	 * interpolate linearly.
	 */
	const auto j = lambda.upper_bound(tau);
	const bool i_valid = (j != lambda.begin());
	auto i = j; if (i_valid) --i;
	const double tau_i = i_valid ? i->first : 0;
	const double lambda_i = i_valid ? i->second.first : 0;
	const double dtau_point = tau - tau_i;

	if (j != lambda.end()) {
		/* interpolate linearly */
		const double tau_j = j->first;
		const double lambda_j = j->second.first;
		const double dtau = tau_j - tau_i;
		const double dlambda = lambda_j - lambda_i;
		return lambda_i + dtau_point * dlambda / dtau;
	} else {
		/* extrapolate constant */
		return lambda_i;
	}
}

double transmission_time_infectiousness::totalhazard(interval_t tau) const
{
    if (tau <= 0)
        return 0.0;

	/* Find i,j with tau_i <= tau < tau_j.
	 * First, find first tau_j > tau. If no such element exists, extrapolate.
	 * Otherwise, set i = j-1 or tau_i = 0 if tau_j is the first element, and
	 * interpolate quadratically.
	 */
	const auto j = lambda.upper_bound(tau);
	const bool i_valid = (j != lambda.begin());
	auto i = j; if (i_valid) --i;
	const double tau_i = i_valid ? i->first : 0;
	const double lambda_i = i_valid ? i->second.first : 0;
	const double Lambda_i = i_valid ? i->second.second : 0;
	const double dtau_point = tau - tau_i;

	if (j != lambda.end()) {
		/* interpolate quadratically */
		const double tau_j = j->first;
		const double lambda_j = j->second.first;
		const double dtau = tau_j - tau_i;
		const double dlambda = lambda_j - lambda_i;
		return Lambda_i + dtau_point * (lambda_i + 0.5 * dlambda * dtau_point / dtau);
	} else {
		/* extrapolate linearly */
		return Lambda_i + (tau - tau_i) * lambda_i;
	}
}

double transmission_time_infectiousness::totalhazard_inverse(interval_t Lambda) const
{
    if (Lambda < 0)
        return NAN;
	
	/* Find i,j with Lambda_i <= Lambda < Lambda_j.
	 * First, find first Lambda_j > Lambda. If no such element exists, extrapolate.
	 * Otherwise, set i = j-1 or Lambda_i = 0 if Lambda_j is the first element.
	 * Finally, compute tau by inverting totalhazard() locally.
	 */
	const auto j = lambda_inverse.upper_bound(Lambda);
	const bool i_valid = (j != lambda_inverse.begin());
	auto i = j; if (i_valid) --i;
	const double Lambda_i = i_valid ? i->first : 0;
	const double tau_i = i_valid ? i->second.first : 0;
	const double lambda_i = i_valid ? i->second.second : 0;
	const double dLambda_point = Lambda - Lambda_i;

	if (j != lambda_inverse.end()) {
		/* interpolate by solving a quadratic equation. We have in totalhazard() that
		 *   dLambda_point = dtau_point * (lambda_i + 0.5 * dlambda * dtau_point / dtau)
         * meaning
		 *   dtau_point^2 * a + dtau_point * b + c
		 * where a = 0.5 * dlambda / dtau, b = lambda_i, c = -dLambda_point. Therefore
		 * dtau_point = (sqrt(b^2 - 4*a*c) - b) / 2a
		 */
		const double tau_j = j->second.first;
		const double lambda_j = j->second.second;
		const double dtau = tau_j - tau_i;
		const double dlambda = lambda_j - lambda_i;
		const double a = 0.5 * dlambda / dtau;
		if (a != 0)
			return tau_i + (std::sqrt(lambda_i*lambda_i + 4*a*dLambda_point) - lambda_i) / (2*a);
		else
			return tau_i + dLambda_point / lambda_i;
	} else {
		/* extrapolate linearly */
		return tau_i + dLambda_point / lambda_i;
	}
}

double transmission_time_infectiousness::hazardbound(interval_t) const
{
    return lambda_max;
}

double transmission_time_infectiousness::survivalprobability(interval_t tau) const
{
    if (tau < 0)
        return 0.0;
    return std::exp(-totalhazard(tau));
}

double transmission_time_infectiousness::survivalquantile(double p) const
{
    if ((p < 0) || (p > 1.0))
        return NAN;
    if (p == 0.0)
        return INFINITY;
    return totalhazard_inverse(-std::log(p));
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------SAMPLE_WITHOUT_REPLACEMENT---------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

sample_without_replacement::sample_without_replacement(
    std::size_t A, std::size_t B, std::size_t k, rng_t &engine)
{
    if (B - A < k)
        throw std::range_error("interval [" + std::to_string(A) + ", " + std::to_string(B) + ") "
                                                                                             "contains too few elements for a sample of size " +
                               std::to_string(k));

    /* ReservoirSample(S[1..n], R[1..k])
     *   for i = 1 to k
     *       R[i] := S[i]
     *   end
     *   W := exp(log(random())/k)
     *   while i <= n
     *       i := i + floor(log(random())/log(1-W)) + 1
     *       if i <= n
     *           R[randomInteger(1,k)] := S[i]  // random index between 1 and k, inclusive
     *           W := W * exp(log(random())/k)
     *       end
     *   end
     * end
     */

    reservoir.reserve(k);
    for (std::size_t i = 0; i < k; ++i)
        reservoir.push_back(A + i);

    double W      = exp(log(std::uniform_real_distribution<double>(0, 1)(engine)));
    std::size_t i = k + 1;
    while (i < (B - A)) {
        const double u = std::uniform_real_distribution<double>(0, 1)(engine);
        i += floor(log(u) / log(1 - W)) + 1;
        if (i < (B - A)) {
            const std::size_t j = std::uniform_int_distribution<std::size_t>(0, k - 1)(engine);
            reservoir[j]        = i;
            const double u      = std::uniform_real_distribution<double>(0, 1)(engine);
            W *= exp(log(u) / k);
        }
    }
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------TRANSMISSION TIME:DETERMINISTIC----------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

double transmission_time_deterministic::sample(rng_t &, interval_t t, double m) const
{
    if (t > value)
        throw std::logic_error("present time cannot be larger than sampled time");
    return value - t;
}

double transmission_time_deterministic::density(interval_t tau) const
{
    return (tau == value) ? INFINITY : 0;
}

double transmission_time_deterministic::hazardrate(interval_t tau) const
{
    return (tau == value) ? INFINITY : 0;
}

double transmission_time_deterministic::hazardbound(interval_t) const
{
    return INFINITY;
}

double transmission_time_deterministic::survivalprobability(interval_t tau) const
{
    return (tau < value) ? 1 : 0;
}

double transmission_time_deterministic::survivalprobability(interval_t tau, interval_t t, double m) const
{
    if ((t < 0) || !std::isfinite(t))
        throw std::range_error("condition t must be non-negative and finite");
    if (m < 1)
        throw std::range_error("m must be positive");
    if (tau < 0)
        return 1.0;
    return (t + tau < value) ? 1 : 0;
}

interval_t transmission_time_deterministic::survivalquantile(double u) const
{
    throw std::runtime_error("not implemented");
}

interval_t transmission_time_deterministic::survivalquantile(double u, interval_t t, double m) const
{
    throw std::runtime_error("not implemented");
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*-----------SUB RNGS---------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

struct sub_seed_adapter
{
    typedef std::uint_least32_t result_type;

    sub_seed_adapter(rng_t &engine_)
        : engine(engine_)
    {
    }

    template <typename It>
    void generate(It begin, It end)
    {
        for (It i = begin; i != end; ++i)
            *i = engine();
    }

    rng_t &engine;
};

sub_rngs::sub_rngs(std::size_t n, rng_t &engine)
{
    // Generate RNGs
    rngs.reserve(n);
    sub_seed_adapter s(engine);
    for (std::size_t i = 0; i < n; ++i)
        rngs.emplace_back(s);
}

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------- VARIOUS HELPER FUNCTIONS ---------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

std::vector<double> rgamma(int n, double a, double b, rng_t &engine)
{
    std::gamma_distribution<double> gam(a, b);

    std::vector<double> vec(n, 0.0); // Initialise empty vector of size n.
    for (int i = 0; i < n; i++)
        vec[i] = gam(engine);

    return vec;
}

double rgamma(double a, double b, rng_t &engine)
{
    std::gamma_distribution<double> gam(a, b);
    return gam(engine);
}

std::vector<double> rand(int n, rng_t &engine)
{
    std::vector<double> vec(n, 0.0);
    std::uniform_real_distribution<> dis(0, 1);
    for (int i = 0; i < n; i++) {
        vec[i] = dis(engine);
    }
    return vec;
}

double rand(double a, double b, rng_t &engine)
{
    std::uniform_real_distribution<> dis(a, b);
    return dis(engine);
}

int poissrnd(double lambda, rng_t &engine)
{
    typedef std::poisson_distribution<int> pois_int_t;
    pois_int_t pois(lambda);
    return (pois(engine));
}

int zipf(int alpha, int n, rng_t &engine)
{
    static int first       = 1; // Static first time flag
    static double harmonic = 0;
    static std::vector<double> c(n + 1, 0); // cumulant function
    double z;                               // Uniform random number (0 < z < 1)
    int zipf_value;                         // Computed exponential value to be returned
    int low  = 1;
    int high = n;
    int mid; // Binary-search bounds

    // Compute normalization constant on first call only
    if (first == 1) {
        for (int i = 1; i <= n; i++) {
            harmonic += 1 / pow(i, alpha);
        }

        c[0] = 0;
        for (int i = 1; i <= n; i++)
            c[i] = c[i - 1] + 1.0 / (pow(i, alpha + 1) * harmonic);

        first = 0;
    }

    // Pull a uniform random number (0 < z < 1)
    do {
        z = rand(0, 1, engine);
    } while ((z == 0) || (z == 1));

    while (true) {
        mid = floor((low + high) / 2);
        if (c[mid - 1] < z && z <= c[mid]) {
            zipf_value = mid;
            break;
        } else if (c[mid] >= z) {
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }

    // Assert that zipf_value is between 1 and N
    assert((zipf_value >= 1) && (zipf_value <= n));

    return (zipf_value);
}

const double PI = 3.1415926535897932384626433832795028842;

double pdf_log_normal(double t, double mean, double variance)
{
    if (t == 0) {
        return 0;
    }
    double mu    = 2 * log(mean) - 0.5 * log(pow(mean, 2.0) + variance);
    double sigma = sqrt(log(1 + variance / pow(mean, 2.0)));

    return 1 / (t * sigma * sqrt(2 * PI)) * exp(-pow(log(t) - mu, 2) / (2 * sigma * sigma));
}
double cdf_log_normal(double t, double mean, double variance)
{
    if (t == 0) {
        return 0;
    }
    double mu    = 2 * log(mean) - 0.5 * log(pow(mean, 2.0) + variance);
    double sigma = sqrt(log(1 + variance / pow(mean, 2.0)));
    return 0.5 * (1 + erf((log(t) - mu) / (sqrt(2) * sigma)));
}
