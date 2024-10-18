#ifndef STATISTICS_H
#define STATISTICS_H

/**
 * @brief Computes the mean over the specified range of values
 * @param begin first element to be averaged
 * @param end one past the last ement to be averaged
 * @return the average value
 */
template<typename Iterator>
double mean(Iterator begin, Iterator end) {
    std::size_t n = 0;
    double v = 0;
    for(Iterator i = begin; i != end; ++i) {
        n += 1;
        v += *i;
    }
    return v / n;
}

/**
 * @brief Computes mean and variance of the specified range of values
 * @param begin first element to be included
 * @param end one past the last ement to be included
 * @return a pair contaning the mean and average
 */
template<typename Iterator>
std::pair<double, double> mean_variance(Iterator begin, Iterator end) {
    const double m = mean(begin, end);
    std::size_t n = 0;
    double v = 0;
    for(Iterator i = begin; i != end; ++i) {
        n += 1;
        v += std::pow((*i - m), 2);
    }
    return std::make_pair(m, v / (n - 1));
}

/**
 * @brief Simple symmetric Z-test (similar to a t-Test but for known variance)
 *
 * @return The symmetric p-value
 */
inline double ztest(double mean_obs, double sd_true, double mean_true) {
	using namespace std;
	const double z = (mean_obs - mean_true) / sd_true;
	return 1 - std::erf(abs(z) / sqrt(2));
}

/**
 * @brief Simple Z-test for the mean of independent observations
 *
 * @param values vector of observed values
 * @param sd_true true standard deviation of observations
 * @param mean_true true mean of observations
 * @return two-sided p-value for the hypothesis that the observations indeed have the specified true mean
 */
inline double ztest_mean(const std::vector<double>& values, double sd_true, double mean_true) {
	using namespace std;
	const std::size_t N = values.size();
	const double m = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
	return ztest(m, sd_true / sqrt(N), mean_true);
}

/**
 * @brief Simple Z-test for the variance of independent observations
 *
 * Note that this test is strictly speaking only correct if the observations are normally distributed
 *
 * @param values vector of observed values
 * @param sd_true true standard deviation of observations
 * @param mean_true true mean of observations
 * @return two-sided p-value for the hypothesis that the observations indeed have the specified true standard deviation
 */
inline double ztest_var(const std::vector<double>& values, double sd_true, double mean_true) {
	using namespace std;
	const std::size_t N = values.size();
	const double var_true = pow(sd_true, 2.0);
	double s = 0.0;
	for(const auto& v: values)
		s += pow(v - mean_true, 2.0);
	s /= (N - 1);
	const double p = ztest(s, 2.0*var_true / sqrt(N - 1), var_true);
	return p;
}

/**
 * @brief One-sample KS test
 *
 * @param values Observed values
 * @param cdf The CDF as a function
 * @return The p-value
 */
inline double kstest(const std::vector<double>& values, std::function<double(double)> cdf) {
	const std::size_t N = values.size();
	double ks = 0.0;
	std::vector<double> vsorted(values.begin(), values.end());
	std::sort(vsorted.begin(), vsorted.end());
	for(std::size_t i=0; i < N; ++i) {
		const double x = vsorted[i];
		const double Ft = cdf(x);
		if ((Ft < 0.0) || (Ft > 1.0))
			throw std::range_error("CDF(" + std::to_string(x) + ") = " + std::to_string(Ft) + " with lies outside of [0,1]");
		const double Fe1 = (double)i / N;
		const double Fe2 = (double)(i+1) / N;
		const double dist = std::max(std::abs(Ft - Fe1), std::abs(Ft - Fe2));
		assert((dist >= 0.0) && (dist <= 1.0));
		ks = std::max(ks, dist);
	}
	boost::math::kolmogorov_smirnov_distribution<double> ks_dist(N);
	return quantile(complement(ks_dist, ks));
}

#endif // STATISTICS_H
