#include "tests/stdafx.h"

#include "random.h"
#include "tests/statistics.h"

#define TEST_DISTRIBUTION_MEAN_VARIANCE(name, m, v, N) \
    TEST_CASE(STRINGIFY(name) " distribution (m=" STRINGIFY(m) ", v=" STRINGIFY(v) ", N=" STRINGIFY(N) ")", "[random]") { \
        std::mt19937 engine; \
        transmission_time_ ## name d(m, v); \
        REQUIRE(d.mean == m); \
        REQUIRE(d.variance == v); \
        std::vector<double> s; \
        s.reserve(N);\
        for(int i=0; i<N; ++i) s.push_back(d.sample(engine, 0, 1)); \
        const double p_mean = ztest_mean(s, sqrt(v), m); \
        REQUIRE(p_mean >= 0.01); \
        const double p_sd = ztest_var(s, sqrt(v), m); \
        REQUIRE(p_sd >= 0.01); \
        const double p_ks = kstest(s, [&d](double x) { return 1.0 - d.survivalprobability(x); }); \
        REQUIRE(p_ks >= 0.01); \
    }

#define TEST_DISTRIBUTION_2PARAM(name, a, b, N) \
    TEST_CASE(STRINGIFY(name) " distribution (a=" STRINGIFY(a) ", b=" STRINGIFY(b) ", N=" STRINGIFY(N) ")", "[random]") { \
        std::mt19937 engine; \
        transmission_time_ ## name d(a, b); \
        std::vector<double> s; \
        s.reserve(N);\
        for(int i=0; i<N; ++i) s.push_back(d.sample(engine, 0, 1)); \
        const double p_mean = ztest_mean(s, sqrt(d.variance), d.mean); \
        REQUIRE(p_mean >= 0.01); \
        const double p_sd = ztest_var(s, sqrt(d.variance), d.mean); \
        REQUIRE(p_sd >= 0.01); \
        const double p_ks = kstest(s, [&d](double x) { return 1.0 - d.survivalprobability(x); }); \
        REQUIRE(p_ks >= 0.01); \
    }

TEST_DISTRIBUTION_MEAN_VARIANCE(lognormal, 2.0, 0.1, 1000)
TEST_DISTRIBUTION_MEAN_VARIANCE(lognormal, 3.0, 1.0, 1000)
TEST_DISTRIBUTION_MEAN_VARIANCE(lognormal, 4.0, 2.0, 1000)
TEST_DISTRIBUTION_MEAN_VARIANCE(lognormal, 5.0, 4.0, 1000)

TEST_DISTRIBUTION_MEAN_VARIANCE(gamma, 2.0, 0.1, 1000)
TEST_DISTRIBUTION_MEAN_VARIANCE(gamma, 3.0, 1.0, 1000)
TEST_DISTRIBUTION_MEAN_VARIANCE(gamma, 4.0, 2.0, 1000)
TEST_DISTRIBUTION_MEAN_VARIANCE(gamma, 5.0, 4.0, 1000)

TEST_DISTRIBUTION_2PARAM(weibull, 1.5, 5.0, 1000)
TEST_DISTRIBUTION_2PARAM(weibull, 2.0, 10.0, 1000)
TEST_DISTRIBUTION_2PARAM(weibull, 4.0, 15.0, 1000)

TEST_CASE("Polynomial rate distribution", "[random]") {
    std::mt19937 engine;

//    transmission_time_polynomial_rate({1.0, 2.0, 3.0});
}

TEST_CASE("Lognormal distribution incremental sampling (pinf=0)", "[random]") {
    std::mt19937 engine;

    transmission_time_lognormal ln(5, 1);

    const int N = 10000;
    const int K = 3;

    /* Sample ordered pairs (s_1, ..., s_k) from the transmission time
     * distribution by taking two independent samples and ordering them.
     */
    std::vector<std::vector<double>> samples_iid(K, std::vector<double> {});
    for(int i=0; i < N; ++i) {
        std::vector<double> s;
        for(int j=0; j < K; ++j) {
            const double t = ln.sample(engine, 0, 1);
            REQUIRE(t);
            s.push_back(t);
        }
        REQUIRE(s.size() == K);
        std::sort(s.begin(), s.end());
        for(int j=0; j < K; ++j)
            samples_iid[j].push_back(s[j]);
    }

    /* Sample ordered pairs (s_1, ..., s_k) from the transmission time
     * distribution by sampling s_i from the conditional distribution
     * with t = s_(i-1) and m = k-i+1. In other words, sample first the
     * minimum, then the second-largest sample conditioned on the value
     * of the minimum and so on. Note that due to the way we defined
     * transmission_time(), the conditional version of the distribtion
     * returns the *difference* between the condition t and the samples
     * value.
     */
    std::vector<std::vector<double>> samples_seq(K, std::vector<double> {});
    for(int i=0; i < N; ++i) {        
        double t = 0;
        for(int j=0; j < K; ++j) {
            const double d = ln.sample(engine, t, K-j);
            REQUIRE(d >= 0);
            t = t + d;
            samples_seq[j].push_back(t);
        }
    }

    /* The two sampling methods should give the same results.
     * Compare mean and variance between the two sampling methods to verify
     */
    for(int j=0; j < K; ++j) {
        REQUIRE(samples_iid[j].size() == N);
        REQUIRE(samples_seq[j].size() == N);
        const auto iid_mv = mean_variance(samples_iid[j].begin(), samples_iid[j].end());
        const auto seq_mv = mean_variance(samples_seq[j].begin(), samples_seq[j].end());
        const double z = ((iid_mv.first - seq_mv.first) * sqrt(N) /
                          (iid_mv.second + seq_mv.second));
        REQUIRE(std::abs(z) <= 3);
    }
}

TEST_CASE("Lognormal distribution incremental sampling (pinf>0)", "[random]") {
    std::mt19937 engine;

    const double PINF=0.2;

    transmission_time_lognormal ln(5, 1, PINF);

    const int N = 10000;
    const int K = 3;

    /* Sample ordered pairs (s_1, ..., s_k) from the transmission time
     * distribution by taking two independent samples and ordering them.
     */
    std::vector<std::vector<double>> samples_iid(K, std::vector<double> {});
    for(int i=0; i < N; ++i) {
        std::vector<double> s;
        for(int j=0; j < K; ++j) {
            const double t = ln.sample(engine, 0, 1);
            REQUIRE(t);
            s.push_back(t);
        }
        REQUIRE(s.size() == K);
        std::sort(s.begin(), s.end());
        for(int j=0; j < K; ++j)
            samples_iid[j].push_back(s[j]);
    }

    /* Sample ordered pairs (s_1, ..., s_k) from the transmission time
     * distribution by sampling s_i from the conditional distribution
     * with t = s_(i-1) and m = k-i+1. In other words, sample first the
     * minimum, then the second-largest sample conditioned on the value
     * of the minimum and so on. Note that due to the way we defined
     * transmission_time(), the conditional version of the distribtion
     * returns the *difference* between the condition t and the samples
     * value.
     */
    std::vector<std::vector<double>> samples_seq(K, std::vector<double> {});
    for(int i=0; i < N; ++i) {        
        double t = 0;
        for(int j=0; j < K; ++j) {
            const double d = std::isfinite(t) ? ln.sample(engine, t, K-j) : INFINITY;
            REQUIRE(d >= 0);
            t = t + d;
            samples_seq[j].push_back(t);
        }
    }

    /* The two sampling methods should give the same results.
     * Compare mean and variance between the two sampling methods to verify
     */
    for(int j=0; j < K; ++j) {
        REQUIRE(samples_iid[j].size() == N);
        REQUIRE(samples_seq[j].size() == N);
        
        /* Check if the number of finite samples agrees with a z-test assuming binomial variance */
        const double iid_fin = std::count_if(samples_iid[j].begin(), samples_iid[j].end(),
                                             [] (auto &e) { return std::isfinite(e); });
        const double seq_fin = std::count_if(samples_seq[j].begin(), samples_seq[j].end(),
                                             [] (auto &e) { return std::isfinite(e); });
        const double z_fin = ((iid_fin - seq_fin) /
                              (sqrt(iid_fin*(N - iid_fin)/N) + sqrt(seq_fin*(N - seq_fin)/N)));
        REQUIRE(std::abs(z_fin) <= 3);
        
        /* Check if the means of the finite samples agree with a z-test */
        auto copy_iid = samples_iid[j];
        copy_iid.erase(std::remove_if(copy_iid.begin(), copy_iid.end(), [] (double t) { return std::isinf(t); }),
                       copy_iid.end());
        const auto iid_mv = mean_variance(copy_iid.begin(), copy_iid.end());
        auto copy_seq = samples_seq[j];
        copy_seq.erase(std::remove_if(copy_seq.begin(), copy_seq.end(), [] (double t) { return std::isinf(t); }),
                       copy_seq.end());
        const auto seq_mv = mean_variance(copy_seq.begin(), copy_seq.end());
        const double z = ((iid_mv.first - seq_mv.first) * sqrt(N) /
                          (iid_mv.second + seq_mv.second));
        REQUIRE(std::abs(z) <= 3);
    }
}