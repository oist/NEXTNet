#include "tests/stdafx.h"

#include "random.h"

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

TEST_CASE("Weibull test", "[random]") {
    std::mt19937 engine;

    double shape = 5;
    double scale = 10;
    transmission_time_weibull wb(shape,scale);

    REQUIRE(abs(wb.mean - 9.1817) /(9.1817) < 0.01);
    REQUIRE(abs(wb.variance - 4.423) /(4.423) < 0.01);
    
    const int N = 100000;
    double m1 = 0;
    double m2c = 0;
    for(int i=0; i < N; ++i) {
        const double t = wb.sample(engine, 0, 1);
        m1 += t;
        m2c += (t - wb.mean) * (t - wb.mean);
    }
    REQUIRE(abs(wb.mean - m1/N) /(m1/N) < 0.01);
    REQUIRE(abs(wb.variance - m2c/N) /(m2c/N) < 0.01);

}

TEST_CASE("Lognormal distribution nonconditional", "[random]") {
    // TODO: Check that t=0, m=1 produces the correct mean and variance
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
