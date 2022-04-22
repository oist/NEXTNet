#include "tests/stdafx.h"

#include "random.h"

#if ENABLE_PLOTTING
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
#endif

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

TEST_CASE("Lognormal distribution nonconditional", "[random]") {
    // TODO: Check that t=0, m=1 produces the correct mean and variance
}

TEST_CASE("Lognormal distribution incremental sampling", "[random]") {
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
