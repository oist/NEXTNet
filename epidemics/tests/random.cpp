#include "tests/stdafx.h"

#include "random.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

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

TEST_CASE("Lognormal distribution sampling", "[random]") {
    std::mt19937 engine;

    transmission_time_lognormal ln(5, 1);
    REQUIRE(ln.survivalprobability(0) == 1.0);
    for(auto t: {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}) {
        for(auto m: {1}) { //, 2, 3}) {
            REQUIRE(ln.survivalprobability(0, t, m) == 1.0);
            const int n = 1000;
            double m1 = 0.0;
            double m2 = 0.0;
            for(int i=0; i < n; ++i) {
                const double v = ln.sample(engine, t, m);
                REQUIRE(v >= 0);
                m1 += v;
                m2 += v * v;
            }
            const double mean = m1 / n;
            const double variance = m2 / n - mean * mean;
            std::cout << "t: " << t << ", m: " << m << ", mean: " << mean << ", variance: " << variance << std::endl;
        }
    }
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
        for(int j=0; j < K; ++j)
            s.push_back(ln.sample(engine, 0, 1));
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
            t = t + ln.sample(engine, t, K-j);
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

#if 0
TEST_CASE("Lognormal distribution plots", "[random]") {
    std::mt19937 engine;

    const double t = 3;
    const int m = 1;
    transmission_time_lognormal ln(5, 1);
    const double n = 1000;
    std::vector<double> q(n, 0);
    std::vector<double> qp(n, 0);
    std::vector<double> u(n, 0);
    std::vector<double> uq(n, 0);
    for(int i=0; i < n; ++i) {
        q[i] = (double(i) / n) * 10.0;
        qp[i] = ln.survivalprobability(q[i], t, m);
        u[i] = (double(i+1) / (n+1)) * 1.0;
        uq[i] = ln.survivalquantile(u[i], t, 1);
    }

//    plt::plot(q, qp);
//    plt::show();

    plt::plot(u, uq);
    plt::show();
}
#endif

