#include "nextnet/tests/stdafx.h"

#include "nextnet/utility.h"
#include "nextnet/tests/statistics.h"

TEST_CASE("integer_set", "[utility]")
{
    typedef integer_set<int> set_t;

    set_t s1;
    s1.insert(1);
    REQUIRE(s1.size() == 1);
    REQUIRE(*s1.begin() == 1);
    auto i = s1.begin();
    REQUIRE(++i == s1.end());
    std::vector<int> v1;
    std::copy(s1.begin(), s1.end(), std::back_inserter(v1));
    REQUIRE(v1 == std::vector{ 1 });
    REQUIRE(s1.erase(1));
    REQUIRE(s1.size() == 0);
    REQUIRE(s1.begin() == s1.end());

    set_t s2;
    s2.insert(1);
    std::vector<int> v2;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v2));
    REQUIRE(v2 == std::vector{ 1 });
    s2.insert(2);
    std::vector<int> v3;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v3));
    REQUIRE(v3 == std::vector{ 1, 2 });
    s2.insert(4);
    std::vector<int> v4;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v4));
    REQUIRE(v4 == std::vector{ 1, 2, 4 });
    s2.insert(3);
    std::vector<int> v5;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v5));
    REQUIRE(v5 == std::vector{ 1, 2, 3, 4 });
    s2.insert(-2);
    std::vector<int> v6;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v6));
    REQUIRE(v6 == std::vector{ -2, 1, 2, 3, 4 });
    s2.insert(6);
    std::vector<int> v7;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v7));
    REQUIRE(v7 == std::vector{ -2, 1, 2, 3, 4, 6 });
    s2.insert(-3);
    std::vector<int> v8;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v8));
    REQUIRE(v8 == std::vector{ -3, -2, 1, 2, 3, 4, 6 });
    s2.insert(0);
    std::vector<int> v9;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v9));
    REQUIRE(v9 == std::vector{ -3, -2, 0, 1, 2, 3, 4, 6 });

    REQUIRE(s2.find(-4) == s2.end());
    REQUIRE(*s2.find(-3) == -3);
    REQUIRE(*s2.find(-2) == -2);
    REQUIRE(s2.find(-1) == s2.end());
    REQUIRE(*s2.find(0) == 0);
    REQUIRE(*s2.find(1) == 1);
    REQUIRE(*s2.find(2) == 2);
    REQUIRE(*s2.find(3) == 3);
    REQUIRE(*s2.find(4) == 4);
    REQUIRE(s2.find(5) == s2.end());
    REQUIRE(*s2.find(6) == 6);
    REQUIRE(s2.find(7) == s2.end());

    s2.erase(2);
    std::vector<int> v10;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v10));
    REQUIRE(v10 == std::vector{ -3, -2, 0, 1, 3, 4, 6 });
    s2.erase(0);
    std::vector<int> v11;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v11));
    REQUIRE(v11 == std::vector{ -3, -2, 1, 3, 4, 6 });
    s2.erase(-3);
    std::vector<int> v12;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v12));
    REQUIRE(v12 == std::vector{ -2, 1, 3, 4, 6 });
    s2.erase(-2);
    std::vector<int> v13;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v13));
    REQUIRE(v13 == std::vector{ 1, 3, 4, 6 });
    s2.erase(4);
    std::vector<int> v14;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v14));
    REQUIRE(v14 == std::vector{ 1, 3, 6 });
    s2.erase(1);
    std::vector<int> v15;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v15));
    REQUIRE(v15 == std::vector{ 3, 6 });
    s2.erase(6);
    std::vector<int> v16;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v16));
    REQUIRE(v16 == std::vector{ 3 });
    s2.erase(3);
    std::vector<int> v17;
    std::copy(s2.begin(), s2.end(), std::back_inserter(v17));
    REQUIRE(v17.begin() == v17.end());

    set_t s3{ 1, 3, 4, 2, 6 };
    std::vector<int> v18;
    std::copy(s3.begin(), s3.end(), std::back_inserter(v18));
    REQUIRE(v18 == std::vector{ 1, 2, 3, 4, 6 });

    const auto i1 = s3.find(2);
    REQUIRE(*i1 == 2);
    const auto i2 = s3.find(4);
    REQUIRE(*i2 == 4);
    const auto i3 = s3.erase(i1, i2);
    REQUIRE(*i3 == 4);
    s3.insert(5);
    std::vector<int> v19;
    std::copy(s3.begin(), s3.end(), std::back_inserter(v19));
    REQUIRE(v19 == std::vector{ 1, 4, 5, 6 });
    const auto i4 = s3.erase(s3.find(1), s3.find(5));
    REQUIRE(*i4 == 5);
    std::vector<int> v20;
    std::copy(s3.begin(), s3.end(), std::back_inserter(v20));
    REQUIRE(v20 == std::vector{ 5, 6 });

    set_t s4(-1, 3);
    REQUIRE_NOTHROW(s4.insert(-1));
    REQUIRE_NOTHROW(s4.insert(3));
    REQUIRE_THROWS_AS(s4.insert(-2), std::range_error);
    REQUIRE_THROWS_AS(s4.insert(4), std::range_error);
}

TEST_CASE("integer_set draw_present", "[utility]")
{
    typedef integer_set<int> set_t;

    set_t s;
    s.insert(1);
    s.insert(2);
    s.insert(3);
    s.insert(5);
    s.insert(6);
    s.insert(8);

    rng_t rng;
    const int N = 10000;
    std::vector<int> counts;
    auto e = s.end();
    --e;
    counts.resize(*e + 1);
    for (int i = 0; i < N; ++i) {
        const int r = s.draw_element(rng);
        counts.at(r) += 1;
    }

    const double M  = N / s.size();
    const double SD = sqrt(M);
    for (std::size_t i = 0; i < counts.size(); ++i) {
        if (s.find((int)i) != s.end())
            REQUIRE(ztest(counts[i], SD, M) >= 0.01);
        else
            REQUIRE(counts[i] == 0);
    }
}

TEST_CASE("integer_set draw_absent", "[utility]")
{
    typedef integer_set<int> set_t;

    set_t s(0, 10);
    s.insert(1);
    s.insert(2);
    s.insert(3);
    s.insert(5);
    s.insert(6);
    s.insert(8);

    rng_t rng;
    const int K = 10;
    const int N = 10000;
    std::vector<int> counts;
    counts.resize(K + 1);
    for (int i = 0; i < N; ++i) {
        const int r = s.draw_complement(rng);
        counts.at(r) += 1;
    }

    const double M  = N / (K + 1 - s.size());
    const double SD = sqrt(M);
    for (std::size_t i = 0; i < counts.size(); ++i) {
        if (s.find((int)i) == s.end())
            REQUIRE(ztest(counts[i], SD, M) >= 0.01);
        else
            REQUIRE(counts[i] == 0);
    }
}
