#include "nextnet/tests/stdafx.h"
#include "nextnet/tests/statistics.h"
#include "nextnet/temporal_network.h"

TEST_CASE("dynamic activity driven graph", "[temporal_network]")
{
    rng_t engine(0);

    const int N = 1000;
    const std::vector<double> activity_rates(N, 1.0);
    const double recovery_rate = 1.0;
    const double eta           = 1.0;
    const double m             = 1;
    activity_driven_network g(activity_rates, eta, m, recovery_rate, engine);

    double TMAX = 60;

    while (g.step(engine, TMAX)) {
        /* do nothing */
    }

    double av_k = 0;
    double k2   = 0;
    for (node_t node = 0; node < N; node++) {
        const double k = (double)g.outdegree(node);
        av_k += k / N;
        k2 += k * k / N;
    }
    double sigma    = sqrt(k2 - av_k * av_k);
    const double SE = sigma / sqrt((double)N);

    const double a = 1.0;
    const double b = recovery_rate;

    double expected_k = m * eta * a / (eta * a + b) * (1 + pow(b / (eta * a + b), 2));
    INFO(av_k);
    INFO(expected_k);

    REQUIRE(abs(av_k - expected_k) < 1.96 * SE);
}

TEST_CASE("dynamic empirical graph", "[temporal_network]")
{
    rng_t engine;

    double dt = 4;
    empirical_contact_network g(TEST_DATA_DIR "/test_empirical_network.txt", empirical_contact_network::finite_duration, dt);
    REQUIRE(g.nodes() == 14);

    /* evolve network up to time 0 */
    while (g.next(engine) <= 0.0)
        g.step(engine, 0.0);

    REQUIRE(g.neighbour(0, 0) == 13);
    REQUIRE(g.outdegree(0) == 2);
    REQUIRE(g.outdegree(10) == 0);
}

TEST_CASE("dynamic Erdös-Reyni", "[temporal_network]")
{
    const int N       = 20;
    const int E       = N * (N - 1) / 2;
    const double P    = 0.1;
    const double TAU  = 2;
    const double TMAX = 1000 * TAU;

    /* Keep a list of pairs (t, f) for every edge which stores the times
     * at which the edge appears (t, 1) or disappears (t, 0)
     */
    std::array<std::vector<std::pair<double, bool>>, E> edge_presence;

    std::mt19937 engine;
    temporal_erdos_reyni g(N, (N - 1) * P, TAU, engine);

    /* Initial edge_presence */
    for (node_t a = 0; a < N; ++a) {
        const int k = g.outdegree(a);
        for (int j = 0; j < k; ++j) {
            const node_t b = g.neighbour(a, j);
            const int e    = (int)edge_index_undirected(a, b);
            if (a > b)
                continue;
            REQUIRE(edge_presence[e].empty());
            edge_presence[e].emplace_back(0, true);
        }
    }
    for (int e = 0; e < E; ++e) {
        if (edge_presence[e].empty())
            edge_presence[e].emplace_back(0, false);
    }

    /* Evolve dynamic network */
    while (auto ev = g.step(engine, TMAX)) {
        const network_event_t event = *ev;
        const int e                 = (int)edge_index_undirected(event.source_node, event.target_node);

        switch (event.kind) {
            case network_event_kind::neighbour_added:
                if (!edge_presence[e].back().second)
                    /* First of fwd/rev pair for this event, add it */
                    edge_presence[e].emplace_back(event.time, true);
                else {
                    /* Second of fwd/rev pair for this event, should be there already */
                    REQUIRE(edge_presence[e].back().first == event.time);
                    REQUIRE(edge_presence[e].back().second == true);
                }
                break;

            case network_event_kind::neighbour_removed:
                if (edge_presence[e].back().second)
                    /* First of fwd/rev pair for this event, add it */
                    edge_presence[e].emplace_back(event.time, false);
                else {
                    /* Second of fwd/rev pair for this event, should be there already */
                    REQUIRE(edge_presence[e].back().first == event.time);
                    REQUIRE(edge_presence[e].back().second == false);
                }
                break;

            default:
                break;
        }
    }

    /* Compute average presence and absence times for every edge */
    for (int e = 0; e < E; ++e) {
        const auto &p = edge_presence[e];
        if (p.empty())
            continue;

        /* Compute present probability and times between events */
        std::vector<double> times_present;
        std::vector<double> times_absent;
        const auto end = p.end();
        auto curr      = p.begin();
        auto prev      = p.end();
        for (; curr != end; prev = curr++) {
            const double dt = curr->first - ((prev != p.end()) ? prev->first : 0.0);
            REQUIRE(((prev == p.end()) || (prev->second != curr->second)));
            if (curr->second) {
                /* appeared */
                times_absent.push_back(dt);
            } else {
                /* disappeared */
                times_present.push_back(dt);
            }
        }

        const double p_present        = std::accumulate(times_present.begin(), times_present.end(), 0.0) / TMAX;
        const double avg_time_present = std::accumulate(times_present.begin(), times_present.end(), 0.0) / times_present.size();
        const double avg_time_absent  = std::accumulate(times_absent.begin(), times_absent.end(), 0.0) / times_absent.size();

        /* Test that we see the exepcted number of events (where an edge first appearing and then disappearing is
         * a single event), mean present/absent times and presence probability. The statistics are the following:
         *   Edge process: Markov process with rate alpha = P/TAU of edge appearing, beta = (1-P)/TAU disappearing,
         *   No of events: Counting process with mean waiting time mu = (1/alpha + 1/beta) / 2,
         *                 variance of waiting time sigma^2 = 1/alpha^2 + 1/beta^2. The mean and
         *                 variance of the no. of events is then TMAX/mu and TMAX*sigma^2 / mu^3
         *                 as TMAX goes to infinity (per counting process theory)
         *   Avg(present): Average of exponentially distributed waiting times with mean TAU/(1-P)
         *   Avg(absent) : Average of exponentially distributed waiting times with mean TAU/P
         *   P(present)  : Scaled sum of TMAX/mu exponentially distributed waiting times with mean TAU/(1-P)
         * We use a simple z-test, i.e. assume normality and known variance of all tested quantities.
         * Since we test each edge separately, the p-value threshold of 0.01 is Bonferroni-corrected by
         * dividing by the number of edges.
         */

        /* Number of events */
        const double events = times_present.size();
        const double alpha = P / TAU, beta = 1 / TAU - P / TAU;
        const double mu          = 1 / alpha + 1 / beta;
        const double sigma2      = 1 / pow(alpha, 2) + 1 / pow(beta, 2);
        const double pval_events = ztest(events, sqrt(TMAX * sigma2 / pow(mu, 3)), TMAX / mu);
        REQUIRE(pval_events >= 0.01 / E);

        /* Present and absent times */
        const double pval_tpresent = ztest(avg_time_present, TAU / ((1 - P) * sqrt(times_absent.size())), TAU / (1 - P));
        REQUIRE(pval_tpresent >= 0.01 / E);
        const double pval_tabsent = ztest(avg_time_absent, TAU / (P * sqrt(times_absent.size())), TAU / P);
        REQUIRE(pval_tabsent >= 0.01 / E);

        /* Presence probability */
        const double pval_ppresent = ztest(p_present, std::sqrt((1 - P) / TMAX), P);
        REQUIRE(pval_ppresent > 0.01 / E);
    }

    std::cout << "---" << std::endl;
}
