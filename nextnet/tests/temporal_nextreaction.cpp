#include "tests/stdafx.h"
#include "tests/simulate.h"
#include "tests/analytical.h"
#include "tests/plot.h"

#include "random.h"
#include "NextReaction.h"
#include "temporal_network.h"
#include "algorithm.h"
#include "statistics.h"

namespace {
struct dynamic_single_edge : virtual network
    , virtual temporal_network
{
    bool edge_present;
    std::vector<absolutetime_t> times;

    template <typename... Args>
    dynamic_single_edge(bool _edge_present, Args &&...args)
        : edge_present(_edge_present)
        , times(std::forward<Args>(args)...)
    {
        std::reverse(times.begin(), times.end());
        if (!times.empty()) {
            double t = times[0];
            for (double tp : times) {
                if (tp > t)
                    throw std::range_error("edges flip times must be passed to dynamic_single_edge in ascending order");
                t = tp;
            }
        }
    }

    virtual node_t nodes() override { return 2; }

    virtual node_t neighbour(node_t node, int neighbour_index) override
    {
        if ((node == 0) && (neighbour_index == 0) && edge_present) return 1;
        return -1;
    }

    virtual int outdegree(node_t node) override
    {
        if (node == 0) return 1;
        return -1;
    }

    virtual absolutetime_t next(rng_t &engine, absolutetime_t) override
    {
        if (times.empty()) return INFINITY;
        return times.back();
    }

    virtual std::optional<network_event_t> step(rng_t &, absolutetime_t max_time = NAN) override
    {
        if (times.empty())
            return std::nullopt;
        const absolutetime_t t = times.back();
        times.pop_back();
        edge_present = !edge_present;
        return network_event_t{
            .kind        = edge_present ? network_event_kind::neighbour_added : network_event_kind::neighbour_removed,
            .source_node = 0,
            .target_node = 1,
            .weight      = 1.0,
            .time        = t
        };
    }
};
} // namespace

TEST_CASE("Effective transmission time distribution", "[temporal_nextreaction]")
{
    rng_t engine;

    using namespace std::string_literals;

    const std::size_t M = 1000;

    // Edge is absent [0.5, 1.0], [1.8, 2.0], [2.6, 3.2], [4.0, 5.0]
    const double INFECTION_TIME                       = 0.2;
    const bool EDGE_STATE_INITIAL                     = true;
    const std::vector<absolutetime_t> EDGE_FLIP_TIMES = { 0.5, 1.0, 1.8, 2.0, 2.6, 3.2, 4.0, 5.0 };

    // Transmission rate is 1.0*tau
    const std::vector<double> COEFFS = { 0.0, 1.0 };

    // Collect M transmission times
    std::vector<double> t_sim, y_sim_new, y_sim_total;
    average_trajectories(engine, [&](rng_t &engine) {
		struct {
			std::unique_ptr<dynamic_single_edge> g;
			std::unique_ptr<transmission_time_polynomial_rate> psi;
			std::unique_ptr<simulate_next_reaction> nr;
			std::unique_ptr<simulate_on_temporal_network> simulator;
		} env;
		env.g = std::make_unique<dynamic_single_edge>(EDGE_STATE_INITIAL, EDGE_FLIP_TIMES);
		env.psi = std::make_unique<transmission_time_polynomial_rate>(COEFFS);
		env.nr = std::make_unique<simulate_next_reaction>(*env.g.get(), *env.psi.get(), nullptr);
		env.nr->add_infections({ std::make_pair(0, INFECTION_TIME)});
		env.simulator = std::make_unique<simulate_on_temporal_network>(*env.nr.get());
		return env; }, [](network_or_epidemic_event_t any_ev) {
		/* Translate event into a pair (time, delta) */
		if (std::holds_alternative<epidemic_event_t>(any_ev)) {
			/* Epidemic event */
			const auto ev = std::get<epidemic_event_t>(any_ev);
			if (ev.kind == epidemic_event_kind::infection)
				return std::make_pair(ev.time, 1.0);
		}
		return std::make_pair((double)NAN, (double)NAN); }, t_sim, y_sim_total, y_sim_new, INFINITY, M, 1);

    // Check distribution
    const auto F = [&](double t) {
        // Compute the CDF of the effective transmission time, i.e. compute
        //   exp( - int_0^t lambda(t - t_infection) * epsilon(t') dt')
        // where t_infection is the infection time of the node,
        // lambda(tau) = sum_i c[i] * tau^i is the hazard rate as used
        // by transmission_time_polynomial_rate, and epsilon(tau) = 1
        // if the edge exists at time epsilon, i.e. if
        //   tau in [0, t_1] u [t_2, t_3] u ...
        bool edge_state = !EDGE_STATE_INITIAL;
        double t1       = 0.0;
        double t2       = 0.0;
        double s        = 0.0;
        for (std::size_t j = 0; j < EDGE_FLIP_TIMES.size(); j += 1, t1 = t2) {
            // update current interval, skip intervals were edge is absent
            t2         = EDGE_FLIP_TIMES[j];
            edge_state = !edge_state;
            if (!edge_state)
                continue;
            // integrate over [0, t] intersect [t_infection, infty] intersect [t1, t2]
            const double tstart = std::max(t1, INFECTION_TIME);
            const double tend   = std::min(t2, t);
            if (tstart >= tend)
                continue;
            for (std::size_t i = 0; i < COEFFS.size(); ++i) {
                s -= std::pow(tstart - INFECTION_TIME, i + 1) * COEFFS[i] / (i + 1);
                s += std::pow(tend - INFECTION_TIME, i + 1) * COEFFS[i] / (i + 1);
            }
        }
        return 1.0 - exp(-s);
    };

    std::vector<double> t_ecdf;
    for (std::size_t i = 0; i < t_sim.size(); ++i)
        t_ecdf.push_back((double)(i + 1) / t_sim.size());
    const double dt = t_sim.back() / 1000;
    std::vector<double> t_grid;
    std::vector<double> t_cdf;
    for (double t = 0; t <= t_sim.back(); t += dt) {
        t_grid.push_back(t);
        t_cdf.push_back(F(t));
    }

#if ENABLE_PLOTTING
    plot("dynamic.transmissiontime.pdf", "Transmission time on a fluctuating edge", [&](auto &gp, auto &p) {
        p.add_plot1d(std::make_pair(t_sim, t_ecdf), "with lines title 'empirical CDF'"s);
        p.add_plot1d(std::make_pair(t_grid, t_cdf), "with lines title 'theoretical CDF'"s);
    });
#endif

    const double pval = kstest(t_sim, F);
    CHECK(pval >= 0.01);
}

TEST_CASE("Epidemic on empirical network nb2 with finite edge durations", "[temporal_nextreaction]")
{
    rng_t engine(0);

    empirical_contact_network g(TEST_DATA_DIR "/college.tab", empirical_contact_network::finite_duration, 3);
    transmission_time_gamma psi(50, 3);
    transmission_time_gamma rho(100, 1);
    simulate_next_reaction nr(g, psi, &rho);
    nr.add_infections({ std::make_pair(0, 0.0) });
    simulate_on_temporal_network sim(nr);

    while (sim.step(engine));
}

TEST_CASE("Epidemic on empirical network nb2 with infitesimal edge durations", "[temporal_nextreaction]")
{
    rng_t engine(0);

    empirical_contact_network g(TEST_DATA_DIR "/college.tab", empirical_contact_network::infitesimal_duration, 3);
    transmission_time_gamma psi(50, 3);
    transmission_time_gamma rho(100, 1);
    simulate_next_reaction nr(g, psi, &rho);
    nr.add_infections({ std::make_pair(0, 0.0) });
    simulate_on_temporal_network sim(nr);

    while (sim.step(engine));
}

TEST_CASE("Epidemic on activity driven graph", "[temporal_nextreaction]")
{
    rng_t engine(0);

    const int N = 10;
    const std::vector<double> activity_rates(N, 0.1);
    const double recovery_rate = 1.0;
    const double eta           = 1.0;
    const double m             = 3;
    activity_driven_network graph(activity_rates, eta, m, recovery_rate, engine);

    graph.advance_time(engine, 30);

    // transmission_time_gamma psi(5,25);
    transmission_time_exponential psi(5);
    // transmission_time_gamma rho(100,1);
    simulate_next_reaction nr(graph, psi);
    nr.add_infections({ std::make_pair(0, 30) });
    simulate_on_temporal_network sim(nr);

    double tmax = 100;
    // while(auto ev = sim.step(engine,tmax)){

    while (true) {

        // std::optional<network_or_epidemic_event_t> any_ev = env.simulator -> step(engine,TMAX);
        std::optional<network_or_epidemic_event_t> any_ev = sim.step(engine, tmax);

        if (any_ev.has_value()) {
            if (std::holds_alternative<epidemic_event_t>(*any_ev)) {
                /* Epidemic event */
                // const auto& ev = std::get<event_t>(*any_ev);
                // std::cerr << "ep ";
                // std::cerr << ev.time << " | " << ev.source_node << "->" << ev.node << std::endl;

            } else if (std::holds_alternative<network_event_t>(*any_ev)) {
                /* Network event */
                // const auto& ev = std::get<network_event_t>(*any_ev);
                // int k = ev.kind == network_event_kind::neighbour_added ? 1 : 0;
                // if ( ev.source_node==8 && ev.target_node==7){
                // 	cnt++;
                // 	std::cerr << "counter :" << cnt << std::endl;
                // }
                // std::cerr << "network " << k << ": ";
                // std::cerr << ev.time << " | " << ev.source_node << "->" << ev.target_node << std::endl;
                // edges_array.push_back(number_of_edges);
            } else {
                throw std::logic_error("unknown event type");
            }

        } else { // no events
            break;
        }
    }
}

#if ENABLE_PLOTTING
TEST_CASE("Plot SIS average trajectories on activity driven network", "[temporal_nextreaction]")
{
    rng_t engine(2);

    using namespace std::string_literals;

    const std::size_t M = 1;

    const double PSI_MEAN     = 3;
    const double RHO_MEAN     = 10;
    const double RHO_VARIANCE = 1;
    const double TMAX         = 2000;

    const int N = 1000;
    std::vector<double> activity_rates(N, 1);
    const double recovery_rate = 1 / 20;
    const double eta           = 1;
    const int m                = 3;

    std::vector<double> t_sim, y_sim_new, y_sim_total;
    average_trajectories(engine, [&](rng_t &engine) {
		struct {
			std::unique_ptr<activity_driven_network> g;
			std::unique_ptr<transmission_time_exponential> psi;
			std::unique_ptr<transmission_time_gamma> rho;
			std::unique_ptr<simulate_next_reaction> nr;
			std::unique_ptr<simulate_on_temporal_network> simulator;
		} env;
		env.g = std::make_unique<activity_driven_network>(activity_rates,eta,m,recovery_rate,engine);
		env.psi = std::make_unique<transmission_time_exponential>(1/PSI_MEAN);
		env.rho = std::make_unique<transmission_time_gamma>(RHO_MEAN, RHO_VARIANCE);
		env.nr = std::make_unique<simulate_next_reaction>(*env.g.get(), *env.psi.get(), env.rho.get());
		env.nr->add_infections({ std::make_pair(0, 0.0)});
		env.simulator = std::make_unique<simulate_on_temporal_network>(*env.nr.get());
		return env; }, [](network_or_epidemic_event_t any_ev) {
		/* Translate event into a pair (time, delta) */
		if (std::holds_alternative<epidemic_event_t>(any_ev)) {
			/* Epidemic event */
			const auto ev = std::get<epidemic_event_t>(any_ev);
			return std::make_pair(ev.time, delta_infected(ev.kind));
		} else if (std::holds_alternative<network_event_t>(any_ev)) {
			/* Network event */
			const auto ev = std::get<network_event_t>(any_ev);
			return std::make_pair(ev.time, 0);
		} else throw std::logic_error("unknown event type"); }, t_sim, y_sim_total, y_sim_new, TMAX, M);

    plot("activity.driven.sis.mean.pdf", "SIS average trajectory on activity driven graph [NextReaction]", [&](auto &gp, auto &p) {
        p.add_plot1d(std::make_pair(t_sim, y_sim_new), "with lines title 'dynamic network'"s);
    });
}
#endif

#if ENABLE_PLOTTING
TEST_CASE("Plot SIS average trajectories on dynamic empirical network", "[temporal_nextreaction]")
{
    rng_t engine(2);

    using namespace std::string_literals;

    const std::size_t M = 1;

    const double PSI_MEAN     = 3;
    const double PSI_VARIANCE = 1;
    const double RHO_MEAN     = 10;
    const double RHO_VARIANCE = 1;
    const double TMAX         = 2000;
    const double DT           = 2;

    auto run = [&](empirical_contact_network::edge_duration_kind duration_kind) -> auto {
        std::vector<double> t, y_tot, y_new;
        average_trajectories(engine, [&](rng_t &engine) {
			struct {
				std::unique_ptr<empirical_contact_network> g;
				std::unique_ptr<transmission_time_gamma> psi;
				std::unique_ptr<transmission_time_gamma> rho;
				std::unique_ptr<simulate_next_reaction> nr;
				std::unique_ptr<simulate_on_temporal_network> simulator;
			} env;
			env.g = std::make_unique<empirical_contact_network>(TEST_DATA_DIR "/college.tab", duration_kind, DT);
			env.psi = std::make_unique<transmission_time_gamma>(PSI_MEAN, PSI_VARIANCE);
			env.rho = std::make_unique<transmission_time_gamma>(RHO_MEAN, RHO_VARIANCE);
			env.nr = std::make_unique<simulate_next_reaction>(*env.g.get(), *env.psi.get(), env.rho.get());
			env.nr->add_infections({ std::make_pair(0, 0.0)});
			env.simulator = std::make_unique<simulate_on_temporal_network>(*env.nr.get());
			return env; }, [](network_or_epidemic_event_t any_ev) {
			/* Translate event into a pair (time, delta) */
			if (std::holds_alternative<epidemic_event_t>(any_ev)) {
				/* Epidemic event */
				const auto ev = std::get<epidemic_event_t>(any_ev);
				return std::make_pair(ev.time, delta_infected(ev.kind));
			} else if (std::holds_alternative<network_event_t>(any_ev)) {
				/* Network event */
				const auto ev = std::get<network_event_t>(any_ev);
				return std::make_pair(ev.time, 0);
			} else throw std::logic_error("unknown event type"); }, t, y_tot, y_new, TMAX, M, 1);
        return std::make_tuple(t, y_tot, y_new);
    };

    auto edge_fin  = run(empirical_contact_network::finite_duration);
    auto edge_infi = run(empirical_contact_network::infitesimal_duration);

    plot("dynamic.empirical.sis.mean.pdf", "SIS average trajectory on dynamic empirical graph [NextReaction]", [&](auto &gp, auto &p) {
        p.add_plot1d(std::make_pair(std::get<0>(edge_fin), std::get<2>(edge_fin)), "with lines title 'finite edge duration'"s);
        p.add_plot1d(std::make_pair(std::get<0>(edge_infi), std::get<2>(edge_infi)), "with lines title 'infitesimal edge duration'"s);
    });
}
#endif

#if ENABLE_PLOTTING
TEST_CASE("Plot SIS average trajectory on dynamic Erdös-Reyni networks", "[temporal_nextreaction]")
{
    using namespace std::string_literals;

    const std::size_t M       = 100;
    const std::size_t N       = 100;
    const std::size_t K       = 3;
    const double TAU          = 10;
    const double PSI_MEAN     = 3;
    const double PSI_VARIANCE = 1;
    const double RHO_MEAN     = 10;
    const double RHO_VARIANCE = 1;
    const double TMAX         = 100;

    rng_t engine;
    std::vector<double> t_sim, y_sim_new, y_sim_total;
    average_trajectories(engine, [&](rng_t &engine) {
		struct {
			std::unique_ptr<temporal_erdos_reyni> g;
			std::unique_ptr<transmission_time_gamma> psi;
			std::unique_ptr<transmission_time_gamma> rho;
			std::unique_ptr<simulate_next_reaction> nr;
			std::unique_ptr<simulate_on_temporal_network> simulator;
		} env;
		env.g = std::make_unique<temporal_erdos_reyni>(N, K, TAU, engine);
		env.psi = std::make_unique<transmission_time_gamma>(PSI_MEAN, PSI_VARIANCE);
		env.rho = std::make_unique<transmission_time_gamma>(RHO_MEAN, RHO_VARIANCE);
		env.nr = std::make_unique<simulate_next_reaction>(*env.g.get(), *env.psi.get(), env.rho.get());
		env.nr->add_infections({ std::make_pair(0, 0.0)});
		env.simulator = std::make_unique<simulate_on_temporal_network>(*env.nr.get());
		return env; }, [](network_or_epidemic_event_t any_ev) {
		/* Translate event into a pair (time, delta) */
		if (std::holds_alternative<epidemic_event_t>(any_ev)) {
			/* Epidemic event */
			const auto ev = std::get<epidemic_event_t>(any_ev);
			return std::make_pair(ev.time, delta_infected(ev.kind));
		} else if (std::holds_alternative<network_event_t>(any_ev)) {
			/* Network event */
			const auto ev = std::get<network_event_t>(any_ev);
			return std::make_pair(ev.time, 0);
		} else throw std::logic_error("unknown event type"); }, t_sim, y_sim_total, y_sim_new, TMAX, M);

    std::vector<double> t_sim_static, y_sim_new_static, y_sim_total_static;
    average_trajectories(engine, [&](rng_t &engine) {
		struct {
			std::unique_ptr<network> g;
			std::unique_ptr<transmission_time_gamma> psi;
			std::unique_ptr<transmission_time_gamma> rho;
			std::unique_ptr<simulation_algorithm> simulator;
		} env;
		env.g = std::make_unique<erdos_reyni>(N, K, engine);
		env.psi = std::make_unique<transmission_time_gamma>(PSI_MEAN, PSI_VARIANCE);
		env.rho = std::make_unique<transmission_time_gamma>(RHO_MEAN, RHO_VARIANCE);
		env.simulator = std::make_unique<simulate_next_reaction>(*env.g.get(), *env.psi.get(), env.rho.get());
		return env; }, t_sim_static, y_sim_total_static, y_sim_new_static, TMAX, M);

    plot("dynamic.er.sis.mean.pdf", "SIS average trajectory on dynamic Erdös-Reyni [NextReaction]", [&](auto &gp, auto &p) {
        p.add_plot1d(std::make_pair(t_sim, y_sim_new), "with lines title 'dynamic network'"s);
        p.add_plot1d(std::make_pair(t_sim_static, y_sim_new_static), "with lines title 'static network'"s);
    });
}
#endif

#if ENABLE_PLOTTING
TEST_CASE("Plot SIR average trajectory on SIRX-Erdös-Reyni network", "[temporal_nextreaction]")
{
    using namespace std::string_literals;

    const std::size_t M       = 50;
    const std::size_t N       = 1000;
    const std::size_t K       = 4;
    const double PSI_MEAN     = 3;
    const double PSI_VARIANCE = 1;
    const double KAPPA0       = 0.01;
    const double KAPPA        = 0.1;
    const double TMAX         = 20;

    rng_t engine;
    std::vector<double> t_sim, y_sim_new, y_sim_total;
    average_trajectories(engine, [&](rng_t &engine) {
		struct {
			std::unique_ptr<erdos_reyni> g_static;
			std::unique_ptr<temporal_sirx_network> g;
			std::unique_ptr<transmission_time_gamma> psi;
			std::unique_ptr<simulate_next_reaction> nr;
			std::unique_ptr<simulate_on_temporal_network> simulator;
		} env;
		env.g_static = std::make_unique<erdos_reyni>(N, K, engine);
		env.g = std::make_unique<temporal_sirx_network>(*env.g_static.get(), KAPPA0, KAPPA);
		env.psi = std::make_unique<transmission_time_gamma>(PSI_MEAN, PSI_VARIANCE);
		env.nr = std::make_unique<simulate_next_reaction>(*env.g.get(), *env.psi.get(), nullptr);
		env.nr->add_infections({ std::make_pair(0, 0.0)});
		env.simulator = std::make_unique<simulate_on_temporal_network>(*env.nr.get());
		return env; }, [](network_or_epidemic_event_t any_ev) {
		/* Translate event into a pair (time, delta) */
		if (std::holds_alternative<epidemic_event_t>(any_ev)) {
			/* Epidemic event */
			const auto ev = std::get<epidemic_event_t>(any_ev);
			return std::make_pair(ev.time, delta_infected(ev.kind));
		} else if (std::holds_alternative<network_event_t>(any_ev)) {
			/* Network event */
			const auto ev = std::get<network_event_t>(any_ev);
			return std::make_pair(ev.time, 0);
		} else throw std::logic_error("unknown event type"); }, t_sim, y_sim_total, y_sim_new, TMAX, M);

    std::vector<double> t_sim_static, y_sim_new_static, y_sim_total_static;
    average_trajectories(engine, [&](rng_t &engine) {
		struct {
			std::unique_ptr<erdos_reyni> g;
			std::unique_ptr<transmission_time_gamma> psi;
			std::unique_ptr<simulation_algorithm> simulator;
		} env;
		env.g = std::make_unique<erdos_reyni>(N, K, engine);
		env.psi = std::make_unique<transmission_time_gamma>(PSI_MEAN, PSI_VARIANCE);
		env.simulator = std::make_unique<simulate_next_reaction>(*env.g.get(), *env.psi.get(), nullptr);
		return env; }, t_sim_static, y_sim_total_static, y_sim_new_static, TMAX, M);

    plot("dynamic.sirx_er.sir.mean.pdf", "Plot SIR average trajectory on SIRX-Erdös-Reyni network [next reaction]", [&](auto &gp, auto &p) {
        p.add_preamble("set logscale y");
        p.add_plot1d(std::make_pair(t_sim, y_sim_new), "with lines title 'dynamic SIRX-Erdös-Reyni'"s);
        p.add_plot1d(std::make_pair(t_sim_static, y_sim_new_static), "with lines title 'static Erdös-Reyni'"s);
    });
}
#endif
