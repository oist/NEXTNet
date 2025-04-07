#include "nextnet/tests/stdafx.h"
#include "nextnet/tests/simulate.h"
#include "nextnet/tests/plot.h"

#include "nextnet/brownian_proximity_network.h"
#include "nextnet/NextReaction.h"
#include "nextnet/algorithm.h"

TEST_CASE("Brownian proximity graph", "[brownian_proximity_graph]")
{
    rng_t engine;

    const std::size_t N = 100;
    const std::size_t K = 5;
    const double r      = 1.0;
    const double D0     = 1.0;
    const double D1     = 1.0;
    const double gamma  = 1.0;
    const double dt     = 0.1;
    brownian_proximity_network g(N, K, r, D0, D1, gamma, dt, engine);

    /* Slowly add infections */
    for (node_t ni = 0; ni <= (int)N; ++ni) {
        /* Check neighbour relationships */
        for (node_t n = 0; n < g.nodes(); ++n) {
            std::vector<double> cn;
            REQUIRE(g.coordinates(n, cn));

            std::unordered_set<node_t> nn1;
            for (node_t m = 0; m < g.nodes(); ++m) {
                std::vector<double> cm;
                REQUIRE(g.coordinates(m, cm));
                const double d = std::sqrt(std::pow(cn[0] - cm[0], 2) + std::pow(cn[1] - cm[1], 2));
                if (d <= 1.0)
                    nn1.insert(m);
            }

            std::unordered_set<node_t> nn2;
            for (index_t i = 0; i < g.outdegree(n); ++i)
                nn2.insert(g.neighbour(n, i));

            REQUIRE(nn1 == nn2);
        }
        if (ni >= (int)N)
            break;

        /* add infection */
        g.notify_epidemic_event(epidemic_event_t{
                                    .kind        = epidemic_event_kind::outside_infection,
                                    .source_node = -1,
                                    .node        = ni,
                                    .time        = (double)ni },
                                engine);

        /* evolve */
        while (g.step(engine, (double)(ni + 1)));
    }
}

TEST_CASE("Epidemic on Brownian proximity graph", "[brownian_proximity_graph]")
{
    rng_t engine;

    const std::size_t N = 100;
    const std::size_t K = 5;
    const double r      = 1.0;
    const double D0     = 0.1;
    const double D1     = 0.1;
    const double gamma  = 1.0;
    const double dt     = 0.1;
    const double Tmax   = 100;
    const double Tstep  = 1;
    brownian_proximity_network g(N, K, r, D0, D1, gamma, dt, engine);
    transmission_time_gamma psi(5, 3);
    // transmission_time_gamma rho(100, 10);
    simulate_next_reaction nr(g, psi, nullptr);
    simulate_on_temporal_network sim(nr);

    nr.add_infections({ std::make_pair(0, 0.0) });
    std::size_t ninfected = 0;
    for (double t = 0; t < Tmax; t += Tstep) {
        while (true) {
            const auto maybe_ev = sim.step(engine, t);
            if (!maybe_ev)
                break;

            const auto ev = *maybe_ev;
            if (std::holds_alternative<network_event_t>(ev)) {
                const auto nw_ev = std::get<network_event_t>(ev);
                _unused(nw_ev);
                /*
                std::cerr << "Network event " << name(nw_ev.kind) << ": "
                        << "time=" << nw_ev.time << ", "
                        << "src=" << nw_ev.source_node <<  ", "
                        << "tgt=" << nw_ev.target_node << std::endl;
                 */
            } else if (std::holds_alternative<epidemic_event_t>(ev)) {
                const auto ep_ev = std::get<epidemic_event_t>(ev);
                ninfected += delta_infected(ep_ev.kind);
                /*
                std::cerr << "Epidemic event " << name(ep_ev.kind) << ": "
                        << "time=" << ep_ev.time << ", "
                        << "src=" << ep_ev.source_node << ", "
                        << "tgt=" << ep_ev.node << std::endl;
                 */
            }
        }
    }

    REQUIRE(ninfected == N);
};

#if ENABLE_PLOTTING
TEST_CASE("Plot SIS average trajectory on Brownian proximity graph", "[brownian_proximity_graph]")
{
    using namespace std::string_literals;

    const std::size_t M       = 1;
    const std::size_t N       = 100;
    const double K            = 3.0;
    const double RADIUS       = 1.0;
    const double D0           = 0.1;
    const double D1           = 0.05;
    const double GAMMA        = 1.0;
    const double TAU          = 10;
    const double PSI_MEAN     = 3;
    const double PSI_VARIANCE = 1;
    const double RHO_MEAN     = 10;
    const double RHO_VARIANCE = 1;
    const double TMAX         = 500;

    rng_t engine;
    std::vector<double> t_sim, y_sim_new, y_sim_total;
    average_trajectories(engine, [&](rng_t &engine) {
		struct {
			std::unique_ptr<brownian_proximity_network> g;
			std::unique_ptr<transmission_time_gamma> psi;
			std::unique_ptr<transmission_time_gamma> rho;
			std::unique_ptr<simulate_next_reaction> nr;
			std::unique_ptr<simulate_on_temporal_network> simulator;
		} env;
		env.g = std::make_unique<brownian_proximity_network>(N, K, RADIUS, D0, D1, GAMMA, engine);
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

    std::vector<double> t_sim_erdos, y_sim_erdos_new, y_sim_erdos_total;
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
		} else throw std::logic_error("unknown event type"); }, t_sim_erdos, y_sim_erdos_total, y_sim_erdos_new, TMAX, M);

    plot("nextreaction_brownian.sis.mean.pdf", "SIS average trajectory on Brownian proximity graph", [&](auto &gp, auto &p) {
        p.add_plot1d(std::make_pair(t_sim, y_sim_new), "with lines title 'Brownian proximity'"s);
        p.add_plot1d(std::make_pair(t_sim_erdos, y_sim_erdos_new), "with lines title 'Dynamic Erdos'"s);
    });
}
#endif
