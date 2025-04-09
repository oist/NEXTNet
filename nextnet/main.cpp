#include <filesystem>

#include "popl.hpp"

#include "nextnet/stdafx.h"
#include "nextnet/random.h"
#include "nextnet/factories/factories.h"

using namespace std;
using namespace std::literals;
using namespace popl;
using namespace factories;

rng_t engine;

/*
 * Errors
 */

struct program_argument_error : public runtime_error
{
    program_argument_error(const string &arg, const string &msg)
        : runtime_error(msg)
        , invalid_argument(arg)
    {
    }

    string invalid_argument;
};

/*
 * Main method
 */

int main(int argc, const char *argv[])
{
	factory<transmission_time>::object_holder_t psi;
	factory<transmission_time>::object_holder_t rho;
	factory<network>::object_holder_t nw;
    vector<algorithm::param_t> alg_params;
    unique_ptr<simulation_algorithm> alg;
    unique_ptr<simulate_on_temporal_network> sotn_alg;
	bool epidemic_events = false;
	bool network_events = false;

	random_engine = &engine;
	
    OptionParser op("options");
    auto help_opt          = op.add<Switch>("h", "help", "produce help message");
    auto psi_opt           = op.add<Value<std::string>>("p", "transmission-time", "transmission time (psi)");
    auto rho_opt           = op.add<Value<std::string>>("r", "recovery-time", "recovery time (rho)");
    auto nw_opt            = op.add<Value<std::string>>("n", "network", "network to simulat on");
    auto alg_opt           = op.add<Value<std::string>>("a", "algorithm", "simulation algorithm to use", "next");
    auto param_opt         = op.add<Value<std::string>>("s", "parameter", "set simulation parameter");
    auto initial_opt       = op.add<Value<node_t>>("i", "initial-infection", "initial infected node");
    auto ev_opt            = op.add<Value<std::string>>("w", "report", "report events (e = epidemic, n = network)", "e");
    auto tmax_opt          = op.add<Value<double>>("t", "stopping-time", "stop simulation at this time");
    auto list_times_opt    = op.add<Switch>("", "list-times", "list distributions");
    auto list_networks_opt = op.add<Switch>("", "list-networks", "list network types");

    try {
        /* Parse arguments */
        op.parse(argc, argv);

        /* Output help if requested, then stop */
        if (help_opt->is_set()) {
            cout << op << endl;
            return 0;
        }

        /* Output list of time distributions / networks if requested, then stop */

        if (list_times_opt->is_set()) {
            cout << "time distributions\n";
            cout << "------------------\n";
            for (const string &s : time_factory.descriptions)
                cout << s << "\n";
            return 0;
        }

        if (list_networks_opt->is_set()) {
            cout << "networks\n";
            cout << "--------\n";
            for (const string &s : network_factory.descriptions)
                cout << s << "\n";
            return 0;
        }

        /* Run simulation */

		if (psi_opt->is_set()) {
			psi = time_factory.make(psi_opt->value());
			cerr << "INFO: Infection time psi = " << time_factory.parse(psi_opt->value()) << std::endl;
		}

		if (rho_opt->is_set()) {
			rho = time_factory.make(rho_opt->value());
			cerr << "INFO: Infection time rho = " <<  time_factory.parse(rho_opt->value()) << std::endl;
		}

		if (nw_opt->is_set()) {
			nw = network_factory.make(nw_opt->value());
			cerr << "INFO: Network nw = " <<  network_factory.parse(nw_opt->value()) << std::endl;
		}
		
		epidemic_events = (ev_opt->value().find("e") != std::string::npos);
		network_events = (ev_opt->value().find("n") != std::string::npos);

        for (size_t i = 0; i < param_opt->count(); ++i) {
            const std::string &p = param_opt->value(i);
            std::size_t j        = p.find('=');
            if (j == p.npos)
                throw program_argument_error("parameter", "invalid parameter setting '"s + p + "', does not contain '='");
            const std::string pname  = p.substr(0, j);
            const std::string pvalue = p.substr(j + 1, p.npos);
            alg_params.push_back({ pname, pvalue });
        }

        {
            auto a_i = algorithms.find(alg_opt->value());
            if (a_i == algorithms.end())
                throw program_argument_error("algorithm", "unknown algorithm "s + alg_opt->value());
            algorithm &alg_factory = a_i->second;

            if (!nw.first)
                throw program_argument_error("network", "no network specified");

            if (!psi.first)
                throw program_argument_error("transmission-time", "no transmission time distribution specified");

            alg = alg_factory.create(*nw.first.get(), *psi.first.get(), rho.first.get(), alg_params);
        }

        for (size_t i = 0; i < initial_opt->count(); ++i) {
            const node_t node = initial_opt->value(i) - 1;
            if ((node < 0) || (node >= nw.first->nodes()))
                throw program_argument_error("initial-infection", "invalid initial node "s + boost::lexical_cast<string>(initial_opt->value(i)));
            alg->add_infections({ { node, 0.0 } });
        }
        if (initial_opt->count() == 0)
            std::cerr << "WARNING: No initially infected node specified with -initial-infection, no epidemic will commence" << std::endl;

		if (dynamic_cast<temporal_network *>(nw.first.get()) != nullptr) {
			std::cerr << "INFO: Simulating on a temporal network" << std::endl;
			sotn_alg.reset(new simulate_on_temporal_network(*alg.get()));
		}
    } catch (program_argument_error &e) {
        cerr << op << endl;
        cerr << "Error in argument " << e.invalid_argument << ": " << e.what() << endl;
        return 1;
    } catch (runtime_error &e) {
        cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (exception &e) {
        cerr << "Internal error: " << e.what() << std::endl;
        return 127;
    }

    try {
        cout << "time" << '\t';
        cout << "epidemic_step" << '\t';
        cout << "network_step" << '\t';
        cout << "kind" << '\t';
        cout << "node" << '\t';
        cout << "neighbour" << '\t';
        cout << "total_infected" << '\t';
        cout << "total_reset" << '\t';
        cout << "infected" << '\n';

        const double tmax          = tmax_opt->is_set() ? tmax_opt->value() : INFINITY;
        size_t epidemic_step       = 0;
        size_t network_step        = 0;
        double infected            = 0;
        double total_infected      = 0;
        double total_reset         = 0;
        while (true) {
            // Execute next event
            const std::optional<network_or_epidemic_event_t> any_ev_opt =
                sotn_alg ? sotn_alg->step(engine, tmax) : alg->step(engine, tmax);

            // Stop if there are no more events
            if (!any_ev_opt)
                break;
            const network_or_epidemic_event_t any_ev = *any_ev_opt;

            // Fill columns
            double time;
            const char *kind;
            node_t node;
            node_t neighbour = -1;
            if (std::holds_alternative<epidemic_event_t>(any_ev)) {
                // Epidemic event
                const auto ev = std::get<epidemic_event_t>(any_ev);
                // Update state
                time = ev.time;
                ++epidemic_step;
                switch (ev.kind) {
                    case epidemic_event_kind::outside_infection:
                    case epidemic_event_kind::infection:
                        ++total_infected;
                        ++infected;
                        break;
                    case epidemic_event_kind::reset:
                        ++total_reset;
                        --infected;
                        break;
                    default:
                        break;
                }
                // Report event?
                if (!epidemic_events)
                    continue;
                // Fill row
                kind = name(ev.kind);
                node = ev.node + 1;
                if (ev.kind == epidemic_event_kind::infection)
                    neighbour = ev.source_node + 1;
            } else if (std::holds_alternative<network_event_t>(any_ev)) {
                // Network event
                const auto ev = std::get<network_event_t>(any_ev);
                // Update state
                time = ev.time;
                ++network_step;
                // Report event?
                if (!network_events)
                    continue;
                // Fill row
                kind      = name(ev.kind);
                node      = ev.source_node + 1;
                neighbour = ev.target_node + 1;
            } else
                throw logic_error("unknown event type");

            // Output
            cout << time << '\t';
            cout << epidemic_step << '\t';
            cout << network_step << '\t';
            cout << kind << '\t';
            cout << node << '\t';
            cout << neighbour << '\t';
            cout << total_infected << '\t';
            cout << total_reset << '\t';
            cout << infected << '\n';
        }

        return 0;
    } catch (exception &e) {
        cerr << "Internal error: " << e.what() << std::endl;
        return 127;
    }
}
