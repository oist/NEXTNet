#include "popl.hpp"

#include "nextnet/stdafx.h"
#include "nextnet/random.h"
#include "nextnet/network_io.h"
#include "nextnet/factories/factories.h"

using namespace std;
using namespace literals;
using namespace popl;
using namespace factories;

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
    rng_t engine;
    factory<transmission_time>::object_holder_t psi;
	factory<transmission_time>::object_holder_t rho;
	factory<network>::object_holder_t nw;
    vector<algorithm::param_t> alg_params;
    unique_ptr<simulation_algorithm> alg;
    unique_ptr<simulate_on_temporal_network> sotn_alg;
    unique_ptr<ostream> nw_out_file;
    unique_ptr<ostream> out_file;
    ostream* nw_out = nullptr;
    ostream* out = nullptr;

	bool epidemic_events = false;
	bool network_events = false;

    OptionParser op("options");
    auto help_opt          = op.add<Switch>("h", "help", "produce help message");
    auto psi_opt           = op.add<Value<string>>("p", "transmission-time", "transmission time (psi)");
    auto rho_opt           = op.add<Value<string>>("r", "recovery-time", "recovery time (rho)");
    auto nw_opt            = op.add<Value<string>>("n", "network", "network to simulat on");
    auto alg_opt           = op.add<Value<string>>("a", "algorithm", "simulation algorithm to use", "next");
    auto param_opt         = op.add<Value<string>>("s", "parameter", "set simulation parameter");
    auto initial_opt       = op.add<Value<node_t>>("i", "initial-infection", "initial infected node");
    auto ev_opt            = op.add<Value<string>>("w", "report", "report events (e = epidemic, n = network)", "e");
    auto tmax_opt          = op.add<Value<double>>("t", "stopping-time", "stop simulation at this time");
    auto seed_opt          = op.add<Value<std::size_t>>("z", "seed", "random number generator seed");
    auto out_nw_opt        = op.add<Value<string>>("g", "output-network", "file to output network to");
    auto out_opt           = op.add<Value<string>>("o", "output", "output file");
    auto list_times_opt    = op.add<Switch>("", "list-times", "list distributions");
    auto list_networks_opt = op.add<Switch>("", "list-networks", "list network types");

    /* Parse arguments and setup things */

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

        /* Whether to output networks and/or epidemic events */

        epidemic_events = (ev_opt->value().find("e") != string::npos);
        network_events = (ev_opt->value().find("n") != string::npos);

        /* Parse algorithm parameters */

        for (size_t i = 0; i < param_opt->count(); ++i) {
            const string &p = param_opt->value(i);
            size_t j        = p.find('=');
            if (j == p.npos)
                throw program_argument_error("parameter", "invalid parameter setting '"s + p + "', does not contain '='");
            const string pname  = p.substr(0, j);
            const string pvalue = p.substr(j + 1, p.npos);
            alg_params.push_back({ pname, pvalue });
        }

        /* Setup rng */
        if (seed_opt->is_set()) {
            std::seed_seq seed { seed_opt->value() };
            engine.seed(seed);
        } else {
            std::random_device rd;
            std::seed_seq seed { rd() };
            engine.seed(seed);
        }
        random_engine = &engine;

        /* Parse times and networks */

        if (psi_opt->is_set())
            psi = time_factory.make(psi_opt->value());

        if (rho_opt->is_set())
            rho = time_factory.make(rho_opt->value());

        if (nw_opt->is_set())
            nw = network_factory.make(nw_opt->value());

        /* Create simulation */

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

        /* Add initial infections */

        for (size_t i = 0; i < initial_opt->count(); ++i) {
            const node_t node = initial_opt->value(i) - 1;
            if ((node < 0) || (node >= nw.first->nodes()))
                throw program_argument_error("initial-infection", "invalid initial node "s + boost::lexical_cast<string>(initial_opt->value(i)));
            alg->add_infections({ { node, 0.0 } });
        }
        if (initial_opt->count() == 0)
            cerr << "WARNING: No initially infected node specified with -initial-infection, no epidemic will commence" << endl;

        /* Create simulate_on_temporal_network algorithm if necessary */

        if (dynamic_cast<temporal_network *>(nw.first.get()) != nullptr)
			sotn_alg.reset(new simulate_on_temporal_network(*alg.get()));
    } catch (program_argument_error &e) {
        cerr << op << endl;
        cerr << "Error in argument " << e.invalid_argument << ": " << e.what() << endl;
        return 1;
    } catch (runtime_error &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    } catch (exception &e) {
        cerr << "Internal error: " << e.what() << endl;
        return 127;
    }

    /* Open outputs */

    try {
        if (out_opt->is_set() && (out_opt->value() != "-"s)) {
            out_file = make_unique<ofstream>(out_opt->value(), ios_base::out | ios_base::trunc);
            if (!*out_file)
                throw runtime_error("failed to create "s + out_opt->value());
            out = out_file.get();
        } else
            out = &cout;

        if (out_nw_opt->is_set() && (out_nw_opt->value() != "-"s)) {
            nw_out_file = make_unique<ofstream>(out_nw_opt->value(), ios_base::out | ios_base::trunc);
            if (!*nw_out_file)
                throw runtime_error("failed to create "s + out_nw_opt->value());
            nw_out = nw_out_file.get();
        }
        else if (out_nw_opt->is_set())
            nw_out = &cout;

    } catch (runtime_error &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    } catch (exception &e) {
        cerr << "Internal error: " << e.what() << endl;
        return 127;
    }

    try {
        *out << "#algorithm = " << alg_opt->value();
        if (sotn_alg)
            *out << " (temporal)";
        *out << endl;

        if (psi_opt->is_set())
            *out << "#psi = " << time_factory.parse(psi_opt->value()) << endl;

        if (rho_opt->is_set())
            *out << "#rho = " <<  time_factory.parse(rho_opt->value()) << endl;

        if (nw_opt->is_set()) {
            *out << "#nw = " <<  network_factory.parse(nw_opt->value()) << endl;
            if (nw_out != out) {
                *out << "#nw ";
                output_network_meta(*out, *nw.first);
            }
        }

        if (nw_out) {
            output_adjacencylist(*nw_out, *nw.first);
            nw_out->flush();
        }

        *out << "time" << '\t';
        *out << "epidemic_step" << '\t';
        *out << "network_step" << '\t';
        *out << "kind" << '\t';
        *out << "node" << '\t';
        *out << "neighbour" << '\t';
        *out << "total_infected" << '\t';
        *out << "total_reset" << '\t';
        *out << "infected" << '\n';

        const double tmax          = tmax_opt->is_set() ? tmax_opt->value() : INFINITY;
        size_t epidemic_step       = 0;
        size_t network_step        = 0;
        double infected            = 0;
        double total_infected      = 0;
        double total_reset         = 0;
        while (true) {
            // Execute next event
            const optional<network_or_epidemic_event_t> any_ev_opt =
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
            if (holds_alternative<epidemic_event_t>(any_ev)) {
                // Epidemic event
                const auto ev = get<epidemic_event_t>(any_ev);
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
            } else if (holds_alternative<network_event_t>(any_ev)) {
                // Network event
                const auto ev = get<network_event_t>(any_ev);
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
            *out << time << '\t';
            *out << epidemic_step << '\t';
            *out << network_step << '\t';
            *out << kind << '\t';
            *out << node << '\t';
            *out << neighbour << '\t';
            *out << total_infected << '\t';
            *out << total_reset << '\t';
            *out << infected << '\n';
        }

        return 0;
    } catch (exception &e) {
        cerr << "Internal error: " << e.what() << endl;
        return 127;
    }
}
