//
//  network_factory.cpp
//  NEXTNet
//
//  Created by Florian G. Pflug on 08.04.25.
//

#include "nextnet/stdafx.h"
#include "nextnet/network.h"
#include "nextnet/weighted_network.h"
#include "nextnet/temporal_network.h"
#include "nextnet/brownian_proximity_network.h"
#include "nextnet/factories/factories.h"

namespace factories {

/* Arguments of networks */

empirical_contact_network::edge_duration_kind parse_contact_kind(const std::string &s)
{
    if (s == "instantaneous")
        return empirical_contact_network::infitesimal_duration;
    else if (s == "finite")
        return empirical_contact_network::finite_duration;
    else
        throw std::range_error("contact_kind must be 'instantaneous' or 'finite'");
}

std::string render_contact_kind(const empirical_contact_network::edge_duration_kind &s)
{
    switch (s) {
        case empirical_contact_network::infitesimal_duration:
            return "instantaneous";
        case empirical_contact_network::finite_duration:
            return "finite";
        default:
            return "?";
    }
}

struct network_ref
{
    static network_ref parse(const std::string &s)
    {
        return { network_factory.parse(s), network_factory.make(s) };
    }

    static std::string render(const network_ref &s)
    {
        std::stringstream buf;
        buf << s.expr;
        return buf.str();
    }

    network_ref(parsed_expression_t expr_, factory<network>::object_holder_t nw_)
        : expr(expr_)
        , nw(std::move(nw_.first))
        , holder(std::move(nw_.second))
    {
    }

    operator network &() const { return *nw.get(); };

    parsed_expression_t expr;
    std::shared_ptr<network> nw;
    std::shared_ptr<std::any> holder;
};

/* Argument declarations */

DECLARE_ARGUMENT_5(file, istream_ref, std::nullopt, istream_ref::parse, istream_ref::render)
DECLARE_ARGUMENT_3(undirected, bool, true)
DECLARE_ARGUMENT_3(simplify, bool, false)
DECLARE_ARGUMENT_3(node_index_base, node_t, 1)
DECLARE_ARGUMENT_3(column_separator, char, ' ')
DECLARE_ARGUMENT_3(weight_separator, char, ':')
DECLARE_ARGUMENT_5(contact_kind, empirical_contact_network::edge_duration_kind,
                   empirical_contact_network::infitesimal_duration, parse_contact_kind, render_contact_kind);
DECLARE_ARGUMENT_3(dt, double, 1.0);
DECLARE_ARGUMENT_3(weight, double, 1.0);
DECLARE_ARGUMENT_3(size, node_t, std::nullopt);
DECLARE_ARGUMENT_3(avg_degree, double, std::nullopt);
DECLARE_ARGUMENT_5(weights, std::vector<double>, std::vector<double>{ 1.0 },
                   parse_vector<double>, render_vector<double>);
DECLARE_ARGUMENT_5(probabilities, std::vector<double>, std::vector<double>{ 1.0 },
                   parse_vector<double>, render_vector<double>);
DECLARE_ARGUMENT_3(timescale, double, 1.0);
DECLARE_ARGUMENT_5(activities, std::vector<double>, std::vector<double>{ 1.0 },
                   parse_vector<double>, render_vector<double>);
DECLARE_ARGUMENT_3(eta_sus, double, std::nullopt);
DECLARE_ARGUMENT_3(eta_inf, double, std::nullopt);
DECLARE_ARGUMENT_3(b_sus, double, std::nullopt);
DECLARE_ARGUMENT_3(b_inf, double, std::nullopt);
DECLARE_ARGUMENT_3(k, int, std::nullopt);
DECLARE_ARGUMENT_3(p, double, std::nullopt);
DECLARE_ARGUMENT_3(m, int, std::nullopt);
DECLARE_ARGUMENT_5(degrees, std::vector<int>, std::nullopt, parse_vector<int>, render_vector<int>);
DECLARE_ARGUMENT_5(triangles, std::vector<int>, std::nullopt, parse_vector<int>, render_vector<int>);
DECLARE_ARGUMENT_3(beta, double, std::nullopt);
DECLARE_ARGUMENT_3(edgelength, int, std::nullopt);
DECLARE_ARGUMENT_5(nw, network_ref, std::nullopt,
                   network_ref::parse, network_ref::render);
DECLARE_ARGUMENT_3(kappa, double, std::nullopt);
DECLARE_ARGUMENT_3(kappa0, double, std::nullopt);
DECLARE_ARGUMENT_3(radius, double, 1.0);
DECLARE_ARGUMENT_3(D0, double, 1.0);
DECLARE_ARGUMENT_3(D1, double, 1.0);
DECLARE_ARGUMENT_3(gamma, double, 1.0);
DECLARE_ARGUMENT_3(timestep, double, 0.1);
DECLARE_ARGUMENT_3(reduced_root_degree, bool, true);

/* Network factory */

factory<network> network_factory = factory<network>("network")
                                       .add<empirical_network, file, undirected, simplify, node_index_base,
                                            column_separator>("empirical")
                                       .add<weighted_empirical_network, file, undirected, simplify, node_index_base,
                                            column_separator, weight_separator>("weighted_empirical")
                                       .add<empirical_contact_network, file, contact_kind, dt, weight>("empirical_contact")
                                       .add<erdos_renyi, size, avg_degree, rng>("erdos_renyi")
                                       .add<weighted_erdos_renyi, size, avg_degree, weights, probabilities, rng>("weighted_erdos_renyi")
                                       .add<temporal_erdos_renyi, size, avg_degree, timescale, rng>("temporal_erdos_renyi")
                                       .add<activity_driven_network, activities, m, eta_sus, eta_inf, b_sus, b_inf, rng>("activity_driven")
                                       .add<watts_strogatz, size, k, p, rng>("watts_strogatz")
                                       .add<barabasi_albert, size, rng, m>("barabasi_albert")
                                       .add<config_model, degrees, rng>("config_model")
                                       .add<config_model_clustered_serrano, degrees, triangles, beta, rng>("config_model_clustered")
                                       .add<cubic_lattice_2d, edgelength>("lattice_2d")
                                       .add<cubic_lattice_3d, edgelength>("lattice_3d")
                                       .add<temporal_sirx_network, nw, kappa0, kappa, rng>("temporal_sirx")
                                       .add<brownian_proximity_network, size, avg_degree, radius, D0, D1, gamma, timestep, rng>("brownian_proximity")
                                       .add<acyclic, avg_degree, reduced_root_degree, rng>("acyclic")
                                       .add<fully_connected, size, rng>("fully_connected");

} /* namespace factories */
