//
//  network_io.cpp
//  NEXTNet
//

#include "nextnet/stdafx.h"
#include "nextnet/types.h"
#include "nextnet/network_io.h"
#include "nextnet/network.h"
#include "nextnet/weighted_network.h"

//---------------------------------------------------
//-----Output networks-------------------------------
//---------------------------------------------------

void output_network_meta(std::ostream& dst, network& nw, char dsep)
{
    weighted_network* wnw = dynamic_cast<weighted_network*>(&nw);
    network_embedding* enw = dynamic_cast<network_embedding*>(&nw);
    std::size_t d = enw ? enw->dimensionality() : 0;

    const node_t N = nw.nodes();
    dst << "N=";
    ((N >= 0) ? dst << nw.nodes() : dst << "infinite" ) << ", is_undirected=" << nw.is_undirected();
    dst << ", is_simple=" << nw.is_simple() << ", is_weighted=" << (wnw != nullptr);
    dst << ", is_embedded=" << (enw != nullptr);
    if (enw != nullptr) {
        std::vector<double> x0(d, 0.0);
        std::vector<double> x1(d, 0.0);
        enw->bounds(x0, x1);
        dst << ", dims=" << enw->dimensionality();
        dst << ", lower=(";
        for(std::size_t i = 0; i < d; ++i)
            (i ? dst << dsep : dst) << x0[i];
        dst << "), upper=(";
        for(std::size_t i = 0; i < d; ++i)
            (i ? dst << dsep : dst) << x1[i];
        dst << ")";
    }
    dst << std::endl;
}

void output_adjacencylist(std::ostream& dst, network& nw, bool include_weights, bool include_coords,
                          bool include_meta, bool include_header, char csep, char dsep, char wsep)
{
    weighted_network* wnw = dynamic_cast<weighted_network*>(&nw);
    network_embedding* enw = dynamic_cast<network_embedding*>(&nw);

    const bool is_undirected = nw.is_undirected();
    include_weights = include_weights && (wnw != nullptr);
    include_coords = include_coords && (enw != nullptr);
    std::size_t d = enw ? enw->dimensionality() : 0;

    // Output meta information
    if (include_meta) {
        dst << '#';
        output_network_meta(dst, nw, dsep);
    }

    // Column headers
    if (include_header) {
        dst << "node";
        if (include_coords)
            for(std::size_t i = 0; i < d; ++i)
                dst << (i ? dsep : csep) << "x" << i;
        dst << csep;
        if (include_weights)
            dst << "neighbour:weight";
        else
            dst << "neighbour";
        dst << csep << "..." << std::endl;
    }

    // Output nodes
    std::vector<double> x(d, 0.0);
    for(node_t n=0, N=nw.nodes(); n < N; ++n) {
        dst << (n+1);

        // Output coordinates
        if (include_coords) {
            enw->coordinates(n, x);
            for(std::size_t i = 0; i < d; ++i)
                dst << (i ? dsep : csep) << x[i];
        }

        const index_t l = nw.outdegree(n);
        for(index_t i=0; i < l; ++i) {
            // Get neighbour (and edge weight if applicable)
            double w = NAN;
            const node_t nn = include_weights ? wnw->neighbour(n, i, &w) : nw.neighbour(n, i);

            // Only output edge if (src <= dist) for undirected networks
            if (is_undirected && (n > nn))
                continue;

            // Output neighbour
            dst << csep << (nn+1);

            // Output weight
            if (include_weights)
                dst << wsep << w;
        }

        dst << '\n';
    }
}
