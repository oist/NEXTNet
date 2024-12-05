#include "stdafx.h"
#include "algorithm.h"
#include "random.h"
#include "types.h"

using namespace std::string_literals;

//--------------------------------------
//----------SIMULATION_ALGORITHM--------
//--------------------------------------

void simulation_algorithm::notify_infected_contact(network_event_t event, rng_t& engine)
{
	throw std::logic_error("instantenous contacts are not implemented for this simulation algorithm");
}

//--------------------------------------
//----SIMULATE_ON_DYNAMIC_NETWORK-------
//--------------------------------------

simulate_on_temporal_network::simulate_on_temporal_network(simulation_algorithm& sim)
	:network(dynamic_cast<temporal_network*>(&sim.get_network())), simulation(sim)
{
	if (network == NULL)
		throw std::runtime_error("underlying simulation must use a dynamic_network");
};

absolutetime_t
simulate_on_temporal_network::next(rng_t& engine, absolutetime_t maxtime)
{
	const double sim_next = simulation.next(engine);
	const double nw_next = network->next(engine, std::min(sim_next, maxtime));
	return std::min(sim_next, nw_next);
}

std::optional<network_or_epidemic_event_t>
simulate_on_temporal_network::step(rng_t& engine, absolutetime_t maxtime)
{
	while (true) {
		const absolutetime_t nexttime = next(engine, maxtime);
		if (std::isinf(nexttime) || (nexttime > maxtime))
			return std::nullopt;

		if (nexttime == network->next(engine, maxtime)) {
			/* Next event is a network event. Since we don't specify an event filter,
			 * we should get an event for exactly nexttime.
			 */
			const network_event_t ev = network->step(engine, nexttime).value();
			assert(ev.time == nexttime);
			switch (ev.kind) {
				case network_event_kind::neighbour_added: {
					/* neighbour added, check if the source node is infected, otherwise nothing to do */
					auto it = infected_neighbour_state.find(ev.source_node);
					if (it != infected_neighbour_state.end()) {
						/* source infected, find the state of the neighbour */
						auto& neighbours = it->second;
						auto it2 = neighbours.find(ev.target_node);
						/* check state */
						if ((it2 == neighbours.end()) || (it2->second == neighbour_state_t::admissible)) {
							/* neighbour state is admissible (default). Activate edge */
							simulation.notify_infected_node_neighbour_added(ev, engine);
						}
						else if (it2->second == neighbour_state_t::masked) {
							/* neighbour is masked (edge was removed previously & hasn't fired). Unmask it. */
							it2->second = neighbour_state_t::admissible;
						}
						/* nothing to do if the neighbour state is transmitted */
					}
					break;
				}
				case network_event_kind::neighbour_removed: {
					/* neighbour removed, check if the source node is infected, otherwise nothing to do */
					auto it = infected_neighbour_state.find(ev.source_node);
					if (it != infected_neighbour_state.end()) {
						/* source infected, find the state of the neighbour */
						auto& neighbours = it->second;
						auto it2 = neighbours.find(ev.target_node);
						/* check state */
						if ((it2 == neighbours.end()) || (it2->second == neighbour_state_t::admissible)) {
							/* neighbour state is admissible (default), set state to masked */
							neighbours[ev.target_node] = neighbour_state_t::masked;
						}
						else if (it2->second == neighbour_state_t::masked) {
							/* neighbour state is already masked, this should not happen */
							throw std::logic_error("duplicate neighbour_removed event: t=" + std::to_string(ev.time)+ ", " + std::to_string(ev.source_node)+"->" + std::to_string(ev.target_node));
						}
						/* nothing to do if the neighbour state is transmitted */
					}
					break;
				}
				case network_event_kind::instantenous_contact: {
					/* instantenous contact, check if the source node is infected, otherwise nothing to do */
					auto it = infected_neighbour_state.find(ev.source_node);
					if (it != infected_neighbour_state.end()) {
						/* source infected, find the state of the neighbour */
						auto& neighbours = it->second;
						auto it2 = neighbours.find(ev.target_node);
						/* check state */
						if ((it2 == neighbours.end()) || (it2->second == neighbour_state_t::admissible)) {
							/* neighbour state is admissible (default), report contact */
							simulation.notify_infected_contact(ev, engine);
						}
						/* nothing to do if the neighbour state is masked or transmitted */
					}
					break;
				}
				default:
					throw std::logic_error("invalid network event kind: "s + name(ev.kind));
			}
			return ev;
		} else if (nexttime == simulation.next(engine)) {
			/* Next event is a simulation event, event cannot be empty
			 * We specify an event filter that ensures that infections traversing
			 * masked edges are blocked. Note that since we specify a filter, there
			 * won't necessarily be an event to return before nexttime. In that case,
			 * we start from the top.
			 */
			std::function<bool(epidemic_event_t)> evf = std::bind(&simulate_on_temporal_network::simulation_event_filter,
														 this, std::placeholders::_1);
			std::optional<epidemic_event_t> maybe_ev = simulation.step(engine, nexttime, evf);
			if (!maybe_ev)
				continue;
			/* Got a simulation event */
			const epidemic_event_t ev = *maybe_ev;
			switch (ev.kind) {
				case epidemic_event_kind::outside_infection:
				case epidemic_event_kind::infection: {
					/* infection event, initialize empty neighbour table, this makes edges admissible (by default) */
					infected_neighbour_state.emplace(ev.node, neighbours_states_t());
					break;
				}
				case epidemic_event_kind::reset: {
					/* recovery/reset event. remove nodes' neighbour table */
					const std::size_t r = infected_neighbour_state.erase(ev.node);
					assert(r == 1);
					_unused(r);
					break;
				}
				default:
					throw std::logic_error("invalid epidemic event kind: "s + name(ev.kind));
			}
			/* Inform the graph in case it wants to react */
			network->notify_epidemic_event(ev, engine);
			return ev;
		} else {
			throw std::logic_error("invalid next event time");
		}
	}
}

bool simulate_on_temporal_network::simulation_event_filter(epidemic_event_t ev)
{
	switch (ev.kind) {
		case epidemic_event_kind::infection: {
			/* Infection event, find the state of the neighbour in the source node */
			auto it = infected_neighbour_state.find(ev.source_node);
			if (it == infected_neighbour_state.end())
				throw std::logic_error("missing neighbour table for infected node");
			auto& neighbours = it->second;
			const auto it2 = neighbours.find(ev.node);
			/* check state */
			if ((it2 == neighbours.end()) || (it2->second == neighbour_state_t::admissible)) {
				/* state is admissible, don't block and set state to transmitted */
				neighbours[ev.node] = neighbour_state_t::transmitted;
				/* /!\ uncomment to allow the nodes to fire multiple times*/
				// neighbours[ev.node] = neighbour_state_t::admissible;
				return true;
			}
			else if (it2->second == neighbour_state_t::masked) {
				/* state is masked, block event and update to admissible */
				it2->second = neighbour_state_t::admissible;
				return false;
			}
			else {
				/* if the state is transmitted, no further infection events should occur */
				throw std::logic_error("infection along already transmitted edge: t=" + std::to_string(ev.time)+ ", " + std::to_string(ev.source_node)+"->" + std::to_string(ev.node));
			}
		}
		
		default:
			/* Don't block event */
			return true;
	}
}
