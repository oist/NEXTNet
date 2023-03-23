#include "stdafx.h"
#include "algorithm.h"
#include "random.h"
#include "types.h"

using namespace std::string_literals;

//------------------------------------------
//--epidemic_on_dynamic_network_simulation--
//------------------------------------------

absolutetime_t
epidemic_on_dynamic_network_simulation::next(rng_t& engine)
{
	if (std::isnan(network_next))
		network_next = network->next(engine);
	if (std::isnan(simulation_next))
		simulation_next = simulation->next(engine);
	
	return std::min(network_next, simulation_next);
}

std::optional<network_or_epidemic_event_t>
epidemic_on_dynamic_network_simulation::step(rng_t& engine, absolutetime_t nexttime)
{
	while (true) {
		if (std::isnan(nexttime))
			nexttime = next(engine);
		
		if (std::isinf(nexttime))
			return std::nullopt;
		
		if (nexttime == network_next) {
			/* Next event is a network event, and the event cannot be empty */
			const network_event_t ev = network->step(engine, nexttime).value();
			assert(ev.time == nexttime);
			nexttime = NAN;
			switch (ev.kind) {
				case network_event_kind::neighbour_added: {
					/* Check if the source node is infected */
					auto it = infected_neighbour_state.find(ev.source_node);
					if (it != infected_neighbour_state.end()) {
						/* Source node is infected. Check if the new neighbour was previously known */
						neighbour_state_t& neighbours = it->second;
						auto it2 = neighbours.find(ev.target_node);
						if (it2 == neighbours.end()) {
							/* Unregistered neighbour. Edge can't be active, just activate it */
							simulation->notify_infected_node_neighbour_added(ev, engine);
						}
						else {
							/* Neighbour previously existed. Edge was removed and is now re-added, unmask it (false) */
							if (it2->second == false)
								throw std::logic_error("duplicate neighbour_added event");
							it2->second = false ;
						}
					}
					break;
				}
				case network_event_kind::neighbour_removed: {
					/* Check if the source node is infected */
					auto it = infected_neighbour_state.find(ev.source_node);
					if (it != infected_neighbour_state.end()) {
						/* Souce node is infected. Mark edge towards removed neighbour as masked (true). */
						neighbour_state_t& neighbours = it->second;
						/* Update neighbour sate (will create an entry if none exists!) */
						neighbours[ev.target_node] = true;
					}
					break;
				}
				default:
					throw std::logic_error("invalid network event kind: "s + name(ev.kind));
			}
			return ev;
		} else if (nexttime == simulation_next) {
			/* Next event is a simulation event, event cannot be empty
			 * We specify an event filter that ensures that infections traversing
			 * masked edges are blocked. Note that since we specify a filter, there
			 * won't necessarily be an event to return before nexttime. In that case,
			 * we start from the top.
			 */
			std::function<bool(event_t)> evf = std::bind(&epidemic_on_dynamic_network_simulation::simulation_event_filter,
														 this, std::placeholders::_1);
			std::optional<event_t> maybe_ev = simulation->step(engine, nexttime, evf);
			nexttime = NAN;
			if (!maybe_ev)
				continue;
			/* Got a simulation event */
			const event_t ev = *maybe_ev;
			switch (ev.kind) {
				case event_kind::outside_infection:
				case event_kind::infection: {
					/* New infection. Create empty neighbour table which implicitly marks all neighours as unmasked */
					infected_neighbour_state.emplace(ev.source_node, neighbour_state_t());
				}
				case event_kind::reset: {
					/* Remove node's neighbour table */
					const std::size_t r = infected_neighbour_state.erase(ev.source_node);
					assert(r == 1);
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

bool epidemic_on_dynamic_network_simulation::simulation_event_filter(event_t ev)
{
	switch (ev.kind) {
		case event_kind::infection: {
			/* Infection event. Find neighbour table entry for source node of infection event */
			auto it = infected_neighbour_state.find(ev.source_node);
			if (it == infected_neighbour_state.end())
				throw std::logic_error("missing neighbour table for infected node");
			const neighbour_state_t& neighbours = it->second;
			/* Block event if there's a neighbour table entry marking the edge as masked (true) */
			const auto it2 = neighbours.find(ev.node);
			return (it2 == neighbours.end()) || !it2->second;
		}
		
		default:
			/* Don't block event */
			return false;
	}
}
