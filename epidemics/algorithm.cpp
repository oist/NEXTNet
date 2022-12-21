#include "stdafx.h"
#include "algorithm.h"
#include "random.h"
#include "types.h"

using namespace std::string_literals;

//------------------------------------------
//--epidemic_on_dynamic_network_simulation--
//------------------------------------------

absolutetime_t
epidemic_on_dynamic_network_simulation::next()
{
	if (std::isnan(network_next))
		network_next = network->next();
	if (std::isnan(simulation_next))
		simulation_next = simulation->next();
	
	return std::min(network_next, simulation_next);
}

std::optional<network_or_epidemic_event_t>
epidemic_on_dynamic_network_simulation::step(rng_t& engine, absolutetime_t nexttime)
{
	if (std::isnan(nexttime))
		nexttime = next();
	
	if (std::isinf(nexttime))
		return std::nullopt;
	
	if (nexttime == network_next) {
		/* Next event is a network event, and the event cannot be empty */
		network_event_t ev = network->step(engine, nexttime).value();
		switch (ev.kind) {
			case network_event_kind::neighbour_added: {
				/* Check if the source node is infected */
				auto it = infected_neighbour_state.find(ev.node);
				if (it != infected_neighbour_state.end()) {
					/* Souce node is infected, insert neighbour table entry marking neighbour as unmasked */
					it->second.emplace(ev.neighbour, true);
					/* TODO: Mark edge as active */
				}
				break;
			}
			case network_event_kind::neighbour_removed: {
				/* Check if the source node is infected */
				auto it = infected_neighbour_state.find(ev.node);
				if (it != infected_neighbour_state.end()) {
					/* Souce node is infected, mark neighbour as masked in neighbour table */
					auto it2 = it->second.find(ev.neighbour);
					if (it2 == it->second.end())
						throw std::logic_error("missing neighbour table entry");
					/* Mark neighbour as currently masked */
					it2->second = 0;
				}
				break;
			}
			default:
				throw std::logic_error("invalid network event kind: "s + name(ev.kind));
		}
		return ev;		
	} else if (nexttime == simulation_next) {
		/* Next event is a simulation event, event cannot be empty */
		event_t ev = simulation->step(engine, nexttime).value();
		switch (ev.kind) {
			case event_kind::infection:
			case event_kind::outside_infection: {
				/* Create neighbour table, mark all currently existing neighbours as unmasked */
				auto r = infected_neighbour_state.emplace(ev.node, neighbour_state_t());
				assert(r.second);
				neighbour_state_t& neighbours = r.first->second;
				for(int i=0; i < network->outdegree(ev.node); ++i)
					neighbours.insert(std::make_pair(network->neighbour(ev.node, i), true));
				break;
			}
			case event_kind::reset: {
				/* Remove node's neighbour table */
				const std::size_t r = infected_neighbour_state.erase(ev.node);
				assert(r == 1);
				break;
			}
			default:
				throw std::logic_error("invalid epidemic event kind: "s + name(ev.kind));
		}
		/* Inform the graph in case it wants to react */
		network->notify_epidemic_event(ev);
		return ev;
	} else {
		throw std::logic_error("invalid next event time");
	}
}
