#ifndef LAMBDA_IDENTITIES_ENTITY_UPDATER_H
#define LAMBDA_IDENTITIES_ENTITY_UPDATER_H

#include <thread>
#include <vector>
#include "entity.h"
#include <functional>
#include <io/parallel_exec.h>

namespace lambda {
	struct EntityUpdate {
		Entity* entity;
		std::function<void(Entity*)> update;
	};

	struct EntityUpdater {
		void add(Entity* entity, std::function<void(Entity*)> update) {
			_updates.push_back({ entity, update });
		}

		void update() {
			parallel_exec<EntityUpdater>(this, parallel_update, _updates.size());
		}

	private:
		std::vector<EntityUpdate> _updates;

		static void parallel_update(const ParallelExecState<EntityUpdater>* state) {
			for (int i = state->start; i < state->self->_updates.size(); i += state->stagger)
				state->self->_updates[i].update(state->self->_updates[i].entity);
		}
	};
}

#endif
