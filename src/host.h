#ifndef LAMBDA_HOST_H
#define LAMBDA_HOST_H

#include <vector>
#include <identities/entity.h>
#include <identities/entity_updater.h>
#include <identities/transformations.h>

namespace lambda {
	struct HostEntity {
		HostEntity(int entityIndex, std::vector<Entity>& collection)
			: entityIndex(entityIndex), collection(collection) {}
		Entity& entity() const { return collection[entityIndex]; }
	private:
		int entityIndex;
		std::vector<Entity>& collection;
	};

	struct Host {
		HostEntity new_entity() {
			_entities.push_back(Entity(transformations));
			return HostEntity((int)(_entities.size() - 1), _entities);
		}

		void frame_cleanup() {
			for (int i = _entities.size() - 1; i >= 0; i--) {
				if (_entities[i].is_destroyed())
					_entities.erase(_entities.begin() + i);
			}
		}

	private:
		EntityUpdater _updater;
		std::vector<Entity> _entities;
		Transformations transformations;
	};
}

#endif
