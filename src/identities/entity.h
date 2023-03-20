#ifndef LAMBDA_IDENTITIES_ENTITY_H
#define LAMBDA_IDENTITIES_ENTITY_H

#include <typeindex>
#include <unordered_map>
#include "entity_data.h"
#include "transformations.h"

namespace lambda {
	struct Entity {
		TransformClaim transform;

		Entity(Transformations& transformations) : transform(transformations) {}

		void add_data(std::unique_ptr<EntityData> data) {
			_database.emplace(typeid(data), data.release());
		}

		void remove_data(std::type_index id) {
			delete _database[id];
			_database.erase(id);
		}

	private:
		std::unordered_map<std::type_index, EntityData*> _database;
	};
}

#endif
