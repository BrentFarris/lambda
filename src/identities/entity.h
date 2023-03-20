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

		template<typename T>
		const T* data() {
			return static_cast<T*>(_database[typeid(T)]);
		}

		template<typename T>
		T* mutable_data() {
			return static_cast<T*>(_database[typeid(T)]);
		}

		bool is_destroyed() const { return _destroyed; }
		void destroy() { _destroyed = true; }

	private:
		bool _destroyed;
		std::unordered_map<std::type_index, EntityData*> _database;
	};
}

#endif
