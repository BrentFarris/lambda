#ifndef LAMBDA_IDENTITIES_TRANSFORMATIONS_H
#define LAMBDA_IDENTITIES_TRANSFORMATIONS_H

#include <thread>
#include <vector>
#include <optional>
#include <math/vector3.h>
#include <math/matrix4x4.h>

enum class TransformRequestType {
	SetPosition,
	SetRotation,
	SetScale,
	ChangePosition,
	ChangeRotation,
	ChangeScale,
};

struct TransformRequest {
	TransformRequestType type;
	union {
		Vector3 translate;
		Vector3 rotate;
		Vector3 scale;
	};
};

struct Transform {
	Matrix4x4 matrix;
	Vector3 position;
	Vector3 rotation;
	Vector3 scale;
	int parent;

	Transform() : parent(-1) {}

	void reset() {
		matrix.reset();
		position = Vector3(0.0F);
		rotation = Vector3(0.0F);
		scale = Vector3(0.0F);
	}

	void transform() {
		for (const TransformRequest& request : _requests) {
			switch (request.type) {
			case TransformRequestType::SetPosition:
				position = request.translate;
				break;
			case TransformRequestType::SetRotation:
				rotation = request.rotate;
				break;
			case TransformRequestType::SetScale:
				scale = request.scale;
				break;
			case TransformRequestType::ChangePosition:
				position += request.translate;
				break;
			case TransformRequestType::ChangeRotation:
				rotation += request.rotate;
				break;
			case TransformRequestType::ChangeScale:
				scale += request.scale;
				break;
			}
		}
		_requests.clear();
		matrix.reset();
		matrix.scale(scale);
		matrix.rotate(rotation);
		matrix.translate(position);
	}

	void set_parent(int newParent) {
		parent = newParent;
	}

	void set_position(const Vector3& position) {
		_requests.push_back({ TransformRequestType::SetPosition, position });
	}

	void set_rotation(const Vector3& rotation) {
		_requests.push_back({ TransformRequestType::SetRotation, rotation });
	}

	void set_scale(const Vector3& scale) {
		_requests.push_back({ TransformRequestType::SetScale, scale });
	}

	void change_position(const Vector3& delta) {
		_requests.push_back({ TransformRequestType::ChangePosition, delta });
	}

	void change_rotation(const Vector3& delta) {
		_requests.push_back({ TransformRequestType::ChangeRotation, delta });
	}

	void change_scale(const Vector3& delta) {
		_requests.push_back({ TransformRequestType::ChangeScale, delta });
	}

private:
	std::vector<TransformRequest> _requests;
};

struct TransformClaim;

/**
 * Rather than transformation matrices being individually tracked by entities,
 * we will collect all transforms into this manager. This allows for queuing up
 * transformation changes and going through them all at once
 */
struct Transformations {
	friend TransformClaim;
private:
	struct ParallelTransformState {
		Transformations* self;
		std::optional<std::thread> thread;
		int start;
		int stagger;
	};

public:
	Transformations() = default;
	Transformations(int size) : _transformations(size) {}

	int claim() {
		if (_available.empty()) {
			_transformations.push_back({});
			return (int)_transformations.size() - 1;
		} else {
			int index = _available.back();
			_available.pop_back();
			_transformations[index].reset();
			return index;
		}
	}

	void release(int index) {
		_available.push_back(index);
	}

	void update() {
		parallel_exec(parallel_transform);
		parallel_exec(parallel_parents);
	}

private:
	std::vector<Transform> _transformations;
	std::vector<int> _available;
	std::vector<ParallelTransformState> _parallels;

	Transform& at(int index) { return _transformations[index]; }

	static void parallel_transform(const ParallelTransformState* state) {
		for (int i = state->start; i < state->self->_transformations.size(); i += state->stagger)
			state->self->_transformations[i].transform();
	}

	static void parallel_parents(const ParallelTransformState* state) {
		for (int i = state->start; i < state->self->_transformations.size(); i += state->stagger) {
			int parent = state->self->_transformations[i].parent;
			while (parent != -1) {
				state->self->_transformations[i].matrix *= state->self->_transformations[parent].matrix;
				parent = state->self->_transformations[parent].parent;
			}
		}
	}

	void parallel_exec(void(*call)(const ParallelTransformState*)) {
		int threadCount = (int)std::thread::hardware_concurrency();
		_parallels.reserve(threadCount);
		for (int i = 0; i < threadCount; ++i) {
			ParallelTransformState state;
			state.self = this;
			state.start = i;
			state.stagger = threadCount;
			_parallels.push_back(state);
			size_t len = _parallels.size();
			_parallels[len - 1].thread = std::thread(call, &_parallels[len - 1]);
		}
		for (auto& p : _parallels)
			if (p.thread->joinable())
				p.thread->join();
		_parallels.clear();
	}
};

/**
 * "Entities" should reference this structure instead of a raw index directly. This
 * has convenience of doing what needs to be done for read/write queuing. Use this
 * in place of where transform pointers would usually go in game engines
 */
struct TransformClaim {
	TransformClaim(Transformations& transformations) : _index(transformations.claim()),
		_transformations(transformations), _transform(transformations.at(_index)) {}
	~TransformClaim() { _transformations.release(_index); }
	const Matrix4x4& matrix() const { return _transform.matrix; }
	void set_position(const Vector3& position) { _transform.set_position(position); }
	void set_rotation(const Vector3& rotation) { _transform.set_position(rotation); }
	void set_scale(const Vector3& scale) { _transform.set_position(scale); }
	void translate(const Vector3& delta) { _transform.change_position(delta); }
	void rotate(const Vector3& delta) { _transform.change_rotation(delta); }
	void scale(const Vector3& delta) { _transform.change_scale(delta); }
	Vector3 position() const { return _transform.position; }
	Vector3 rotation() const { return _transform.rotation; }
	Vector3 scale() const { return _transform.scale; }
	void set_parent(const TransformClaim& parent) { set_parent_by_index(parent._index); }
	void set_parent_by_index(int parent) { _transform.set_parent(parent); }

private:
	int _index;
	// This is an aggregation since this pointer will never outlive Transformations
	Transform& _transform;
	Transformations& _transformations;
};

#endif
