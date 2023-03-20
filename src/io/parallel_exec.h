#ifndef LAMBDA_IO_PARALLEL_EXEC_H
#define LAMBDA_IO_PARALLEL_EXEC_H

namespace lambda {
	template <class T>
	struct ParallelExecState {
		T* self;
		std::optional<std::thread> thread;
		int start;
		int stagger;
	};

	template <class T>
	static void parallel_exec(T* self, void(*call)(const ParallelExecState<T>*), size_t collectionSize) {
		size_t threadCount = std::thread::hardware_concurrency();
		int stagger = (int)std::min(threadCount, collectionSize);
		std::vector<ParallelExecState<T>> parallels{};
		parallels.reserve(stagger);
		for (int i = 0; i < stagger; ++i) {
			parallels.push_back({});
			ParallelExecState<T>& state = parallels.back();
			state.self = self;
			state.start = i;
			state.stagger = stagger;
			state.thread = std::thread(call, &parallels[parallels.size() - 1]);
		}
		for (auto& p : parallels)
			if (p.thread->joinable())
				p.thread->join();
		parallels.clear();
	}
}

#endif
