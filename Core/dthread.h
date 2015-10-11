#ifndef __DTHREAD_H
#define __DTHREAD_H
//--------------------------------------------------------------------------
static size_t load;		// workload
//--------------------------------------------------------------------------
namespace DThread {
//--------------------------------------------------------------------------
#ifdef _OPENMP
static int NUM_THREADS = omp_get_max_threads();
#else
static int NUM_THREADS = thread::hardware_concurrency();
#endif
//--------------------------------------------------------------------------
static inline void CalcWorkload(size_t total)
{
	load = (total / NUM_THREADS) + ((total % NUM_THREADS) > 0);
}
//--------------------------------------------------------------------------
template<class Functor>
static inline void Start(Functor action)
{
#ifdef _OPENMP
	#pragma omp parallel
	{
		action(omp_get_thread_num());
	}
#else
	static vector<thread> workers(NUM_THREADS);
	int __thread_id = 0;
	for (auto &w : workers) {
		w = thread(action, __thread_id);
		__thread_id++;
	}
	for (auto &w : workers) {
		w.join();
	}
#endif
}
//--------------------------------------------------------------------------
} // namespace DThread
//--------------------------------------------------------------------------
#define NUM_THREADS	DThread::NUM_THREADS
//--------------------------------------------------------------------------
#define calc_working_range(start, stop)					\
	size_t start = __thread_id * load;					\
	size_t stop = start + load;							\
	if (stop > VS) stop = VS;
//--------------------------------------------------------------------------
#endif
