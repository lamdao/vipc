//--------------------------------------------------------------------------
// dthread.h - Multi-thread control, core of high-performance computing
//--------------------------------------------------------------------------
// Author: Lam H. Dao <daohailam(at)yahoo(dot)com>
//--------------------------------------------------------------------------
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//--------------------------------------------------------------------------
#ifndef __DTHREAD_H
#define __DTHREAD_H
//--------------------------------------------------------------------------
#include <mutex>
#include <condition_variable>
//--------------------------------------------------------------------------
using namespace std;
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
static const int MAX_NUM_THREADS = thread::hardware_concurrency();
//--------------------------------------------------------------------------
static inline void CalcWorkload(size_t total)
{
	load = (total / NUM_THREADS) + ((total % NUM_THREADS) > 0);
}
//--------------------------------------------------------------------------
class Permission
{
private:
	bool ok = false;
	condition_variable cv;
	mutex mtx;
public:
	void wait() {
		unique_lock<mutex> lock(mtx);
		cv.wait(lock, [this] { return ok; });
	}
	void allow() {
		ok = true;
		cv.notify_all();
	}
	void clear() {
		ok = false;
	}
};
//--------------------------------------------------------------------------
inline void SetAffinity(thread &thr, int n)
{
#ifdef USE_AFFINITY
#ifdef _WIN32
	SetThreadAffinityMask((HANDLE)thr.native_handle(), 1<<(n%MAX_NUM_THREADS));
#else
	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);
	CPU_SET(n % MAX_NUM_THREADS, &cpuset);
	pthread_setaffinity_np(thr.native_handle(), sizeof(cpu_set_t), &cpuset);
#endif
#endif
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
	static Permission permission;
	int __thread_id = 0;

	permission.clear();
	for (auto &w : workers) {
		w = thread([&](int __thread_id) {
				permission.wait();
				action(__thread_id);
			}, __thread_id);
		SetAffinity(w, __thread_id);
		__thread_id++;
	}

	permission.allow();
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
