#ifndef __STDTHREAD_BARRIER_H
#define __STDTHREAD_BARRIER_H
//--------------------------------------------------------------------------
#include <omp.h>
#include <map>
#include <chrono>
#include <thread>
#include <algorithm>
//--------------------------------------------------------------------------
// This is ugly but kinda effective in reseting barrier (?!)
//--------------------------------------------------------------------------
class Barrier
{
private:
	int nthreads;
#ifndef _OPENMP
	std::map<std::thread::id,bool> ready;

	inline bool all_ready()
	{
		if (ready.size() < nthreads)
			return false;

		for (auto &r : ready) {
			if (!r.second) return false;
		}
		return true;
	}
#endif
public:
	explicit Barrier(int nthreads): nthreads(nthreads) {}

	inline void wait()
	{
	#ifdef _OPENMP
		#pragma omp barrier
	#else
		std::thread::id n = this_thread::get_id();
		ready[n] = true;	// not busy anymore
		while (true) {		// wait for all others to finish
			if (all_ready())
				break;
			this_thread::sleep_for(std::chrono::milliseconds(5));
		}
		this_thread::sleep_for(std::chrono::milliseconds(10*ready.size()));
		ready[n] = false;	// go back to busy state
	#endif
	}
};
//--------------------------------------------------------------------------
#endif
