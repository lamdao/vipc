#ifndef __STDTHREAD_BARRIER_H
#define __STDTHREAD_BARRIER_H
//--------------------------------------------------------------------------
#include <vector>
#include <chrono>
#include <thread>
#include <algorithm>
//--------------------------------------------------------------------------
// This is ugly but kinda effective in reseting barrier (?!)
//--------------------------------------------------------------------------
class Barrier
{
private:
	std::vector<bool> ready;

	inline bool all_ready()
	{
		for (auto r : ready) {
			if (!r) return false;
		}
		return true;
	}
public:
	explicit Barrier(int nthreads): ready(nthreads) {}

	inline void wait(const int id)
	{
		ready[id] = true;	// not busy anymore
		while (true) {		// wait for all others to finish
			if (all_ready())
				break;
			this_thread::sleep_for(std::chrono::milliseconds(5));
		}
		this_thread::sleep_for(std::chrono::milliseconds(10*ready.size()));
		ready[id] = false;	// go back to busy state
	}
};
//--------------------------------------------------------------------------
#endif
