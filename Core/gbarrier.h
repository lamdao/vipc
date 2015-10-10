#ifndef __STDTHREAD_BARRIER_H
#define __STDTHREAD_BARRIER_H
//--------------------------------------------------------------------------
#include <mutex>
#include <atomic>
#include <chrono>
#include <condition_variable>
//--------------------------------------------------------------------------
// This is good, however it's hard to implement reseting (?!).
// The barrier.h has ugly implementation but it "works"
//--------------------------------------------------------------------------
class Barrier
{
private:
	std::mutex mtx;
	std::condition_variable cv;
	std::atomic<size_t> count, nthreads;
public:
	explicit Barrier(std::size_t nthr): count(0), nthreads(nthr) {
	}
	void wait() {
		std::unique_lock<std::mutex> lock(mtx);
		if (++count == nthreads) {
			cv.notify_all();
		} else {
			cv.wait(lock, [this] { return count == nthreads; });
		}
	}
};
//--------------------------------------------------------------------------
#endif
