#ifndef __STDTHREAD_BARRIER_H
#define __STDTHREAD_BARRIER_H
//--------------------------------------------------------------------------
#include <mutex>
#include <condition_variable>
//--------------------------------------------------------------------------
class Barrier
{
	enum State {Up, Down};
private:
#ifndef _OPENMP
	std::mutex mtx;
	std::condition_variable cv;
#endif
	size_t count, nthreads;
	State state;
public:
	explicit Barrier(std::size_t nthr):
				count(nthr), nthreads(nthr),
				state(State::Down) {}
	void wait() {
	#ifdef _OPENMP
		#pragma omp barrier
	#else
		std::unique_lock<std::mutex> lock(mtx);
		if (state == State::Up) {
			if (++count == nthreads) {
				state = State::Down;
				cv.notify_all();
			} else {
				cv.wait(lock, [this] { return state == State::Down; });
			}
		} else {
			if (--count == 0) {
				state = State::Up;
				cv.notify_all();
			} else {
				cv.wait(lock, [this] { return state == State::Up; });
			}
		}
	#endif
	}
};
//--------------------------------------------------------------------------
#endif
