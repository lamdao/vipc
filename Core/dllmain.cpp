//------------------------------------------------------------------------
#ifdef _WIN32
#define STRICT
#define VC_EXTRALEAN
#include <windows.h>
#include <winuser.h>
#endif
//--------------------------------------------------------------------------
#include <memory.h>
#include <math.h>
#include <omp.h>
//--------------------------------------------------------------------------
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>
//--------------------------------------------------------------------------
#ifdef _WIN32
#define DLAPI	extern "C" __declspec(dllexport)
#else
#define DLAPI	extern "C"
#endif
//--------------------------------------------------------------------------
#define uchar	unsigned char
#define ushort	unsigned short
#define uint	unsigned int
//------------------------------------------------------------------------
#define	UC		uchar
#define	US		ushort
#define	UI		uint
#define	FL		float
#define	DB		double
//--------------------------------------------------------------------------
#ifndef TRUE
#define TRUE	1
#define FALSE	0
#endif
//--------------------------------------------------------------------------
enum DataTypes {
	DT_NONE = 0,
	DT_BYTE = 1,
	DT_SHORT = 2,
	DT_INT = 3,
	DT_FLOAT = 4,
	DT_DOUBLE = 5,
	DT_USHORT = 12,
	DT_UINT = 13
};
//--------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------
// Library entry
//------------------------------------------------------------------------
#ifdef _WIN32
BOOL APIENTRY DllMain(HINSTANCE hModule, DWORD dwReason, LPVOID lpReserved)
{
	switch (dwReason){
		case DLL_PROCESS_ATTACH:
			DisableThreadLibraryCalls(hModule);
			break;
		case DLL_PROCESS_DETACH:
			break;
	}

	return TRUE;
}
#endif
//--------------------------------------------------------------------------
// Volume info/utils
//--------------------------------------------------------------------------
static int VX, VY, VZ, VT;	// volume dimensions (X,Y,Z) and type
//--------------------------------------------------------------------------
static size_t VP, VS;		// volume plane size (VP), volume size (VS)
//------------------------------------------------------------------------
static size_t load;		// workload
//--------------------------------------------------------------------------
static int NUM_THREADS = thread::hardware_concurrency();
//--------------------------------------------------------------------------
static vector<thread> workers(NUM_THREADS);
//--------------------------------------------------------------------------
static inline size_t pos2idx(int x, int y, int z)
{
	return (size_t)z * VP + (size_t)y * VX + (size_t)x;
}
//--------------------------------------------------------------------------
static inline void idx2pos(size_t n, int &x, int &y, int &z)
{
	z = (int)(n / VP); n = n % VP;
	y = (int)(n / VX);
	x = (int)(n % VX);
}
//------------------------------------------------------------------------
static inline void calc_workload(size_t total)
{
	load = (total / NUM_THREADS) + ((total % NUM_THREADS) > 0);
}
//--------------------------------------------------------------------------
template<class Functor>
static inline void thread_start(Functor action) {
#ifdef _OPENMP
	#pragma omp parallel
	{
		action(omp_get_thread_num());
	}
#else
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
#define calc_working_range(start, stop)					\
	size_t start = __thread_id * load;					\
	size_t stop = start + load;							\
	if (stop > VS) stop = VS;
//------------------------------------------------------------------------
static inline void vdim_setup(char *cdim, int &vx, int &vy, int &vz,
							int &vt, size_t &vp, size_t &vs)
{
	int *dim = (int *)cdim;
	vx = dim[1],
	vy = dim[2];
	vz = dim[3];
	vt = dim[4];
	vp = (size_t)vx * (size_t)vy;
	vs = (size_t)vp * (size_t)vz;
}
//------------------------------------------------------------------------
static inline void vdim_setup(char *cdim)
{
	vdim_setup(cdim, VX, VY, VZ, VT, VP, VS);
	calc_workload(VS);
}
//------------------------------------------------------------------------
// Kernel Voxel management
//------------------------------------------------------------------------
class kernel_voxel {
public:
	bool cube;
	int w, h, d;
	int radius, cmax;
	size_t size, count;
	double sigma, vmax;
private:
	typedef struct {
		int dx, dy, dz;
	} koffset;
	vector<koffset> data;
public:
	kernel_voxel() {
		w = h = d = 0;
		size = count = 0;
		cmax = radius = 0;
		cube = false;
	}
	void setup(int width) {
		cube = true;
		w = h = d = width;
		size = count = w * w * w;
		radius = w / 2;
	}
	void setup(char *cdim, void *mask) {
		int *dim = (int *)cdim;
		if (dim[0] == 1) {
			setup(dim[1]);
			return;
		}
		w = dim[1], h = dim[2], d = dim[3];
		int t = dim[4], n = dim[5];
		switch (t) {
		case  1: setup(w, h, d, n, (UC *)mask); break;
		case  2: setup(w, h, d, n, (US *)mask); break;
		case  3: setup(w, h, d, n, (UI *)mask); break;
		case  4: setup(w, h, d, n, (FL *)mask); break;
		case  5: setup(w, h, d, n, (DB *)mask); break;
		case 12: setup(w, h, d, n, (US *)mask); break;
		case 13: setup(w, h, d, n, (UI *)mask); break;
		default: setup(w, h, d, n, (UC *)mask); break;
		}
		size = n;
	}
	template<class ktype>
	void setup(int w, int h, int d, int n, ktype *mask)
	{
		int p = w * h;
		int hd = d / 2, hh = h / 2, hw = w / 2;

		data.clear();
		for (int i = 0; i < n; i++) {
			if (!mask[i]) continue;

			int k = i % p;
			koffset ko = {
				k % w - hw,
				k / w - hh,
				i / p - hd
			};
			data.push_back(ko);
		}
		count = data.size();
	}
	template<class ktype>
	void setup(int n, ktype *mask)
	{
		w = h = d = n;
		setup(n, n, n, n*n*n, mask);
	}
	inline bool operator ()(int n, int vx, int vy, int vz, int &kx, int &ky, int &kz) {
		kx = vx + data[n].dx;
		ky = vy + data[n].dy;
		kz = vz + data[n].dz;
		return kx >= 0 && kx < VX && ky >= 0 && ky < VY && kz >= 0 && kz < VZ;
	}
};
//------------------------------------------------------------------------
static kernel_voxel kvm;
//------------------------------------------------------------------------
#define for_each_kernel_voxel(kx,ky,kz,action)			\
	int x, y, z;										\
	idx2pos(idx, x, y, z);								\
	for (int nk = 0; nk < kvm.count; nk++)	{			\
		int kx, ky, kz;									\
		if (!kvm(nk, x, y, z, kx, ky, kz)) continue;	\
		action;											\
	}
//--------------------------------------------------------------------------
// Dummy function used for release DLL from IDL
//--------------------------------------------------------------------------
DLAPI int done(int argc, char *argv[])
{
	return TRUE;
}
//------------------------------------------------------------------------
