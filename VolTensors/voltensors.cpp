//--------------------------------------------------------------------------
#include "dllmain.cpp"
#include "barrier.h"
#include "eig3.h"
//--------------------------------------------------------------------------
static int VX1, VY1, VZ1;
//--------------------------------------------------------------------------
static float *vol = NULL;	// volume
static float *edv = NULL;	// eigen values of each voxel in volume
static float *evx = NULL;	// eigen vectors of each voxel in volume
//--------------------------------------------------------------------------
static void SetupParameters(int argc, char **argv)
{
	vdim_setup(argv[0]);
	vol = (float *)argv[1];
	edv = (float *)argv[2];
	evx = (float *)argv[3];	// can be NULL if don't need

	VX1 = VX - 1;
	VY1 = VY - 1;
	VZ1 = VZ - 1;
}
//--------------------------------------------------------------------------
static inline void vxdist(double *v, float *evx)
{
	static int t[] = {0,3,6,1,4,7,2,5,8};
	for (int i = 0; i < 9; i++) *evx++ = (float)v[t[i]];
}
//--------------------------------------------------------------------------
static inline void dvdist(double *d, float *edv)
{
	edv[0] = (float)d[0];
	edv[1] = (float)d[1];
	edv[2] = (float)d[2];
}
//--------------------------------------------------------------------------
#define calc_dx(vol, dx, start, stop) {									\
	for(size_t idx = start; idx < stop; idx++) {						\
		int x, y, z;													\
		idx2pos(idx, x, y, z);											\
		if (x > 0 && y > 0 && z > 0 && x < VX1 && y < VY1 && z < VZ1)	\
			dx[idx] = (vol[idx+1] - vol[idx-1]) / 2.0;					\
	}																	\
}
//--------------------------------------------------------------------------
#define calc_dy(vol, dy, start, stop) {									\
	for(size_t idx = start; idx < stop; idx++) {						\
		int x, y, z;													\
		idx2pos(idx, x, y, z);											\
		if (x > 0 && y > 0 && z > 0 && x < VX1 && y < VY1 && z < VZ1)	\
			dy[idx] = (vol[idx+VX] - vol[idx-VX]) / 2.0;				\
	}																	\
}
//--------------------------------------------------------------------------
#define calc_dz(vol, dz, start, stop) {									\
	for(size_t idx = start; idx < stop; idx++) {						\
		int x, y, z;													\
		idx2pos(idx, x, y, z);											\
		if (x > 0 && y > 0 && z > 0 && x < VX1 && y < VY1 && z < VZ1)	\
			dz[idx] = (vol[idx+VP] - vol[idx-VP]) / 2.0;				\
	}																	\
}
//--------------------------------------------------------------------------
void CalcTensors()
{
	vector<float> tmp(VS * 7);
	Barrier barrier(NUM_THREADS);

	DThread::Start([&barrier, &tmp](const int __thread_id) {
		float *dxx = &tmp[VS * 1];
		float *dxy = &tmp[VS * 2];
		float *dxz = &tmp[VS * 3];
		float *dyy = &tmp[VS * 4], *dyx = dxy;
		float *dyz = &tmp[VS * 5];
		float *dzz = &tmp[VS * 6], *dzy = dyz, *dzx = dxz;

		calc_working_range(start, stop);

		calc_dx(vol, tmp, start, stop); barrier.wait();
		calc_dx(tmp, dxx, start, stop);
		calc_dy(tmp, dxy, start, stop);
		calc_dz(tmp, dxz, start, stop);

		calc_dy(vol, tmp, start, stop); barrier.wait();
		calc_dy(tmp, dyy, start, stop);
		calc_dz(tmp, dyz, start, stop);

		calc_dz(vol, tmp, start, stop); barrier.wait();
		calc_dz(tmp, dzz, start, stop);

		barrier.wait();
		for (size_t idx = start; idx < stop; idx++) {
			double hx[9], vx[9], d[3];
			hx[0] = dxx[idx];
			hx[1] = dxy[idx];
			hx[2] = dxz[idx];
			hx[3] = dyx[idx];
			hx[4] = dyy[idx];
			hx[5] = dyz[idx];
			hx[6] = dzx[idx];
			hx[7] = dzy[idx];
			hx[8] = dzz[idx];
			eigen_decomposition(hx, vx, d);
			dvdist(d, &edv[3*idx]);
			if (evx != NULL)
				vxdist(vx, &evx[9*idx]);
		}
	});
}
//--------------------------------------------------------------------------
// A faster way is using approximate method to calculate
// Hessian matrix (H) from Jacobian vector (J)
// J=[dx,dy,dz]  H = J'J
//--------------------------------------------------------------------------
void CalcApproximateTensors()
{
	DThread::Start([&](const int __thread_id) {
		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int x, y, z;
			idx2pos(idx, x, y, z);
			if (x <= 0 || y <= 0 || z <= 0 || x >= VX1 || y >= VY1 || z >= VZ1)
				continue;

			double dx = (vol[idx+ 1] - vol[idx- 1]) / 2.0;
			double dy = (vol[idx+VX] - vol[idx-VX]) / 2.0;
			double dz = (vol[idx+VP] - vol[idx-VP]) / 2.0;
			double hx[9] = {dx * dx, dx * dy, dx * dz,
							dy * dx, dy * dy, dy * dz,
							dz * dx, dz * dy, dz * dz};
			double vx[9], d[3];
			eigen_decomposition(hx, vx, d);
			dvdist(d, &edv[3*idx]);
			if (evx != NULL) {
				vxdist(vx, &evx[9*idx]);
			}
		}
	});
}
//--------------------------------------------------------------------------
DLAPI int vol_calc_tensors(int argc, char *argv[])
{
	SetupParameters(argc, argv);
	CalcApproximateTensors();
	return 0;
}
//--------------------------------------------------------------------------
