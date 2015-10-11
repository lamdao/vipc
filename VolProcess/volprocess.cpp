//------------------------------------------------------------------------
#include "dllmain.cpp"
//------------------------------------------------------------------------
#include <limits.h>
#include <float.h>
//------------------------------------------------------------------------
static void *vol = NULL;	// volume
static void *tmp = NULL;	// temporary buffer for volume (vol)
//--------------------------------------------------------------------------
// Segmented/Gradient Volume detection
//--------------------------------------------------------------------------
template<class dtype>
dtype vol_find_binvalue(dtype *vol)
{
	dtype *last = &vol[VS];
	while (*vol == 0 && vol < last) vol++;
	return vol < last ? *vol : 1;
}
//--------------------------------------------------------------------------
template<class dtype>
bool gradient_detect(dtype *vol)
{
	bool found = false;
	dtype bvalue = vol_find_binvalue(vol);

	auto thr = [&] (const int __thread_id) {
		calc_working_range(start, stop);
		for (size_t idx = start; !found && idx < stop; idx++) {
			if (vol[idx] != 0 && vol[idx] != bvalue) {
				found = true;
				break;
			}
		}
	};
	DThread::Start(thr);
	return found;
}
//--------------------------------------------------------------------------
DLAPI int vol_is_binary(int argc, char *argv[])
{
	vdim_setup(argv[0]);
	vol = argv[1];

	switch (VT % 10) {
	case 1: return gradient_detect((UC *)vol);
	case 2:	return gradient_detect((US *)vol);
	case 3: return gradient_detect((UI *)vol);
	case 4: return gradient_detect((FL *)vol);
	case 5: return gradient_detect((DB *)vol);
	}

	return 0;
}
//--------------------------------------------------------------------------
// Volume utility functions
//--------------------------------------------------------------------------
int RW, RH, RD, RT;
size_t RP, RS;
//--------------------------------------------------------------------------
// Volume Rescale
//--------------------------------------------------------------------------
static double SX = 1.0, SY = 1.0, SZ = 1.0;
//--------------------------------------------------------------------------
#define collect_voxels_in_range(rstart, rstop, vox)					\
	double rx = SX * x; bx = (int)rx;								\
	double ry = SY * y; by = (int)ry;								\
	double rz = SZ * z; bz = (int)rz;								\
	if (cx != bx || cy != by || cz != bz) {							\
		cx = bx; cy = by; cz = bz;									\
		memset(vox, 0, sizeof(vox));								\
		for (z = rstart; z <= rstop; z++) {							\
			int kz = bz + z;										\
			if (kz < 0 || kz >= RD) continue;						\
			size_t zz = (size_t)kz * RP;							\
			for (y = rstart; y <= rstop; y++) {						\
				int ky = by + y;									\
				if (ky < 0 || ky >= RH) continue;					\
				size_t yy = (size_t)ky * RW;						\
				for (x = rstart; x <= rstop; x++) {					\
					int kx = bx + x;								\
					if (kx < 0 || kx >= RW) continue;				\
					v[z-(rstart)][y-(rstart)][x-(rstart)] =			\
						src[zz+yy+kx];								\
				}													\
			}														\
		}															\
	}
//--------------------------------------------------------------------------
double DMAX[6] = {
	0,
	UCHAR_MAX,
	USHRT_MAX,
	UINT_MAX,
	FLT_MAX-1,
	DBL_MAX-1,
};
//--------------------------------------------------------------------------
#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#define M_PI2	9.86960440108935799230
//--------------------------------------------------------------------------
// Lanczos(3) interpolation 
//--------------------------------------------------------------------------
static inline double lanczos(double t)
{
	if (t >= -0.01 && t <= 0.01)
		return 1.0;
	double d = t * t * M_PI2;
	return 3 * sin(M_PI * t) * sin(M_PI * t / 3) / d;
}
//--------------------------------------------------------------------------
static inline double lanczos(double *v, double x, double y, double z)
{
	double vlx[6], vly[6];
	double *lx = vlx, *ly = vly;
	for (int i = -2; i <= 3; i++) {
		*lx++ = lanczos(x - i);
		*ly++ = lanczos(y - i);
	}

	double r = 0.0;
	for (int k = -2; k <= 3; k++) {
		double lz = lanczos(z - k);
		for (int j = 0; j < 6; j++) {
			double yz = vly[j] * lz;
			for (int i = 0; i < 6; i++) {
				r += (*v++) * vlx[i] * yz;
			}
		}
	}

	return r;
}
//--------------------------------------------------------------------------
template<class dtype>
void vol_resize_lanczos()
{
	const dtype dmax = (dtype)DMAX[VT % 10];
	dtype *src = (dtype *)vol;
	dtype *dst = (dtype *)tmp;

	auto thr = [=](const int __thread_id) {
		double v[6][6][6];
		int cx = -1, cy = -1, cz = -1;

		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int x, y, z, bx, by, bz;
			idx2pos(idx, x, y, z);
			collect_voxels_in_range(-2, 3, v);
			double dv = lanczos(&v[0][0][0], rx - bx, ry - by, rz - bz);
			dst[idx] = (dtype)(dv < 0 ? 0 : (dv > dmax ? dmax : dv));
		}
	};
	DThread::Start(thr);
}
//--------------------------------------------------------------------------
// Akima interpolation
// !!! Works nicely with 1-D data, but not with N-D at the moment
//--------------------------------------------------------------------------
#define N	6
static inline double akima(double *v, double t)
{
	double dm[N+4], w[N], b[N], c[N], d[N], *f1 = &dm[2], *f2 = dm;
	double mv[N+4], *md = &mv[2], *mm = mv, *mp = &mv[N+1];

	for (int i = 0; i < N-1; i++) {
		md[i] = v[i+1] - v[i];
	}
	mm[1] = 2.0 * md[0] - md[1];
	mm[0] = 2.0 * mm[1] - md[0];
	mp[0] = 2.0 * md[N-2] - md[N-3];
	mp[1] = 2.0 * mp[0] - md[N-2];
	for (int i = 0; i < N+2; i++) {
		dm[i] = fabs(mv[i+1] - mv[i]);
	}

	for (int i = 0; i < N; i++) {
		w[i] = f1[i] + f2[i];
	}
	for (int i = 0; i < N; i++) {
		b[i] = (f1[i] * mv[i+1] + f2[i] * mv[i+2]) / w[i];
	}
	for (int i = 0; i < N-1; i++) {
		c[i] = 3.0 * md[i] - 2.0 * b[i] - b[i+1];
		d[i] = b[i] + b[i+1] - 2.0 * md[i];
	}
	return ((t * d[2] + c[2]) * t + b[2]) * t + v[2];
}
//--------------------------------------------------------------------------
template<class dtype>
void vol_resize_akima()
{
	const dtype dmax = (dtype)DMAX[VT % 10];
	dtype *src = (dtype *)vol;
	dtype *dst = (dtype *)tmp;

	auto thr = [=](const int __thread_id) {
		double v[N][N][N];
		int cx = -1, cy = -1, cz = -1;

		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int x, y, z, bx, by, bz;
			idx2pos(idx, x, y, z);

			collect_voxels_in_range(-N/2+1, N/2, v);

			double vtx[N];
			for (x = 0; x < N; x++) {
				double vty[N];
				for (y = 0; y < N; y++) {
					double vtz[N];
					for (z = 0; z < N; z++) {
						vtz[z] = v[z][y][x];
					}
					vty[y] = akima(vtz, rz - bz);
				}
				vtx[x] = akima(vty, ry - by);
			}
			double dv = akima(vtx, rx - bx);
			dst[idx] = (dtype)(dv < 0 ? 0 : (dv > dmax ? dmax : dv));
		}
	};
	DThread::Start(thr);
}
#undef N
//--------------------------------------------------------------------------
// Cubic interpolation
//--------------------------------------------------------------------------
static inline double catmull_spline(double p[4], double t)
{
	double t2 = t * t, t3 = t2 * t;
	return 0.5 * (
		(2 * p[1]) + (-p[0] + p[2]) * t +
		(2 * p[0] - 5 * p[1] + 4 * p[2] - p[3]) * t2 +
		(-p[0] + 3 * p[1] - 3 * p[2] + p[3]) * t3
	);
}
//--------------------------------------------------------------------------
static inline double cubic_spline(double p[4], double t)
{
	return p[1] + 0.5 * t * (
		p[2] - p[0] + t * (
			2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] +
			t * (3.0 * (p[1] - p[2]) + p[3] - p[0])
		)
	);
}
//--------------------------------------------------------------------------
template<class dtype>
void vol_resize_tricubic()
{
	const dtype dmax = (dtype)DMAX[VT % 10];
	dtype *src = (dtype *)vol;
	dtype *dst = (dtype *)tmp;

	auto thr = [=](const int __thread_id) {
		double v[4][4][4];
		int cx = -1, cy = -1, cz = -1;

		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int x, y, z, bx, by, bz;
			idx2pos(idx, x, y, z);

			collect_voxels_in_range(-1, 2, v);

			double *vp = &v[0][0][0];
			for (x = 0; x < 4; x++) {
				double vty[4];
				for (y = 0; y < 4; y++) {
					double vtz[4];
					for (z = 0; z < 4; z++) {
						vtz[z] = v[z][y][x];
					}
					vty[y] = catmull_spline(vtz, rz - bz);
				}
				vp[x] = catmull_spline(vty, ry - by);
			}
			double dv = catmull_spline(vp, rx - bx);
			dst[idx] = (dtype)(dv < 0 ? 0 : (dv > dmax ? dmax : dv));
		}
	};
	DThread::Start(thr);
}
//--------------------------------------------------------------------------
// Linear interpolation
//--------------------------------------------------------------------------
static inline double calc_fragment(double r, int &n, int max)
{
	n = (int)r;
	if (n < max)
		return r - n;
	n = n - 1;
	return 1.0;
}
//--------------------------------------------------------------------------
template<class dtype>
void vol_resize_trilinear()
{
	dtype *src = (dtype *)vol;
	dtype *dst = (dtype *)tmp;

	const int RWS1 = RW-1, RHS1 = RH-1, RDS1 = RD-1;

	auto thr = [=](const int __thread_id) {
		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int x, y, z;
			idx2pos(idx, x, y, z);

			int px, py, pz;
			double fx = calc_fragment(SX * x, px, RWS1);
			double fy = calc_fragment(SY * y, py, RHS1);
			double fz = calc_fragment(SZ * z, pz, RDS1);
			double fr = 1.0 - fz;
			dtype *dp = &src[pz * RP + py * RW + px];
			double i1 = fr * dp[0] + fz * dp[RP]; dp += RW;
			double i2 = fr * dp[0] + fz * dp[RP]; dp +=  1;
			double j2 = fr * dp[0] + fz * dp[RP]; dp -= RW;
			double j1 = fr * dp[0] + fz * dp[RP];
			dst[idx] = (dtype)(
				(i1 * (1. - fy) + i2 * fy) * (1. - fx) +
				(j1 * (1. - fy) + j2 * fy) * fx
			);
		}
	};
	DThread::Start(thr);
}
//--------------------------------------------------------------------------
// Nearest neighbor interpolation
//--------------------------------------------------------------------------
template<class dtype>
void vol_resize_nearest()
{
	dtype *src = (dtype *)vol;
	dtype *dst = (dtype *)tmp;

	auto thr = [=](const int __thread_id) {
		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int x, y, z;
			idx2pos(idx, x, y, z);
			x = (int)(SX * x + 0.5); if (x >= RW) x = RW-1;
			y = (int)(SY * y + 0.5); if (y >= RH) y = RH-1;
			z = (int)(SZ * z + 0.5); if (z >= RD) z = RD-1;
			dst[idx] = src[z * RP + y * RW + x];
		}
	};
	DThread::Start(thr);
}
//--------------------------------------------------------------------------
#define	def_resize_func(name)							\
void apply_##name()										\
{														\
	switch (VT % 10) {									\
	case 1: vol_resize_##name<UC>(); break;				\
	case 2:	vol_resize_##name<US>(); break;				\
	case 3: vol_resize_##name<UI>(); break;				\
	case 4: vol_resize_##name<FL>(); break;				\
	case 5: vol_resize_##name<DB>(); break;				\
	}													\
}
//--------------------------------------------------------------------------
def_resize_func(nearest);
def_resize_func(trilinear);
def_resize_func(tricubic);
def_resize_func(lanczos);
def_resize_func(akima);
//--------------------------------------------------------------------------
DLAPI int vol_resize(int argc, char *argv[])
{
	vdim_setup(argv[0]); tmp = argv[1];
	vdim_setup(argv[2], RW, RH, RD, RT, RP, RS); vol = argv[3];

	SX = (double)RW / VX;
	SY = (double)RH / VY;
	SZ = (double)RD / VZ;

	int mode = *(int *)argv[4];
	switch (mode) {
		case 0: apply_nearest(); break;
		case 1: apply_trilinear(); break;
		case 2: apply_tricubic(); break;
		case 3: apply_lanczos(); break;
		case 4: apply_akima(); break;
	}
	return 0;
}
//--------------------------------------------------------------------------
// Volume Mirroring
//--------------------------------------------------------------------------
static inline int get_mirror(int r, int limit)
{
	if (r >= limit) r = limit - (r - limit + 1);
	if (r < 0) r = -r;
	return r;
}
//--------------------------------------------------------------------------
template<class dtype>
void vol_mirror(dtype *src, dtype *dst)
{
	const size_t SVP = RW * RH;
	const int RW2 = (VX - RW) / 2;
	const int RH2 = (VY - RH) / 2;
	const int RD2 = (VZ - RD) / 2;

	auto thr = [=](const int __thread_id) {
		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int x, y, z;
			idx2pos(idx, x, y, z);
			int vx = get_mirror(x - RW2, RW);
			int vy = get_mirror(y - RH2, RH);
			int vz = get_mirror(z - RD2, RD);
			dst[idx] = src[vz * SVP + vy * RW + vx];
		}
	};
	DThread::Start(thr);
}
//--------------------------------------------------------------------------
DLAPI int vol_mirror(int argc, char *argv[])
{
	int *sdim = (int *)argv[0];	// src vol dimensions
	int *rdim = (int *)argv[1];	// dst vol dimensions (result vol)
	vol = argv[2];	// src vol
	tmp = argv[3];	// dst vol
	VX = rdim[1], VY = rdim[2],	VZ = rdim[3], VT = rdim[4];
	VP = (size_t)VX * (size_t)VY;
	VS = (size_t)VP * (size_t)VZ;

	RW = sdim[1], RH = sdim[2], RD = sdim[3];

	switch (VT % 10) {
	case 1: vol_mirror((UC *)vol, (UC *)tmp); break;
	case 2:	vol_mirror((US *)vol, (US *)tmp); break;
	case 3: vol_mirror((UI *)vol, (UI *)tmp); break;
	case 4: vol_mirror((FL *)vol, (FL *)tmp); break;
	case 5: vol_mirror((DB *)vol, (DB *)tmp); break;
	}

	return 0;
}
//--------------------------------------------------------------------------
// Generate distance to center map
//--------------------------------------------------------------------------
bool centerized = true, nosqrt = false;
//--------------------------------------------------------------------------
template<class dtype>
void vol_cdist_center(dtype *vol)
{
	const int cx = VX/2, cy = VY/2, cz = VZ/2;

	auto thr = [=](const int __thread_id) {
		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int x, y, z;
			idx2pos(idx, x, y, z);
			x -= cx; y -= cy; z -= cz;
			double d = x*x + y*y + z*z;
			vol[idx] = nosqrt ? (dtype)d : (dtype)sqrt(d);
		}
	};
	DThread::Start(thr);
}
//--------------------------------------------------------------------------
template<class dtype>
void vol_cdist_shift(dtype *vol)
{
	const int cx = VX/2, cy = VY/2, cz = VZ/2;

	auto thr = [=](const int __thread_id) {
		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int x, y, z;
			idx2pos(idx, x, y, z);
			if (x >= cx) x -= VX;
			if (y >= cy) y -= VY;
			if (z >= cz) z -= VZ;
			double d = x*x + y*y + z*z;
			vol[idx] = nosqrt ? (dtype)d : (dtype)sqrt(d);
		}
	};
	DThread::Start(thr);
}
//--------------------------------------------------------------------------
void vol_cdistance()
{
	if (centerized) {
		if (VT == DT_DOUBLE)
			vol_cdist_center((DB *)vol);
		else
			vol_cdist_center((FL *)vol);
	} else {
		if (VT == DT_DOUBLE)
			vol_cdist_shift((DB *)vol);
		else
			vol_cdist_shift((DB *)vol);
	}
}
//--------------------------------------------------------------------------
DLAPI int vol_cdistance(int argc, char *argv[])
{
	vdim_setup(argv[0]);
	vol = (void *)argv[1];
	nosqrt = ((int *)argv[2])[0] > 0;
	centerized = ((int *)argv[2])[1] > 0;

	vol_cdistance();
	return 0;
}
//--------------------------------------------------------------------------
template<class dtype>
void vol_hmcalc(dtype *vol, double *p)
{
	double D0 = p[0], H = p[1], L = p[2], C = -p[3];

	H -= L;
	D0 *= D0;

	auto thr = [=](const int __thread_id) {
		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			vol[idx] = (dtype)(H * (1 - exp(C * vol[idx] / D0)) + L);
		}
	};
	DThread::Start(thr);
}
//--------------------------------------------------------------------------
DLAPI int vol_homomorphic(int argc, char *argv[])
{
	vdim_setup(argv[0]);
	vol = (void *)argv[1];
	nosqrt = ((int *)argv[2])[0] > 0;
	centerized = ((int *)argv[2])[1] > 0;

	vol_cdistance();
	if (VT == DT_DOUBLE)
		vol_hmcalc((DB *)vol, (double *)argv[3]);
	else
		vol_hmcalc((FL *)vol, (double *)argv[3]);

	return 0;
}
//--------------------------------------------------------------------------
