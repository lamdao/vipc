//------------------------------------------------------------------------
#include "dllmain.cpp"
#include "barrier.h"
//------------------------------------------------------------------------
static void *vol = NULL;	// volume
static void *tmp = NULL;	// temporary buffer for volume (vol)
static void *ext = NULL;
//--------------------------------------------------------------------------
#include "convolution.cpp"
//--------------------------------------------------------------------------
// Volume kernel filters
//--------------------------------------------------------------------------
#define	calc_box_start(idx)								\
		int x, y, z;									\
		size_t pdx = idx;								\
		idx2pos(pdx, x, y, z);							\
		int nx, ny, nz;									\
		int bz = z - kvm.radius;						\
		int by = y - kvm.radius;						\
		int bx = x - kvm.radius;						\
		nx = ny = nz = kvm.w;							\
		if (bz < 0) {									\
			nz += bz;									\
			bz = 0;										\
		}												\
		if (bz < z) {									\
			pdx -= VP * (z - bz);						\
		}												\
		if (by < 0) {									\
			ny += by;									\
			by = 0;										\
		}												\
		if (by < y) {									\
			pdx -= VX * (y - by);						\
		}												\
		if (bx < 0) {									\
			nx += bx;									\
			bx = 0;										\
		}												\
		if (bx < x) {									\
			pdx -= (x - bx);							\
		}
//--------------------------------------------------------------------------
#define foreach_box_position(action)					\
	while (bz < VZ && nz > 0) {							\
		register size_t cdx = pdx;						\
		int cy = by, dy = 0;							\
		for (; dy < ny && cy < VY; dy++, cy++) {		\
			int cx = bx, dx = 0;						\
			for (; dx < nx && cx < VX; dx++, cx++)		\
				action;									\
			cdx += VX;									\
		}												\
		pdx += VP;										\
		bz += 1;										\
		nz -= 1;										\
	}
//--------------------------------------------------------------------------
#define get_current_values(bsrc)	(bsrc)[cdx+dx]
//--------------------------------------------------------------------------
#define execute(fx,...)									\
	switch (VT % 10) {									\
	case  1: fx<UC>(__VA_ARGS__); break;				\
	case  2: fx<US>(__VA_ARGS__); break;				\
	case  3: fx<UI>(__VA_ARGS__); break;				\
	case  4: fx<FL>(__VA_ARGS__); break;				\
	case  5: fx<DB>(__VA_ARGS__); break;				\
	}
//--------------------------------------------------------------------------
// Kernel filter
static void vol_kfilter_setup(int argc, char **argv)
{
	vdim_setup(argv[0]);
	vol = (void *)argv[1];
	tmp = (void *)argv[2];
	kvm.setup(argv[3], argv[4]);
	ext = (void *)argv[5];
}
//--------------------------------------------------------------------------
// Iterative filter
static void vol_ifilter_setup(int argc, char **argv)
{
	vdim_setup(argv[0]);
	vol = (void *)argv[1];
	ext = (void *)argv[2];
}
//--------------------------------------------------------------------------
// Volume Median filter
//--------------------------------------------------------------------------
template<class dtype>
dtype qselect(register dtype *x, int n)
{
	#define swap(a, b) {dtype t = x[a];x[a] = x[b];x[b] = t;}
	int k = n / 2;
	int l = 0;
	int r = n - 1;
	while (l < r) {
		dtype v = x[k];
		swap(k, r);
		int p = l;
		for (int i = l; i < r; i++) {
			if (x[i] < v) {
				swap(i, p);
				p += 1;
			}
		}
		swap(r, p);
		if (p == k) break;
		if (p < k)
			l = p + 1;
		else
			r = p - 1;
	}
	return x[k];
}
//--------------------------------------------------------------------------
template<class dtype>
void median_kernel_filter()
{
	dtype *dres = (dtype *)tmp;
	dtype *dsrc = (dtype *)vol;
	thread_start([=](const int __thread_id) {
		vector<dtype> dtmp(kvm.size);
		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int count = 0;
			for_each_kernel_voxel(kx, ky, kz, {
				dtmp[count++] = dsrc[pos2idx(kx, ky, kz)];
			});
			dres[idx] = qselect(&dtmp[0], count);
		}
	});
}
//--------------------------------------------------------------------------
template<class dtype>
void median_box_filter()
{
	dtype *dres = (dtype *)tmp;
	dtype *dsrc = (dtype *)vol;
	thread_start([=](const int __thread_id){
		vector<dtype> dtmp(kvm.size);
		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int count = 0;
			calc_box_start(idx);
			foreach_box_position({
				dtmp[count++] = get_current_values(dsrc);
			});
			dres[idx] = qselect(&dtmp[0], count);
		};
	});
}
//--------------------------------------------------------------------------
DLAPI int vol_median_filter(int argc, char *argv[])
{
	vol_kfilter_setup(argc, argv);
	if (kvm.cube) {
		execute(median_box_filter);
	}
	else {
		execute(median_kernel_filter);
	}
	return 0;
}
//--------------------------------------------------------------------------
// Calculate median using cdf
//--------------------------------------------------------------------------
static inline int cdfappend(int *cdf, int value, int &count, int prev=-1)
{
	int last = 3 * count;
	cdf[last  ] = value;
	cdf[last+1] = prev;
	cdf[last+2] = 1;
	count += 1;
	return last;
}
//--------------------------------------------------------------------------
void cdfadd(int *cdf, int value, int &tail, int &count)
{
	if (tail < 0) {
		tail = cdfappend(cdf, value, count);
		return;
	}

	int last;
	int n = tail;
	int p = -1;
	while (cdf[n] > value) {
		if (cdf[n+1] < 0) {	// reach to head
			cdf[n+1] = cdfappend(cdf, value, count);
			return;
		}
		p = n;
		n = cdf[n+1];
	}
	if (cdf[n] == value) {
		cdf[n+2] += 1;
		count++;
		return;
	}
	last = cdfappend(cdf, value, count, n);
	if (p >= 0)
		cdf[p+1] = last;
	else
		tail = last;
}
//--------------------------------------------------------------------------
static inline int cdfsel(int *cdf, int tail, int N)
{
	int p = 0;
	N = N / 2 + (N & 1);
	while (tail >= 0) {
		p += cdf[tail+2];
		if (p >= N)
			break;
		tail = cdf[tail+1];
	}
	return cdf[tail];
}
//--------------------------------------------------------------------------
template<class dtype>
void mediancdf_box_filter()
{
	dtype *dres = (dtype *)tmp;
	dtype *dsrc = (dtype *)vol;
	auto thr = [=](const int __thread_id){
		vector<int> dtmp(3 * kvm.size);
		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			calc_box_start(idx);

			int count = 0, tail = -1;
			foreach_box_position({
				cdfadd(&dtmp[0], (int)get_current_values(dsrc), tail, count);
			});
			dres[idx] = (dtype)cdfsel(&dtmp[0], tail, count);
		}
	};
	thread_start(thr);
}
//--------------------------------------------------------------------------
DLAPI int vol_mediancdf_filter(int argc, char *argv[])
{
	vol_kfilter_setup(argc, argv);
	execute(mediancdf_box_filter);
	return 0;
}
//--------------------------------------------------------------------------
// Volume mean filter
//--------------------------------------------------------------------------
template<class dtype>
void mean_kernel_filter()
{
	dtype *dres = (dtype *)tmp;
	dtype *dsrc = (dtype *)vol;

	thread_start([=](const int __thread_id) {
		calc_working_range(start, stop);
		for (size_t idx=start; idx < stop; idx++) {
			int count = 0;
			double d = 0.0;
			for_each_kernel_voxel(kx, ky, kz, {
				d += dsrc[pos2idx(kx, ky, kz)];
				count++;
			});
			dres[idx] = (dtype)(d / count);
		}
	});
}
//--------------------------------------------------------------------------
template<class dtype>
void mean_box_filter()
{
	dtype *dres = (dtype *)tmp;
	dtype *dsrc = (dtype *)vol;

	thread_start([=](const int __thread_id){
		calc_working_range(start, stop);
		for (size_t idx = start; idx < stop; idx++) {
			int count = 0;
			double d = 0.0;
			calc_box_start(idx);
			foreach_box_position({
				d += get_current_values(dsrc);
				count++;
			});
			dres[idx] = (dtype)(d / count);
		}
	});
}
//--------------------------------------------------------------------------
DLAPI int vol_mean_filter(int argc, char *argv[])
{
	vol_kfilter_setup(argc, argv);
	if (kvm.cube) {
		execute(mean_box_filter);
	}
	else {
		execute(mean_kernel_filter);
	}
	return 0;
}
//--------------------------------------------------------------------------
// Sigma Filter
//--------------------------------------------------------------------------
#define sigma_filter_main(collector)									\
{																		\
	dtype *dsrc = (dtype *)vol;											\
	dtype *dtmp = (dtype *)tmp;											\
	vector<double> data(kvm.count);										\
	calc_working_range(start, stop);									\
	for (size_t idx = start; idx < stop; idx++) {						\
		int n, count = 0;												\
		double total = 0.0;												\
		double value = (double)dsrc[idx];								\
		double range = collector(dsrc, idx, total, &data[0], count);	\
		double vmax = value + range;									\
		double vmin = value - range;									\
		double pass = 0.0;												\
		for (n = count-1, count = 0; n >= 0; n--) {						\
			double d = data[n];											\
			if (d < vmin || d > vmax) continue;							\
			pass += d;													\
			count++;													\
		}																\
		value = (count >= kvm.cmax ?									\
						(pass/count) :									\
						(total-value)/(kvm.count-1)) + 0.5;				\
		dtmp[idx] = (dtype)(value > kvm.vmax ? kvm.vmax : value);		\
	}																	\
}
//--------------------------------------------------------------------------
template<class dtype>
double calc_kernel_var_range(dtype *vol, size_t idx, double &sum, double *data, int &count)
{
	double var = 0.0;
	for_each_kernel_voxel(kx, ky, kz, {
		double d = (double)vol[pos2idx(kx, ky, kz)];
		data[count++] = d;
		var += d * d;
		sum += d;
	});
	double mean = sum / count;
	return kvm.sigma * sqrt(var / count - mean * mean);
}
//--------------------------------------------------------------------------
template<class dtype>
void sigma_kernel_filter()
{
	thread_start([](const int __thread_id){
		sigma_filter_main(calc_kernel_var_range);
	});
}
//--------------------------------------------------------------------------
template<class dtype>
double calc_box_var_range(dtype *vol, size_t idx, double &sum, double *data, int &count)
{
	double var = 0.0;
	calc_box_start(idx);
	foreach_box_position({
		double d = (double)get_current_values(vol);
		data[count++] = d;
		var += d * d;
		sum += d;
	});
	double mean = sum / count;
	return kvm.sigma * sqrt(var / count - mean * mean);
}
//--------------------------------------------------------------------------
template<class dtype>
void sigma_box_filter()
{
	thread_start([](const int __thread_id){
		sigma_filter_main(calc_box_var_range);
	});
}
//--------------------------------------------------------------------------
void sigma_filter_init()
{
	double *sfp = (double *)ext;
	kvm.sigma = sfp[0];
	kvm.cmax = (int)(sfp[1] * kvm.count + 0.99999);
	kvm.vmax = sfp[2];
}
//--------------------------------------------------------------------------
DLAPI int vol_sigma_filter(int argc, char *argv[])
{
	vol_kfilter_setup(argc, argv);

	sigma_filter_init();
	if (kvm.cube) {
		execute(sigma_box_filter);
	}
	else {
		execute(sigma_kernel_filter);
	}
	return 0;
}
//--------------------------------------------------------------------------
// Utility functions for iterative filters (anisotropic, laplace filter)
//--------------------------------------------------------------------------
double VK = 180.0, DT = 1.0 / 14;
//--------------------------------------------------------------------------
template<class dtype>
void ifilter_malloc()
{
	tmp = malloc(VS * sizeof(dtype));
}
//--------------------------------------------------------------------------
int ifilter_prepare()
{
	int niters = 1;
	if (ext) {
		double *p = (double *)ext;
		VK = p[0];
		DT = p[1];
		niters = (int)p[2];
		if (niters < 1) niters = 1;
	}

	switch (VT % 10) {
	case  1: ifilter_malloc<UC>(); break;
	case  2: ifilter_malloc<US>(); break;
	case  3: ifilter_malloc<UI>(); break;
	case  4: ifilter_malloc<FL>(); break;
	case  5: ifilter_malloc<DB>(); break;
	default: ifilter_malloc<DB>(); break;
	}
	return niters;
}
//--------------------------------------------------------------------------
void ifilter_cleanup(int niters)
{
	if ((niters & 1) == 1) {
		int vs[] = { 0,
			sizeof(UC),
			sizeof(US),
			sizeof(UI),
			sizeof(FL),
			sizeof(DB),
		};
		memcpy(vol, tmp, VS * vs[VT % 10]);
	}
	free(tmp);
}
//--------------------------------------------------------------------------
// Anisotropic filter
//--------------------------------------------------------------------------
#define gflow(v)	((v) * exp(-pow((v)/VK,2)))
//--------------------------------------------------------------------------
template<class dtype>
void anisotropic_filter(int niters)
{
	Barrier barrier(NUM_THREADS);

	auto thr = [&barrier, niters](const int __thread_id){
		dtype *dres = (dtype *)tmp;
		dtype *dsrc = (dtype *)vol;

		calc_working_range(start, stop);
		for (int ni = 1; ni <= niters; ni++) {
			for (size_t idx = start; idx < stop; idx++) {
				int x, y, z;
				idx2pos(idx, x, y, z);

				double v = dsrc[idx];
				double e = x >    0 ? dsrc[idx - 1] - v : -v;
				double w = x < VX - 1 ? v - dsrc[idx + 1] : v;
				double n = y >    0 ? dsrc[idx - VX] - v : -v;
				double s = y < VY - 1 ? v - dsrc[idx + VX] : v;
				double u = z >    0 ? dsrc[idx - VP] - v : -v;
				double d = z < VZ - 1 ? v - dsrc[idx + VP] : v;
				dres[idx] = (dtype)(v + DT*(gflow(e) - gflow(w) + gflow(n) - gflow(s) + gflow(u) - gflow(d)));
			}
			barrier.wait(__thread_id);
			if (ni < niters) {
				dtype *t = dsrc;
				dsrc = dres;
				dres = t;
			}
		}
	};
	thread_start(thr);
}
//--------------------------------------------------------------------------
DLAPI int vol_anisotropic_filter(int argc, char *argv[])
{
	vol_ifilter_setup(argc, argv);
	int niters = ifilter_prepare();
	execute(anisotropic_filter, niters);
	ifilter_cleanup(niters);
	return 0;
}
//--------------------------------------------------------------------------
// Laplace filter
//--------------------------------------------------------------------------
#define SQRT2	1.414213562373095145474621858739
#define SQRT3	1.732050807568877193176604123437
//--------------------------------------------------------------------------
static double wlaplace[] = {
	1.0 / SQRT3, 1.0 / SQRT2, 1.0 / SQRT3,
	1.0 / SQRT2, 1.0, 1.0 / SQRT2,
	1.0 / SQRT3, 1.0 / SQRT2, 1.0 / SQRT3,

	1.0 / SQRT2, 1.0, 1.0 / SQRT2,
	1.0, 0.0, 1.0,
	1.0 / SQRT2, 1.0, 1.0 / SQRT2,

	1.0 / SQRT3, 1.0 / SQRT2, 1.0 / SQRT3,
	1.0 / SQRT2, 1.0, 1.0 / SQRT2,
	1.0 / SQRT3, 1.0 / SQRT2, 1.0 / SQRT3,
};
//--------------------------------------------------------------------------
template<class dtype>
void laplace_filter(int niters)
{
	Barrier barrier(NUM_THREADS);
	auto thr = [&barrier, niters](const int __thread_id) {
		dtype *dres = (dtype *)tmp;
		dtype *dsrc = (dtype *)vol;

		calc_working_range(start, stop);
		for (int ni = 1; ni <= niters; ni++) {
			for (size_t idx = start; idx < stop; idx++) {
				int x, y, z;
				idx2pos(idx, x, y, z);

				int n = 0;
				double v = 0.0;
				double *wl = wlaplace;
				for (int dz = -1; dz <= 1; dz++, wl += 9) {
					int kz = z + dz;
					if (kz < 0 || kz >= VZ)
						continue;
					double *wt = wl;
					for (int dy = -1; dy <= 1; dy++, wt += 3) {
						int ky = y + dy;
						if (ky < 0 || ky >= VY)
							continue;
						double *w = wt;
						for (int dx = -1; dx <= 1; dx++, w++) {
							int kx = x + dx;
							if (*w != 0.0 && kx >= 0 && kx < VX) {
								v += *w * dsrc[pos2idx(kx, ky, kz)];
								n++;
							}
						}
					}
				}
				dres[idx] = (dtype)(dsrc[idx] + DT * (v - (double)n * dsrc[idx]) / 2.0);
			}
			barrier.wait(__thread_id);

			if (ni < niters) {
				dtype *t = dsrc;
				dsrc = dres;
				dres = t;
			}
		}
	};

	thread_start(thr);
}
//--------------------------------------------------------------------------
DLAPI int vol_laplace_filter(int argc, char *argv[])
{
	vol_ifilter_setup(argc, argv);
	int niters = ifilter_prepare();
	execute(laplace_filter, niters);
	ifilter_cleanup(niters);
	return 0;
}
//--------------------------------------------------------------------------
