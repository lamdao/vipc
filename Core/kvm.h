#ifndef __KERNEL_VOXEL_H
#define __KERNEL_VOXEL_H
//--------------------------------------------------------------------------
// Kernel Voxel management
//--------------------------------------------------------------------------
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
//--------------------------------------------------------------------------
static kernel_voxel kvm;
//--------------------------------------------------------------------------
#define for_each_kernel_voxel(kx,ky,kz,action)			\
	int x, y, z;										\
	idx2pos(idx, x, y, z);								\
	for (int nk = 0; nk < kvm.count; nk++)	{			\
		int kx, ky, kz;									\
		if (!kvm(nk, x, y, z, kx, ky, kz)) continue;	\
		action;											\
	}
//--------------------------------------------------------------------------
#endif
