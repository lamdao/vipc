//--------------------------------------------------------------------------
#include "dllmain.cpp"
//--------------------------------------------------------------------------
#define VOL(x,y,z)	vol[pos2idx((x),(y),(z))]
//--------------------------------------------------------------------------
// Morphology
//--------------------------------------------------------------------------
namespace Morphology
{
	bool multiop = false;
	//----------------------------------------------------------------------
	vector<ushort> buffer;
	//----------------------------------------------------------------------
	ushort *vol = NULL;	// volume
	ushort *tmp = NULL;	// temporary buffer for volume (vol)
	//----------------------------------------------------------------------
	// Utility
	//----------------------------------------------------------------------
	void SetupParameters(char **params, bool multiop = false)
	{
		vol = (ushort *)params[0];
		ushort *kmask = (ushort *)params[1];

		int *dim = (int *)params[2];
		VX = dim[0];
		VY = dim[1];
		VZ = dim[2];
		VP = (size_t)VX * (size_t)VY;
		VS = (size_t)VP * (size_t)VZ;
		DThread::CalcWorkload(VS);

		if (buffer.size() != VS) {
			vector<ushort> buf(VS);
			buffer.swap(buf);
		}
		tmp = &buffer[0];

		kvm.setup(dim[3], kmask);

		Morphology::multiop = multiop;
	}
	//----------------------------------------------------------------------
	void UpdateBuffer()
	{
		if (multiop) {
			ushort *t = tmp;
			tmp = vol;
			vol = t;
		} else {
			memcpy(vol, tmp, VS * sizeof(ushort));
		}
	}
	//----------------------------------------------------------------------
	// Gray morphological operations
	//----------------------------------------------------------------------
	namespace GrayScale {
		ushort vmax = 0;
		void FindMaxValue()
		{
			vector<ushort> vx(NUM_THREADS);
			DThread::Start([&vx](const int __thread_id) {
				ushort &m = vx[__thread_id];
				calc_working_range(start, stop);
				for (size_t idx = start; idx < stop; idx++) {
					if (m < vol[idx])
						m = vol[idx];
				}
			});

			vmax = 0;
			for (auto v : vx) {
				if (v > vmax)
					vmax = v;
			}
		}
		void Dilate()
		{
			DThread::Start([](const int __thread_id) {
				calc_working_range(start, stop);
				for (size_t idx = start; idx < stop; idx++) {
					int smax = 0;
					for_each_kernel_voxel(kx, ky, kz, {
						int d = (int)VOL(kx, ky, kz);
						if (d >= vmax) {
							smax = vmax;
							goto set_gdilate_value;
						} else if (d > smax) {
							smax = d;
						}
					});
				set_gdilate_value:
					tmp[idx] = smax;
				}
			});

			UpdateBuffer();
		}
		//------------------------------------------------------------------
		void Erode()
		{
			DThread::Start([](const int __thread_id) {
				calc_working_range(start, stop);
				for (size_t idx = start; idx < stop; idx++) {
					int smin = vmax;
					for_each_kernel_voxel(kx, ky, kz, {
						int d = (int)VOL(kx, ky, kz);
						if (d <= 0) {
							smin = 0;
							goto set_gerode_value;
						} else if (d < smin) {
							smin = d;
						}
					});
				set_gerode_value:
					tmp[idx] = smin;
				}
			});
			UpdateBuffer();
		}
	};
	//----------------------------------------------------------------------
	// Bin morphological operations
	//----------------------------------------------------------------------
	namespace Binary {
		void Dilate()
		{
			DThread::Start([](const int __thread_id) {
				calc_working_range(start, stop);
				for (size_t idx = start; idx < stop; idx++) {
					for_each_kernel_voxel(kx, ky, kz, {
						int d = (int)VOL(kx, ky, kz);
						if (d > 0) {
							tmp[idx] = d;
							goto check_next_voxel;
						}
					});
					tmp[idx] = 0;
				check_next_voxel:;
				}
			});
			UpdateBuffer();
		}
		//------------------------------------------------------------------
		void Erode() {
			DThread::Start([](const int __thread_id) {
				calc_working_range(start, stop);
				for (size_t idx = start; idx < stop; idx++) {
					for_each_kernel_voxel(kx, ky, kz, {
						int d = (int)VOL(kx, ky, kz);
						if (d == 0) {
							tmp[idx] = 0;
							goto check_next_voxel;
						}
					});
					tmp[idx] = vol[idx];
				check_next_voxel:;
				}
			});
			UpdateBuffer();
		}
	};
};
//--------------------------------------------------------------------------
// Library interfaces
//--------------------------------------------------------------------------
DLAPI int p_morph_gdilate(int argc, char *argv[])
{
	Morphology::SetupParameters(argv);
	Morphology::GrayScale::FindMaxValue();
	Morphology::GrayScale::Dilate();
	return TRUE;
}
//--------------------------------------------------------------------------
DLAPI int p_morph_gerode(int argc, char *argv[])
{
	Morphology::SetupParameters(argv);
	Morphology::GrayScale::FindMaxValue();
	Morphology::GrayScale::Erode();
	return TRUE;
}
//--------------------------------------------------------------------------
DLAPI int p_morph_gopen(int argc, char *argv[])
{
	Morphology::SetupParameters(argv, true);
	Morphology::GrayScale::FindMaxValue();
	Morphology::GrayScale::Erode();
	Morphology::GrayScale::Dilate();
	return TRUE;
}
//--------------------------------------------------------------------------
DLAPI int p_morph_gclose(int argc, char *argv[])
{
	Morphology::SetupParameters(argv, true);
	Morphology::GrayScale::FindMaxValue();
	Morphology::GrayScale::Dilate();
	Morphology::GrayScale::Erode();
	return TRUE;
}
//--------------------------------------------------------------------------
DLAPI int p_morph_bdilate(int argc, char *argv[])
{
	Morphology::SetupParameters(argv);
	Morphology::Binary::Dilate();
	return TRUE;
}
//--------------------------------------------------------------------------
DLAPI int p_morph_berode(int argc, char *argv[])
{
	Morphology::SetupParameters(argv);
	Morphology::Binary::Erode();
	return TRUE;
}
//--------------------------------------------------------------------------
DLAPI int p_morph_bopen(int argc, char *argv[])
{
	Morphology::SetupParameters(argv, true);
	Morphology::Binary::Erode();
	Morphology::Binary::Dilate();
	return TRUE;
}
//--------------------------------------------------------------------------
DLAPI int p_morph_bclose(int argc, char *argv[])
{
	Morphology::SetupParameters(argv, true);
	Morphology::Binary::Dilate();
	Morphology::Binary::Erode();
	return TRUE;
}
//--------------------------------------------------------------------------
// Multiplex functions
//--------------------------------------------------------------------------
#define mfx(op)															\
DLAPI int p_morph_##op(int argc, char *argv[])							\
{																		\
	bool g = ((int *)argv[2])[5] != 0;									\
	return g ? p_morph_g##op(argc, argv) : p_morph_b##op(argc, argv);	\
}
//--------------------------------------------------------------------------
mfx(dilate);
mfx(erode);
mfx(open);
mfx(close);
//--------------------------------------------------------------------------
