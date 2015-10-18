//--------------------------------------------------------------------------
// convolution.cpp - Volume convolution
//--------------------------------------------------------------------------
// Author: Lam H. Dao <daohailam(at)yahoo(dot)com>
//--------------------------------------------------------------------------
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//--------------------------------------------------------------------------
#include "fftw3.h"
//------------------------------------------------------------------------
#define DC	fftw_complex
#define FC	fftwf_complex
//------------------------------------------------------------------------
#define FFTW_DInitThreads()								\
if (fftw_init_threads()) {								\
	fftw_plan_with_nthreads(NUM_THREADS);				\
}
#define FFTW_FInitThreads()								\
if (fftwf_init_threads()) {								\
	fftwf_plan_with_nthreads(NUM_THREADS);				\
}
//--------------------------------------------------------------------------
// Volume Convolution
//--------------------------------------------------------------------------
template<class ctype, class dtype>
void Convolve(ctype *v, ctype *k, size_t total)
{
	for (size_t n = 0; n < total; n++) {
		dtype *a = v[n];
		dtype *b = k[n];
		dtype r = a[0] * b[0] - a[1] * b[1];
		dtype i = a[0] * b[1] + a[1] * b[0];
		a[0] = r;
		a[1] = i;
	}
}
//--------------------------------------------------------------------------
template<class dtype>
void Normalize(dtype *vol, size_t total)
{
	for (size_t n = 0; n < total; n++) {
		vol[n] /= VS;
	}
}
//--------------------------------------------------------------------------
void vol_fconvol(float *vol, float *ker, int vc)
{
	FFTW_FInitThreads();
	fftwf_complex *v = fftwf_alloc_complex(vc);
	fftwf_complex *k = fftwf_alloc_complex(vc);
	fftwf_plan fwd_v2c = fftwf_plan_dft_r2c_3d(VZ, VY, VX, vol, v, FFTW_ESTIMATE);
	fftwf_plan fwd_k2c = fftwf_plan_dft_r2c_3d(VZ, VY, VX, ker, k, FFTW_ESTIMATE);
	fftwf_plan fwd_c2r = fftwf_plan_dft_c2r_3d(VZ, VY, VX, v, vol, FFTW_ESTIMATE);
	fftwf_execute(fwd_v2c); fftwf_destroy_plan(fwd_v2c);
	fftwf_execute(fwd_k2c); fftwf_destroy_plan(fwd_k2c);
	Convolve<FC, FL>(v, k, vc);
	fftwf_execute(fwd_c2r); fftwf_destroy_plan(fwd_c2r);
	Normalize(vol, VS);
	fftw_free(k);
	fftw_free(v);
}
//--------------------------------------------------------------------------
void vol_dconvol(double *vol, double *ker, int vc)
{
	FFTW_DInitThreads();
	fftw_complex *v = fftw_alloc_complex(vc);
	fftw_complex *k = fftw_alloc_complex(vc);
	fftw_plan fwd_v2c = fftw_plan_dft_r2c_3d(VZ, VY, VX, vol, v, FFTW_ESTIMATE);
	fftw_plan fwd_k2c = fftw_plan_dft_r2c_3d(VZ, VY, VX, ker, k, FFTW_ESTIMATE);
	fftw_plan fwd_c2r = fftw_plan_dft_c2r_3d(VZ, VY, VX, v, vol, FFTW_ESTIMATE);
	fftw_execute(fwd_v2c); fftw_destroy_plan(fwd_v2c);
	fftw_execute(fwd_k2c); fftw_destroy_plan(fwd_k2c);
	Convolve<DC, DB>(v, k, vc);
	fftw_execute(fwd_c2r); fftw_destroy_plan(fwd_c2r);
	Normalize(vol, VS);
	fftw_free(k);
	fftw_free(v);
}
//--------------------------------------------------------------------------
DLAPI int vol_convol(int argc, char *argv[])
{
	vdim_setup(argv[0]);
	int vc = (VX / 2 + 1) * VY * VZ;
	if (VT == DT_DOUBLE) {
		vol_dconvol((DB *)argv[1], (DB *)argv[2], vc);
	} else if (VT == DT_FLOAT) {
		vol_fconvol((FL *)argv[1], (FL *)argv[2], vc);
	} else {
		fprintf(stderr, "Data type is not support at the moment.\n");
		fflush(stderr);
	}
	return 0;
}
//--------------------------------------------------------------------------
// Apply mask/weight in frequency domain
//--------------------------------------------------------------------------
template<class ctype, class dtype>
void ApplyMask(ctype *v, dtype *k, size_t total)
{
	for (size_t n = 0; n < total; n++) {
		dtype *a = v[n];
		a[0] *= *k;
		a[1] *= *k;
		k++;
	}
}
//--------------------------------------------------------------------------
void vol_flt_maskfreq(float *vol, float *ker, int vc)
{
	FFTW_FInitThreads();
	fftwf_complex *v = fftwf_alloc_complex(vc);
	fftwf_plan fwd_v2c = fftwf_plan_dft_r2c_3d(VZ, VY, VX, vol, v, FFTW_ESTIMATE);
	fftwf_plan fwd_c2r = fftwf_plan_dft_c2r_3d(VZ, VY, VX, v, vol, FFTW_ESTIMATE);
	fftwf_execute(fwd_v2c); fftwf_destroy_plan(fwd_v2c);
	ApplyMask(v, ker, vc);
	fftwf_execute(fwd_c2r); fftwf_destroy_plan(fwd_c2r);
	fftw_free(v);

	Normalize(vol, VS);
}
//--------------------------------------------------------------------------
void vol_dbl_maskfreq(double *vol, double *ker, int vc)
{
	FFTW_DInitThreads();
	fftw_complex *v = fftw_alloc_complex(vc);
	fftw_plan fwd_v2c = fftw_plan_dft_r2c_3d(VZ, VY, VX, vol, v, FFTW_ESTIMATE);
	fftw_plan fwd_c2r = fftw_plan_dft_c2r_3d(VZ, VY, VX, v, vol, FFTW_ESTIMATE);
	fftw_execute(fwd_v2c); fftw_destroy_plan(fwd_v2c);
	ApplyMask(v, ker, vc);
	fftw_execute(fwd_c2r); fftw_destroy_plan(fwd_c2r);
	fftw_free(v);

	Normalize(vol, VS);
}
//--------------------------------------------------------------------------
DLAPI int vol_maskfreq(int argc, char *argv[])
{
	vdim_setup(argv[0]);
	int vc = (VX / 2 + 1) * VY * VZ;
	if (VT == DT_DOUBLE) {
		vol_dbl_maskfreq((DB *)argv[1], (DB *)argv[2], vc);
	}
	else if (VT == DT_FLOAT) {
		vol_flt_maskfreq((FL *)argv[1], (FL *)argv[2], vc);
	} else {
		fprintf(stderr, "Data type is not support at the moment.\n");
		fflush(stderr);
	}
	return 0;
}
//--------------------------------------------------------------------------
