//--------------------------------------------------------------------------
// vdim.h - Volume dimensions initialization/conversion utilities
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
#ifndef __VOLUME_UTILS_H
#define __VOLUME_UTILS_H
//--------------------------------------------------------------------------
// Volume info and utility functions
//--------------------------------------------------------------------------
static int VX, VY, VZ, VT;	// volume dimensions (X,Y,Z) and type
//--------------------------------------------------------------------------
static size_t VP, VS;		// volume plane size (VP), volume size (VS)
//--------------------------------------------------------------------------
// Convert from a position [x,y,z] to a linear index in volume(VX,VY,VZ)
//--------------------------------------------------------------------------
static inline size_t pos2idx(int x, int y, int z)
{
	return (size_t)z * VP + (size_t)y * VX + (size_t)x;
}
//--------------------------------------------------------------------------
// Convert from a linear index to a position [x,y,z] in volume(VX,VY,VZ)
//--------------------------------------------------------------------------
static inline void idx2pos(size_t n, int &x, int &y, int &z)
{
	z = (int)(n / VP); n = n % VP;
	y = (int)(n / VX);
	x = (int)(n % VX);
}
//--------------------------------------------------------------------------
// Volume dimensions initialize from outside caller
//--------------------------------------------------------------------------
static inline void vdim_setup(void *cdim, int &vx, int &vy, int &vz,
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
//--------------------------------------------------------------------------
static inline void vdim_setup(void *cdim)
{
	vdim_setup(cdim, VX, VY, VZ, VT, VP, VS);
	DThread::CalcWorkload(VS);
}
//--------------------------------------------------------------------------
#endif
