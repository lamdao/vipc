#ifndef __TYPEDEFS_H
#define __TYPEDEFS_H
//--------------------------------------------------------------------------
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
//--------------------------------------------------------------------------
#endif
