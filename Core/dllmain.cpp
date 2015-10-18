//--------------------------------------------------------------------------
// dllmain.cpp - Common entry for shared library
//--------------------------------------------------------------------------
// Author: Lam H. Dao <daohailam(at)yahoo(dot)com>
//--------------------------------------------------------------------------
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//------------------------------------------------------------------------
#include "typedefs.h"
#include "dthread.h"
#include "vdim.h"
#include "kvm.h"
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
// Dummy function used for release DLL from IDL
//--------------------------------------------------------------------------
DLAPI int done(int argc, char *argv[])
{
	return TRUE;
}
//------------------------------------------------------------------------
