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
