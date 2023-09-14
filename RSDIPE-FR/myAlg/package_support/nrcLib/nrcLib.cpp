// nrcLib.cpp : Defines the entry point for the DLL application.
//

#include "stdafx.h"
#include "nrcLib.h"

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    switch (ul_reason_for_call)
	{
		case DLL_PROCESS_ATTACH:
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH:
			break;
    }
    return TRUE;
}


// This is an example of an exported variable
NRCLIB_API int nNrcLib=0;

// This is an example of an exported function.
NRCLIB_API int fnNrcLib(void)
{
	return 42;
}

// This is the constructor of a class that has been exported.
// see nrcLib.h for the class definition
CNrcLib::CNrcLib()
{ 
	return; 
}

