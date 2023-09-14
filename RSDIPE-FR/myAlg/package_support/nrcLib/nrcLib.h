
// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the NRCLIB_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// NRCLIB_API functions as being imported from a DLL, wheras this DLL sees symbols
// defined with this macro as being exported.
#ifdef NRCLIB_EXPORTS
#define NRCLIB_API __declspec(dllexport)
#else
#define NRCLIB_API __declspec(dllimport)
#endif

// This class is exported from the nrcLib.dll
class NRCLIB_API CNrcLib {
public:
	CNrcLib(void);
	// TODO: add your methods here.
};

extern NRCLIB_API int nNrcLib;

NRCLIB_API int fnNrcLib(void);

NRCLIB_API void jacobi(float **a, int n, float d[], float **v, int *nrot);