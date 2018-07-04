#ifndef XOP_STANDARD_HEADERS						// Skip if XOP standard headers have already been included
	#define XOP_STANDARD_HEADERS 1
	#define _IGORXOP_ 1
	
	// Create platform-identifying macros
	#if defined(TARGET_OS_MAC)
		#define MACIGOR			// Defined for 32-bit and 64-bit Igor on Macintosh
		#ifdef __LP64__
			#define IGOR64		// Defined for 64-bit Igor only
		#else
			#define IGOR32		// Defined for 32-bit Igor only
		#endif
	#endif
	#if (defined(_WIN32) || defined(_WIN64))
		#define WINIGOR			// Defined for 32-bit and 64-bit Igor on Macintosh
		#ifdef _WIN64
			#define IGOR64		// Defined for 64-bit Igor only
		#else
			#define IGOR32		// Defined for 32-bit Igor only
		#endif
	#endif
	
	#ifdef MACIGOR				// Compiling for Macintosh [
		#include <ctype.h>
		#include <string.h>
		#include <stdlib.h>
		#include <stdio.h>
		#include <stddef.h>
		// #include <math.h>		// This is included via fp.h in the precompiled headers.
		
		#ifdef __LITTLE_ENDIAN__
			#define XOP_LITTLE_ENDIAN	// We are compiling for Intel Macintosh.
		#endif

		typedef struct RECT {			// Windows RECT required for WinRectToMacRect and MacRectToWinRect in XOPSupport.c.
			long left;
			long top;
			long right;
			long bottom;
		} RECT;

		#ifdef __cplusplus
			#define HOST_EXPORT extern "C"		// Declares function or variable as being exported from IGOR to XOPs.
			#define HOST_IMPORT extern "C"		// Declares function or variable as being imported by IGOR from an XOP.
		#else
			#define HOST_EXPORT					// Declares function or variable as being exported from IGOR to XOPs.
			#define HOST_IMPORT					// Declares function or variable as being imported by IGOR from an XOP.
		#endif		
	#endif						// Compiling for Macintosh ]

	#ifdef WINIGOR				// Compiling for Windows [
		#include <Windows.h>		// This creates the _WINDOWS_ symbol.
		
		#ifdef SetPort				// SetPort is defined in WinSpool.h
			#undef SetPort			// But we use SetPort in the Macintosh sense.
		#endif
		
		#define XOP_LITTLE_ENDIAN	// HR, 030130: Was LITTLE_ENDIAN but Apple stole that from us (see endian.h).

		#undef DUPLICATE		// This is defined in XOP.h but also in NB30 for use with NetBios 3.0.

		/*	HOST_EXPORT declares a routine or variable exported from the host (IGOR) to one or more DLLs (XOPs)
			In IGOR, HOST_EXPORT is defined as __declspec(dllexport). However, in an XOP, HOST_EXPORT is
			defined as __declspec(dllimport). This allows IGOR and an XOP to share a prototype file
			where routines exported from IGOR are marked as HOST_EXPORT.
			
			HOST_IMPORT marks a routine imported by IGOR from a DLL. Currently, this is not
			actually used.
		*/
		#ifdef __cplusplus
			#define HOST_EXPORT extern "C" __declspec(dllimport)	// Declares function or variable as being exported from IGOR to XOPs.
			#define HOST_IMPORT extern "C" __declspec(dllexport)	// Declares function or variable as being imported by IGOR from an XOP.
		#else
			#define HOST_EXPORT __declspec(dllimport)				// Declares function or variable as being exported from IGOR to XOPs.
			#define HOST_IMPORT __declspec(dllexport)				// Declares function or variable as being imported by IGOR from an XOP.
		#endif		
		
		#include <ctype.h>
		#include <string.h>
		#include <stdlib.h>
		#include <stdio.h>
		#include <stddef.h>
		#include <math.h>
	#endif						// Compiling for Windows ]

	#include "WMTypes.h"		// Data types defined by WaveMetrics
	
	#ifdef WINIGOR
		#ifdef WM_WINMAC_SUPPORT
			/*	This provides support used by elaborate WaveMetrics XOPs such as the
				the Data Browser and Gizmo. This support is not available outside WaveMetrics because
				XOPs that use it might need to be recompiled each time we revised Igor itself.
			*/
			#include "XOPWMWinMacSupport.h"
		#else
			// This provides support for a small number of Macintosh routines that are needed by Windows XOPs.
			#include "XOPWinMacSupport.h"
		#endif
	#endif
	
	#include "IgorXOP.h"
	#include "IgorErrors.h"
	#include "XOP.h"
	#include "XOPSupport.h"
	#ifdef MACIGOR
		#include "XOPSupportMac.h"
	#endif
	#ifdef WINIGOR
		#include "XOPSupportWin.h"
	#endif
#endif				// XOP_STANDARD_HEADERS
