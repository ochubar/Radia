/*	These data types are used by WaveMetrics to support code that will be
	compiled for both 32-bit and 64-bit executables.
	
	This header is included by XOPStandardHeaders.h.
	
	The long type should be avoided because it is 32 bits on WIN64 and 64 bits on MAC64.
*/

#ifndef WM_TYPES_H

#define WM_TYPES_H

/*	These types are defined on Macintosh in MacTypes.h. We define them here so that
	we can use the same types regardless of platform.	
*/
#ifdef WINIGOR
	typedef signed char SInt8;
	typedef unsigned char UInt8;
	typedef signed short SInt16;
	typedef unsigned short UInt16;
	typedef signed long SInt32;
	typedef unsigned long UInt32;
	typedef signed __int64 SInt64;
	typedef unsigned __int64 UInt64;
#endif

#define LegacyLong long			// A long that must remain a long because, for example, it is used with a system API that requires a long.

typedef unsigned int TickCountInt;	// Holds a tick count. A tick is roughly 1/60th of a second.

/*	The remaining types are 32-bit values when compiling a 32-bit executable
	and 64-bit values when compiling a 64-bit executable.
*/

#ifdef MACIGOR
	#ifdef IGOR64					// 64-bits on Macintosh
		#error("64-bits is not yet supported on Macintosh")
	#else							// 32-bits on Macintosh
		typedef long PSInt;				// Pointer-sized int - an int that may hold a pointer
		typedef long BCInt;				// A byte count. Can exceed 2^31-1 when running 64-bits.
		typedef long IndexInt;			// An index used to index off of a pointer. Also used to hold a wave point number. Can exceed 2^31-1 when running 64-bits.
		typedef long CountInt;			// A count that should be 32 bits in a 32-bit app and 64-bits in a 64-bit app such as the size of a wave dimension.
	#endif
#endif

#ifdef WINIGOR
	#ifdef IGOR64					// 64-bits on Windows
		typedef __int64 PSInt;			// Pointer-sized int - an int that may hold a pointer
		typedef __int64 BCInt;			// A byte count. Can exceed 2^31-1 when running 64-bits.
		typedef __int64 IndexInt;		// An index used to index off of a pointer. Also used to hold a wave point number. Can exceed 2^31-1 when running 64-bits.
		typedef __int64 CountInt;		// A count that should be 32 bits in a 32-bit app and 64-bits in a 64-bit app such as the size of a wave dimension.
	#else							// 32-bits on Windows
		typedef long PSInt;				// Pointer-sized int - an int that may hold a pointer
		typedef long BCInt;				// A byte count. Can exceed 2^31-1 when running 64-bits.
		typedef long IndexInt;			// An index used to index off of a pointer. Also used to hold a wave point number. Can exceed 2^31-1 when running 64-bits.
		typedef long CountInt;			// A count that should be 32 bits in a 32-bit app and 64-bits in a 64-bit app such as the size of a wave dimension.
	#endif
#endif

#endif		// WM_TYPES_H
