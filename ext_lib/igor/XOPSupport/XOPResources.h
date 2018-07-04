/*	XOPResources.h -- equates used in creating XOP-specific resources.

	HR, 10/6/96
		Created this file from parts of XOPTypes.r.
		XOPTypes.r is needed for Macintosh only.
		This file is needed for Macintosh and Windows.
*/

/*	XOP protocol version Code

	Use this in first field of the XOPI resource.
	See discussion in XOP.h for details.
*/
#define XOP_VERSION 4

// XOP Toolkit version. Use this in the fifth field of the XOPI resource.
#define XOP_TOOLKIT_VERSION 600			// 600 means XOP Toolkit 6.00.


// Development System Codes - Use DEV_SYS_CODE in the second field of the XOPI resource.
#define DEV_SYS_CODE_ALL_OTHERS 0		// Use this for any type of development not listed below.
#define DEV_SYS_CODE_MAC_MPW 1			// Obsolete.
#define DEV_SYS_CODE_MAC_MACH 2			// Use this for Macintosh Mach-O XOPs.
#ifdef MACIGOR
	#define DEV_SYS_CODE DEV_SYS_CODE_MAC_MACH
#else
	#define DEV_SYS_CODE DEV_SYS_CODE_ALL_OTHERS
#endif



/* Operation Categories - Used in the XOPC resource. */
#define displayOp 1						// Other Windows, (ie not graph, table, or layout)
#define waveOp 2						// About Waves
#define dataOp 4						// Wave Analysis
#define cursorOp 8						// Controls and cursors
#define utilOp 16						// Programming and utilities
#define XOPOp 32						// External operation
#define compilableOp 64					// XOP supports Operation Handler and can be used in user function. Supported in Igor Pro 5 or later. Ignored by earlier versions.
#define graphOp 128						// Graph window operation
#define tableOp 256						// Table window operation
#define layoutOp 512					// Layout window operation
#define allWinOp 1024					// All windows operation
#define drawOp 2048						// Drawing tools operation
#define ioOp 4096						// Input/output operation
#define printOp 8192					// Printing operation
#define statsOp 16384					// Statistics operation		// Recognized by Igor Pro 6.02 or later
#define threadSafeOp 32768				// Thread safe operation		// Recognized by Igor Pro 6.02 or later

/* Flags used in XMN1 menu flags field (copied from XOP.h) */
#define SHOW_MENU_AT_LAUNCH 1				// For XMN1 resource menuFlags field */
#define SHOW_MENU_WHEN_ACTIVE 2				// For XMN1 resource menuFlags field */


/* Flags used in XMI1 menu item flags field (copied from XOP.h) */
#define ITEM_REQUIRES_WAVES 1				// For XMI1 resource itemFlags field
#define ITEM_REQUIRES_GRAPH 2				// For XMI1 resource itemFlags field
#define ITEM_REQUIRES_TABLE 4				// For XMI1 resource itemFlags field
#define ITEM_REQUIRES_LAYOUT 8				// For XMI1 resource itemFlags field
#define ITEM_REQUIRES_GRAPH_OR_TABLE 16		// For XMI1 resource itemFlags field
#define ITEM_REQUIRES_TARGET 32				// For XMI1 resource itemFlags field
#define ITEM_REQUIRES_PANEL 64				// For XMI1 resource itemFlags field
#define ITEM_REQUIRES_NOTEBOOK 128			// For XMI1 resource itemFlags field
#define ITEM_REQUIRES_GRAPH_OR_PANEL 256	// For XMI1 resource itemFlags field
#define ITEM_REQUIRES_DRAW_WIN 512			// For XMI1 resource itemFlags field
#define ITEM_REQUIRES_PROC_WIN 1024			// For XMI1 resource itemFlags field


/* Function categories - Used in the XOPF resource. */
#define F_TRIG 1							// trig functions
#define F_EXP 2								// exponential functions
#define F_SPEC 4							// special functions
#define F_CMPLX 8							// complex functions
#define F_TIMEDATE 0x10			/* 16 */	// time and date functions (was F_COMPARE before 5/28/93)
#define F_ROUND 0x20			/* 32 */	// rounding functions
#define F_CONV 0x40				/* 64 */	// conversion functions
#define F_WAVE 0x80				/* 128 */	// wave functions
#define F_UTIL 0x100			/* 256 */	// Programming and utilities functions
#define F_NUMB 0x200			/* 512 */	// number functions (e.g. Pi)
#define F_ANLYZWAVES 0x400		/* 1024 */	// wave analysis (CurveFit, Convolve, etc)
#define F_IO 0x800				/* 2048 */	// input/output (files, movies, FIFOs, Paths)
#define F_WINDOWS 0x1000		/* 4096 */	// graphs and other windows
#define F_EXTERNAL 0x2000		/* 8192 */	// external functions (set automatically by Igor for XFUNCs)
#define F_THREADSAFE 0x4000		/* 16384 */	// Marks external function as thread safe for Igor Pro 6 or later.
#define F_STR 0x8000			/* 32768 */	// strings functions


/*	Number types for use in the result or parameter fields of an XOPF resource.
	For example:
		NT_FP64,			// A double-precision floating point result or parameter
		NT_FP64 | NT_CMPLX,	// A double-precision, complex floating point result or parameter
		NT_FP32,			// OBSOLETE. Do not use for result or parameter.
		HSTRING_TYPE		// A string result or parameter.
		WAVE_TYPE			// A wave parameter
*/
#define NT_CMPLX 1							// Complex numbers
#define NT_FP32 2							// 32 bit fp numbers (OBSOLETE)
#define NT_FP64 4							// 64 bit fp numbers
#define WAVE_TYPE 0x4000					// Special type that signifies that a parameter is a wave.
#define HSTRING_TYPE 0x2000					// Special type that signifies that a parameter is a string.
#define FV_REF_TYPE 0x1000					// Signifies pass-by-reference. Used for numeric and string parameters only.
#define FV_STRUCT_TYPE 0x0200				// Signifies structure parameter. Must be used with FV_REF_TYPE.
#define DATAFOLDER_TYPE	0x100				// Data folder reference (requires Igor Pro 6.20 or later)


