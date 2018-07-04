// IgorXOP.h -- Miscellaneous equates for interfacing XOP to Igor.

/*	These equates come from .h files used in compiling Igor and include
	various information that an XOP might need in communicating with Igor.
*/

#pragma pack(2)	// All structures passed between Igor and XOP are two-byte aligned.

// Miscellaneous
#define TUStuffHandle Handle
#define waveHndl Handle
#define DataFolderHandle Handle


// From DataFolderBits.h. Used by Data Browser to determine events of interest to it.
#define kDFCB_NewChildFolder	1
#define kDFCB_KillChildFolder	2
#define kDFCB_RenameChildFolder	4
#define kDFCB_NewWave			8
#define kDFCB_KillWave			0x10
#define kDFCB_RenameWave		0x20
#define kDFCB_NewVariable		0x40
#define kDFCB_KillVariable		0x80
#define kDFCB_RenameVariable	0x100
#define kDFCB_NewString			0x200
#define kDFCB_KillString		0x400
#define kDFCB_RenameString		0x800
#define kDFCD_LockWave			0x1000					// AG27MAY03, Igor Pro 5


// From WinGenMacs.c. Used for DoWindowRecreationDialog callback.
enum CloseWinAction {
	kCloseWinCancel = 0,
	kCloseWinSave = 1,
	kCloseWinReplace = 2,
	kCloseWinNoSave = 3
};


// From Igor.h
#ifndef NIL
	#define NIL 0L
#endif
#ifndef FALSE				// Conditional compilation is needed because Metrowerks
	#define FALSE 0			// defines TRUE and FALSE in their MacHeaders.c.
#endif
#ifndef TRUE
	#define TRUE -1
#endif

#define MAXCMDLEN 400		// HR, 10/2/93 -- changed from 200 to 400 for Igor 2.0.

#define WAVE_OBJECT 1
#define WAVEARRAY_OBJECT 2
#define VAR_OBJECT 3
#define STR_OBJECT 4
#define STRUCT_OBJECT 5
#define XOPTARGWIN_OBJECT 5		// HR, 980714. Igor Pro 3.13B03.
#define GRAPH_OBJECT 6
#define TABLE_OBJECT 7
#define LAYOUT_OBJECT 8
#define PANEL_OBJECT 9			// HR, 10/2/93. Igor Pro 2.0.
#define NOTEBOOK_OBJECT 10
#define DATAFOLDER_OBJECT 11	// HR, 7/7/95. Igor Pro 3.0.
#define PATH_OBJECT 12			// HR, 7/28/95. Igor Pro 3.0. Symbolic path.
#define PICT_OBJECT 13			// HR, 7/28/95. Igor Pro 3.0. Picture.
#define ANNOTATION_OBJECT 14	// JP, 4/24/98. Igor Pro 3.13.
#define CONTROL_OBJECT 15		// JP, 4/24/98. Igor Pro 3.13.

// #defines for identifying windows.
#define GRAF_MASK 1
#define SS_MASK 2
#define PL_MASK 4
#define PICT_MASK 8
#define MW_MASK 16
#define TEXT_MASK 32
#define PANEL_MASK 64
#define PROC_MASK 128
#define MOVIE_MASK 256
#define HELP_MASK 512
#define XOP_TARGET_MASK 4096			// HR, 980706: Added for XOP target windows.
#define ALL_MASK -1

#define FIRSTWINCODE CMDWIN
#define CMDWIN 1
#define WMDIALOGWIN 2				// HR, 11/19/92 -- was 10.
#define OLD_PROCWIN 2				// HR, 11/19/92 -- PROCWIN is now 10.
#define GRAFWIN 3
#define SSWIN 4
#define PLWIN 5
#define PICTWIN 6
#define MWWIN 7						// Notebook window (may be plain or formatted text).
#define TEXTWIN 8
#define PANELWIN 9
#define PROCWIN 10					// HR, 11/19/92 -- was 2.
#define MOVIEWIN 11
#define HELPWIN 12
#define HELPDLOGWIN 13
#define XOPWIN 14
#define XOPTARGWIN 15				// To group all XOP target windows together in the windows menu.
#define LASTWINCODE	XOPTARGWIN


/*	Name space codes.
	For use with UniqueName2 callback (Igor Pro 2.0 or later)
	and CheckName callback (Igor Pro 3.0 or later).
*/
#define MAIN_NAME_SPACE 1			// Igor's main name space (waves, variables, windows).
#define PATHS_NAME_SPACE 2			// Igor's symbolic path name space.
#define PICTS_NAME_SPACE 3			// Igor's picture name space.
#define WINDOWS_NAME_SPACE 4		// Igor's windows name space. Added in Igor Pro 3.0.
#define DATAFOLDERS_NAME_SPACE 5	// Igor's data folders name space. Added in Igor Pro 3.0.

// From IgorMenus.h
// These are the menu IDs that you can use in XMI1 resources to attach menu items to Igor menus.
#define APPLEID 1
#define FILEID 2
#define EDITID 3
#define WAVEFORMID 4
#define DATAID 4					// HR, 10/2/93 -- old "Waves" menu is now called "Data"
#define ANALYSISID 5
#define MACROID 6
#define WINDOWSID 7
#define MISCID 8
#define LAYOUTID 10					// HR, 10/2/93 -- this was 9 prior to Igor 2.0
#define GRAPHID 12					// HR, 10/2/93 -- added for Igor Pro 2.0
#define PANELID 13					// HR, 10/2/93 -- added for Igor Pro 2.0
#define TABLEID 14					// HR, 10/2/93 -- added for Igor Pro 2.0
#define PROCEDUREID 15				// HR, 10/2/93 -- added for Igor Pro 2.0
#define NOTEBOOKID 16				// HR, 10/2/93 -- added for Igor Pro 2.0
#define LOAD_SUB_ID 50
#define SAVE_SUB_ID 51
// #define SAVEGRAPHICS_SUB_ID 52	// HR, 981105: The Save Graphics submenu was removed in Igor Pro 3.1.
#define OPEN_FILE_SUB_ID 55
#define CONTROL_WIN_SUB_ID 56
#define NEW_WIN_SUB_ID 58
#define MISC_OPS_SUB_ID 59
#define APPEND_TO_GRAPH_SUBID 89	// HR, 3/2/96 -- added for Igor Pro 3.0

#define T_COMMA		1			// terminator codes
#define T_RPAREN 	2
#define T_SEMI		4
#define T_RBRACK	8
#define T_RCBRACE	16
#define T_NORM		(T_COMMA | T_SEMI)
#define T_CRP		(T_COMMA | T_RPAREN)
#define T_CRB		(T_COMMA | T_RCBRACE)
#define T_CRBRACK	(T_COMMA | T_RBRACK)	// LH, 3/18/90


// From IgorMath.h
#define NT_CMPLX 1				// complex numbers
#define NT_FP32 2				// 32 bit fp numbers
#define NT_FP64 4				// 64 bit fp numbers
#define NT_I8 8					// 8 bit signed integer (changed 1/21/91)
#define NT_I16 	0x10			// 16 bit integer numbers
#define NT_I32 	0x20			// 32 bit integer numbers
#define NT_UNSIGNED 0x40		// Makes above signed integers unsigned. NOTE: Requires Igor Pro 3.0 or later.
#define DATAFOLDER_TYPE 0x100	// Data type for DFREF waves (waves containing data folder references)
#define WAVE_TYPE 0x4000		// Data type for wave-reference waves (waves containing wave references)


// From wave.h
#define MAX_WAVE_NAME 31		// maximum length of wave name -- not including the null
								//	NOTE: Prior to Igor 3.0, this was 18 and we recommended that you use MAX_OBJ_NAME (31) instead of MAX_WAVE_NAME.
							
#define OLD_MAX_WAVE_NAME 18	// MAX_WAVE_NAME prior to Igor Pro 3.0.

#define TEXT_WAVE_TYPE 0		// The wave type code for text waves. Added in Igor Pro 3.0.

#define MAX_DIMENSIONS 10		// Maximum number of dimensions in a multi-dimension object.
								// In Igor 3.0, the max is actually 4 but this may increase in the future.
#define ROWS 0					// Dimension 0 is the row dimension.
#define COLUMNS 1				// Dimension 1 is the column dimension.
#define LAYERS 2				// Dimension 2 is the layer dimension.
#define CHUNKS 3				// Dimension 3 is the chunk dimension.

#define MAX_UNIT_CHARS 49		// Max number of characters in a units string, not including trailing null, in Igor Pro 3.0 or later.
								// Prior to Igor Pro 3.0, the maximum was 3 characters.

#define MAX_DIM_LABEL_CHARS 31	// Max chars in a dimension label, not including trailing null.
								
#define kMDWaveAccessMode0 0	// Access code for MDAccessNumericWaveData. Used by Igor for future compatibility check.

// From WM.h
#define UNKNOWNCURSOR 0
#define ARROWCURSOR 1
#define WATCHCURSOR 2
#define IBEAMCURSOR 3
#define HANDCURSOR 4
#define SPINNINGCURSOR 5
#define CROSSHAIRCURSOR 6
#define MAX_OBJ_NAME 31			// maximum length of: variables,macros,annotations.
#define MAX_LONG_NAME 255		// Added in 6.00D00. Used for double names.

#ifdef MACIGOR
	// HR, 080728: XOP Toolkit 5.09 - Support long file names on Macintosh.
	#define MAX_VOLUMENAME_LEN 255			// Maximum length of volume name
	#define MAX_DIRNAME_LEN 255				// Maximum length of directory name
	#define MAX_FILENAME_LEN 255			// Maximum length of file name
	#define MAX_PATH_LEN 511				// Maximum length of path name. This was 511 in Igor 3.0 so I am leaving it as 511 even though it is not clear whether the Mac OS support more than 255.
#endif
#ifdef WINIGOR
	#define MAX_VOLUMENAME_LEN 255			// Maximum length of volume name (e.g., "C:")
	#define MAX_DIRNAME_LEN 255				// Maximum length of directory name
	#define MAX_FILENAME_LEN 255			// maximum length of file name
	#define MAX_PATH_LEN 259				// maximum length of path name
#endif


/*	This is used to select one of two ways of doing something or to allow both.
	For an example, see VolumeNameLength() in WMFileUtils.c.
*/
typedef enum PlatformCode {
	kMacPlatform=1,				// This is stored on disk. The value must not be changed.
	kWinPlatform,				// This is stored on disk. The value must not be changed.
	kMacOrWinPlatform,
	#ifdef MACIGOR
		kCurrentPlatform=kMacPlatform
	#endif
	#ifdef WINIGOR
		kCurrentPlatform=kWinPlatform
	#endif
} PlatformCode;

#define CR_STR "\015"			// Can be used as follows: XOPNotice("Test"CR_STR);
#define CR_CHAR '\015'
#define LF_STR "\012"
#define LF_CHAR '\012'

// From Functions.h

// These are used to identify parameters to external functions as waves, data folder references, strings or names
#define WAVE_TYPE 0x4000		// Added to number types above to signify parameter is wave
#define DATAFOLDER_TYPE 0x100	// Signifies parameter is a DFREF
#define HSTRING_TYPE 0x2000		// Signifies parameter is a handle to a string

// These are used to test parameter types returned by GetUserFunctionInfo.
#define FV_REF_TYPE 0x1000		// Signifies pass-by-reference
#define FV_FUNC_TYPE 0x0400		// Signifies a function reference
#define WAVE_Z_TYPE	0x8000		// Identifies WAVE/Z argument
#define FV_STRUCT_TYPE 0x0200	// Requires Igor Pro 5.03 or later

struct NVARRec {				// Used for NVAR structure fields.
	PSInt urH;					// 32 bits in Igor32, 64 bits in Igor64
	IndexInt index;				// 32 bits in Igor32, 64 bits in Igor64
};
typedef struct NVARRec NVARRec;

struct SVARRec {				// Used for SVAR structure fields.
	PSInt urH;					// 32 bits in Igor32, 64 bits in Igor64
	IndexInt index;				// 32 bits in Igor32, 64 bits in Igor64
};
typedef struct SVARRec SVARRec;

typedef PSInt FUNCREF;			// Used for FUNCREF structure fields.

struct FunctionInfo {			// Used by GetUserFunctionInfo.
	char name[MAX_OBJ_NAME+1];
	int compilationIndex;
	int functionID;
	int subType;
	int isExternalFunction;
	int returnType;
	int reserved[25];					// Do not use. Reserved for future use.
	int numOptionalParameters;
	int numRequiredParameters;
	int totalNumParameters;
	int parameterTypes[100];
};
typedef struct FunctionInfo FunctionInfo;
typedef FunctionInfo* FunctionInfoPtr;

typedef void* UserFunctionThreadInfoPtr;

// From TextUtils.h

// structure for getting info about text utility document
struct TUDocInfo {			// 10/23/93: added for Igor Pro 2.0D83
	short version;						// version number of this structure
	short permission;					// 0 = read only, 1 = read/write
	short fileType;						// for future use
	int paragraphs;						// total number of paragraphs in document
	char reserved[256];					// for future use
};
typedef struct TUDocInfo TUDocInfo;
typedef struct TUDocInfo *TUDocInfoPtr;
#define TUDOCINFO_VERSION 1

struct TULoc {							// identifies a location in a text utility document
	int paragraph;						// location's paragraph
	unsigned short pos;					// character offset in paragraph for text paragraph
};
typedef struct TULoc TULoc;
typedef struct TULoc *TULocPtr;

/*	When to erase message in status area.
	This is a bitwise parameter used with the TUSetStatus() callback.
	The status is always changed if a new message comes along.
	This controls if and when it will be erased before a new message comes.
*/
#define TU_ERASE_STATUS_NEVER 0
#define TU_ERASE_STATUS_WHEN_SELECTION_CHANGES 1
#define TU_ERASE_STATUS_WHEN_WINDOW_ACTIVATED 2
#define TU_ERASE_STATUS_WHEN_WINDOW_DEACTIVATED 4
#define TU_ERASE_STATUS_WHEN_DOC_MODIFIED 8
#define TU_ERASE_STATUS_WHEN_ANYTHING_HAPPENS -1


// From CmdWin.h
// modes for PutCmdLine()
#define INSERTCMD 1					// insert text at current insertion point
#define FIRSTCMD 2					// insert text in front of cmd buffer
#define FIRSTCMDCRHIT 3				// insert text in front of cmd buffer and set crHit
#define REPLACEFIRSTCMD 4			// replace first line of cmd buffer with text
#define REPLACEALLCMDSCRHIT 5 		// replace all lines of cmd buffer with text and set crHit
#define REPLACEALLCMDS 6			// replace all lines of cmd buffer with text


// From ColorTables.h
typedef Handle IgorColorTableHandle;
typedef struct IgorColorSpec {
	UInt32 value;					// index or other value
	RGBColor rgb;					// true color
} IgorColorSpec;


// From CommandUtils.h
// structure for getting and passing wave range information
struct WaveRangeRec {
	// the following fields are set by GetWaveRange based on command line
	double x1, x2;	// *** 4/21/90 -- changed from float to double
	int rangeMode;					// bit 0 set if start specified, bit 1 set if end specified
	int isBracket;					// true if range specified by [] instead of ()
	int gotRange;					// true if /R=range was present

	// next, you setup these fields
	Handle waveHandle;
	CountInt minPoints;				// min number of points in acceptable range
	
	// Then, following fields are setup by CalcWaveRange
	IndexInt p1, p2;
	int wasBackwards;				// truth p1 > p2 before being swapped
};
typedef struct WaveRangeRec WaveRangeRec;
typedef struct WaveRangeRec *WaveRangeRecPtr;


// From Variables.h
#define VAR_GLOBAL	0x4000				// bit flag for type parameter of Variable XOP callback

struct NumVarValue{						// Used in Get/SetDataFolderObject call.
	int numType;			// NT_FP64 possibly ORed with NT_CMPLX (if variable is complex).
	int spare;				// For future use - set to zero.
	double realValue;
	double imagValue;
};
typedef struct NumVarValue NumVarValue;
typedef NumVarValue* NumVarValuePtr;

union DataObjectValue {
	waveHndl waveH;						// Use this if the object is a wave.
	NumVarValue nv;						// Use this if the object is a numeric variable.
	Handle strH;						// Use this if the object is a string variable.
	DataFolderHandle dfH;				// Use this if the object is a data folder.
	char spare[64];						// For possible future use.
};
typedef union DataObjectValue DataObjectValue;
typedef DataObjectValue* DataObjectValuePtr;


// From Save.h
#define SAVE_TYPE_SAVE 1				// experiment save type codes
#define SAVE_TYPE_SAVEAS 2
#define SAVE_TYPE_SAVEACOPY 3
#define SAVE_TYPE_STATIONERY 4

#define LOAD_TYPE_NEW 1					// experiment load type codes
#define LOAD_TYPE_OPEN 2
#define LOAD_TYPE_REVERT 3
#define LOAD_TYPE_STATIONERY 4
#define LOAD_TYPE_MERGE 5				// Added in Igor Pro 5.

#define EXP_UNPACKED 0					// experiment file type codes
#define EXP_PACKED 1


// From OperationHandler.h
#define kUseCMDMessageForInterpreting 1
#define kOperationIsThreadSafe 2		// HR, 070507: Pass to RegisterOperation if operation is thread safe.

struct DataFolderAndName {
	DataFolderHandle dfH;
	char name[MAX_OBJ_NAME+1];
};
typedef struct DataFolderAndName DataFolderAndName;
typedef struct DataFolderAndName *DataFolderAndNamePtr;

// Options used with GetOperationDestWave.
enum {
	kOpDestWaveOverwriteOK = 1,
	kOpDestWaveChangeExistingWave = 2,
	kOpDestWaveOverwriteExistingWave = 4,
	kOpDestWaveMakeFreeWave = 8,
	kOpDestWaveMustAlreadyExist = 16
};

struct WaveRange {
	waveHndl waveH;
	double startCoord;					// Start point number or x value
	double endCoord;					// End point number or x value
	int rangeSpecified;					// 0 if user specified no range. 1 if user specified range.
	int isPoint;						// 0 if user specified range using X values. 1 if user specified range using points.
};
typedef struct WaveRange WaveRange;
typedef struct WaveRange *WaveRangePtr;

// This is what is in the runtime parameter structure for a mode 1 structure parameter.
// For a mode 0 structure parameter, the runtime parameter structure contains just a pointer to the structure.
struct IgorStructInfo {
	void* structPtr;					// Pointer to the structure.
	unsigned int structSize;			// Size of structure in bytes.
	char structTypeName[MAX_OBJ_NAME+1];
	int moduleSerialNumber;				// Used by Igor to determine the procedure file in which structure is defined.
	unsigned char reserved[32];			// Reserved for future use.
};
typedef struct IgorStructInfo IgorStructInfo;
typedef struct IgorStructInfo *IgorStructInfoPtr;

#pragma pack()	// Restore default structure packing
