/*	This file provides support for a small number of Macintosh routines that
	are needed by Windows XOPs.
*/

#ifdef WINIGOR				// Compiling for Windows [
	// These are for Macintosh emulation on Windows so that we can write portable XOPs.

	#pragma pack(2)		// All structures passed between Igor and XOP are two-byte aligned.

	#ifdef __cplusplus
	extern "C" {
	#endif

	typedef UInt32 OSType;
	typedef char* Ptr;
	typedef Ptr* Handle;

	typedef struct Point {
		short v;
		short h;
	} Point;

	typedef struct Rect {
		short top;
		short left;
		short bottom;
		short right;
	} Rect;

	typedef struct RGBColor {
		unsigned short red;
		unsigned short green;
		unsigned short blue;
	} RGBColor;

	struct EventRecord {
		SInt32 what;
		SInt32 message;
		UInt32 when;
		Point where;
		SInt32 modifiers;
	};
	typedef struct EventRecord EventRecord;

	// Menu-related emulation.
	typedef Handle MenuHandle;
	
	// These menu-related routines are used by XOPSupport but should not be needed by an XOP.
	HOST_EXPORT void WMDrawMenuBar(void);
	HOST_EXPORT void WMDeleteMenu(short menuID);
	HOST_EXPORT void WMInsertMenu(MenuHandle menuH, short beforeID);
	HOST_EXPORT MenuHandle GetMenuHandle(short menuID);
	HOST_EXPORT MenuHandle WMGetMenuFromModule(HMODULE hModule, short resourceID);
	HOST_EXPORT int GetMenuID(MenuHandle theMenu);
	HOST_EXPORT void SetMenuID(MenuHandle mH, int newMenuID);
	
	// These menu-related routines may be needed by an XOP.
	HOST_EXPORT short CountMItems(MenuHandle theMenu);
	HOST_EXPORT void CheckItem(MenuHandle theMenu,short item,int checked);
	HOST_EXPORT void DisableItem(MenuHandle menuH, short item);
	HOST_EXPORT void EnableItem(MenuHandle menuH, short item);
	HOST_EXPORT void getmenuitemtext(MenuHandle theMenu, short item, char* itemString);
	HOST_EXPORT void setmenuitemtext(MenuHandle menuH, short item, char* itemString);
	HOST_EXPORT void DeleteMenuItem(MenuHandle theMenu, short item);
	HOST_EXPORT void insertmenuitem(MenuHandle menuH, char *itemString, short afterItemOneBased);
	HOST_EXPORT void appendmenu(MenuHandle menuH, char *itemString);

	// Memory-related emulation.
	HOST_EXPORT int MemError(void);
	HOST_EXPORT Ptr NewPtr(BCInt size);
	HOST_EXPORT BCInt GetPtrSize(Ptr p);
	HOST_EXPORT void SetPtrSize(Ptr p, BCInt size);
	HOST_EXPORT int PtrToHand(Ptr srcPtr, Handle* destHandlePtr, BCInt size);
	HOST_EXPORT int PtrAndHand(void* p, Handle h, BCInt size);
	HOST_EXPORT void DisposePtr(Ptr p);
	HOST_EXPORT Handle NewHandle(BCInt size);
	HOST_EXPORT int HandToHand(Handle* destHandlePtr);
	HOST_EXPORT int HandAndHand(Handle hand1, Handle hand2);
	HOST_EXPORT BCInt GetHandleSize(Handle h);
	HOST_EXPORT void SetHandleSize(Handle h, BCInt newSize);
	HOST_EXPORT void DisposeHandle(Handle h);

	// These are obsolete and were removed in XOP Toolkit 6.
	// HOST_EXPORT int HGetState(Handle h);
	// HOST_EXPORT void HSetState(Handle h,int state);
	// HOST_EXPORT void HLock(Handle h);
	// HOST_EXPORT void HUnlock(Handle h);
	// HOST_EXPORT void MoveHHi(Handle h);

	// Miscellaneous emulation.
	HOST_EXPORT TickCountInt TickCount(void);
	
	// Other support from Igor.lib.
	HOST_EXPORT int WindowsErrorToIgorError(DWORD winErrCode);
	HOST_EXPORT DWORD WMGetLastError(void);
	HOST_EXPORT int SetScrapData(void* theWindow, int length, UInt32 theType, void *source);
	int GetStandardFileWinPath(char* fullPath);
	int SetStandardFileWinPath(const char* fullPath);
	void GetWinIndString(HMODULE hModule, char* theString, int strListID, int index);
	HMODULE GetXOPModule(struct IORec** ioRecHandle);
	HMODULE GetIgorModule(void);
	HWND GetIgorClientHWND(void);
	int HandleXOPWinMessage(HMODULE hXOPModule, HWND hwnd, UINT iMsg, WPARAM wParam, LPARAM lParam, int beforeOrAfter);

	#ifdef __cplusplus
	}
	#endif

	#pragma pack()		// Reset structure alignment to default.
	
#endif						// Compiling for Windows ]
