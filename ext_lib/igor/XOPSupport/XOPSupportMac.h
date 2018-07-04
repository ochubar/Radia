// This file contains equates and prototypes that are needed on Macintosh only.


#ifdef __cplusplus
extern "C" {						/* This allows C++ to call the XOPSupport routines */
#endif

// Misc utilities.
// IsMacOSX has been removed. Use #ifdef MACIGOR.
void debugstr(const char* text);
void paramtext(const char* param0, const char* param1, const char* param2, const char* param3);


/* Resource routines (in XOPSupportMac.c) */
int XOPRefNum(void);
Handle GetXOPResource(int resType, int resID);
Handle GetXOPNamedResource(int resType, const char* name);


/* Mac-specific dialog utilities (in XOPSupportMac.c) */
int XOPGetDialogItemAsControl(DialogPtr theDialog, int itemNumber, ControlHandle* controlHPtr);
DialogPtr GetXOPDialog(int dialogID);
void DisposeXOPDialog(DialogPtr theDialog);
void XOPDialog(ModalFilterUPP filterProc, short *itemPtr);
void DoXOPDialog(short* itemHitPtr);


/* Mac-specific window utilities (in XOPSupportMac.c) */
WindowPtr GetXOPWindow(int windowID, char* wStorage, WindowPtr behind);


// Mac-specific menu utilities (in XOPSupportMac.c)
short CountMItems(MenuHandle theMenu);
void CheckItem(MenuHandle theMenu,short itemNumber,int checked);
void DisableItem(MenuHandle menuH, short itemNumber);
void EnableItem(MenuHandle menuH, short itemNumber);
void getmenuitemtext(MenuHandle theMenu, short itemNumber, char* itemString);
void setmenuitemtext(MenuHandle menuH, short itemNumber, const char* itemString);
void insertmenuitem(MenuHandle menuH, const char* itemString, short afterItemOneBased);
void appendmenu(MenuHandle menuH, const char* itemString);

#ifdef __cplusplus
}
#endif
