/*	This file contains routines that are Macintosh-specific.
	This file is used only when compiling for Macintosh.
*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h

// IsMacOSX has been removed. Use #ifdef MACIGOR.

void
debugstr(const char* message)			// Sends debug message to low level debugger (e.g., Macsbug).
{
	char ctemp[256];
	unsigned char ptemp[256];
	int len;
	
	len = strlen(message);
	if (len >= sizeof(ctemp))
		len = sizeof(ctemp)-1;
	strncpy(ctemp, message, len);
	ctemp[len] = 0;
	CopyCStringToPascal(ctemp, ptemp);
	DebugStr(ptemp);
}

/*	paramtext(param0, param1, param2, param3)
	
	Thread Safety: paramtext is not thread-safe.
*/
void
paramtext(const char* param0, const char* param1, const char* param2, const char* param3)
{
	unsigned char p0[256];
	unsigned char p1[256];
	unsigned char p2[256];
	unsigned char p3[256];

	CopyCStringToPascal(param0, p0);
	CopyCStringToPascal(param1, p1);
	CopyCStringToPascal(param2, p2);
	CopyCStringToPascal(param3, p3);
	ParamText(p0, p1, p2, p3);
}

/*	Resource Routines -- for dealing with resources

	These routines are for accessing Macintosh resource forks.
	There are no Windows equivalent routines, so these routines can't be used in
	Windows or cross-platform XOPs.
*/

/*	XOPRefNum()

	Returns XOP's resource file reference number.
	
	Thread Safety: XOPRefNum is thread-safe but there is little that you can do with it that is thread-safe.
*/
int
XOPRefNum(void)
{
	return((*(*XOPRecHandle)->stuffHandle)->xopRefNum);
}

/*	GetXOPResource(resType, resID)

	Tries to get specified handle from XOP's resource fork.
	Does not search any other resource forks and does not change curResFile.
	
	Thread Safety: GetXOPResource is not thread-safe.
*/
Handle
GetXOPResource(int resType, int resID)
{
	Handle rHandle;
	int curResNum;
	
	curResNum = CurResFile();
	UseResFile(XOPRefNum());
	rHandle = Get1Resource(resType, resID);
	UseResFile(curResNum);
	return(rHandle);
}

/*	GetXOPNamedResource(resType, name)

	Tries to get specified handle from XOP's resource fork.
	Does not search any other resource forks and does not change curResFile.
	
	Thread Safety: GetXOPNamedResource is not thread-safe.
*/
Handle
GetXOPNamedResource(int resType, const char* name)
{
	Handle rHandle;
	unsigned char pName[256];
	int curResNum;
	
	curResNum = CurResFile();
	UseResFile(XOPRefNum());
	CopyCStringToPascal(name, pName);
	rHandle = Get1NamedResource(resType, pName);
	UseResFile(curResNum);
	return rHandle;
}


/*	Dialog Routines

	These are Macintosh dialog routines that have no Windows equivalent.
*/

static void
ShowHideContextualHelp(int code)		// -1 means toggle, 0 means hide, 1 means show contextual help window.
{
	CallBack1(SHOW_HIDE_CONTEXTUAL_HELP, (void*)code);
}

/*	XOPGetDialogItemAsControl(theDialog, itemNumber, controlHPtr)

	Returns via controlHPtr the control handle for the specified dialog item.
	
	This routine displays an error message in Igor's history area if the
	item in question does not have an associated control handle. This would be
	a programming bug that the XOP programmer must correct.
	
	Long ago, the Mac OS provided a set of Dialog Manager calls that you could call
	to manipulate dialog items. Then Apple introduced new capabilities implemented
	by a part of the OS called the Appearance Manager. The method for manipulating
	a dialog item (e.g., setting the text of an EditText field) depends on whether
	you use Appearance Manager features or not. All modern code uses Appearance
	Manager techniques.
	
	To avoid the need to support two ways of doing the same thing, the Carbon XOP
	Toolkit dialog support routines ASSUME that the dialog is Appearance-Manager savvy.
	You make your dialog Appearance savvy by:
	
	1. Including a dlgx resource with the same resource ID as your DLOG resource.
	2. Including the kDialogFlagsUseControlHierarchy bit in the dlgx resource.
	3. Using Appearance Manager programming methods for implementing your dialog.
	
	Including the kDialogFlagsUseControlHierarchy bit in the dlgx resource causes
	the Mac OS to create a control handle for each item in your dialog, if the item
	is not already defined as a CNTL item.
	
	If this routine prints an error message for a valid dialog item number,
	this indicates that you have not made your dialog resources and/or code
	Appearance-Manager savvy.
	
	Thread Safety: XOPGetDialogItemAsControl is not thread-safe.
*/
int
XOPGetDialogItemAsControl(DialogPtr theDialog, int itemNumber, ControlHandle* controlHPtr)
{
	int err;
	
	err = GetDialogItemAsControl(theDialog, itemNumber, controlHPtr);
	if (err==0 && *controlHPtr==NULL)			// Should never happen.
		err = -1;
	if (err) {
		char message[256];
		sprintf(message, "XOP Bug: XOPGetDialogItemAsControl got error %d for item number %d."CR_STR, err, itemNumber);
		XOPNotice(message);
	}
	return err;
}

/*	XOPDialogFilter(dialog, eventPtr, itemHitPtr)
	
	Handles the mapping of the enter key to the default button, the mapping
	of the escape key to the cancel button, and setting the cursor depending on
	the control under it. This requires that you previously called SetDialogDefaultItem,
	SetDialogCancelItem and SetDialogTracksCursor as shown in the sample XOPs.
	
	Also handles updating Igor windows if the dialog is moved. However, it does not
	update XOP windows.
	
	Thread Safety: XOPDialogFilter is not thread-safe.
*/
static pascal Boolean
XOPDialogFilter(DialogPtr theDialog, EventRecord *eventPtr, short *itemHitPtr)
{
	int callStdFilter;
	int result;

    result = 0;					// Means ModalDialog should handle the event.
    callStdFilter = 1;
    
    switch (eventPtr->what) {
    	case updateEvt:
    		if ((WindowPtr)eventPtr->message != GetDialogWindow(theDialog)) {
				// Mac OS X does the updating automatically, apparently using offscreen bitmaps.
    		}
    		break;
    	
    	case keyDown:
    		{
    			int keyCode = eventPtr->message & charCodeMask;		// Low order byte is character code.
    			switch(keyCode) {
    				case kHelpCharCode:
    					if (eventPtr->modifiers & optionKey) {
    						// option-help shows or hides Igor's contextual help window.
    						ShowHideContextualHelp(-1);
    						result = -1;							// We handled the event.
    						*itemHitPtr = 0;						// This is not a hit in any item.
    						callStdFilter = 0;
    					}
    					break;
    			}
    		}
    		break;
    }
    
    if (callStdFilter)
		result = StdFilterProc(theDialog, eventPtr, itemHitPtr);
	
	return result;
}

static ModalFilterUPP gXOPDialogFilterUPP = NULL;	// Created by GetXOPDialog, disposed by DisposeXOPDialog.

/*	GetXOPDialog(dialogID)

	This routine is implemented on Macintosh only.
	
	Thread Safety: GetXOPDialog is not thread-safe.
*/
DialogPtr
GetXOPDialog(int dialogID)
{
	DialogPtr theDialog;
	int saveResFile;
	
	saveResFile = CurResFile();
	UseResFile(XOPRefNum());
	theDialog = GetNewDialog(dialogID, NULL, (WindowPtr)-1);
	UseResFile(saveResFile);

	if (theDialog == NULL)
		return NULL;

	gXOPDialogFilterUPP = NewModalFilterUPP(XOPDialogFilter);
	
	return theDialog;
}

/*	DisposeXOPDialog(dialogID)

	This routine is implemented on Macintosh only.
	
	Thread Safety: DisposeXOPDialog is not thread-safe.
*/
void
DisposeXOPDialog(DialogPtr theDialog)
{
	DisposeDialog(theDialog);

	if (gXOPDialogFilterUPP != NULL) {
		DisposeModalFilterUPP(gXOPDialogFilterUPP);
		gXOPDialogFilterUPP = NULL;
	}
	
	SetDialogBalloonHelpID(-1);		// Tell Igor's contextual help that the dialog is finished.
}

/*	XOPDialog(filterProc, itemPtr)

	This routine merely calls the Mac toolbox ModalDialog routine. In the days
	of 68K Macintoshes, it did some additional work that is no longer needed.
	
	NOTE:	If you are compiling a PowerMac native XOP, the filterProc parameter
			must be the address of a routine descriptor or NULL, not the direct address
			of your filter routine.

	This routine is implemented on Macintosh only.
	
	Thread Safety: XOPDialog is not thread-safe.
*/
void
XOPDialog(ModalFilterUPP filterProc, short *itemPtr)
{
	ModalDialog(filterProc, itemPtr);
}

/*	DoXOPDialog(itemHitPtr)

	This routine is implemented on Macintosh only.
	
	Thread Safety: DoXOPDialog is not thread-safe.
*/
void
DoXOPDialog(short* itemHitPtr)
{
	XOPDialog(gXOPDialogFilterUPP, itemHitPtr);
}

/*	Window Routines

	These are Macintosh window routines that have no Windows equivalent.
*/

/*	GetXOPWindow(windowID, wStorage, behind)

	This routine is for Macintosh only.

	Creates a window using the WIND resource with resource ID=windowID.
	This resource must be in the XOPs resource fork.
	wStorage and behind work as for the toolbox trap GetNewWindow.
	
	Returns a window pointer if everything OK or NULL if couldn't get window.
	Sets the windowKind field of the window equal to the XOPs refNum.
	YOU MUST NOT CHANGE THE windowKind FIELD.
	
	HR, 980317: Changed GetNewWindow to GetNewCWindow.
	
	Thread Safety: GetXOPWindow is not thread-safe.
*/
WindowPtr
GetXOPWindow(int windowID, char* wStorage, WindowPtr behind)
{
	int curResNum;
	Handle wHandle;
	WindowPtr wPtr = NULL;
	
	curResNum = CurResFile();
	UseResFile(XOPRefNum());
	wHandle = Get1Resource('WIND', windowID);	/* make sure resource exists in XOP file */
	if (wHandle) {
		wPtr = GetNewCWindow(windowID, wStorage, behind);
		if (wPtr)
			SetWindowKind(wPtr, XOPRefNum());
		ReleaseResource(wHandle);
	}
	UseResFile(curResNum);
	return(wPtr);
}


/*	Menu Routines

	These are Classic Mac routines that are not supported by the Carbon SDK.
	However, we need to support them in the XOP Toolkit because they are part
	of a cross-platform (Mac and Win) set of routines and we don't want to
	force Windows programmers to change their code just because Apple decided
	not to support some old routines.
*/

/*	CountMItems(menuH)

	Thread Safety: CountMItems is not thread-safe.
*/
short
CountMItems(MenuHandle menuH)
{
	return CountMenuItems(menuH);
}

/*	CheckItem(menuH, itemNumber, checked)

	Thread Safety: CheckItem is not thread-safe.
*/
void
CheckItem(MenuHandle menuH, short itemNumber, int checked)
{
	CheckMenuItem(menuH, itemNumber, checked);
}

/*	DisableItem(menuH, itemNumber)

	Thread Safety: DisableItem is not thread-safe.
*/
void
DisableItem(MenuHandle menuH, short itemNumber)
{
	DisableMenuItem(menuH, itemNumber);
}

/*	EnableItem(menuH, itemNumber)

	Thread Safety: EnableItem is not thread-safe.
*/
void
EnableItem(MenuHandle menuH, short itemNumber)
{
	EnableMenuItem(menuH, itemNumber);
}

/*	getmenuitemtext(theMenu, itemNumber, itemString)

	Thread Safety: getmenuitemtext is not thread-safe.
*/
void
getmenuitemtext(MenuHandle theMenu, short itemNumber, char* itemString)
{
	unsigned char pItemString[256];
	
	GetMenuItemText(theMenu, itemNumber, pItemString);
	CopyPascalStringToC(pItemString, itemString);
}

/*	setmenuitemtext(menuH, itemNumber, itemString)

	Thread Safety: setmenuitemtext is not thread-safe.
*/
void
setmenuitemtext(MenuHandle menuH, short itemNumber, const char* itemString)
{
	unsigned char pItemString[256];
	
	CopyCStringToPascal(itemString, pItemString);
	SetMenuItemText(menuH, itemNumber, pItemString);
}

/*	insertmenuitem(menuH, itemString, afterItemOneBased)

	Thread Safety: insertmenuitem is not thread-safe.
*/
void
insertmenuitem(MenuHandle menuH, const char* itemString, short afterItemOneBased)
{
	unsigned char pItemString[256];
	
	CopyCStringToPascal(itemString, pItemString);
	InsertMenuItem(menuH, pItemString, afterItemOneBased);
}

/*	appendmenu(menuH, itemString)

	Thread Safety: appendmenu is not thread-safe.
*/
void
appendmenu(MenuHandle menuH, const char* itemString)
{
	unsigned char pItemString[256];
	
	CopyCStringToPascal(itemString, pItemString);
	AppendMenuItemText(menuH, pItemString);
}
