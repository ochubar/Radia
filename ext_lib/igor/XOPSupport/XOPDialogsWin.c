/*	This file contains utilities for Windows XOPs that create dialogs.
	Most of these routines are also available on Macintosh.  Comments
	identify those that are Windows-specific.
*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h

/*	XOPEmergencyAlert(message)
	
	This routine used by the XOP Toolkit for dire emergencies only.
	You should not need it. Use XOPOKAlert instead.
	
	Thread Safety: XOPEmergencyAlert is not thread-safe.
*/
void
XOPEmergencyAlert(const char* message)
{
	MessageBox(NULL, message, "Emergency", MB_OK);
}

/*	XOPOKAlert(title, message)
	
	Thread Safety: XOPOKAlert is not thread-safe.
*/
void
XOPOKAlert(const char* title, const char* message)
{
	if (!CheckRunningInMainThread("XOPOKAlert"))
		return;

	MessageBox(IgorClientHWND(), (char*)message, (char*)title, MB_OK);
}

/*	XOPOKCancelAlert(title, message)

	Returns 1 for OK, -1 for cancel.
	
	Thread Safety: XOPOKCancelAlert is not thread-safe.
*/
int
XOPOKCancelAlert(const char* title, const char* message)
{
	int result;

	if (!CheckRunningInMainThread("XOPOKCancelAlert"))
		return -1;

	
	result = MessageBox(IgorClientHWND(), (char*)message, (char*)title, MB_OKCANCEL);
	if (result == IDOK)
		return 1;
	return -1;
}

/*	XOPYesNoAlert(title, message)

	Returns 1 for yes, 2 for no.
	
	Thread Safety: XOPYesNoAlert is not thread-safe.
*/
int
XOPYesNoAlert(const char* title, const char* message)
{
	int result;
	
	if (!CheckRunningInMainThread("XOPYesNoAlert"))
		return 2;

	result = MessageBox(IgorClientHWND(), (char*)message, (char*)title, MB_YESNO);
	if (result == IDYES)
		return 1;
	return 2;
}

/*	XOPYesNoCancelAlert(title, message)

	Returns 1 for yes, 2 for no, -1 for cancel.
	
	Thread Safety: XOPYesNoCancelAlert is not thread-safe.
*/
int
XOPYesNoCancelAlert(const char* title, const char* message)
{
	int result;
	
	if (!CheckRunningInMainThread("XOPYesNoCancelAlert"))
		return -1;

	result = MessageBox(IgorClientHWND(), (char*)message, (char*)title, MB_YESNOCANCEL);
	if (result == IDYES)
		return 1;
	if (result == IDNO)
		return 2;
	return -1;
}

/*	SetDialogBalloonHelpID(balloonHelpID)

	On the Macintosh, sets the resource ID for the 'hdlg' resource to be used
	for contextual help. If balloonHelpID is -1, this indicates that no balloon
	help is to be used.
	
	Currently does nothing on Windows.
	
	Thread Safety: SetDialogBalloonHelpID is not thread-safe.
*/
void
SetDialogBalloonHelpID(int balloonHelpID)
{
	if (!CheckRunningInMainThread("SetDialogBalloonHelpID"))
		return;
}

/*	GetDBox(theDialog, itemID, box)
	
	Thread Safety: GetDBox is not thread-safe.
*/
void
GetDBox(HWND theDialog, int itemID, Rect *boxPtr)
{
	HWND hwnd;
	RECT wr;
	
	hwnd = GetDlgItem(theDialog, itemID);
	if (hwnd == NULL) {
		SetRectEmpty(&wr);
	}
	else {
		GetWindowRect(hwnd, &wr);						// GetWindowRect returns the position in screen coordinates.
		ScreenToClient(theDialog, (POINT*)&wr.left);
		ScreenToClient(theDialog, (POINT*)&wr.right);
	}
	WinRectToMacRect(&wr, boxPtr);
}

/*	HiliteDControl(theDialog, itemNumber, enable)

	Thread Safety: HiliteDControl is not thread-safe.
*/
void
HiliteDControl(HWND theDialog, int itemID, int enable)
{
	if (enable)
		EnableDControl(theDialog, itemID);
	else
		DisableDControl(theDialog, itemID);
}

/*	DisableDControl(theDialog, itemID)

 	Disable checkbox or radio button.
	
	Thread Safety: DisableDControl is not thread-safe.
*/
void
DisableDControl(HWND theDialog, int itemID)
{
	HWND hwnd;
	
	hwnd = GetDlgItem(theDialog, itemID);
	if (hwnd == NULL)
		return;
	EnableWindow(hwnd, 0);
}

/*	EnableDControl(theDialog, itemID)

 	Enables of radio button.
	
	Thread Safety: EnableDControl is not thread-safe.
*/
void
EnableDControl(HWND theDialog, int itemID)
{
	HWND hwnd;
	
	hwnd = GetDlgItem(theDialog, itemID);
	if (hwnd == NULL)
		return;
	EnableWindow(hwnd, 1);
}

/*	SetRadBut(theDialog, firstID, lastID, itemID)

	Turns the radio button identified by itemID on and all others off.
	It is assumed that all items between firstID and lastID are radio buttons
	and that itemID is between firstID and lastID.
	
	Thread Safety: SetRadBut is not thread-safe.
*/
void
SetRadBut(HWND theDialog, int firstID, int lastID, int itemID)
{
	int i;
	
	for(i=firstID; i<=lastID; i++)
		SetCheckBox(theDialog, i, i==itemID);
}

/*	GetRadBut(theDialog, itemID)

 	Returns state of radio button.
	
	Thread Safety: GetRadBut is not thread-safe.
*/
int
GetRadBut(HWND theDialog, int itemID)
{
	return SendDlgItemMessage(theDialog, itemID, BM_GETCHECK, 0, 0) == BST_CHECKED; 
}

/*	GetDText(theDialog, itemID, text)

	Gets text from text item in dialog.
	Returns the number of characters in the text.
	You must allocate 256 bytes for the text parameter.
	
	Thread Safety: GetDText is not thread-safe.
*/
int
GetDText(HWND theDialog, int itemID, char text[256])
{
	GetDlgItemText(theDialog, itemID, text, 256);
	return (int)strlen(text);
}

/*	SetDText(theDialog, itemID, text)

	Sets text in text item in dialog.
	
	A side-effect of this routine is that the Windows OS will send a WM_COMMAND
	message to your dialog procedure, to pass it notifications such as EN_CHANGE.
	This causes problems because you can not distinguish an EN_CHANGE notification
	that arises from SetDText from an EN_CHANGE notification that arises from the
	user changing the value.
	
	To solve this problem, SetDText sets the gSetDTextInProgress global and
	IsWinDialogItemHitMessage tests it. IsWinDialogItemHitMessage returns true only
	if gSetDTextInProgress is 0.
	
	Thread Safety: SetDText is not thread-safe.
*/
int gSetDTextInProgress = 0;			// Tested by IsWinDialogItemHitMessage in XOPSupportWin.c.
void
SetDText(HWND theDialog, int itemID, const char* text)
{
	int temp;
	
	temp = gSetDTextInProgress;
	gSetDTextInProgress = 1;
	
	SetDlgItemText(theDialog, itemID, text);
	
	gSetDTextInProgress = temp;
}

/*	GetDInt(theDialog, itemID, theInt)

	Gets integer from text item in dialog.
	Returns zero if a number was read from text, non-zero if no number read.
	
	Thread Safety: GetDInt is not thread-safe.
*/
int
GetDInt(HWND theDialog, int itemID, int *theInt)
{
	char temp[256];

	*theInt = 0;		// In case sscanf finds no number.
	GetDText(theDialog, itemID, temp);
	return(sscanf(temp, "%d", theInt) != 1);
}

/*	SetDInt(theDialog, itemID, theInt)

	Sets text item in dialog to integer.
	
	Thread Safety: SetDInt is not thread-safe.
*/
void
SetDInt(HWND theDialog, int itemID, int theInt)
{
	char temp[32];
	
	sprintf(temp, "%d", theInt);
	SetDText(theDialog, itemID, temp);
}

/*	GetDInt64(theDialog, theItem, valPtr)

	Gets 64-bit integer from text item in dialog.
	Returns zero if a number was read from text, non-zero if no number read.
	
	Thread Safety: GetDInt64 is not thread-safe.
*/
int
GetDInt64(HWND theDialog, int theItem, SInt64* valPtr)
{
	char temp[256];
	int result;
	
	*valPtr = 0;		// In case sscanf finds no number.
	GetDText(theDialog, theItem, temp);
	result = sscanf(temp, "%lld", valPtr) != 1;
	return result;
}

/*	SetDInt64(theDialog, theItem, val)

	Sets text item in dialog to specified value.
	
	Thread Safety: SetDInt64 is not thread-safe.
*/
void
SetDInt64(HWND theDialog, int theItem, SInt64 val)
{
	char temp[32];
	
	sprintf(temp, "%lld", val);
	SetDText(theDialog, theItem, temp);
}

/*	GetDDouble(theDialog, itemID, theDouble)

	Gets 64 bit float from text item in dialog.
	Returns zero if a number was read from text, non-zero if no number read.
	
	Thread Safety: GetDDouble is not thread-safe.
*/
int
GetDDouble(HWND theDialog, int itemID, double *theDouble)
{
	char temp[256];
	double d;
	int result;

	*theDouble = 0;		// In case sscanf finds no number.
	GetDText(theDialog, itemID, temp);
	result = sscanf(temp, "%lf", &d) != 1;
	*theDouble = d;
	return result;
}

/*	SetDDouble(theDialog, itemID, theDouble)

	Sets text item in dialog to 64 bit float.
	
	Thread Safety: SetDDouble is not thread-safe.
*/
void
SetDDouble(HWND theDialog, int itemID, double *theDouble)
{
	char temp[32];
	
	sprintf(temp, "%g", *theDouble);			// HR, 980402: Changed from "%lg" to "%g".
	SetDText(theDialog, itemID, temp);
}

/*	ToggleCheckBox(theDialog, itemID)

	Toggles state of check box and returns new value.
	
	Thread Safety: ToggleCheckBox is not thread-safe.
*/
int
ToggleCheckBox(HWND theDialog, int itemID)
{
	int state;
	
	state = GetCheckBox(theDialog, itemID);
	state = !state;
	SetCheckBox(theDialog, itemID, state);
	return state;
}

/*	GetCheckBox(theDialog, itemID)

 	Returns state of check box.
	
	Thread Safety: GetCheckBox is not thread-safe.
*/
int
GetCheckBox(HWND theDialog, int itemID)
{
	return SendDlgItemMessage(theDialog, itemID, BM_GETCHECK, 0, 0) == BST_CHECKED; 
}

/*	SetCheckBox(theDialog, itemID, state)

	Sets state of check box.
	Returns the new value.
	
	Thread Safety: SetCheckBox is not thread-safe.
*/
int
SetCheckBox(HWND theDialog, int itemID, int state)
{
	SendDlgItemMessage(theDialog, itemID, BM_SETCHECK, state, 0); 
	return state;
}

/*	SelEditItem(theDialog, itemID)

	Selects the entire text of the edit item specified by itemID and
	sets the focus to that item.
	
	On Macintosh, if itemID is 0, it selects the entire text for the
	current editText item. However, on Windows, if itemID is 0, it
	does nothing.
	
	Thread Safety: SelEditItem is not thread-safe.
*/
void
SelEditItem(HWND theDialog, int itemID)
{
	char text[256];
	int start, end;

	if (itemID == 0)
		return;
	
	start = 0;
	end = GetDText(theDialog, itemID, text);
	SendDlgItemMessage(theDialog, itemID, EM_SETSEL, start, end); 
	SetFocus(GetDlgItem(theDialog, itemID));
}

/*	SelMacEditItem(theDialog, itemID)

	This is a NOP on Windows.

	On Macintosh, selects the entire text of the edit item specified by itemID and
	sets the focus to that item. If itemID is 0, it selects the entire text for the
	current editText item.

	In Macintosh dialogs, usually only edit items can have the focus. Thus, there
	are times when it is convenient for the program to set the focus to a particular
	edit item. On Windows, any item can have the focus and it is usually best to
	leave this in the control of the user. This routine achieves this by selecting
	a particular edit item on Macintosh while doing nothing on Windows.
	
	Thread Safety: SelMacEditItem is not thread-safe.
*/
void
SelMacEditItem(HWND theDialog, int itemID)
{
}

/*	DisplayDialogCmd(theDialog, itemID, cmd)

	Displays the command in an IGOR-style dialog. See GBLoadWaveDialog.c
	for an example.
	
	itemID is the item number of the dialog item in which the command
	is to be displayed. On the Macintosh, this must be a user item. On Windows,
	it must be an EDITTEXT item. See GBLoadWaveDialog.rc for the styles
	to use for this item.
	
	Thread Safety: DisplayDialogCmd is not thread-safe.
*/
void
DisplayDialogCmd(HWND theDialog, int itemID, const char* cmd)
{
	SetDText(theDialog, itemID, cmd);
}

/*	SetDialogPort(theDialog)
	
	SetDialogPort sets the current GrafPort on Macintosh and does nothing on
	Windows.
	
	On Macintosh, it returns the current GrafPort before SetDialogPort was
	called. On Windows, it returns theDialog. This routine exists solely to avoid
	the need for an ifdef when you need to deal with the Macintosh current GrafPort
	in cross-platform dialog code. 
	
	Thread Safety: SetDialogPort is not thread-safe.
*/
HWND
SetDialogPort(HWND theDialog)
{
	return theDialog;
}

/*	ShowDialogWindow(theDialog)

	Makes the dialog window visible.
	
	Thread Safety: ShowDialogWindow is not thread-safe.
*/
void
ShowDialogWindow(HWND theDialog)
{
	ShowWindow(theDialog, SW_SHOW);
}

/* PopMenu Routines */

/*	InitPopMenus(theDialog)

	This routine is currently a NOP on Windows.
	
	Thread Safety: InitPopMenus is not thread-safe.
*/
void
InitPopMenus(HWND theDialog)
{
}

/*	ItemIsPopMenu(theDialog, itemID)

	Returns the truth that the item is a popup menu (a drop-down combo box).
	
	Thread Safety: ItemIsPopMenu is not thread-safe.
*/
int
ItemIsPopMenu(HWND theDialog, int itemID)
{
	HWND hwnd;
	char className[32];
	DWORD style;
	
	hwnd = GetDlgItem(theDialog, itemID);
	if (hwnd == NULL)
		return 0;
	
	if (GetClassName(hwnd, className, sizeof(className)) == 0)
		return 0;					// Should not happen.
	
	if (CmpStr(className, "COMBOBOX") == 0) {
		style = GetWindowLong(hwnd, GWL_STYLE);
		if (style & CBS_DROPDOWNLIST)
			return 1;
		return 0;
	}
	
	return 0;
}

/*	CreatePopMenu(theDialog, popupItemID, titleItemID, itemList, initialItem)

	Initializes a popup menu, implemented as a combo box, in the specified dialog.
	
	popupItemID is the dialog item ID for the popup menu. This must be the ID
	of a drop-down list box item, created in your resource script file using
	the COMBOBOX keyword and the style (CBS_DROPDOWNLIST | CBS_HASSTRINGS).
	
	titleItemID is the dialog item ID for the static text title for the popup menu.
	Prior to XOP Toolkit 5, on Macintosh, this item was highlighted when the user
	clicked on the popup menu. As of XOP Toolkit 5, it is no longer used by must be
	present for backward compatibility.
	
	itemList is a semicolon-separated list of items to insert into the menu.
	For example, "Red;Green;Blue".
	
	initialItem is the 1-based number of the item in the popup menu that should
	be initially selected.
	
	Thread Safety: CreatePopMenu is not thread-safe.
*/
int
CreatePopMenu(HWND theDialog, int popupItemID, int titleItemID, const char* itemList, int initialItem)
{
	AddPopMenuItems(theDialog, popupItemID, itemList);
	SetPopItem(theDialog, popupItemID, initialItem);
	return 0;
}

/*	GetPopMenu(theDialog, itemID, selItem, selStr)

 	Returns currently selected item ID and string.
	Pass NULL for selStr if you don't care about it.
	
	NOTE: For consistency with the Macintosh version of this routine,
	selItem is 1-based, not 0-based as is usual for combo box controls.
	
	Thread Safety: GetPopMenu is not thread-safe.
*/
void
GetPopMenu(HWND theDialog, int itemID, int *selItem, char *selStr)
{
	HWND hwnd;
	LRESULT index;
	
	hwnd = GetDlgItem(theDialog, itemID);
	if (hwnd == NULL)
		return;								// Should not happen.
	
	index = SendMessage(hwnd, CB_GETCURSEL, 0, 0);
	if (index == CB_ERR)
		*selItem = 1;
	else
		*selItem = (int)(index + 1);
	
	if (selStr != NULL) {
		if (index == CB_ERR)
			*selStr = 0;
		else
			SendMessage(hwnd, CB_GETLBTEXT, index, (LPARAM)selStr);
	}
}

/*	SetPopMatch(theDialog, itemID, selStr)

 	Sets currently selected item to that which matches string (case insensitive).
	Returns item number or zero if there is no match.
	
	NOTE: For consistency with the Macintosh version of this routine,
	the returned item number is 1-based, not 0-based as is usual for
	combo box controls.
	
	Thread Safety: SetPopMatch is not thread-safe.
*/
int
SetPopMatch(HWND theDialog, int itemID, const char *selStr)
{
	HWND hwnd;
	LRESULT index;
	
	hwnd = GetDlgItem(theDialog, itemID);
	if (hwnd == NULL)
		return 0;							// Should not happen.
	
	index = SendMessage(hwnd, CB_FINDSTRINGEXACT, -1, (LPARAM)selStr);
	if (index == CB_ERR)
		return 0;

	SendMessage(hwnd, CB_SETCURSEL, index, 0);
	return (int)(index+1);
}

/*	SetPopItem(theDialog, itemID, itemNumber)

 	Makes itemNumber the currently selected item.

	NOTE: For consistency with the Macintosh version of this routine,
	itemNumber is 1-based, not 0-based as is usual for ombo box controls.
	
	Thread Safety: SetPopItem is not thread-safe.
*/
void
SetPopItem(HWND theDialog, int itemID, int itemNumber)
{
	SendDlgItemMessage(theDialog, itemID, CB_SETCURSEL, itemNumber-1, 0);
}

/*	KillPopMenus(theDialog)

 	This routine is a NOP on Windows.
	
	Thread Safety: KillPopMenus is not thread-safe.
*/
void
KillPopMenus(HWND theDialog)
{
}

/*	AddPopMenuItems(theDialog, itemID, itemList)

	Adds the contents of itemList to the combo box control in the dialog specified
	by theDialog whose identifier is itemID.
	
	itemList can be a single item ("Red") or a semicolon-separated list of items
	("Red;Green;Blue;"). The trailing semicolon is optional.
	
	Thread Safety: AddPopMenuItems is not thread-safe.
*/
void
AddPopMenuItems(HWND theDialog, int itemID, const char* itemList)
{
	FillPopMenu(theDialog, itemID, itemList, (int)strlen(itemList), 10000);
}

/*	FillPopMenu(theDialog, itemID, itemList, itemListLen, afterItem)

	Sets the contents of the combo box control in the dialog specified
	by theDialog whose identifier is itemID.
	
	itemList is a semicolon separated list of items to be put into the popup menu.
	itemListLen is the total number of characters in itemList.
	afterItem specifies where the items in itemList are to appear in the menu.
		afterItem = 0			new items appear at beginning of popup menu
		afterItem = 10000		new items appear at end of popup menu
		afterItem = item number	new items appear after specified existing item number.
		
	NOTE: For consistency with the Macintosh version of this routine, afterItem is 1-based,
	not 0-based as is usual for combo box controls.
	
	Thread Safety: FillPopMenu is not thread-safe.
*/
void
FillPopMenu(HWND theDialog, int itemID, const char* itemList, int itemListLen, int afterItem)
{
	HWND hwnd;
	const char *p1;
	const char *p2;
	char item[256];
	int len, itemLen;
	
	hwnd = GetDlgItem(theDialog, itemID);
	if (hwnd == NULL)
		return;								// Should not happen.

	p1 = itemList;
	while (itemListLen > 0) {
		if (p2 = strchr(p1, ';'))
			len = (int)(p2 - p1);
		else
			len = itemListLen;				// The last one.
		itemLen = len;
		if (itemLen > 255)
			itemLen = 255;
		strncpy(item, p1, itemLen);
		item[itemLen] = 0;
		if (afterItem >= 10000)
			SendMessage(hwnd, CB_ADDSTRING, 0, (LPARAM)item);
		else
			SendMessage(hwnd, CB_INSERTSTRING, afterItem, (LPARAM)item);
		p1 += len+1;
		itemListLen -= len+1;
		afterItem++;
	}
}

/*	FillWavePopMenu(theDialog, itemID, match, options, afterItem)

	Sets the contents of the combo box control in the dialog specified
	by theDialog whose identifier is itemID.
	
	Puts names of waves into the combo box.
	match and options are as for the Igor WaveList() function:
		match = "*" for all waves
		options = "" for all waves
		options = "WIN: Graph0" for waves in graph0 only.
	
	afterItem is as for FillPopMenu, described above. For consistency with the
	Macintosh version of this routine, afterItem is 1-based, not 0-based as is
	usual for combo box controls.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: FillWavePopMenu is not thread-safe.
*/
int
FillWavePopMenu(HWND theDialog, int itemID, const char *match, const char *options, int afterItem)
{
	Handle listHandle;
	int result;	
	
	listHandle = NewHandle(0L);
	result = WaveList(listHandle, match, ";", options);
	FillPopMenu(theDialog, itemID, *listHandle, (int)GetHandleSize(listHandle), afterItem);
	DisposeHandle(listHandle);
	return result;
}

/*	FillPathPopMenu(theDialog, itemID, match, options, afterItem)

	Sets the contents of the combo box control in the dialog specified
	by theDialog whose identifier is itemID.
	
	Puts names of Igor paths into the combo box.
	match and options are as for the Igor PathList() function:
		match = "*" for all paths
		options = "" for all paths
	
	afterItem is as for FillPopMenu, described above. For consistency with the
	Macintosh version of this routine, afterItem is 1-based, not 0-based as is
	usual for combo box controls.
	
	Thread Safety: FillPathPopMenu is not thread-safe.
*/
int
FillPathPopMenu(HWND theDialog, int itemID, const char *match, const char *options, int afterItem)
{
	Handle listHandle;
	int result;	
	
	listHandle = NewHandle(0L);
	result = PathList(listHandle, match, ";", options);
	FillPopMenu(theDialog, itemID, *listHandle, (int)GetHandleSize(listHandle), afterItem);
	DisposeHandle(listHandle);
	return result;
}

/*	FillWindowPopMenu(theDialog, itemID, char *match, char *options, int afterItem)

	Sets the contents of the combo box control in the dialog specified
	by theDialog whose identifier is itemID.
	
	Puts names of Igor windows into the combo box.
	match and options are as for the Igor WinList() function:
		match = "*" for all windows
		options = "" for all windows
		options = "WIN: 1" for all graphs		( bit 0 selects graphs)
		options = "WIN: 2" for all tables		( bit 1 selects graphs)
		options = "WIN: 4" for all layouts		( bit 2 selects graphs)
		options = "WIN: 3" for all graphs and tables
	
	afterItem is as for FillPopMenu, described above. For consistency with the
	Macintosh version of this routine, afterItem is 1-based, not 0-based as is
	usual for combo box controls.
	
	Thread Safety: FillWindowPopMenu is not thread-safe.
*/
int
FillWindowPopMenu(HWND theDialog, int itemID, const char *match, const char *options, int afterItem)
{
	Handle listHandle;
	int result;	
	
	listHandle = NewHandle(0L);
	result = WinList(listHandle, match, ";", options);
	FillPopMenu(theDialog, itemID, *listHandle, (int)GetHandleSize(listHandle), afterItem);
	DisposeHandle(listHandle);
	return result;
}

/*	DeletePopMenuItems(theDialog, itemID, afterItem)

	Deletes the contents of the existing dialog combo box.
	
	afterItem is 1-based. Pass 0 to delete all items.
	
	NOTE: For consistency with the Macintosh version of this routine, afterItem is 1-based,
	not 0-based as is usual for combo box controls. If afterItem is zero, all items are deleted.
	
	Thread Safety: DeletePopMenuItems is not thread-safe.
*/
void
DeletePopMenuItems(HWND theDialog, int itemID, int afterItem)
{
	LRESULT index, numItems, numItemsToDelete;
	HWND hwnd;
	
	hwnd = GetDlgItem(theDialog, itemID);
	if (hwnd == NULL)							// Should not happen.
		return;

	numItems = SendMessage(hwnd, CB_GETCOUNT, 0, 0);
	if (numItems == CB_ERR)
		return;									// Should not happen.
	numItemsToDelete = numItems - afterItem;
	
	index = numItems - 1;						// Delete from the end.
	while (numItemsToDelete > 0) {
		SendMessage(hwnd, CB_DELETESTRING, index, 0);
		index -= 1;
		numItemsToDelete -= 1;
	}	
}
