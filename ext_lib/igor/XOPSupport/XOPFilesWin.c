/*	Contains platform-specific file-related routines.
	Platform-independent file-related routines are in XOPFiles.c
*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h

/*	XOPCreateFile(fullFilePath, overwrite, macCreator, macFileType)

	Creates a file with the location and name specified by fullFilePath.
	
	fullFilePath must be an HFS path (using colon separators) on Macintosh and a Windows path
	(using backslashes) on Windows.
	
	On Macintosh, the elements of fullFilePath may exceed the normal HFS 31 character limit.

	If overwrite is true and a file by that name already exists, it first
	deletes the conflicting file. If overwrite is false and a file by that
	name exists, it returns an error.
	
	macFileType is ignored on Windows. On Macintosh, it is used to set
	the new file's type. For example, use 'TEXT' for a text file.
	
	macCreator is ignored on Windows. On Macintosh, it is used to set
	the new file's creator code. For example, use 'IGR0' (last character is zero)
	for an file.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPCreateFile is thread-safe with Igor Pro 6.20 or later.
*/
int
XOPCreateFile(const char* fullFilePath, int overwrite, int macCreator, int macFileType)
{
	HANDLE fileH;
	DWORD accessMode, shareMode;
	int err;
	
	if (FullPathPointsToFile(fullFilePath)) {
		if (overwrite) {
			if (err = XOPDeleteFile(fullFilePath))
				return err;
		}
		else {
			return FILE_CREATE_ERROR;
		}
	}
	
	err = 0;
	accessMode = GENERIC_READ | GENERIC_WRITE;
	shareMode = 0;
	fileH = CreateFile(fullFilePath, accessMode, shareMode, NULL, CREATE_NEW, FILE_ATTRIBUTE_NORMAL, NULL);
	if (fileH == INVALID_HANDLE_VALUE)
		err = WMGetLastError();
	else
		CloseHandle(fileH);
	return err;
}

/*	XOPDeleteFile(fullFilePath)

	Deletes the file specified by fullFilePath.
	
	fullFilePath must be an HFS path (using colon separators) on Macintosh and a Windows path
	(using backslashes) on Windows.
	
	On Macintosh, the elements of fullFilePath may exceed the normal HFS 31 character limit.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPDeleteFile is thread-safe with Igor Pro 6.20 or later.
*/
int
XOPDeleteFile(const char* fullFilePath)
{
	int err;

	err = 0;
	if (DeleteFile(fullFilePath) == 0)
		err = WMGetLastError();
	return err;
}

/*	XOPOpenFile(fullFilePath, readOrWrite, fileRefPtr)

	If readOrWrite is zero, opens an existing file for reading and returns a file reference
	via fileRefPtr.

	If readOrWrite is non-zero, opens an existing file for writing or creates a new
	file if none exists and returns a file reference via fileRefPtr.

	fullFilePath must be an HFS path (using colon separators) on Macintosh and a Windows path
	(using backslashes) on Windows.
	
	On Macintosh, the elements of fullFilePath may exceed the normal HFS 31 character limit.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPOpenFile is thread-safe.
*/
int
XOPOpenFile(const char* fullFilePath, int readOrWrite, XOP_FILE_REF* fileRefPtr)
{
	char path[MAX_PATH_LEN+1];
	
	if (strlen(fullFilePath) > MAX_PATH_LEN)
		return PATH_TOO_LONG;
	strcpy(path, fullFilePath);

	*fileRefPtr = fopen(path, readOrWrite ? "wb" : "rb");
	if (*fileRefPtr == NULL)
		return FILE_OPEN_ERROR;
	return 0;
}

/*	FullPathPointsToFile(fullPath)

	Returns 1 if the path points to an existing file, 0 if it points to a folder or
	does not point to anything.
	
	fullPath may be a Macintosh or a Windows path.

	This routine is typically used by a file-loader XOP when it decides if it has
	enough information to load the file or needs to display an open file dialog.	
	
	Thread Safety: FullPathPointsToFile is thread-safe with Igor Pro 6.20 or later.
*/
int
FullPathPointsToFile(const char* fullPath)
{
	char nativePath[MAX_PATH_LEN+1];
	DWORD attributes;
	int err;
	
	if (err = GetNativePath(fullPath, nativePath))
		return err;

	attributes = GetFileAttributes(nativePath);
	if (attributes == 0xFFFFFFFF)					// Error?
		return 0;

	if ((attributes & FILE_ATTRIBUTE_DIRECTORY) != 0)
		return 0;									// Points to a folder.
	
	return 1;
}

/*	FullPathPointsToFolder(fullPath)

	Returns true if the path points to an existing folder, false if it points to a file or
	does not point to anything.
	
	fullPath may be a Macintosh or a Windows path.
	
	Thread Safety: FullPathPointsToFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
FullPathPointsToFolder(const char* fullPath)
{
	char nativePath[MAX_PATH_LEN+1];
	DWORD attributes;
	int err;
	
	if (err = GetNativePath(fullPath, nativePath))
		return err;

	attributes = GetFileAttributes(nativePath);
	if (attributes == 0xFFFFFFFF)					// Error?
		return 0;

	if ((attributes & FILE_ATTRIBUTE_DIRECTORY) == 0)
		return 0;									// Points to a file.
	
	return 1;
}

/*	GetNativePath(filePathIn, filePathOut)
	
	Call this to make sure that a file path uses the conventions regarding
	colons and backslashes of the current platform.
	
	It copies filePathIn to filePathOut. If filePathIn does not use the conventions
	of the current platform, it converts filePathOut to use those conventions.
	
	filePathOut can point to the same memory as filePathIn or it can
	point to different memory.
	
	filePathOut must be declared to hold MAX_PATH_LEN+1 characters.
	
	Function result is 0 if OK or an error code (e.g., PATH_TOO_LONG).
	
	Thread Safety: GetNativePath is thread-safe with Igor Pro 6.20 or later.
*/
int
GetNativePath(const char* filePathIn, char filePathOut[MAX_PATH_LEN+1])
{
	int err;
	
	if (strlen(filePathIn) > MAX_PATH_LEN)
		return PATH_TOO_LONG;
		
	if (filePathOut != filePathIn)
		strcpy(filePathOut, filePathIn);
	
	err = MacToWinPath(filePathOut);

	return err;
}

static UINT_PTR APIENTRY	// Hook for open or save file dialogs.
OpenOrSaveFileNameHook(HWND hdlg, UINT uiMsg, WPARAM wParam, LPARAM lParam)
{	
	HWND hMainDlg;
	
	/*	Because we use the OFN_EXPLORER flag and we specify a hook function, Windows
		creates a child dialog for us and the hdlg parameter to this hook is the child
		dialog.
	*/
	hMainDlg = GetParent(hdlg);
	if (hMainDlg == NULL)
		return 0;

	switch(uiMsg) {
		case WM_INITDIALOG:
			/*	HR, 090121, XOPSupport 5.09: Because we now use OFN_ENABLESIZING, the OS positions
				sizes the dialog. However, without this hook, the dialog initially comes up in the top/left
				corner of the frame window. Therefore I decided to leave the hook in. After this hook
				positions the window, the OS repositions and resizes it which may cause a brief flash.
			*/
			PositionWinDialogWindow(hMainDlg, NULL);
			break;
	}
	return 0;			// Let default dialog box procedure process the message.
}

/*	XOPOpenFileDialog(prompt, fileFilterStr, fileIndexPtr, initialDir, filePath)

	Displays the open file dialog.
	
	Returns 0 if the user chooses a file or -1 if the user cancels or another
	non-zero number in the event of an error. Returns the full path to the
	file via filePath. In the event of a cancel, filePath is unmodified.
	filePath is a native path (using colons on Macintosh, backslashes on Windows).
	
	prompt sets the dialog window title.
	
	If fileFilterStr is "", then the open file dialog displays all types
	of files, both on Macintosh and Windows. If fileFilterStr is not "",
	it identifies the type of files to display.

	fileFilterStr provides control over the Enable popup menu which the Macintosh Navigation
	Manager displays in the Open File dialog.  For example, the string
		"Text Files:TEXT,IGTX:.txt,.itx;All Files:****:;"
	results in two items in the Enable popup menu. The first says "Text Files"
	and displays any file whose Macintosh file type is TEXT or IGTX as well as any
	file whose file name extension is ".txt" or ".itx". The second item says "All Files"
	and displays all files.
	
	For further details on the fileFilterStr on Macintosh, see the comments in XOPNavOpenFileDialog.
	
	On Windows, fileFilterStr is constructed as for the lpstrFilter field of
	the OPENFILENAME structure for the Windows GetOpenFileName routine. For
	example, to allow the user to select text files and Igor Text files, use
	"Text Files (*.txt)\0*.txt\0Igor Text Files (*.itx)\0*.itx\0All Files (*.*)\0*.*\0\0".
	Note that the string ends with two null characters (\0\0).
	
	fileIndexPtr is ignored if it is NULL. If it is not NULL, then
	*fileIndexPtr is the one-based index of the file type filter to be initially
	selected. In the example given above, setting *fileIndexPtr to 2 would select
	the Igor Text file filter on entry to the dialog. On exit from the dialog,
	*fileIndexPtr is set to the index of the file filter string that the user last
	selected.  
	
	initialDir can be "" or it can point to a full path to a directory. It
	determines the directory that will be initially displayed in the open file
	dialog. If "", the directory will be the last directory that was seen
	in the open or save file dialogs. If initialDir points to a valid path to a directory,
	then this directory will be initially displayed in the dialog. On Macintosh,
	initialDir is a Macintosh HFS path. On Windows, it is a Windows path.
	
	Returns via filePath the full path to the file that the user chose. filePath
	is unchanged if the user cancels. filePath is a Macintosh HFS path on Macintosh
	and a Windows path on Windows. filePath must point to a buffer of
	at least MAX_PATH_LEN+1 bytes.
	
	On Windows, the initial value of filePath sets the initial contents of
	the File Name edit control in the open file dialog. The following values
	are valid:
		""									If there is no initial file name
		a file name
		a full Mac or Win path to a file	Allowed as of XOP Toolkit 5.04
	
	On Macintosh, the initial value of filePath is not currently used. It should be set
	the same as for Windows because it may be used in the future.
	
	In the event of an error other than a cancel, XOPOpenFileDialog displays
	an error dialog. This should never or rarely happen.
	
	WINDOWS NOTES
	
	The dialog will appear in the upper left corner of the screen. This is
	because Windows provides no straight-forward way to set the position of
	the dialog.
	
	Thread Safety: XOPOpenFileDialog is not thread-safe.
*/
int
XOPOpenFileDialog(
	const char* prompt,
	const char* fileFilterStr, int* fileIndexPtr,
	const char* initialDir,
	char filePath[MAX_PATH_LEN+1])
{
	OPENFILENAME ofn;
	char filePath2[MAX_PATH_LEN+1];
	char initialDir2[MAX_PATH_LEN+1];

	if (!CheckRunningInMainThread("XOPOpenFileDialog"))
		return NOT_IN_THREADSAFE;
	
	if (*fileFilterStr == 0)
		fileFilterStr = "All Files (*.*)\0*.*\0\0";
		
	if (*initialDir == 0) {
		GetStandardFileWinPath(initialDir2);	// Get Igor's open file dialog directory.
	}
	else {
		strcpy(initialDir2, initialDir);
		SetStandardFileWinPath(initialDir);		// Sets initial directory for next open file dialog. This will be overridden below, but not if the user cancels.
	}
		
	/*	HR, 040928, XOP Toolkit 5.04
		Previously this copied filePath to filePath2. This was correct because the filePath parameter
		was supposed to be either "" or just the proposed file name. However, I incorrectly passed
		a full path for the filePath parameter in all of the sample XOPs. This mistake undoubtedly
		leaked into users' XOPs. Therefore I now allow filePath to be either "", just a file name,
		or a full path.
	*/
	// strcpy(filePath2, filePath);				// HR, 010815: Previously filePath2 was set to "" which prevented the File Name item in the Windows Open File dialog from being preset as the comment above says it should be.
	GetLeafName(filePath, filePath2);

	MemClear(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = IgorClientHWND();
	ofn.lpstrFile = filePath2;
	ofn.nMaxFile = MAX_PATH_LEN+1;
	ofn.lpstrFilter = fileFilterStr;
	ofn.nFilterIndex = fileIndexPtr==NULL ? 1 : *fileIndexPtr;
	ofn.lpstrTitle = prompt;
	ofn.lpstrFileTitle = NULL;
	ofn.lpstrInitialDir = initialDir2;
	ofn.lpfnHook = OpenOrSaveFileNameHook;		// Needed to set position of the dialog. Otherwise, it is in top/left corner of screen.
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
	ofn.Flags |= OFN_EXPLORER;
	ofn.Flags |= OFN_ENABLEHOOK;				// Needed so that hook will be called.
	ofn.Flags |= OFN_ENABLESIZING;				// HR, 090121: Added this to get resizeable dialog.
	ofn.Flags |= OFN_HIDEREADONLY;
	ofn.Flags |= OFN_NOCHANGEDIR;				// Changing the current directory causes problems. e.g., if set to floppy disk and the floppy is removed, the system slows down.

	if (GetOpenFileName(&ofn) == 0) {
		int err;
		err = CommDlgExtendedError();			// err will be zero if cancel.
		if (err == 0)
			return -1;

		// We got an error other than cancel.
		*filePath2 = 0;							// HR, 021114: Clear possible bad fields
		*initialDir2 = 0;						// and try again.
		if (GetOpenFileName(&ofn) != 0) {		// Succeeded this time?
			err = 0;
		}
		else {
			if (CommDlgExtendedError() == 0)
				return -1;						// User canceled.
			
			// Report the original error.
			err = WindowsErrorToIgorError(err);
			IgorError("XOPOpenFileDialog", err);
			return err;
		}
	}
	
	if (fileIndexPtr != NULL)
		*fileIndexPtr = ofn.nFilterIndex;
	
	strcpy(filePath, filePath2);
	SetStandardFileWinPath(filePath);			// Update Igor's open file dialog directory.

	return 0;
}

/*	XOPSaveFileDialog(prompt, fileFilterStr, fileIndexPtr, initialDir, defaultExtensionStr, filePath)

	Displays the save file dialog.
	
	Returns 0 if the user provides a file name or -1 if the user cancels or another
	non-zero number in the event of an error.
	
	Returns the full path to the file via filePath. filePath is both an input and an
	output as explained below. In the event of a cancel, filePath is unmodified.
	filePath is a Macintosh HFS path on Macintosh and a Windows path on Windows.
	
	On Windows, prompt sets the dialog caption. On Macintosh, it sets a prompt
	string in the dialog.
	
	fileFilterStr is now used to control the contents of the Format popup menu
	in the Save File dialog.
	
	On Macintosh, if there is only one format in which you can save the file,
	pass "" for fileFilterStr. This will cause the Format menu to be hidden.
	If you can save the file in more than one format, pass a string like this:
		"Plain Text:TEXT:.txt;Igor Text:IGTX:.itx;"
		
	This would give you a Format menu like this:
		Plain Text
		Igor Text
	
	fileFilterStr on Macintosh

		fileFilterStr consists of sections terminated by a semicolon. For example,
		here is one section:
			"Data Files:TEXT:.dat;"
			
		Each section consists of three components: a menu item string (e.g., Data Files)
		to be displayed in the Format popup menu, a Macintosh file type (e.g., TEXT),
		and an extension (e.g., .dat).

		At present, only the menu item string and extension are used.

		The Macintosh file type is currently not used. If there is no meaningful Macintosh
		file type, leave the file type component empty.

		If there is no meaningful extension, leave the extension component empty.

	fileFilterStr on Windows
	
		On Windows, fileFilterStr identifies the types of files to display and the types
		of files that can be created. It is constructed as for the lpstrFilter
		field of the OPENFILENAME structure for the Windows GetSaveFileName routine.
		For example, to allow the user to save as a text file or as an Igor Text file,
		use "Text Files (*.txt)\0*.txt\0Igor Text Files (*.itx)\0*.itx\0\0". Note that
		the string ends with two null characters (\0\0). If fileFilterStr is "", this
		behaves the same as "Text Files (*.txt)\0*.txt\0\0". 

	fileIndexPtr it is ignored if it is NULL. If it is not NULL, then *fileIndexPtr
	is the one-based index of the file type filter to be initially selected.
	In the example given above, setting *fileIndexPtr to 2 would select the Igor
	Text file type on entry to the dialog. On exit from the dialog, *fileIndexPtr
	is set to the index of the file type string that the user last selected.
	
	initialDir can be "" or it can point to a full path to a directory. It
	determines the directory that will be initially displayed in the save file
	dialog. If "", the directory will be the last directory that was seen in the
	open or save file dialogs. If initialDir points to a valid path to a directory,
	then this directory will be initially displayed in the dialog. On Macintosh,
	initialDir is a Macintosh HFS path. On Windows, it is a Windows path. 
	
	defaultExtensionStr points to the extension to be added to the
	file name if the user does not enter an extension. For example, pass "txt"
	to have ".txt" appended if the user does not enter an extension. If you don't
	want any extension to be added in this case, pass NULL.
	
	Prior to XOP Toolkit 6.00, defaultExtensionStr was ignored on Macintosh.
	
	Returns via filePath the full path to the file that the user chose
	or "" if the user cancelled. The path is a Macintosh HFS path on Macintosh
	and a Windows path on Windows. filePath must point to a buffer of
	at least MAX_PATH_LEN+1 bytes.
	
	On Windows and Macintosh, the initial value of filePath sets the initial contents of
	the File Name edit control in the save file dialog. The following values
	are valid:
		""									If there is no initial file name
		a file name
		a full Mac or Win path to a file
	
	In the event of an error other than a cancel, XOPSaveFileDialog displays
	an error dialog. This should never or rarely happen.
	
	WINDOWS NOTES
	
	The dialog will appear in the upper left corner of the screen. This is
	because Windows provides no straight-forward way to set the position of
	the dialog.
	
	Thread Safety: XOPSaveFileDialog is not thread-safe.
*/
int
XOPSaveFileDialog(
	const char* prompt,
	const char* fileFilterStr, int* fileIndexPtr,
	const char* initialDir,
	const char* defaultExtensionStr,
	char filePath[MAX_PATH_LEN+1])
{
	OPENFILENAME ofn;
	char filePath2[MAX_PATH_LEN+1];
	char initialDir2[MAX_PATH_LEN+1];
	
	if (!CheckRunningInMainThread("XOPSaveFileDialog"))
		return NOT_IN_THREADSAFE;
	
	if (*fileFilterStr == 0)
		fileFilterStr = "Text Files (*.txt)\0*.txt\0\0";
		
	if (*initialDir == 0) {
		GetStandardFileWinPath(initialDir2);	// Get Igor's save file dialog directory.
	}
	else {
		strcpy(initialDir2, initialDir);
		SetStandardFileWinPath(initialDir);		// Sets initial directory for next save file dialog. This will be overridden below, but not if the user cancels.
	}
		
	/*	HR, 040928, XOP Toolkit 5.04
		Previously this copied filePath to filePath2. This was correct because the filePath parameter
		was supposed to be either "" or just the proposed file name. However, I incorrectly passed
		a full path for the filePath parameter in all of the sample XOPs. This mistake undoubtedly
		leaked into users' XOPs. Therefore I now allow filePath to be either "", just a file name,
		or a full path.
	*/
	// strcpy(filePath2, filePath);				// HR, 010815: Previously filePath2 was set to "" which prevented the File Name item in the Windows Open File dialog from being preset as the comment above says it should be.
	GetLeafName(filePath, filePath2);

	MemClear(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = IgorClientHWND();
	ofn.lpstrFile = filePath2;
	ofn.nMaxFile = MAX_PATH_LEN+1;
	ofn.lpstrFilter = fileFilterStr;
	ofn.nFilterIndex = fileIndexPtr==NULL ? 1 : *fileIndexPtr;
	ofn.lpstrDefExt = defaultExtensionStr;
	ofn.lpstrTitle = prompt;
	ofn.lpstrFileTitle = NULL;
	ofn.lpstrInitialDir = initialDir2;
	ofn.lpfnHook = OpenOrSaveFileNameHook;		// Needed to set position of the dialog. Otherwise, it is in top/left corner of screen.
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_OVERWRITEPROMPT;
	ofn.Flags |= OFN_EXPLORER;
	ofn.Flags |= OFN_ENABLEHOOK;				// Needed so that hook will be called.
	ofn.Flags |= OFN_ENABLESIZING;				// HR, 090121: Added this to get resizeable dialog.
	ofn.Flags |= OFN_HIDEREADONLY;
	ofn.Flags |= OFN_NOCHANGEDIR;				// Changing the current directory causes problems. e.g., if set to floppy disk and the floppy is removed, the system slows down.

	if (GetSaveFileName(&ofn) == 0) {
		int err;
		err = CommDlgExtendedError();			// err will be zero if cancel.
		if (err == 0)
			return -1;

		// We got an error other than cancel.
		*filePath2 = 0;							// HR, 021114: Clear possible bad fields
		*initialDir2 = 0;						// and try again.
		if (GetSaveFileName(&ofn) != 0) {		// Succeeded this time?
			err = 0;
		}
		else {
			if (CommDlgExtendedError() == 0)
				return -1;						// User canceled.
			
			// Report the original error.
			err = WindowsErrorToIgorError(err);
			IgorError("XOPSaveFileDialog", err);
			return err;
		}
	}
	
	if (fileIndexPtr != NULL)
		*fileIndexPtr = ofn.nFilterIndex;
	
	strcpy(filePath, filePath2);
	SetStandardFileWinPath(filePath);			// Update Igor's open file dialog directory.

	return 0;
}
