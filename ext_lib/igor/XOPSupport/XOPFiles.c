/*	Contains independent-specific file-related routines.
	Platform-specific file-related routines are in XOPFilesMac.c and XOPFilesWin.c.
*/

/*	This file contains utilities for XOPs that open, read and write files.
	This includes file-loader specific routines.
	
	HR, 10/8/96: Split these routines out from XOPSupport.c.
*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h

// In XOPFilesMac.c and XOPFilesWin.c.
// int XOPCreateFile(const char* fullFilePath, int overwrite, int macCreator, int macFileType);

// In XOPFilesMac.c and XOPFilesWin.c.
// int XOPDeleteFile(const char* fullFilePath);

// In XOPFilesMac.c and XOPFilesWin.c.
// int XOPOpenFile(const char* fullFilePath, int readOrWrite, XOP_FILE_REF* fileRefPtr);

/*	XOPCloseFile(fileRef)

	Closes the referenced file.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPCloseFile is thread-safe.
*/
int
XOPCloseFile(XOP_FILE_REF fileRef)
{
	if (fclose(fileRef))
		return FILE_CLOSE_ERROR;
	return 0;
}

/*	XOPReadFile(fileRef, count, buffer, numBytesReadPtr)

	Reads count bytes from the referenced file into the buffer.
	
	If numBytesReadPtr is not NULL, stores the number of bytes read in
	*numBytesReadPtr.
	
	Returns 0 if OK or an error code.
	
	If bytes remain to be read in the file and you ask to read more bytes
	than remain, the remaining bytes are returned and the function result is
	zero. If no bytes remain to be read in the file and you ask to read bytes,
	no bytes are returned and the function result is FILE_EOF_ERROR.
	
	XOPReadFile is appropriate when you are reading data of variable size, in
	which case you do not want to consider it an error if the end of file is reached
	before reading all of the bytes that you requested. If you are reading a
	record of fixed size, use use XOPReadFile2 instead of XOPReadFile.
	
	Thread Safety: XOPReadFile is thread-safe.
*/
int
XOPReadFile(XOP_FILE_REF fileRef, UInt32 count, void* buffer, UInt32* numBytesReadPtr)
{
	UInt32 numBytesRead;
	
	if (count == 0) {
		if (numBytesReadPtr != NULL)
			*numBytesReadPtr = 0;
		return 0;
	}
	
	clearerr(fileRef);
	numBytesRead = (UInt32)fread(buffer, 1, count, fileRef);
	if (numBytesReadPtr != NULL)
		*numBytesReadPtr = numBytesRead;
	if (ferror(fileRef))
		return FILE_READ_ERROR;
	if (numBytesRead==0 && XOPAtEndOfFile(fileRef))
		return FILE_EOF_ERROR;			// We were at the end of file when asked to read some bytes.
	return 0;
}

/*	XOPReadFile2(fileRef, count, buffer, numBytesReadPtr)

	Reads count bytes from the referenced file into the buffer.
	
	If numBytesReadPtr is not NULL, stores the number of bytes read in
	*numBytesReadPtr.
	
	Returns 0 if OK or an error code.
	
	If bytes remain to be read in the file and you ask to read more bytes
	than remain, the remaining bytes are returned and the function result is
	FILE_EOF_ERROR.
	
	XOPReadFile2 is appropriate when you are reading a record of fixed size, in
	which case you want to consider it an error if the end of file is reached
	before reading all of the bytes in the record. If you are reading a record
	of variable size then you should use XOPReadFile instead of XOPReadFile2.
	
	Thread Safety: XOPReadFile2 is thread-safe.
*/
int
XOPReadFile2(XOP_FILE_REF fileRef, UInt32 count, void* buffer, UInt32* numBytesReadPtr)
{
	UInt32 numBytesRead;
	
	if (count == 0) {
		if (numBytesReadPtr != NULL)
			*numBytesReadPtr = 0;
		return 0;
	}
	
	clearerr(fileRef);
	numBytesRead = (UInt32)fread(buffer, 1, count, fileRef);
	if (numBytesReadPtr != NULL)
		*numBytesReadPtr = numBytesRead;
	if (ferror(fileRef))
		return FILE_READ_ERROR;
	if (numBytesRead < count) {				// We did not read all of the bytes requested.
		if (XOPAtEndOfFile(fileRef))
			return FILE_EOF_ERROR;			// We hit the end of file.
		return FILE_READ_ERROR;				// Some other occurred but ferror did not reflect it.
	}
	return 0;
}

/*	XOPReadFile64(fileRef, count, buffer, numBytesReadPtr)

	Reads count bytes from the referenced file into the buffer.

	Use XOPReadFile64 when you need to read potentially greater than 4GB
	in one read call. When called from a 32-bit XOP, XOPReadFile64 is limited
	to 4GB and returns an error if count is greater than 4GB.
	
	If numBytesReadPtr is not NULL, stores the number of bytes read in
	*numBytesReadPtr.
	
	Returns 0 if OK or an error code.
	
	If bytes remain to be read in the file and you ask to read more bytes
	than remain, the remaining bytes are returned and the function result is
	FILE_EOF_ERROR.
	
	Thread Safety: XOPReadFile64 is thread-safe.
*/
int
XOPReadFile64(XOP_FILE_REF fileRef, SInt64 count, void* buffer, SInt64* numBytesReadPtr)
{
	if (count == 0) {
		if (numBytesReadPtr != NULL)
			*numBytesReadPtr = 0;
		return 0;
	}

	#ifdef IGOR32
	{
		UInt32 numBytesRead32;
		int err;

		if (count > UINT_MAX)
			return FILE_READ_ERROR;

		err = XOPReadFile2(fileRef, (UInt32)count, buffer, &numBytesRead32);
		if (numBytesReadPtr != NULL)
			*numBytesReadPtr = numBytesRead32;
		return err;
	}
	#endif

	#ifdef IGOR64
	{
		SInt64 numBytesRead;
		clearerr(fileRef);
		numBytesRead = fread(buffer, 1, count, fileRef);
		if (numBytesReadPtr != NULL)
			*numBytesReadPtr = numBytesRead;
		if (ferror(fileRef))
			return FILE_READ_ERROR;
		if (numBytesRead < count) {				// We did not read all of the bytes requested.
			if (XOPAtEndOfFile(fileRef))
				return FILE_EOF_ERROR;			// We hit the end of file.
			return FILE_READ_ERROR;				// Some other occurred but ferror did not reflect it.
		}
	}
	#endif

	return 0;
}

/*	XOPWriteFile(fileRef, count, buffer, numBytesWrittenPtr)

	Writes count bytes from the buffer to the referenced file.
	
	If numBytesWrittenPtr is not NULL, stores the number of bytes written in
	*numBytesWrittenPtr.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPWriteFile is thread-safe.
*/
int
XOPWriteFile(XOP_FILE_REF fileRef, UInt32 count, const void* buffer, UInt32* numBytesWrittenPtr)
{
	UInt32 numBytesWritten;
	
	if (count == 0) {
		if (numBytesWrittenPtr != NULL)
			*numBytesWrittenPtr = 0;
		return 0;
	}
	
	numBytesWritten = (UInt32)fwrite(buffer, 1, count, fileRef);
	if (numBytesWrittenPtr != NULL)
		*numBytesWrittenPtr = numBytesWritten;
	if (numBytesWritten != count)
		return FILE_WRITE_ERROR;
	return 0;
}

/*	XOPWriteFile64(fileRef, count, buffer, numBytesWrittenPtr)

	Writes count bytes from the buffer to the referenced file.

	Use XOPWriteFile64 when you need to write potentially greater than 4GB
	in one write call. When called from a 32-bit XOP, XOPWriteFile64 is limited
	to 4GB and returns an error if count is greater than 4GB.
	
	If numBytesWrittenPtr is not NULL, stores the number of bytes written in
	*numBytesWrittenPtr.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPWriteFile64 is thread-safe.
*/
int
XOPWriteFile64(XOP_FILE_REF fileRef, SInt64 count, const void* buffer, SInt64* numBytesWrittenPtr)
{
	if (count == 0) {
		if (numBytesWrittenPtr != NULL)
			*numBytesWrittenPtr = 0;
		return 0;
	}

	#ifdef IGOR32
	{
		UInt32 numBytesWritten32;
		int err;

		if (count > UINT_MAX)
			return FILE_READ_ERROR;

		err = XOPWriteFile(fileRef, (UInt32)count, buffer, &numBytesWritten32);
		if (numBytesWrittenPtr != NULL)
			*numBytesWrittenPtr = numBytesWritten32;
		return err;
	}
	#endif

	#ifdef IGOR64
	{
		SInt64 numBytesWritten;
		numBytesWritten = fwrite(buffer, 1, count, fileRef);
		if (numBytesWrittenPtr != NULL)
			*numBytesWrittenPtr = numBytesWritten;
		if (numBytesWritten != count)
			return FILE_WRITE_ERROR;
	}
	#endif

	return 0;
}

/*	XOPGetFilePosition(fileRef, filePosPtr)

	Returns via filePosPtr the current file position of the referenced file.
	
	This routine does not work with files greater than 4GB. For very large files
	use XOPGetFilePosition2 instead.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPGetFilePosition is thread-safe.
*/
int
XOPGetFilePosition(XOP_FILE_REF fileRef, UInt32* filePosPtr)
{
	SInt32 pos;
	
	pos = ftell(fileRef);
	if (pos == -1L)
		return FILE_POS_ERROR;
	*filePosPtr = pos;
	return 0;
}

/*	XOPSetFilePosition(fileRef, filePos, mode)

	Sets the current file position in the referenced file.
	
	If mode is -1, then filePos is relative to the start of the file.
	If mode is 0, then filePos is relative to the current file position.
	If mode is 1, then filePos is relative to the end of the file.
	
	This routine does not work with files greater than 2GB. For very large files
	use XOPSetFilePosition2 instead.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPSetFilePosition is thread-safe.
*/
int
XOPSetFilePosition(XOP_FILE_REF fileRef, SInt32 filePos, int mode)
{
	int seekMode;
	
	switch(mode) {
		case -1:
			seekMode = SEEK_SET;
			break;
		case 0:
			seekMode = SEEK_CUR;
			break;
		case 1:
			seekMode = SEEK_END;
			break;
		default:
			return FILE_POS_ERROR;
	}
	
	if (fseek(fileRef, filePos, seekMode) != 0)
		return FILE_POS_ERROR;
	return 0;
}

/*	XOPGetFilePosition2(fileRef, filePosPtr)

	Returns via filePosPtr the current file position of the referenced file.
	
	XOPGetFilePosition2 is the same as XOPGetFilePosition except that the file position
	parameter is an SInt64* rather than a UInt32* and it works with files larger than 2 GB.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPGetFilePosition2 is thread-safe.
*/
int
XOPGetFilePosition2(XOP_FILE_REF fileRef, SInt64* filePosPtr)
{
	fpos_t pos;
	
	if (fgetpos(fileRef,&pos))
		return FILE_POS_ERROR;
	*filePosPtr = pos;
	return 0;
}

/*	XOPSetFilePosition2(fileRef, filePos)

	Sets the current file position in the referenced file.
	
	XOPSetFilePosition2 is the same as XOPSetFilePosition except that the file position
	parameter is an SInt64 rather than an SInt32 and it works with files larger than 2 GB
	and also it is lacking the mode parameter.
	
	dFilePos is relative to the start of the file.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPSetFilePosition2 is thread-safe.
*/
int
XOPSetFilePosition2(XOP_FILE_REF fileRef, SInt64 filePos)
{
	fpos_t pos;
	
	pos = (fpos_t)filePos;
	
	if (fsetpos(fileRef, &pos))
		return FILE_POS_ERROR;

	return 0;
}

/*	XOPAtEndOfFile(fileRef)

	Returns 1 if the current file position is at the end of file, 0 if not.
	
	Thread Safety: XOPAtEndOfFile is thread-safe.
*/
int
XOPAtEndOfFile(XOP_FILE_REF fileRef)
{
	if (feof(fileRef))				// Hit end of file?
		return 1;
	return 0;
}

/*	XOPNumberOfBytesInFile(fileRef, numBytesPtr)

	Returns via numBytesPtr the total number of bytes in the referenced file.
	
	This routine does not work with files greater than 4GB. For very large files
	use XOPNumberOfBytesInFile2 instead.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPNumberOfBytesInFile is thread-safe.
*/
int
XOPNumberOfBytesInFile(XOP_FILE_REF fileRef, UInt32* numBytesPtr)
{
	SInt32 originalPos;

	originalPos = ftell(fileRef);
	if (fseek(fileRef, 0, SEEK_END) != 0)
		return FILE_POS_ERROR;
	*numBytesPtr = ftell(fileRef);
	if (*numBytesPtr == -1L)
		return FILE_POS_ERROR;
	if (fseek(fileRef, originalPos, SEEK_SET) != 0)
		return FILE_POS_ERROR;
	return 0;
}

/*	XOPNumberOfBytesInFile2(fileRef, numBytesPtr)

	Returns via numBytesPtr the total number of bytes in the referenced file.
	
	XOPNumberOfBytesInFile2 is the same as XOPNumberOfBytesInFile except that the numBytesPtr
	parameter is an SInt64* rather than a UInt32* and it works with files larger than 2 GB.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: XOPNumberOfBytesInFile2 is thread-safe.
*/
#ifdef MACIGOR
int
XOPNumberOfBytesInFile2(XOP_FILE_REF fileRef, SInt64* numBytesPtr)
{
	fpos_t originalPos;
	fpos_t endPos;
	
	if (fgetpos(fileRef,&originalPos))
		return FILE_POS_ERROR;

	if (fseeko(fileRef, 0, SEEK_END) != 0)
		return FILE_POS_ERROR;

	if (fgetpos(fileRef,&endPos)) {
		fsetpos(fileRef,&originalPos);
		return FILE_POS_ERROR;
	}

	*numBytesPtr = endPos;

	if (fsetpos(fileRef,&originalPos))
		return FILE_POS_ERROR;

	return 0;
}
#endif
#ifdef WINIGOR
/*	HR, 080926, 1.63:

	I previously tried to implement this like the Macintosh version except using fseek
	instead of fseek0 because Windows does not have fseek0. This seemed to work but then
	inexplicably no longer worked.

	I have changed the implementation to use _filelengthi64.
*/
_CRTIMP __int64 __cdecl _filelengthi64(int);
int
XOPNumberOfBytesInFile2(XOP_FILE_REF fileRef, SInt64* numBytesPtr)
{
	int fd;			// C runtime file descriptor
	__int64 fileLength;

	fd = _fileno(fileRef);

	fileLength = _filelengthi64(fd);
	if (fileLength < 0)
		return FILE_POS_ERROR;

	*numBytesPtr = fileLength;

	return 0;
}
#endif

/*	XOPReadLine(fileRef, buffer, bufferLength)

	buffer points to a buffer into which the line of data is to be read.
	
	bufferLength is the size of the buffer. The buffer can hold bufferLength-1
	characters, plus the terminating null character.
	
	A line in the file may end with:
		<end-of-file>
		CR
		LF
		CRLF
	
	XOPReadLine reads the next line of text into the buffer and null-terminates it.
	The terminating CR, LF, or CRLF is not stored in the buffer.
	
	If numBytesReadPtr is not NULL, stores the number of bytes read in
	*numBytesReadPtr.
	
	Returns 0 if OK or a non-zero error code.
	
	The function result will be LINE_TOO_LONG_IN_FILE if there is not enough room in the
	buffer to read the entire line. It will be FILE_EOF_ERROR if we hit the end-of-file
	before reading any characters. It will be zero if we read any characters
	(even just a CR or LF) before hitting the end of the file.
	
	This routine was designed for simplicity of use. For applications that require
	blazing speed (e.g., reading files containing tens of thousands of lines or more),
	a more complex buffering scheme can improve performance considerably.
	
	Thread Safety: XOPReadLine is thread-safe.
*/
int
XOPReadLine(XOP_FILE_REF fileRef, char* buffer, UInt32 bufferLength, UInt32* numBytesReadPtr)
{
	char* bufPtr;
	char ch;
	int done;
	UInt32 count, numBytesRead;
	int numBytesToPutBack;
	int err;
	
	/*	The nominalLineLength variable is used to attempt to tune this routine
		to the actual line length encountered.
	
		Ideally, we would read the exact right number of characters each time
		we are called to read a line. This is not possible. The next best thing
		is to read slightly more characters than we need so that we go through
		the outer loop below only once.
	*/
	static UInt32 nominalLineLength = 80;
	
	bufPtr = buffer;
	
	done = 0;
	numBytesToPutBack = 0;
	numBytesRead = 0;
	while(!done) {
		count = bufferLength - numBytesRead - 1;		// Maximum that could be read.
		if (count <= 0) {
			err = LINE_TOO_LONG_IN_FILE;
			break;
		}
		if (count > nominalLineLength)
			count = nominalLineLength;					// Maximum that we want to read at one time.
		
		err = XOPReadFile(fileRef, count, bufPtr, &count);
		if (err != 0) {
			if (numBytesRead>0 && XOPAtEndOfFile(fileRef))
				err = 0;
			break;
		}
		
		while(count > 0) {
			if (*bufPtr==CR_CHAR) {
				if (count > 1) {
					// Check for LF following CR.
					numBytesToPutBack = count-1;
					if (bufPtr[1] == LF_CHAR)
						numBytesToPutBack -= 1;
				}
				else {
					// Need to read another character in order to check for LF following CR.
					XOPReadFile(fileRef, 1, &ch, &count);
					if (count>0 && ch!=LF_CHAR)
						numBytesToPutBack = 1;
				}
				done = 1;
				break;
			}
			else {
				if (*bufPtr==LF_CHAR) {
					numBytesToPutBack = count-1;
					done = 1;
					break;
				}
			}
			bufPtr += 1;
			numBytesRead += 1;
			count -= 1;
		}
	}
	
	if (numBytesToPutBack > 0)
		XOPSetFilePosition(fileRef, -numBytesToPutBack, 0);
	
	if (err == 0)
		nominalLineLength = numBytesRead + 10;
	
	buffer[numBytesRead] = 0;

	if (numBytesReadPtr != NULL)
		*numBytesReadPtr = numBytesRead;

	return err;
}

// In XOPFilesMac.c and XOPFilesWin.c.
// int FullPathPointsToFile(const char* fullPath);

// In XOPFilesMac.c and XOPFilesWin.c.
// int FullPathPointsToFolder(const char* fullPath);

/*	WinToMacPath(path)

	This routine converts a Windows path into a Macintosh path by replacing
	':\' with ':' and  '\' with ':'. HOWEVER, it does not change a UNC volume name.
	Thus,
		C:\A\B\C				=>		C:A:B:C
		\\server\share\A\B\C	=>		\\server\share:A:B:C
	The volume name is "\\server\share" whether the path is a Mac path or
	a Windows path.
	
	Also, leading periods are changed to colons. For example, "..\FolderA\FileB"
	is changed to "::FolderA:FileB".
	
	If the path is already a Mac path, it will do nothing.
	
	NOTE:	The path may be shorter on output than in was on input
			(':\' or '.\' changed to ':').
			
	This routine is Asian-character-set aware.
	
	Function result is 0 if OK or error code.
	
	Thread Safety: WinToMacPath is thread-safe with Igor Pro 6.20 or later.
*/
int
WinToMacPath(char path[MAX_PATH_LEN+1])
{
	return (int)CallBack1(WIN_TO_MAC_PATH, path);
}

/*	MacToWinPath(path)

	This routine converts a Macintosh path into a Windows path by replacing
	':' with ':\' at the start of a full path and replacing ':' with '\' elsewhere.
	HOWEVER, it does not change a UNC volume name. Thus,
		C:A:B:C					=>		C:\A\B\C
		\\server\share:A:B:C	=>		\\server\share\A\B\C
	The volume name is "\\server\share" whether the path is a Mac path or
	a Windows path.

	NOTE: If a Mac path contains a '\' character, the resulting path will
		  not work as a Windows path. Therefore, '\' characters must not be
		  used in Mac paths.
	
	Also, leading colons are changed to periods. For example, "::FolderA:FileB"
	is changed to "..\FolderA\FileB".
	
	If the path is already a Windows path, it will do nothing.
	
	NOTE:	The path may be longer on output than in was on input ('C:' changed to
			'C:\' or ':' changed to '.\'). The buffer is assumed to be MAX_PATH_LEN+1
			characters long. MacToWinPath will not overwrite the buffer. It will
			generate an error if the output path can not fit in MAX_PATH_LEN characters.
	
	Function result is 0 if OK or an error code (e.g., PATH_TOO_LONG).
	
	Thread Safety: MacToWinPath is thread-safe with Igor Pro 6.20 or later.
*/
int
MacToWinPath(char path[MAX_PATH_LEN+1])
{
	return (int)CallBack1(MAC_TO_WIN_PATH, path);
}

// In XOPFilesMac.c and XOPFilesWin.c.
// int GetNativePath(const char* filePathIn, char filePathOut[MAX_PATH_LEN+1]);

/*	EscapeBackslashesInUNCVolumeName(macFilePath)
	
	This routine is used when we are generating a literal string containing a file path
	that may refer to a Windows server. For a Windows server, the volume name is a
	Universal Name Convention name, something like \\server\share. This is true even
	if the path is a Macintosh path. If we are generating a command using this path,
	as in the LoadWave dialog, we need to escape the backslashes because Igor interprets
	backslashes as escape characters.
	
	Returns the number of characters added to the string or -1 if the string is too long.
	Note that macFilePath is assumed to be MAX_PATH_LEN characters long.
	
	Technical Note: This routine will not work on a UNC name that uses Asian
	characters if the second byte of any character is equivalent to a backslash
	character. This is because the strchr function below does not know about
	Asian characters. I would be possible but difficult to fix this. The situation
	should be very rare if it happens at all.
	
	Thread Safety: EscapeBackslashesInUNCVolumeName is thread-safe.
*/
int
EscapeBackslashesInUNCVolumeName(char macFilePath[MAX_PATH_LEN+1])
{	
	int maxCharsThatCanBeAdded;
	char macFilePathIn[MAX_PATH_LEN+1];
	char* pIn;
	char* pOut;
	char* pBackslash;
	char backslashChar;

	backslashChar = '\\';
	if (macFilePath[0]!=backslashChar || macFilePath[1]!=backslashChar)
		return 0;									// This is not a Windows server name.

	maxCharsThatCanBeAdded = MAX_PATH_LEN - (int)strlen(macFilePath);
	if (maxCharsThatCanBeAdded < 3)
		return -1;									// We need to add three backslashes.
	
	strcpy(macFilePathIn, macFilePath);
	macFilePath[2] = backslashChar;					// Add first extra backslash.
	macFilePath[3] = backslashChar;					// Add second extra backslash.
	
	pIn = &macFilePathIn[2];						// Points to the original input data.
	pOut = &macFilePath[4];							// Points to the output data.
	pBackslash = strchr(pIn, backslashChar);		// Points to the third backslash in the input data. NULL if there is no third backslash - this would be a bogus server name.

	while(*pIn != 0) {
		if (pIn == pBackslash)
			*pOut++ = backslashChar;				// Add third extra backslash.
		*pOut++ = *pIn++;		
	}
	*pOut = 0;
	
	return 3;
}

/*	GetDirectoryAndFileNameFromFullPath(fullFilePath, dirPath, fileName)

	fullFilePath is a full path to a file. It may be a Macintosh path (using colons)
	or a Windows path (using backslashes).
	
	On output, dirPath is the full native path to the folder containing the file.
	This path includes a trailing colon (Macintosh) or backslash (Windows).
	
	On output, fileName contains just the name of the file.
	
	Returns 0 if OK or an error code.
	
	Note that GetDirectoryAndFileNameFromFullPath does not know or care if
	the file exists or if the directories referenced in the input path exist.
	It merely separates the file name part from the full path.
	
	A simple implementation of this routine would merely search for colon or
	backslash characters. However, this simple approach causes problems on Asian
	systems that use two-byte characters. The problem is that the second byte of
	a two-byte character may have the same code as a colon or backslash, causing
	the simple implementation to mistakenly take it for a path separator.
	
	GetDirectoryAndFileNameFromFullPath takes a more complex approach to avoid
	this problem. To achieve this, the routine has to know the character encoding
	governing the fullFilePath parameter. GetDirectoryAndFileNameFromFullPath assumes
	that the system default character encoding governs the fullFilePath parameter.
	
	Thread Safety: GetDirectoryAndFileNameFromFullPath is thread-safe with Igor Pro 6.20 or later.
*/
int
GetDirectoryAndFileNameFromFullPath(const char* fullFilePath, char dirPath[MAX_PATH_LEN+1], char fileName[MAX_FILENAME_LEN+1])
{
	return (int)CallBack3(GET_DIR_AND_FILE_FROM_FULL_PATH, (void*)fullFilePath, dirPath, fileName);
}

/*	IsFullPathAndFileName(filePath)

	filePath can be a Windows or Macintosh path.
	
	Returns 1 if filePath has the syntax of a full path and file name, 0 if not.
	
	Note that this does not tell you if filePath actually points to an existing file.
	For that, use FullPathPointsToFile.
	
	Thread Safety: IsFullPathAndFileName is thread-safe with Igor Pro 6.20 or later.
*/
static int
IsFullPathAndFileName(const char* filePath)
{
	char macFilePath[MAX_PATH_LEN+1];
	const char* p;
	int err;
	
	// WinToMacPath does nothing if it is already a Macintosh path.
	strcpy(macFilePath, filePath);
	if (err = WinToMacPath(macFilePath))
		return err;
	
	if (*macFilePath == ':')						// Starts with colon ?
		return 0;									// It's a partial path.
	
	p = strrchr2(macFilePath, ':');					// Find last colon.
	if (p == NULL)
		return 0;									// Has no colon.
		
	p += 1;											// Make sure it doesn't end with colon.
	if (*p)											// Next char is not null ?
		return 1;

	return 0;
}

/*	GetLeafName(filePath, name)

	filePath is either "", a valid file name or a valid path to a file.
	
	name must be able to hold MAX_FILENAME_LEN+1 bytes.
	
	Returns via name the leaf part of the path, if it is a path or the contents of filePath
	if is not a path.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: GetLeafName is thread-safe with Igor Pro 6.20 or later.
*/
int
GetLeafName(const char* filePath, char name[MAX_FILENAME_LEN+1])
{
	char macFilePath[MAX_PATH_LEN+1];
	const char* p;
	int err;
	
	*name = 0;
	
	// WinToMacPath does nothing if it is already a Macintosh path.
	strcpy(macFilePath, filePath);
	if (err = WinToMacPath(macFilePath))
		return err;
	
	p = strrchr2(macFilePath, ':');					// Find last colon.
	if (p == NULL)									// Has no colon?
		p = filePath;								// filePath is either "" or just a file name.
	else
		p += 1;										// filePath is either "" or a path and file name.
	
	if (strlen(p) > MAX_FILENAME_LEN)
		return STR_TOO_LONG;
	strcpy(name, p);

	return 0;
}

/*	GetFullPathFromSymbolicPathAndFilePath(symbolicPathName, filePath, fullFilePath)
	
	symbolicPathName is the name of an Igor symbolic path or "" if no symbolic
	path is to be used.
	
	filePath is either a full path, a partial path, or a simple file name.
	It may use Macintosh or Windows path conventions.
	
	fullFilePath is an output and will contain the full path to the file referenced
	by symbolicPathName and filePath. The returned path will use Macintosh path
	conventions on Macintosh and Windows path conventions on Windows.
	
	This routine is used by file loader XOPs to get a full path to a file based
	on the typical inputs to a file loader, namely an optional Igor symbolic path
	and an optional full or partial file path or file name.
	
	The two most common cases are:
		LoadWave <full path to file>
		LoadWave/P=<symbolic path name> <file name>
		
	where <file name> conotes a simple file name.
	
	Less common cases that this routine also handles are:
		LoadWave/P=<symbolic path name> <full path to file>		// Symbolic path is ignored.
		LoadWave/P=<symbolic path name> <partial path to file>
	
	In the following cases, the full path to the file can not be determined,
	so GetFullPathFromSymbolicPathAndFilePath returns an error. This would cause
	a file loader to display an open file dialog:
		LoadWave <file name>
		LoadWave <partial path to file>
	
	filePath and fullFilePath may point to the same storage, in which case the
	output string will overwrite the input string.
	
	This routine does not check that the output path is valid or points to an existing
	file. This makes the routine useable for applications in which you are creating the
	file as well as applications in which you are reading the file. If you want to verify
	that the output path points to an existing file, use the FullPathPointsToFile routine.
	
	Returns 0 if it was able to create the full path or an error code if not.
	
	Thread Safety: GetFullPathFromSymbolicPathAndFilePath is thread-safe with Igor Pro 6.20 or later.
*/
int
GetFullPathFromSymbolicPathAndFilePath(const char* symbolicPathName, const char filePath[MAX_PATH_LEN+1], char fullFilePath[MAX_PATH_LEN+1])
{
	char symbolicPathPath[MAX_PATH_LEN+1];
	char filePath2[MAX_PATH_LEN+1];
	int err;
	
	if (err = GetNativePath(filePath, filePath2))
		return err;

	if (*symbolicPathName==0) {
		strcpy(fullFilePath, filePath2);
	}
	else {
		if (IsFullPathAndFileName(filePath2)) {		// HR, 020806: Used IsFullPathAndFileName instead of FullPathPointsToFile.
			// Although we have a symbolic path, we don't use it because the input path was already a full path.
			strcpy(fullFilePath, filePath2);
		}
		else {
			if (err = GetPathInfo2(symbolicPathName, symbolicPathPath))		// symbolicPathPath is native.
				return err;
		
			if (err = ConcatenatePaths(symbolicPathPath, filePath2, fullFilePath))
				return err;
		}
	}
	
	// HR, 020806: Added this.
	if (!IsFullPathAndFileName(fullFilePath))
		return FILE_OPEN_ERROR;
	
	return 0;
}

/*	ConcatenatePaths(pathIn1, nameOrPathIn2, pathOut[MAX_PATH_LEN+1])

	Concatenates pathIn1 and nameOrPathIn2 into pathOut. pathOut will be
	a native path. The input paths may use Macintosh or Windows conventions.
	
	pathIn1 is a full path to a directory. It can end with zero or one separator.
	
	nameOrPathIn2 is either a file name, a folder name or a partial path to a file
	or folder. It can end with zero or one separator.
	
	pathOut can point to the same memory as either of the input parameters.

	The target file or folder does not need to already exist.
	
	For pathIn1, any of the following are legal.
		"hd:FolderA:FolderB"
		"hd:FolderA:FolderB:"
		"C:\FolderA\FolderB"
		"C:\FolderA\FolderB\"
	
	For nameOrPathIn2, any of the following are legal.
		"FileA"
		"FolderC"
		":FolderC"
		"\FolderC"
		".FolderC"					// Legal in a Windows path only.
		"::FolderC"
		"\\FolderC"
		"..FolderC"					// Legal in a Windows path only.
		"FolderC:FileA"
		"FolderC\FileA"
		"\FolderC:FileA"
		"\FolderC\FileA"
	
	Here are some examples.
		"hd:FolderA:FolderB:"				+	"FolderC"		=>		"hd:FolderA:FolderB:FolderC"
		"hd:FolderA:FolderB:"				+	":FolderC"		=>		"hd:FolderA:FolderB:FolderC"
		"hd:FolderA:FolderB:"				+	"::FolderC"		=>		"hd:FolderA:FolderC"

		"C:\FolderA\FolderB\"				+	"FolderC"		=>		"C:\FolderA\FolderB\FolderC"
		"C:\FolderA\FolderB\"				+	"\FolderC"		=>		"C:\FolderA\FolderB\FolderC"
		"C:\FolderA\FolderB\"				+	"\\FolderC"		=>		"C:\FolderA\FolderC"

		"\\server\share\FolderA\FolderB\"	+	"FolderC"		=>		"\\server\share\FolderA\FolderB\FolderC"
		"\\server\share\FolderA\FolderB\"	+	"\FolderC"		=>		"\\server\share\FolderA\FolderB\FolderC"
		"\\server\share\FolderA\FolderB\"	+	"\\FolderC"		=>		"\\server\share\FolderA\FolderC"
	
	Multiple colons or backslashes in nameOrPathIn2 mean that we want to back up,
	starting from the folder specified by pathIn1.
	
	Returns 0 or error code. In case of error, the contents of pathOut is undefined.
	
	Thread Safety: ConcatenatePaths is thread-safe with Igor Pro 6.20 or later.
*/
int
ConcatenatePaths(const char* pathIn1, const char* nameOrPathIn2, char pathOut[MAX_PATH_LEN+1])
{
	return (int)CallBack3(CONCATENATE_PATHS, (void*)pathIn1, (void*)nameOrPathIn2, (void*)pathOut);
}

/*	ParseFilePath(mode, pathIn, separator, whichEnd, whichElement, pathOut)

	The ParseFilePath function provides the ability to manipulate file paths and to extract
	sections of file paths.
	
	This XOPSupport routine works the same as the Igor built-in ParseFilePath function.
	For further explanation, see the documentation for ParseFilePath in the Igor Reference help file.
	
	The output is returned via the pathOut parameter which must be able to hold MAX_PATH_LEN+1 bytes.
	
	The function result is 0 for success or an Igor error code.

	Added in Igor Pro 6.20B03. If you call this when running with an earlier version,
	it will return IGOR_OBSOLETE.
	
	Thread Safety: ParseFilePath is thread-safe with Igor Pro 6.20 or later.
*/
int
ParseFilePath(int mode, const char* pathIn, const char* separator, int whichEnd, int whichElement, char pathOut[MAX_PATH_LEN+1])
{
	return (int)CallBack6(PARSE_FILE_PATH, (void*)mode, (void*)pathIn, (void*)separator, (void*)whichEnd, (void*)whichElement, pathOut);
}

/*	SpecialDirPath(pathID, domain, flags, createDir, pathOut)

	The SpecialDirPath function returns a full path to a file system directory specified
	by pathID and domain. It provides a programmer with a way to access directories
	of special interest, such as the preferences directory, the Igor Pro User Files
	directory and the desktop directory.
	
	This XOPSupport routine works the same as the Igor built-in SpecialDirPath function.
	For further explanation, see the documentation for SpecialDirPath in the Igor Reference help file.
	
	The output is returned via the pathOut parameter which must be able to hold MAX_PATH_LEN+1 bytes.
	
	The function result is 0 for success or an Igor error code.

	Added in Igor Pro 6.20B03. If you call this when running with an earlier version,
	it will return IGOR_OBSOLETE.
	
	Thread Safety: SpecialDirPath is thread-safe with Igor Pro 6.20 or later.
*/
int
SpecialDirPath(const char* pathID, int domain, int flags, int createDir, char pathOut[MAX_PATH_LEN+1])
{
	return (int)CallBack5(SPECIAL_DIR_PATH, (void*)pathID, (void*)domain, (void*)flags, (void*)createDir, pathOut);
}

// In XOPFilesMac.c and XOPFilesWin.c.
// int XOPOpenFileDialog(const char* prompt, const char* fileFilterStr, int* fileIndexPtr, const char* initialDir, char filePath[MAX_PATH_LEN+1]);

// In XOPFilesMac.c and XOPFilesWin.c.
// int XOPSaveFileDialog(const char* prompt, const char* fileFilterStr, int* fileIndexPtr, const char* initialDir, const char* defaultExtensionStr, char filePath[MAX_PATH_LEN+1]);

/*	FileLoaderMakeWave(column, waveName, numPoints, fileLoaderFlags, waveHandlePtr)

	FileLoaderMakeWave makes a wave with numPoints points and numeric type as specified
	by fileLoaderFlags.
	
	The function result is 0 or an error code.
	It returns a handle to the wave via waveHandlePtr.
	
	fileLoaderFlags is interpreted using the standard file loader flag bit definitions.
	See XOP.h.
	
	NOTE:	In the event of a name conflict, FileLoaderMakeWave can change the contents of
			waveName. waveName must be able to hold MAX_OBJ_NAME characters.
	
	Thread Safety: FileLoaderMakeWave is thread-safe with Igor Pro 6.20 or later.
*/
int
FileLoaderMakeWave(int column, char *waveName, CountInt numPoints, int fileLoaderFlags, waveHndl *waveHandlePtr)
{	
	int type, overwrite;
	char newName[MAX_OBJ_NAME+2];
	char temp[128];
	int result;

	type = (fileLoaderFlags & FILE_LOADER_DOUBLE_PRECISION) ? NT_FP64:NT_FP32;
	overwrite = fileLoaderFlags & FILE_LOADER_OVERWRITE;
	result = MakeWave(waveHandlePtr,waveName,numPoints,type,overwrite);
	
	/*	If error other than NOMEM or NAME_WAV_CONFLICT, it's probably a conflict
		with an operation or function. Try again using a different wave name.
	*/
	if (result && result != NOMEM && result != NAME_WAV_CONFLICT) {
		sprintf(newName, "X_%s", waveName);
		SanitizeWaveName(newName, column);
		if (!(fileLoaderFlags&FILE_LOADER_QUIET)) {
			sprintf(temp, "Name conflict making %s, name changed to %s\015", waveName, newName);
			XOPNotice(temp);
		}
		result = MakeWave(waveHandlePtr,newName,numPoints,type,overwrite);
		if (result == 0)
			strcpy(waveName, newName);
	}

	if (result && !(fileLoaderFlags&FILE_LOADER_QUIET)) {
		sprintf(temp, "Error making %s\015", waveName);
		XOPNotice(temp);
	}
	return(result);
}

/*	DoSetFileLoaderOutputVariables(runningInUserFunction, fileNameOrPath, numWavesLoaded, waveNames)

	See SetFileLoaderOutputVariables and SetRuntimeFileLoaderOutputVariables for details.
	
	Returns 0 or error code.
	
	Thread Safety: DoSetFileLoaderOutputVariables is thread-safe with Igor Pro 6.20 or later.
*/
static int
DoSetFileLoaderOutputVariables(int runningInUserFunction, const char* fileNameOrPath, int numWavesLoaded, const char* waveNames)
{
	const char* fullPathPtr;
	const char* fileNamePtr;
	int result;
	
	fullPathPtr = fileNameOrPath;						// Assume full path.
	
	// Find start of the leaf name.
	fileNamePtr = strrchr2(fileNameOrPath, ':');
	if (fileNamePtr!=NULL && fileNamePtr[1]!='\\') {	// Found colon other than the one in "C:\" ?
		fileNamePtr += 1;								// Point to character after the last colon.
	}
	else {
		fileNamePtr = strrchr2(fileNameOrPath, '\\');
		if (fileNamePtr != NULL) {
			fileNamePtr += 1;							// Point to character after the last backslash.
		}
		else {
			// No colon and no backslash.
			fileNamePtr = fileNameOrPath;
			fullPathPtr = "";							// No full path.
		}
	}
	
	// HR, 981103: Added setting of S_path for XOP Toolkit 3.1.
	// For backward compatibility, we skip creating S_path if the XOP did not pass in a full path.
	if (*fullPathPtr != 0) {							// We have a full path?
		char dirPath[MAX_PATH_LEN+1];
		char fileName2[MAX_FILENAME_LEN+1];
		
		if (GetDirectoryAndFileNameFromFullPath(fullPathPtr, dirPath, fileName2) == 0) {
			if (FullPathPointsToFolder(dirPath)) {
				WinToMacPath(dirPath);
				if (runningInUserFunction) {
					if (result = SetOperationStrVar("S_path", dirPath))
						return result;
				}
				else {
					if (result = SetIgorStringVar("S_path", dirPath, 0))
						return result;
				}
			}
		}
	}
	
	if (runningInUserFunction) {
		if (result = SetOperationStrVar("S_fileName", fileNamePtr))
			return result;
		
		if (result = SetOperationStrVar("S_waveNames", waveNames))
			return result;
		
		if (result = SetOperationNumVar("V_flag", numWavesLoaded))
			return result;
	
	}
	else {
		if (result = SetIgorStringVar("S_fileName", fileNamePtr, 0))
			return result;
		
		if (result = SetIgorStringVar("S_waveNames", waveNames, 0))
			return result;
		
		if (result = SetIgorIntVar("V_flag", numWavesLoaded, 0))
			return result;
	}
	
	return 0;
}

/*	SetFileLoaderOutputVariables(fileNameOrPath, numWavesLoaded, waveNames)

	If your external operation uses Operation Handler, use SetFileLoaderOperationOutputVariables
	instead of this routine.
	
	Call SetFileLoaderOutputVariables after a file load to set the "standard"
	file loader output globals:
		S_fileName			The name of the file loaded.
		S_path				The full path to the folder containing the file. See note below.
		V_flag				The number of waves loaded.
		S_waveNames			Semicolon-separate list of wave names
							(e.g. "wave0;wave1;wave2;").
							
	fileNameOrPath can be either just the file name (e.g., "Data File") or a full
	path including a file name (e.g., "hd:Data Folder:Data File"). Passing a full
	path is recommended. If it is a full path, it can use either Macintosh or Windows syntax.
	
	If fileNameOrPath is a full path, SetFileLoaderOutputVariables stores the path
	to the folder containing the file in S_path and stores the simple file name in
	S_fileName. In this case, the path uses Macintosh path syntax and includes
	a trailing colon. 
	
	If fileNameOrPath is a simple file name, SetFileLoaderOutputVariables does not
	set or create S_path and stores the simple file name in S_fileName. You should
	pass a full path so that S_path will be set.
	
	Returns 0 or error code.
	
	Thread Safety: SetFileLoaderOutputVariables is thread-safe with Igor Pro 6.20 or later.
*/
int
SetFileLoaderOutputVariables(const char* fileNameOrPath, int numWavesLoaded, const char* waveNames)
{
	return DoSetFileLoaderOutputVariables(0, fileNameOrPath, numWavesLoaded, waveNames);
}

/*	SetFileLoaderOperationOutputVariables(runningInUserFunction, fileNameOrPath, numWavesLoaded, waveNames)

	Call SetFileLoaderOperationOutputVariables at the end of an Operation Handler file
	load operation to set the standard file loader output globals.
	
	When called from the command line it creates global variables.
	
	When called from a user-defined function or from a macro it creates local variables.
	
	You obtain the value to pass for the runningInUserFunction parameter from the
	calledFromFunction field of your operation’s runtime parameter structure.

	When you register your operation via RegisterOperation, you must specify that your
	operation sets the numeric variable V_flag and the string variables S_fileName, S_path,
	and S_waveNames. See SimpleLoadWaveOperation.c for an example.

	SetFileLoaderOperationOutputVariables sets the following standard file loader output globals:
		S_fileName			The name of the file loaded.
		S_path				The full path to the folder containing the file. See note below.
		V_flag				The number of waves loaded.
		S_waveNames			Semicolon-separate list of wave names
							(e.g. "wave0;wave1;wave2;").
							
	fileNameOrPath can be either just the file name (e.g., "Data File") or a full
	path including a file name (e.g., "hd:Data Folder:Data File"). Passing a full
	path is recommended. If it is a full path, it can use either Macintosh or Windows syntax.
	
	If fileNameOrPath is a full path, SetFileLoaderOperationOutputVariables stores
	the path to the folder containing the file in S_path and stores the simple file
	name in S_fileName. In this case, the path uses Macintosh path syntax and includes
	a trailing colon. 
	
	If fileNameOrPath is a simple file name, SetFileLoaderOperationOutputVariables
	does not set or create S_path and stores the simple file name in S_fileName.
	You should pass a full path so that S_path will be set.
	
	Returns 0 or error code.
	
	Thread Safety: SetFileLoaderOperationOutputVariables is thread-safe with Igor Pro 6.20 or later.
*/
int
SetFileLoaderOperationOutputVariables(int runningInUserFunction, const char* fileNameOrPath, int numWavesLoaded, const char* waveNames)
{
	return DoSetFileLoaderOutputVariables(runningInUserFunction, fileNameOrPath, numWavesLoaded, waveNames);
}
	
/*	PrepareLoadIgorData(ldiPtr, refNumPtr, topFIHPtr)

	This routine is for use by the WaveMetrics Data Browser only.
	
	ldiPtr is a pointer to a LoadDataInfo structure. This structure is private to WaveMetrics.
	topFIHPtr is a pointer to a LoadFileInfo structure handle. This structure is private to WaveMetrics.
	
	Thread Safety: PrepareLoadIgorData is not thread-safe.
*/
int
PrepareLoadIgorData(struct LoadDataInfo* ldiPtr, int* refNumPtr, struct LoadFileInfo*** topFIHPtr)
{
	if (!CheckRunningInMainThread("PrepareLoadIgorData"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(PREPARE_LOAD_IGOR_DATA, ldiPtr, refNumPtr, topFIHPtr);
}
	
/*	LoadIgorData(ldiPtr, refNum, topFIH, destDataFolderH)

	This routine is for use by the WaveMetrics Data Browser only.
	
	ldiPtr is a pointer to a LoadDataInfo structure. This structure is private to WaveMetrics.
	topFIH is a handle to a LoadFileInfo structure. This structure is private to WaveMetrics.
	
	Thread Safety: LoadIgorData is not thread-safe.
*/
int
LoadIgorData(struct LoadDataInfo* ldiPtr, int refNum, struct LoadFileInfo** topFIH, DataFolderHandle destDataFolderH)
{
	if (!CheckRunningInMainThread("LoadIgorData"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(DO_LOAD_IGOR_DATA, ldiPtr, (void*)refNum, topFIH, destDataFolderH);
}
	
/*	EndLoadIgorData(ldiPtr, refNum, topFIH)

	This routine is for use by the WaveMetrics Data Browser only.
	
	ldiPtr is a pointer to a LoadDataInfo structure. This structure is private to WaveMetrics.
	topFIH is a handle to a LoadFileInfo structure. This structure is private to WaveMetrics.
	
	Thread Safety: EndLoadIgorData is not thread-safe.
*/
int
EndLoadIgorData(struct LoadDataInfo* ldiPtr, int refNum, struct LoadFileInfo** topFIH)
{
	if (!CheckRunningInMainThread("EndLoadIgorData"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(END_LOAD_IGOR_DATA, ldiPtr, (void*)refNum, topFIH);
}

/*	SaveIgorData(sdiPtr, topDataFolderH)

	This routine is for use by the WaveMetrics Data Browser only.
	
	sdiPtr is a pointer to a SaveDataInfo structure. This structure is private to WaveMetrics.
	
	Thread Safety: SaveIgorData is not thread-safe.
*/
int
SaveIgorData(struct SaveDataInfo* sdiPtr, DataFolderHandle topDataFolderH)
{
	if (!CheckRunningInMainThread("SaveIgorData"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack2(DO_SAVE_IGOR_DATA, sdiPtr, topDataFolderH);
}
