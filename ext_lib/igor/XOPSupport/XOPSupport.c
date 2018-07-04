// XOPSupport.c - Support routines for Igor XOPs

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
#include "NamedFIFO.h"

// Global Variables.
IORecHandle XOPRecHandle;						// The XOP's ioRecHandle.
int XOPModified = 0;
int igorVersion = 0;							// Set by XOPInit. Example: 128 for version 1.28.
static PSInt gMainThreadID = 0;					// HR, 091123: Used in RunningInMainThread.
static int gUseThreadsafeCallbackMethod = 0;	// HR, 091118: Igor now supports threadsafe callback method.


// *** Utility Routines ***

// Capitalize removed from XOP Toolkit 6. Use CmpStr instead of calling Capitalize and strcmp.

/*	CmpStr(str1, str2)

	Compares two C strings case insensitively.
	
	Result is the same as for the strcmp function.
	
	Thread Safety: CmpStr is thread-safe.
*/
int
CmpStr(const char *s1, const char *s2)
{
	unsigned const char *str1= (unsigned const char *) s1;
	unsigned const char *str2= (unsigned const char *) s2;
	unsigned int c1,c2;
	int result= 0;						/* assume s1 == s2 */
	
	do {
		c1= *str1++;
		if( c1 >= 'a' && c1 <= 'z' )	/* if (islower(c1))     */
			c1 -= 'a' - 'A';			/*     c1= toupper(c1); */
	
		c2= *str2++;
		if( c2 >= 'a' && c2 <= 'z' )	/* if (islower(c2))     */
			c2 -= 'a' - 'A';			/*     c2= toupper(c2); */

		if( c1 != c2 ) {
			if( c1 < c2 )
				result= -1;				/* s1 < s2 */
			else
				result= 1;				/* s1 > s2 */
			break;
		}

	} while ( c1 );	/* falls through if c1 == 0 (and c2 == 0) because of above if(); s1 == s2 */
	
	return result;
}

/*	strchr2(str, ch)

	strchr2 is like the standard C strchr function except that it is
	Asian-language-aware. 
	
	Returns a pointer to the first occurrence of ch in the null-terminated string
	str or NULL if there is no such occurrence. 
	
	On a system that uses an Asian script system as the default script, strchr2
	knows about two-byte characters. For example, if you are searching for a
	backslash in a full path and if the path contains Asian characters, and if the
	second byte of an Asian character has the same code as the backslash character,
	strchr will mistakenly find this second byte while strchr2 will not. 
	
	On a system that does not use an Asian script system as the default script,
	strchr2 is just like strchr. 
	
	Thread Safety: strchr2 is thread-safe with Igor Pro 6.20 or later.
*/
const char*
strchr2(const char* str, int ch)
{
	return (const char*)CallBack2(STRCHR2, (void*)str, (void*)ch);
}

/*	strrchr2(str, ch)

	strrchr2 is like the standard C strrchr function except that it is
	Asian-language-aware. 
	
	Returns a pointer to the last occurrence of ch in the null-terminated string
	str or NULL if there is no such occurrence. 
	
	On a system that uses an Asian script system as the default script, strrchr2
	knows about two-byte characters. For example, if you are searching for a
	backslash in a full path and if the path contains Asian characters, and if the
	second byte of an Asian character has the same code as the backslash character,
	strrchr will mistakenly find this second byte while strrchr2 will not.
	
	On a system that does not use an Asian script system as the default script,
	strrchr2 is just like strrchr. 
	
	Thread Safety: strrchr2 is thread-safe with Igor Pro 6.20 or later.
*/
const char*
strrchr2(const char* str, int ch)
{
	return (const char*)CallBack2(STRRCHR2, (void*)str, (void*)ch);
}

/*	EscapeSpecialCharacters(input, inputLength, output, outputBufferSize, numCharsOutPtr)

	Converts CR, LF, tab, double-quote and backslash to escape sequences. This is used
	mostly when you are generating a command to be executed on Igor's command line or
	via Igor's Execute operation. Conversion is necessary because Igor interprets backslash
	as an escape character.
	
	For example, if you are generating a file loading command and the file uses a Windows
	UNC path of the form \\server\share, the backslashes in the command must be escaped to
	give \\\\server\\share because Igor interprets \\ as an escape sequence that evaluates
	to a single backslash.
	
	input is a pointer to a string to be escaped. Input may be null-terminated or not.
	
	inputLength is the number of characters to be escaped.

	NOTE: If input is null-terminated and you want output to also be null-terminated,
	then inputLength must include the null terminator character.
	
	output is a pointer to a buffer in which escaped output text is returned.
	
	outputBufferSize is the total size of the output buffer.
	
	numCharsOutPtr is set to the number of characters stored in output. This includes
	the null-terminator if it was included in inputLength.

	input and output can point to the same memory.
	
	Returns 0 if OK or an error code.
	
	Example:
		char str[32];
		int length;
		// Assume str is set to contain a C string containing characters that may need to be escaped.
		length = strlen(str) + 1;		// Include null-terminator
		err = EscapeSpecialCharacters(str, length, str, sizeof(str), &length);
	
	Thread Safety: EscapeSpecialCharacters is thread-safe.
*/
int
EscapeSpecialCharacters(const char* input, int inputLength, char* output, int outputBufferSize, int* numCharsOutPtr)
{
	const char* pIn;
	char ch;
	char* pOut;
	int i;
	int outputCharsAvailable;
	int tempOutputBufferAllocated;
	int err;
	
	err = 0;
	*numCharsOutPtr = 0;
	
	tempOutputBufferAllocated = 0;
	if (output == input) {
		output = (char*)NewPtr(outputBufferSize);
		if (output == NULL)
			return NOMEM;	
		tempOutputBufferAllocated = 1;
	}

	pIn = input;
	pOut = output;

	outputCharsAvailable = outputBufferSize;
	for(i=0; i<inputLength; i+=1) {
		ch = *pIn++;
		
		if (outputCharsAvailable <= 0) {
			err = STR_TOO_LONG;
			break;
		}
		
		if (outputCharsAvailable < 2) {
			switch(ch) {
				case '\r':
				case '\n':
				case '\t':
				case '"':
				case '\\':
					err = STR_TOO_LONG;
					break;
			}
			if (err != 0)
				break;
		}

		switch(ch) {
			case '\r':
				*pOut++ = '\\';
				*pOut++ = 'r';
				*numCharsOutPtr += 2;
				outputCharsAvailable -= 2;
				break;

			case '\n':
				*pOut++ = '\\';
				*pOut++ = 'n';
				*numCharsOutPtr += 2;
				outputCharsAvailable -= 2;
				break;

			case '\t':
				*pOut++ = '\\';
				*pOut++ = 't';
				*numCharsOutPtr += 2;
				outputCharsAvailable -= 2;
				break;

			case '\"':
				*pOut++ = '\\';
				*pOut++ = '"';
				*numCharsOutPtr += 2;
				outputCharsAvailable -= 2;
				break;

			case '\\':
				*pOut++ = '\\';
				*pOut++ = '\\';
				*numCharsOutPtr += 2;
				outputCharsAvailable -= 2;
				break;

			default:
				*pOut++ = ch;
				*numCharsOutPtr += 1;
				outputCharsAvailable -= 1;
				break;
		}
	}
	
	if (tempOutputBufferAllocated) {
		memcpy((char*)input, output, *numCharsOutPtr);
		DisposePtr(output);
	}
	
	return err;
}

/*	UnEscapeSpecialCharacters(input, inputLength, output, outputBufferSize, numCharsOutPtr)

	Converts escape sequences for CR, LF, tab, double-quote and backslash to CR, LF, tab,
	double-quote and backslash. You would call this if you receive text that you know is
	escaped. See EscapeSpecialCharacters for a discussion of why text might be escaped.
	
	input is a pointer to a string to be unescaped. Input may be null-terminated or not.
	
	inputLength is the number of characters to be unescaped.

	NOTE: If input is null-terminated and you want output to also be null-terminated,
	then inputLength must include the null terminator character.
	
	output is a pointer to a buffer in which unescaped output text is returned.
	
	outputBufferSize is the total size of the output buffer.
	
	numCharsOutPtr is set to the number of characters stored in output. This includes
	the null-terminator if it was included in inputLength.

	input and output can point to the same memory.
	
	Returns 0 if OK or an error code.
	
	Example:
		char str[32];
		int length;
		// Assume str is set to contain a C string containing escape sequences.
		length = strlen(str) + 1;		// Include null-terminator
		err = UnEscapeSpecialCharacters(str, length, str, sizeof(str), &length);
	
	Thread Safety: UnEscapeSpecialCharacters is thread-safe.
*/
int
UnEscapeSpecialCharacters(const char* input, int inputLength, char* output, int outputBufferSize, int* numCharsOutPtr)
{
	const char* pIn;
	char ch;
	char* pOut;
	int i;
	int outputCharsAvailable;
	int tempOutputBufferAllocated;
	int err;
	
	err = 0;
	*numCharsOutPtr = 0;
	
	tempOutputBufferAllocated = 0;
	if (output == input) {
		output = (char*)NewPtr(outputBufferSize);
		if (output == NULL)
			return NOMEM;	
		tempOutputBufferAllocated = 1;
	}

	pIn = input;
	pOut = output;

	outputCharsAvailable = outputBufferSize;
	for(i=0; i<inputLength; i+=1) {
		if (outputCharsAvailable < 1) {
			err = STR_TOO_LONG;
			break;
		}

		ch = *pIn++;

		if (ch!='\\' || i==(inputLength-1)) {
			// This is not a backslash or it is the last character of the input.
			*pOut++ = ch;
			*numCharsOutPtr += 1;
			outputCharsAvailable -= 1;
		}
		else {
			// We hit a backslash. Check the next character.
			ch = *pIn++;
			i += 1;
			switch(ch) {
				case 'r':
					*pOut++ = '\r';			// Replace \r with CR.
					*numCharsOutPtr += 1;
					outputCharsAvailable -= 1;
					break;

				case 'n':
					*pOut++ = '\n';			// Replace \n with LF.
					*numCharsOutPtr += 1;
					outputCharsAvailable -= 1;
					break;

				case 't':
					*pOut++ = '\t';			// Replace \t with tab.
					*numCharsOutPtr += 1;
					outputCharsAvailable -= 1;
					break;

				case '"':
					*pOut++ = '"';			// Replace \" with double-quote.
					*numCharsOutPtr += 1;
					outputCharsAvailable -= 1;
					break;

				case '\\':
					*pOut++ = '\\';			// Replace \\ with backslash.
					*numCharsOutPtr += 1;
					outputCharsAvailable -= 1;
					break;

				default:				// This is not an sequence that we escape so we copy it to the output unchanged.
					if (outputCharsAvailable < 2) {
						err = STR_TOO_LONG;
						break;
					}
					*pOut++ = '\\';
					*pOut++ = ch;
					*numCharsOutPtr += 2;
					outputCharsAvailable -= 2;
					break;
			}
			if (err != 0)
				break;			// err was set in default case above.
		}
	}
	
	if (tempOutputBufferAllocated) {
		memcpy((char*)input, output, *numCharsOutPtr);
		DisposePtr(output);
	}
	
	return err;
}

/*	MemClear(p, n)

	p points to start of memory to clear.
	n is number of bytes to clear.
	
	Thread Safety: MemClear is thread-safe.
*/
void
MemClear(void *p, BCInt n)
{
	memset(p, 0, n);
}

/*	MoveLockHandle was removed in XOP Toolkit 6. Remove all calls to MoveLockHandle,
	HLock, HGetState, HSetState and HUnlock from your source code. These calls are
	all obsolete.
*/

/*	GetCStringFromHandle(h, str, maxChars)

	h is a handle containing a string.

	str is a C string (null-terminated character array).

	maxChars is the maximum number of characters that str can hold, not including the
	null terminator byte.
	
	GetCStringFromHandle transfers the characters from h to str.
	
	If h is NULL, it returns USING_NULL_STRVAR. This is typically a programmer error.
	
	If the characters in h will not fit in str, it returns STR_TOO_LONG.
	
	If the characters fit, it returns 0.
	
	Thread Safety: GetCStringFromHandle is thread-safe.
*/
int
GetCStringFromHandle(Handle h, char* str, int maxChars)
{
	int numBytesInString;
	int err;
	
	err = 0;
	
	*str = 0;

	if (h == NULL)
		return USING_NULL_STRVAR;

	numBytesInString = (int)GetHandleSize(h);
	if (numBytesInString > maxChars) {
		numBytesInString = maxChars;
		err = STR_TOO_LONG;
	}
	
	memcpy(str, *h, numBytesInString);
	str[numBytesInString] = 0;
	
	return err;
}

/*	PutCStringInHandle(str, h)

	str is a C string (null-terminated character array).

	h is a handle in which the C string data is to be stored.
	
	PutCStringInHandle transfers the characters from str to h. Note that
	the trailing null from the C string is not stored in the handle.
	
	If h is NULL, it returns USING_NULL_STRVAR. This is typically a programmer error.
	
	If an out-of-memory occurs when resizing the handle, it returns NOMEM.
	
	If the operation succeeds, it returns 0.
	
	Thread Safety: PutCStringInHandle is thread-safe.
*/
int
PutCStringInHandle(const char* str, Handle h)
{
	int numBytesInString;

	if (h == NULL)
		return USING_NULL_STRVAR;
		
	numBytesInString = (int)strlen(str);
	SetHandleSize(h, numBytesInString);
	if (MemError())
		return NOMEM;
	
	memcpy(*h, str, numBytesInString);
	
	return 0;
}

/*	CheckAbort(timeoutTicks)

	Returns -1 if user is now pressing cmd-dot (Macintosh) or Ctrl-Break (Windows).
	Returns 1 if TickCount > timeoutTicks.
	Returns 0 otherwise.
	However, if timeoutTicks == 0, does not check for timeout.
	
	Actually does check only every .1 second.
	
	Thread Safety: CheckAbort is thread-safe.
*/
int
CheckAbort(TickCountInt timeOutTicks)
{
	TickCountInt ticks;
	static TickCountInt lastTicks = 0;
	
	ticks = TickCount();
	if (ticks < lastTicks+6)
		return(0);
	lastTicks = ticks;	
	if (timeOutTicks && ticks>timeOutTicks)
		return 1;					// Timeout.

	#ifdef MACIGOR					// Test for cmd-dot.
	{
		KeyMap keys;
		UInt32 theKeys[4];
	
		/*	The KeyMap data type is defined weird on MacIntel which makes it hard to use.
			We copy the data to an array of 4 UInt32s here. The data will be big-endian
			even when running on MacIntel.
			
			If running on MacIntel we convert to little-endian because we are treating
			the data as an array of UInt32s. If we treated it as an array of bytes then
			we would not want to byte swap it.
		*/
		GetKeys(keys);
		memcpy(theKeys, &keys, sizeof(theKeys));
		#ifdef __LITTLE_ENDIAN__
			FixByteOrder(theKeys, 4, 4);
		#endif
	
		// HR, 12/5/93. Added check so that user procedures can use option-cmd-dot for abort signal.
		// HR, 2/22/94. Changed to allow cmd-dot or shift-cmd-dot because of European keyboards.
		if (theKeys[1]==0x00808000 || theKeys[1]==0x00808001) {		// Cmd-dot or shift-cmd-dot ?
			if (theKeys[0]==0 && theKeys[2]==0 && theKeys[3]==0)	// No other keys pressed ?
				return -1;
		}
	}
	#endif	
	#ifdef WINIGOR						// Test for Ctrl-Break.
	{
		/*	Control-break is treated like cmd-dot on Mac.
			Alt-control-break is reserved for user procedures.
			
			HR, 9/25/97: For some reason, GetAsyncKeyState does not always return
			a value with the high bit set if I press and hold Ctrl-Break. Thus,
			I have changed the mask from 0x8000 (key down now) to 0x8001 (key down now
			or key was pressed since last call to GetAsyncKeyState).
		*/
		if ((GetAsyncKeyState(VK_CANCEL) & 0x8001) != 0) {		// This tests for Ctrl-Break.
			if ((GetAsyncKeyState(VK_MENU) & 0x8001) == 0)		// Alt key must not be pressed.
				return -1;
		}
	}
	#endif
	return 0;
}

#define F32_INF 0x7f800000L		// 32 bit infinity

/*	IsNaN32(floatPtr)

	Returns truth that floatPtr points to a NaN.
	
	Thread Safety: IsNaN32 is thread-safe.
*/
int
IsNaN32(float *floatPtr)
{
	return(((*(UInt32*)floatPtr) & 0x7fffffff) > F32_INF);
}

#define F64_INF 0x7ff00000L		// 64 bit infinity

/*	IsNaN64(doublePtr)

	Returns truth that doublePtr points to a NaN.
	
	Thread Safety: IsNaN64 is thread-safe.
*/
int
IsNaN64(double *doublePtr)		// Returns truth that doublePtr points to a NaN.
{
	#ifdef XOP_LITTLE_ENDIAN
		return(((*((UInt32*)doublePtr+1)) & 0x7fffffff) > F64_INF);
	#else
		return(((*(UInt32*)doublePtr) & 0x7fffffff) > F64_INF);
	#endif
}

/*	SetNaN32(fPtr)

	Sets *fPtr to NaN.
	
	Thread Safety: SetNaN32 is thread-safe.
*/
void
SetNaN32(float* fPtr)
{
	*fPtr = SINGLE_NAN;
}

/*	SetNaN32(dPtr)

	Sets *dPtr to NaN.
	
	Thread Safety: SetNaN64 is thread-safe.
*/
void
SetNaN64(double* dPtr)
{
	*dPtr = DOUBLE_NAN;
}

/*	IsINF32(floatPtr)

	Returns truth that floatPtr points to a +/- INF.
	
	Thread Safety: IsINF32 is thread-safe.
*/
int
IsINF32(float *floatPtr)
{
	return ((*(UInt32*)floatPtr) & 0x7FFFFFFF) == F32_INF;
}

/*	IsINF64(doublePtr)

	Returns truth that doublePtr points to a +/- INF.
	
	Thread Safety: IsINF64 is thread-safe.
*/
int
IsINF64(double *doublePtr)		// Returns truth that doublePtr points to a +/- INF.
{
	#ifdef XOP_LITTLE_ENDIAN
		return ((*((UInt32*)doublePtr+1)) & 0x7FFFFFFF) == F64_INF;
	#else
		return ((*(UInt32*)doublePtr) & 0x7FFFFFFF) == F64_INF;
	#endif
}

/*	IgorVersion()

	Returns Igor version number times 100 (e.g. 6.11 is returned as 611).
	
	HR, 091118: Made static in XOP Toolkit 6. Use the igorVersion global instead.
	
	Thread Safety: IgorVersion is not thread-safe.
*/
static int
IgorVersion(void)
{
	#ifdef WINIGOR
		HMODULE igorModule;
		char igorName[MAX_FILENAME_LEN+1];
		char* versionBuffer;
		DWORD versionInfoSize;
		DWORD dwHandle;							// A dummy variable for GetFileVersionInfoSize.
		VS_FIXEDFILEINFO* vsp;
		UINT vsLen;
		int units, tenths, hundredths;
		int version;
		
		igorModule = IgorModule();
		if (igorModule == NULL)					// Should not happen.
			return 0;
		if (GetModuleFileName(igorModule, igorName, MAX_FILENAME_LEN+1) == 0)
			return 0;
		versionInfoSize = GetFileVersionInfoSize(igorName, &dwHandle);
		if (versionInfoSize <= 0)
			return 0;
		versionBuffer = (char*)NewPtr(versionInfoSize);
		if (versionBuffer == NULL)
			return 0;
		if (GetFileVersionInfo(igorName, 0L, versionInfoSize, versionBuffer) == 0) {
			DisposePtr(versionBuffer);
			return 0;
		}
		if (VerQueryValue(versionBuffer, "\\", (void**)&vsp, &vsLen) == 0) {
			DisposePtr(versionBuffer);
			return 0;
		}
	
		units = vsp->dwFileVersionMS >> 16;
		tenths = vsp->dwFileVersionMS & 0xFFFF;
		hundredths = vsp->dwFileVersionLS >> 16;
		version = 100*units + 10*tenths + hundredths;
		
		DisposePtr(versionBuffer);
		return version;
	#endif
	
	#ifdef MACIGOR
		int curResFile;
		Handle vHandle;
		char *vPtr;
		int tens, units, tenths, hundredths;
		int version = 0;
		
		curResFile = CurResFile();
		UseResFile((*(*XOPRecHandle)->stuffHandle)->oldApRefNum);	// Use Igor's resource fork.
		vHandle = Get1Resource('vers', 1);
		UseResFile(curResFile);
		if (vHandle) {							// Pick version out of BCD code.
			vPtr = *vHandle;
			#ifdef __LITTLE_ENDIAN__
				tens = (vPtr[3] & 0xF0) >> 4;
				units = vPtr[3] & 0x0F;
				tenths = (vPtr[2] & 0xF0) >> 4;
				hundredths = vPtr[2] & 0x0F;
			#else
				tens = (vPtr[0] & 0xF0) >> 4;
				units = vPtr[0] & 0x0F;
				tenths = (vPtr[1] & 0xF0) >> 4;
				hundredths = vPtr[1] & 0x0F;
			#endif
			version = 1000*tens + 100*units + 10*tenths + hundredths;
		}
		return(version);
	#endif
}

/*	XOPInit(ioRecHandle)

	Does initialization common to all XOPs. ioRecHandle is the parameter passed by host
	application to XOP.
	
	NOTE: XOPInit must be called before any other call to XOPSupport.
	
	Thread Safety: XOPInit is not thread-safe.
*/
void
XOPInit(IORecHandle ioRecHandle)
{
	XOPRecHandle = ioRecHandle;		// Set global rec handle.
	
	// gMainThreadID is used by RunningInMainThread
	if (gMainThreadID == 0) {
		#ifdef MACIGOR
			gMainThreadID = (PSInt)MPCurrentTaskID();
		#endif
		#ifdef WINIGOR
			gMainThreadID = (PSInt)GetCurrentThreadId();
		#endif
	}

	CheckRunningInMainThread("XOPInit");

	igorVersion = IgorVersion();	// igorVersion global added 10/9/93.
	
	if (igorVersion >= 620) {		// HR, 091118: Igor now supports threadsafe callback method.
		// Switch to the threadsafe callback method.
		if (CallBack1(SET_IGOR_CALLBACK_METHOD, (void*)1) == 0)	// This should always succeed.
			gUseThreadsafeCallbackMethod = 1;					// Callbacks will use the threadsafe method.
	}
}

/*	RunningInMainThread()

	Returns 1 if you are running in the main thread, 0 if you are running in another thread.
	
	If you are not running in the main thread then certain calls are illegal.
	For example, an external operation or function can not do a callback to Igor
	if it is running with an Igor Pro version prior to 6.20. Here is code that enforces this:
	
		if (igorVersion < 620) {
			if (!RunningInMainThread())
				return NOT_IN_THREADSAFE;	// NOT_IN_THREADSAFE is an Igor error code.
		}
	
	Thread Safety: RunningInMainThread is thread-safe.
*/
int
RunningInMainThread(void)
{
	PSInt currentThreadID;
	
	#ifdef MACIGOR
		currentThreadID = (PSInt)MPCurrentTaskID();
	#endif
	#ifdef WINIGOR
		currentThreadID = (PSInt)GetCurrentThreadId();
	#endif
	
	if (currentThreadID == gMainThreadID)
		return 1;
	return 0;
}

/*	CheckRunningInMainThread(routineName)

	CheckRunningInMainThread is called by non-thread-safe XOPSupport routines to
	provide the XOP programmer with feedback if the XOP calls a non-threadsafe
	XOPSupport routine from a thread other than the main thread. In this event,
	it displays an error dialog.
	
	To avoid a cascade of such dialogs, the dialog is presented only once per minute.
	
	The return value is the truth that you are running in the main thread.
	
	You can call CheckRunningInMainThread from your own code to provide feedback
	if your non-threadsafe routine is called from a pre-emptive thread. It is
	typically used as follows:
	
	if (!CheckRunningInMainThread("SomeRoutineName"))
		return NOT_IN_THREADSAFE;		// NOT_IN_THREADSAFE is an Igor error code.
	
	Thread Safety: CheckRunningInMainThread is thread-safe.
*/
int
CheckRunningInMainThread(const char* routineName)
{
	int runningInMainThread = RunningInMainThread();
	static TickCountInt lastMessageTicks = 0;
	
	runningInMainThread = RunningInMainThread();
	if (!runningInMainThread) {
		TickCountInt currentTicks = TickCount();
		if (currentTicks >= lastMessageTicks+3600) {	// Prevent a cascade of dialogs
			char temp[256];
			sprintf(temp, "BUG: %s is not threadsafe."CR_STR, routineName);
			XOPEmergencyAlert(temp);
			lastMessageTicks = currentTicks;
		}
	}
	
	return runningInMainThread;
}

/*	SetXOPType(type)

	Sets XOPType field of XOP's ioRecHandle. This informs host of capabilities and/or mode of
	XOP.
	
	Thread Safety: SetXOPType is not thread-safe.
*/
void
SetXOPType(int type)
{
	CheckRunningInMainThread("SetXOPType");
	(*XOPRecHandle)->XOPType = type;
}

/*	SetXOPEntry(entryPoint)

	Sets XOPEntry field of XOP's ioRecHandle. This informs host of routine to call to pass
	messages to XOP after the INIT message.
	
	Thread Safety: SetXOPEntry is not thread-safe.
*/
void
SetXOPEntry(void (*entryPoint)(void))
{
	CheckRunningInMainThread("SetXOPEntry");
	(*XOPRecHandle)->XOPEntry = entryPoint;
}

/*	SetXOPResult(result)

	Sets the result field of the XOP's ioRecHandle.
	
	Thread Safety: SetXOPResult is not thread-safe.
*/
void
SetXOPResult(XOPIORecResult result)
{
	CheckRunningInMainThread("SetXOPResult");
	(*XOPRecHandle)->result = result;
}

/*	GetXOPResult()

	Returns the result field of the XOP's ioRecHandle.
	
	Thread Safety: GetXOPResult is not thread-safe.
*/
XOPIORecResult
GetXOPResult(void)
{
	CheckRunningInMainThread("GetXOPResult");
	return((*XOPRecHandle)->result);
}

/*	SetXOPMessage(message)

	For use by WaveMetrics only.

	Thread Safety: SetXOPMessage is not thread-safe.
*/
void
SetXOPMessage(int message)
{
	CheckRunningInMainThread("SetXOPMessage");
	(*XOPRecHandle)->message = message;
}

/*	GetXOPMessage()

	Returns the message field of the XOP's ioRecHandle.

	Thread Safety: GetXOPMessage is not thread-safe.
*/
int
GetXOPMessage(void)
{
	CheckRunningInMainThread("GetXOPMessage");
	return((*XOPRecHandle)->message);
}

/*	GetXOPStatus()

	Returns status field of XOP's ioRecHandle.

	Thread Safety: GetXOPStatus is not thread-safe.
*/
int
GetXOPStatus(void)
{
	CheckRunningInMainThread("GetXOPStatus");
	return((*XOPRecHandle)->status);
}

/*	GetXOPItem(itemNumber)

	Returns an item from the XOP's ioRecHandle.
	itemNumber is the number of the item to return starting from zero.

	Thread Safety: GetXOPItem is not thread-safe.
*/
XOPIORecParam
GetXOPItem(int itemNumber)
{
	XOPIORecParam item = 0;
	
	CheckRunningInMainThread("GetXOPItem");

	if ((*XOPRecHandle)->numItems > itemNumber)				// Make sure item exists.
		item = (*XOPRecHandle)->items[itemNumber];
	return(item);
}

/*	SetRecHandle(numItems)

	For use by WaveMetrics only.

	Given a handle to an IORec, SetRecHandle sets the number of items in the IORec to numItems.
	
	This is called by CallBack0, CallBack1, . . . You do not need to call it.

	Thread Safety: SetRecHandle is not thread-safe.
*/
static void
SetRecHandle(int numItems)
{
	/*	HR, 091203: This routine no longer sets the size of XOPRecHandle.
		That is done by Igor when the XOP is created. The maximum number
		of items is NUM_IOREC_ITEMS.
	*/
	(*XOPRecHandle)->numItems = numItems;
}


// *** Callback Routines --- for requesting service from Igor ***

/*	CallBack(message)

	Does a callback to Igor passing it the XOPRecHandle after clearing the result field.
	
	message specifies the requested operation from Igor.
	
	Returns result field from XOPRecHandle after callback.
	
	This routine is not threadsafe and is used only when running with Igor 6.12
	or earlier. With Igor 6.20 or later we use a different threadsafe mechanism.

	For use by WaveMetrics only.

	Thread Safety: CallBack is not thread-safe.
*/
static XOPIORecResult
CallBack(int message)
{
	#ifdef MACIGOR
		pascal void (*CallBack)(void*);			// Address of Igor's callback routine.
	#endif
	#ifdef WINIGOR
		void (*CallBack)(void*);				// Address of Igor's callback routine.
	#endif

	SetXOPResult(0L);							// Start with no error.
	SetXOPMessage(message);
	CallBack = (*XOPRecHandle)->callBackProc;
	(*CallBack)(XOPRecHandle);
	return(GetXOPResult());
}

/*	CallBack0(message)

	Does a callback with the specified message which has 0 parameters.
	
	Returns result field from XOPRecHandle after callback.

	For use by WaveMetrics only.

	Thread Safety: CallBack0 is thread-safe with Igor Pro 6.20 or later.
*/
XOPIORecResult
CallBack0(int message)
{
	if (gUseThreadsafeCallbackMethod) {
		XOPCallRec xcr;
		xcr.version = kXOPCallRecVersion;
		xcr.ioRecHandle = XOPRecHandle;
		xcr.message = message;
		xcr.numParameters = 0;
		xcr.parametersP = NULL;
		xcr.reserved1 = xcr.reserved2 = xcr.reserved3 = xcr.reserved4 = 0;
		(*XOPRecHandle)->callBackProc(&xcr);
		return xcr.result;
	}
	else {
		SetRecHandle(0);
		return(CallBack(message));
	}
}

/*	CallBack1(message, item)

	Does a callback with the specified message passing the item as the only
	parameter in the thingList.
	
	Returns result field from XOPRecHandle after callback.

	For use by WaveMetrics only.

	Thread Safety: CallBack1 is thread-safe with Igor Pro 6.20 or later.
*/
XOPIORecResult
CallBack1(int message, void *item)
{
	if (gUseThreadsafeCallbackMethod) {
		XOPCallRec xcr;
		XOPIORecParam parameters[1];
		parameters[0] = (XOPIORecParam)item;
		xcr.version = kXOPCallRecVersion;
		xcr.ioRecHandle = XOPRecHandle;
		xcr.message = message;
		xcr.numParameters = 1;
		xcr.parametersP = parameters;
		xcr.reserved1 = xcr.reserved2 = xcr.reserved3 = xcr.reserved4 = 0;
		(*XOPRecHandle)->callBackProc(&xcr);
		return xcr.result;
	}
	else {
		SetRecHandle(1);
		(*XOPRecHandle)->items[0] = (XOPIORecParam)item;
		return(CallBack(message));
	}
}

/*	CallBack2(message, item0, item1)

	Does a callback with the specified message passing the item0 and item1 as
	parameters in the thingList.
	
	Returns result field from XOPRecHandle after callback.

	For use by WaveMetrics only.

	Thread Safety: CallBack2 is thread-safe with Igor Pro 6.20 or later.
*/
XOPIORecResult
CallBack2(int message, void *item0, void *item1)
{
	if (gUseThreadsafeCallbackMethod) {
		XOPCallRec xcr;
		XOPIORecParam parameters[2];
		parameters[0] = (XOPIORecParam)item0;
		parameters[1] = (XOPIORecParam)item1;
		xcr.version = kXOPCallRecVersion;
		xcr.ioRecHandle = XOPRecHandle;
		xcr.message = message;
		xcr.numParameters = 2;
		xcr.parametersP = parameters;
		xcr.reserved1 = xcr.reserved2 = xcr.reserved3 = xcr.reserved4 = 0;
		(*XOPRecHandle)->callBackProc(&xcr);
		return xcr.result;
	}
	else {
		SetRecHandle(2);
		(*XOPRecHandle)->items[0] = (XOPIORecParam)item0;
		(*XOPRecHandle)->items[1] = (XOPIORecParam)item1;
		return(CallBack(message));
	}
}

/*	CallBack3(message, item0, item1, item2)

	Does a callback with the specified message passing item0, item1 and item1 as
	parameters in the thingList.
	
	Returns result field from XOPRecHandle after callback.

	For use by WaveMetrics only.

	Thread Safety: CallBack3 is thread-safe with Igor Pro 6.20 or later.
*/
XOPIORecResult
CallBack3(int message, void *item0, void *item1, void *item2)
{
	if (gUseThreadsafeCallbackMethod) {
		XOPCallRec xcr;
		XOPIORecParam parameters[3];
		parameters[0] = (XOPIORecParam)item0;
		parameters[1] = (XOPIORecParam)item1;
		parameters[2] = (XOPIORecParam)item2;
		xcr.version = kXOPCallRecVersion;
		xcr.ioRecHandle = XOPRecHandle;
		xcr.message = message;
		xcr.numParameters = 3;
		xcr.parametersP = parameters;
		xcr.reserved1 = xcr.reserved2 = xcr.reserved3 = xcr.reserved4 = 0;
		(*XOPRecHandle)->callBackProc(&xcr);
		return xcr.result;
	}
	else {
		SetRecHandle(3);
		(*XOPRecHandle)->items[0] = (XOPIORecParam)item0;
		(*XOPRecHandle)->items[1] = (XOPIORecParam)item1;
		(*XOPRecHandle)->items[2] = (XOPIORecParam)item2;
		return(CallBack(message));
	}
}

/*	CallBack4(message, item0, item1, item2, item3)

	Does a callback with the specified message passing item0, item1, item2 and item3 as
	parameters in the thingList.
	
	Returns result field from XOPRecHandle after callback.

	For use by WaveMetrics only.

	Thread Safety: CallBack4 is thread-safe with Igor Pro 6.20 or later.
*/
XOPIORecResult
CallBack4(int message, void *item0, void *item1, void *item2, void *item3)
{
	if (gUseThreadsafeCallbackMethod) {
		XOPCallRec xcr;
		XOPIORecParam parameters[4];
		parameters[0] = (XOPIORecParam)item0;
		parameters[1] = (XOPIORecParam)item1;
		parameters[2] = (XOPIORecParam)item2;
		parameters[3] = (XOPIORecParam)item3;
		xcr.version = kXOPCallRecVersion;
		xcr.ioRecHandle = XOPRecHandle;
		xcr.message = message;
		xcr.numParameters = 4;
		xcr.parametersP = parameters;
		xcr.reserved1 = xcr.reserved2 = xcr.reserved3 = xcr.reserved4 = 0;
		(*XOPRecHandle)->callBackProc(&xcr);
		return xcr.result;
	}
	else {
		SetRecHandle(4);
		(*XOPRecHandle)->items[0] = (XOPIORecParam)item0;
		(*XOPRecHandle)->items[1] = (XOPIORecParam)item1;
		(*XOPRecHandle)->items[2] = (XOPIORecParam)item2;
		(*XOPRecHandle)->items[3] = (XOPIORecParam)item3;
		return(CallBack(message));
	}
}

/*	CallBack5(message, item0, item1, item2, item3, item4)

	Does a callback with the specified message passing item0 through item4 as
	parameters in the thingList.
	
	Returns result field from XOPRecHandle after callback.

	For use by WaveMetrics only.

	Thread Safety: CallBack5 is thread-safe with Igor Pro 6.20 or later.
*/
XOPIORecResult
CallBack5(int message, void *item0, void *item1, void *item2, void *item3, void *item4)
{
	if (gUseThreadsafeCallbackMethod) {
		XOPCallRec xcr;
		XOPIORecParam parameters[5];
		parameters[0] = (XOPIORecParam)item0;
		parameters[1] = (XOPIORecParam)item1;
		parameters[2] = (XOPIORecParam)item2;
		parameters[3] = (XOPIORecParam)item3;
		parameters[4] = (XOPIORecParam)item4;
		xcr.version = kXOPCallRecVersion;
		xcr.ioRecHandle = XOPRecHandle;
		xcr.message = message;
		xcr.numParameters = 5;
		xcr.parametersP = parameters;
		xcr.reserved1 = xcr.reserved2 = xcr.reserved3 = xcr.reserved4 = 0;
		(*XOPRecHandle)->callBackProc(&xcr);
		return xcr.result;
	}
	else {
		SetRecHandle(5);
		(*XOPRecHandle)->items[0] = (XOPIORecParam)item0;
		(*XOPRecHandle)->items[1] = (XOPIORecParam)item1;
		(*XOPRecHandle)->items[2] = (XOPIORecParam)item2;
		(*XOPRecHandle)->items[3] = (XOPIORecParam)item3;
		(*XOPRecHandle)->items[4] = (XOPIORecParam)item4;
		return(CallBack(message));
	}
}

/*	CallBack6(message, item0, item1, item2, item3, item4, item5)

	Does a callback with the specified message passing item0 through item5 as
	parameters in the thingList.
	
	Returns result field from XOPRecHandle after callback.

	For use by WaveMetrics only.

	Thread Safety: CallBack6 is thread-safe with Igor Pro 6.20 or later.
*/
XOPIORecResult
CallBack6(int message, void *item0, void *item1, void *item2, void *item3, void *item4, void *item5)
{
	if (gUseThreadsafeCallbackMethod) {
		XOPCallRec xcr;
		XOPIORecParam parameters[6];
		parameters[0] = (XOPIORecParam)item0;
		parameters[1] = (XOPIORecParam)item1;
		parameters[2] = (XOPIORecParam)item2;
		parameters[3] = (XOPIORecParam)item3;
		parameters[4] = (XOPIORecParam)item4;
		parameters[5] = (XOPIORecParam)item5;
		xcr.version = kXOPCallRecVersion;
		xcr.ioRecHandle = XOPRecHandle;
		xcr.message = message;
		xcr.numParameters = 6;
		xcr.parametersP = parameters;
		xcr.reserved1 = xcr.reserved2 = xcr.reserved3 = xcr.reserved4 = 0;
		(*XOPRecHandle)->callBackProc(&xcr);
		return xcr.result;
	}
	else {
		SetRecHandle(6);
		(*XOPRecHandle)->items[0] = (XOPIORecParam)item0;
		(*XOPRecHandle)->items[1] = (XOPIORecParam)item1;
		(*XOPRecHandle)->items[2] = (XOPIORecParam)item2;
		(*XOPRecHandle)->items[3] = (XOPIORecParam)item3;
		(*XOPRecHandle)->items[4] = (XOPIORecParam)item4;
		(*XOPRecHandle)->items[5] = (XOPIORecParam)item5;
		return(CallBack(message));
	}
}

/*	CallBack8(message, item0, item1, item2, item3, item4, item5, item6, item7)

	Does a callback with the specified message passing item0 through item7 as
	parameters in the thingList.
	
	Returns result field from XOPRecHandle after callback.

	For use by WaveMetrics only.

	Thread Safety: CallBack8 is thread-safe with Igor Pro 6.20 or later.
*/
XOPIORecResult
CallBack8(int message, void *item0, void *item1, void *item2, void *item3, void *item4, void *item5, void *item6, void *item7)
{
	if (gUseThreadsafeCallbackMethod) {
		XOPCallRec xcr;
		XOPIORecParam parameters[8];
		parameters[0] = (XOPIORecParam)item0;
		parameters[1] = (XOPIORecParam)item1;
		parameters[2] = (XOPIORecParam)item2;
		parameters[3] = (XOPIORecParam)item3;
		parameters[4] = (XOPIORecParam)item4;
		parameters[5] = (XOPIORecParam)item5;
		parameters[6] = (XOPIORecParam)item6;
		parameters[7] = (XOPIORecParam)item7;
		xcr.version = kXOPCallRecVersion;
		xcr.ioRecHandle = XOPRecHandle;
		xcr.message = message;
		xcr.numParameters = 8;
		xcr.parametersP = parameters;
		xcr.reserved1 = xcr.reserved2 = xcr.reserved3 = xcr.reserved4 = 0;
		(*XOPRecHandle)->callBackProc(&xcr);
		return xcr.result;
	}
	else {
		SetRecHandle(8);
		(*XOPRecHandle)->items[0] = (XOPIORecParam)item0;
		(*XOPRecHandle)->items[1] = (XOPIORecParam)item1;
		(*XOPRecHandle)->items[2] = (XOPIORecParam)item2;
		(*XOPRecHandle)->items[3] = (XOPIORecParam)item3;
		(*XOPRecHandle)->items[4] = (XOPIORecParam)item4;
		(*XOPRecHandle)->items[5] = (XOPIORecParam)item5;
		(*XOPRecHandle)->items[6] = (XOPIORecParam)item6;
		(*XOPRecHandle)->items[7] = (XOPIORecParam)item7;
		return(CallBack(message));
	}
}

/*	IgorError(title, errCode)

	Displays an error alert appropriate for the specified error code.
	
	Title is a short string that identifies what generated the error.
	
	errCode may be an Igor error code (defined in IgorXOP.h), an XOP-defined error code,
	or, when running on Macintosh, a Mac OS error code. To display a message for a Windows
	OS error, convert the code to an Igor code by calling WindowsErrorToIgorError.

	Thread Safety: IgorError is not thread-safe.
*/
void
IgorError(const char *title, int errCode)
{
	if (!CheckRunningInMainThread("IgorError"))
		return;
	CallBack2(IGORERROR, (void*)title, (void *)errCode);
}

/*	GetIgorErrorMessage(errCode, errorMessage)

	Returns via errorMessage the message corresponding to the specified error code.
	
	errCode may be an Igor error code (defined in IgorXOP.h), an XOP-defined error code,
	or, when running on Macintosh, a Mac OS error code. To display a message for a Windows
	OS error, convert the code to an Igor code by calling WindowsErrorToIgorError.

	Do not pass 0 for errCode. There is no error message corresponding to 0.

	The function result is 0 if OK, IGOR_OBSOLETE if the current version of Igor
	is earlier than 3.14, or another non-zero error code if the errCode parameter
	is invalid. If GetIgorErrorMessage fails to get a message, it sets *errorMessage to 0.

	Thread Safety: GetIgorErrorMessage is not thread-safe.
*/
int
GetIgorErrorMessage(int errCode, char errorMessage[256])
{
	int err;

	if (!CheckRunningInMainThread("GetIgorErrorMessage")) {
		*errorMessage = 0;
		return NOT_IN_THREADSAFE;
	}
	
	err = (int)CallBack2(GET_IGOR_ERROR_MESSAGE, (void*)errCode, errorMessage);
	if (err != 0)
		*errorMessage = 0;
	return err;
}

/*	WinInfo(index, typeMask, name, windowRefPtr)

	Returns information about an Igor target window.

	Thread Safety: WinInfo is not thread-safe.
*/
int
WinInfo(int index, int typeMask, char *name, XOP_WINDOW_REF* windowRefPtr)
{
	if (!CheckRunningInMainThread("WinInfo"))
		return 0;	// Return value is an Igor window type.

	return (int)CallBack4(WININFO, (void *)index, (void *)typeMask, name, windowRefPtr);
}


// *** Notice Routines -- for displaying messages in the history area ***

/*	XOPNotice(noticePtr)

	XOPNotice does a callback to Igor to get the notice identified by noticePtr displayed.
	
	noticePtr is a C string (null byte at the end)

	Typically notices should be terminated by carriage return (use CR_STR at end of string).
	
	Thread Safety: XOPNotice is thread-safe with Igor Pro 6.20 or later.
*/
void
XOPNotice(const char *noticePtr)
{
	CallBack1(NOTICE, (void*)noticePtr);
}

/*	XOPNotice2(noticePtr, options)

	XOPNotice2 does a callback to Igor to get the notice identified by noticePtr displayed.
	
	options is defined as follows:
		Bit 0:
			If cleared, the modification state of the history area is not modified. Use this
			if you want to display a notice but don't want that to cause the current experiment
			to be modified.
			
			If set, the history area is marked as modified after the notice is displayed.
		
	All other bits are reserved and must be passed as zero.
	
	noticePtr is a C string (null byte at the end)

	Typically notices should be terminated by carriage return (use CR_STR at end of string).
	
	Thread Safety: XOPNotice2 is thread-safe with Igor Pro 6.20 or later.
*/
void
XOPNotice2(const char *noticePtr, UInt32 options)
{
	CallBack2(NOTICE2, (void*)noticePtr, (void*)options);
}

/*	XOPNotice3(noticePtr, rulerName, options)

	XOPNotice3 does a callback to Igor to display text in the history area. Through the rulerName
	and options parameters, XOPNotice3 provides control over the formatting of the text in the History
	Carbon Copy notebook. If you are not using the History Carbon Copy feature (added in Igor Pro
	6.10B04), use XOPNotice instead of XOPNotice3.
		
	noticePtr is a C string (null byte at the end).
	Typically notices should be terminated by carriage return (use CR_STR at end of string as shown below).
	
	rulerName is the name of a ruler in the History Carbon Copy notebook.
	XOPNotice3 applies the ruler to the current selection in the History Carbon Copy notebook
	before inserting the text. If rulerName is "" then the ruler is not changed.
	
	options is defined as follows:
		Bit 0:
			If set, the history carbon copy ruler is restored after the text is inserted.
			If cleared, the history carbon copy ruler is not restored.
		
		All other bits are reserved and must be passed as zero.
	
	Typically you would call this to use a ruler that you previously created in the history
	carbon copy notebook to control formatting for specific kinds of notices. For example, if
	you created a ruler named Warning, you would send a warning to the history like this:
	
		XOPNotice3("You have been warned!"CR_STR, "Warning", 1);
		
	Bit 0 of options should be set except in the case when you expect to make many calls to
	XOPNotice3 in a row. In this case, you can make it run slightly faster like this:
	
		XOPNotice3("Starting", "Progress", 0);		// Starting a long process
		XOPNotice3(" .", "", 0);					// Add a dot to indicate progress
		. . .
		XOPNotice3(" .", "", 0);					// Add a dot to indicate progress
		XOPNotice3("Done!"CR_STR, "", 0);			// The long process is done
		XOPNotice3("", "Normal", 0);				// Restore normal ruler
		
	In this example, we set the history carbon copy ruler to "Progress" which is the name
	of a ruler that we created for paragraphs containing progress messages. We then called
	XOPNotice3 repeatedly to add a dot to the progress message. In these calls, we specified
	no ruler so the ruler is not changed. At the end, we add a carriage return (CR_STR) to
	make sure we are on a new paragraph and then set the ruler to Normal so that normal
	history carbon copy processing will occur in the future.
	
	Added in Igor Pro 6.10B04. If you call this when running with an earlier version,
	nothing will be sent to the history.
	
	Thread Safety: XOPNotice3 is thread-safe with Igor Pro 6.20 or later.
*/
void
XOPNotice3(const char *noticePtr, const char* rulerName, UInt32 options)
{
	CallBack3(NOTICE3, (void*)noticePtr, (void*)rulerName, (void*)options);
}

/*	XOPResNotice(strListID, index)

	XOPResNotice does a callback to Igor to get a notice displayed. The notice is fetched
	from a resource which must be in the resource fork of the XOP.
	
	The resource must be of type 'STR#'. The resource ID should be between 1100 and 1199.
	These resource IDs are reserved for XOPs.
	
	strListID is the resource ID of the STR# containing the string.
	index is the number of the string in the STR# resource.
	
	The strings in the STR# resource must be 255 characters or less.
	
	Thread Safety: XOPResNotice is not thread-safe.
*/
void
XOPResNotice(int strListID, int index)
{
	char theString[256];
	
	if (!CheckRunningInMainThread("XOPResNotice"))	// GetXOPIndString is not threadsafe
		return;

	GetXOPIndString(theString, strListID, index);
	XOPNotice(theString);
}


// *** Command Routines -- for executing Igor commands from an XOP ***

/*	XOPCommand(cmdPtr)

	XOPCommand does a callback to Igor to execute the command identified by cmdPtr.
	The command appears in the command line while Igor is executing it.
	
	cmdPtr is a C string (null byte at the end)

	XOPCommand returns the result from the command execution.
	
	Thread Safety: XOPCommand is not thread-safe.
*/
int
XOPCommand(const char *cmdPtr)
{
	if (!CheckRunningInMainThread("XOPCommand"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack1(COMMAND, (void*)cmdPtr);
}

/*	XOPSilentCommand(cmdPtr)

	XOPSilentCommand does a callback to Igor to execute the command identified
	by cmdPtr. The command does not appear in the command line while Igor is
	executing it.
	
	cmdPtr is a C string (null byte at the end)

	XOPSilentCommand returns the result from the command execution.
	
	Thread Safety: XOPSilentCommand is not thread-safe.
*/
int
XOPSilentCommand(const char *cmdPtr)
{
	if (!CheckRunningInMainThread("XOPSilentCommand"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack1(SILENT_COMMAND, (void*)cmdPtr);
}

/*	XOPCommand2(cmdPtr, silent, sendToHistory)

	XOPCommand2 does a callback to Igor to execute the command identified by cmdPtr.
	
	cmdPtr is a C string (null byte at the end). The string must consist of one
	line of text not longer than MAXCMDLEN and with no carriage return characters.
	
	If silent is non-zero, the command appears in the command line while Igor is
	executing it.
	
	If sendToHistory is non-zero and if the result from command execution is 0 (no error),
	Igor appends the command to the history with a bullet character before it.

	XOPCommand2 returns the result from the command execution.
	
	Thread Safety: XOPCommand2 is not thread-safe.
*/
int
XOPCommand2(const char *cmdPtr, int silent, int sendToHistory)
{
	if (!CheckRunningInMainThread("XOPCommand2"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(COMMAND2, (void*)cmdPtr, (void*)silent, (void*)sendToHistory);
}

/*	XOPCommand3(cmdPtr, silent, sendToHistory, historyTextHPtr)

	XOPCommand3 does a callback to Igor to execute the command identified by cmdPtr.
	
	cmdPtr is a C string (null byte at the end). The string must consist of one
	line of text not longer than MAXCMDLEN and with no carriage return characters.
	
	If silent is non-zero, the command appears in the command line while Igor is
	executing it.
	
	If sendToHistory is non-zero and if the result from command execution is 0 (no error),
	Igor appends the command to the history with a bullet character before it.
	
	XOPCommand3 returns via *historyTextHPtr any result text inserted into the history
	area as a result of the command. In the event of an error, *historyTextHPtr is set
	to NULL. If no error occurs, *historyTextHPtr is set to a new handle which you own
	and must dispose using DisposeHandle when you no longer need it. *historyTextHPtr
	is not null terminated.

	XOPCommand3 returns the result from the command execution.
	
	Thread Safety: XOPCommand3 is not thread-safe.
*/
int
XOPCommand3(const char *cmdPtr, int silent, int sendToHistory, Handle* historyTextHPtr)
{
	if (!CheckRunningInMainThread("XOPCommand3"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(COMMAND3, (void*)cmdPtr, (void*)silent, (void*)sendToHistory, historyTextHPtr);
}

/*	PutCmdLine(cmd, mode)

	PutCmdLine puts the specified text into Igor's command line using the specified mode.
	See IgorXOP.h for list and description of modes.
	
	Thread Safety: PutCmdLine is not thread-safe.
*/
void
PutCmdLine(const char *cmd, int mode)
{
	if (!CheckRunningInMainThread("PutCmdLine"))
		return;
	CallBack2(PUTCMDLINE, (void*)cmd, (void *)mode);
}

/*	FinishDialogCmd(char *cmd, int mode)

	Called at the end of an Igor style dialog.
	cmd is a C string containing the command generated by the Igor style dialog.
	If mode is 1, puts the command in Igor's command line and starts execution.
	If mode is 2, puts the command in Igor's command line but does not start execution.
	If mode is 3, puts the command in the clipboard.
	
	Thread Safety: FinishDialogCmd is not thread-safe.
*/
void
FinishDialogCmd(const char *cmd, int mode)
{
	if (!CheckRunningInMainThread("FinishDialogCmd"))
		return;

	switch (mode) {
		case 1:
			PutCmdLine(cmd, FIRSTCMDCRHIT);
			break;
		case 2:
			PutCmdLine(cmd, INSERTCMD);
			break;
		case 3:
			#ifdef MACIGOR
			{
				ScrapRef scrap;
				int err;
				
				if (err = ClearCurrentScrap())
					return;
				if (err = GetCurrentScrap(&scrap))
					return;
				if (err = PutScrapFlavor(scrap, 'TEXT', kScrapFlavorMaskNone, strlen(cmd), cmd))
					return;
			}
			#endif
			#ifdef WINIGOR
				SetScrapData(NULL, (int)strlen(cmd), 'TEXT', (char*)cmd);	// This routine is in IGOR.
			#endif
			break;
	}
}

//  *** Variable Access Routines ***

/*	FetchNumVar(varName, doublePtr1, doublePtr2)

	FetchNumVar returns the value of a named variable via the double pointers doublePtr1
	and doublePtr2. The real part is returned via doublePtr1 and the imaginary part via
	doublePtr2. If the numeric variable is not complex then the *doublePtr2 is meaningless.
	However, doublePtr2 is required anyway.
	
	It returns -1 if the variable does not exist or the numeric type of the variable if it does.
	
	Thread Safety: FetchNumVar is thread-safe with Igor Pro 6.20 or later.
*/
int
FetchNumVar(const char *varName, double *doublePtr1, double *doublePtr2)
{
	return (int)CallBack3(FETCHNUMVAR, (void*)varName, doublePtr1, doublePtr2);
}

/*	StoreNumVar(varName, doublePtr1, doublePtr2)

	NOTE: We recommend that you use SetIgorIntVar, SetIgorFloatingVar or SetIgorComplexVar
	instead of StoreNumVar as they are easier to use.

	StoreNumVar stores the value in double1 and double2 into a variable whose name is specified
	by varName. double1 contains the real part and double2 contains the imaginary part.
	double2 is required even if the variable is not complex.
	
	It returns -1 if the variable does not exist or the numeric type of the variable if it does.
	
	This routine was changed to take double pointer arguments for release 2 of the XOP Toolkit.
	See Igor XOP Tookit Tech Note #1 for details.
	
	Thread Safety: StoreNumVar is thread-safe with Igor Pro 6.20 or later.
*/
int
StoreNumVar(const char *varName, const double *doublePtr1, const double *doublePtr2)
{
	return (int)CallBack3(STORENUMVAR, (void*)varName, (void*)doublePtr1, (void*)doublePtr2);
}

/*	FetchStrVar(varName, stringPtr)

	FetchStrVar returns the value of a named string variable via the pointer stringPtr.
	It returns 0 if it was able to fetch the string or an error code if it was not able.
	
	stringPtr should be big enough to hold a 255 byte string.
	
	Thread Safety: FetchStrVar is thread-safe with Igor Pro 6.20 or later.
*/
int
FetchStrVar(const char *varName, char *stringPtr)
{
	return (int)CallBack2(FETCHSTRVAR, (void*)varName, stringPtr);
}

/*	FetchStrHandle(varName)

	FetchStrHandle returns the handle containing the text for the named string
	variable or NULL if no such string variable exists.
	
	The text is not null terminated. Use GetHandleSize to determine the
	number of characters in the string.
	
	You should not dispose of or otherwise modify this handle since it belongs
	to Igor.
	
	Thread Safety: FetchStrHandle is thread-safe with Igor Pro 6.20 or later.
*/
Handle
FetchStrHandle(const char *varName)
{
	return (Handle)CallBack1(FETCHSTRHANDLE, (void*)varName);
}

/*	StoreStrVar(varName, stringPtr)

	NOTE: We recommend that you use SetIgorStringVar instead of StoreStrVar
	as it is easier to use.

	StoreStrVar stores the value in stringPtr into a string variable whose name is specified by
	varName.

	It returns 0 if it was able to store the string or an error code if it was not able.
	
	Thread Safety: StoreStrVar is thread-safe with Igor Pro 6.20 or later.
*/
int
StoreStrVar(const char *varName, const char *stringPtr)
{
	return (int)CallBack2(STORESTRVAR, (void*)varName, (void*)stringPtr);
}

/*	Variable(varName, varType)

	NOTE: We recommend that you use SetIgorIntVar, SetIgorFloatingVar, SetIgorComplexVar
	or SetIgorStringVar instead of Variable as they are easier to use.

	Variable creates an Igor variable with the specified name and type.
	varType is
		0 for a string variable

		NT_FP64 for a double precision floating point variable
		
		(NT_FP64 | NT_CMPLX) for a complex double precision floating point variable

		Use VAR_GLOBAL in type to force the variable or string to be global.
		Example: (NT_FP64 | VAR_GLOBAL)
		
		If you don't use VAR_GLOBAL then the variable or string will be
		global if executed from the command line or local if executed from
		a macro.
	
	It returns 0 if it was able to make the variable or an error code
	if it was not able.
	
	Thread Safety: Variable is thread-safe with Igor Pro 6.20 or later.
*/
int
Variable(const char *varName, int varType)
{
	return (int)CallBack2(VARIABLE, (void*)varName, (void *)varType);
}

/*	VariableList(listHandle, match, sep, varTypeCode)

	Thread Safety: VariableList is not thread-safe.
*/
int
VariableList(Handle listHandle, const char *match, const char *sep, int varTypeCode)
{
	if (!CheckRunningInMainThread("VariableList"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(VARIABLELIST, listHandle, (void*)match, (void*)sep, (void*)varTypeCode);
}

/*	StringList(listHandle, match, sep, varTypeCode)

	Thread Safety: StringList is not thread-safe.
*/
int
StringList(Handle listHandle, const char *match, const char *sep)
{
	if (!CheckRunningInMainThread("StringList"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(STRINGLIST, listHandle, (void*)match, (void*)sep);
}

/*	SetIgorIntVar(numVarName, value, forceGlobal)

	Creates an Igor numeric variable if it does not already exist. Sets the variable to an integer value.

	When called from the command line or from a user-defined function, SetIgorIntVar always creates
	and/or sets a global variable in the current data folder.

	When called from a macro, SetIgorIntVar creates and/or sets a global variable in the
	current data folder if forceGlobal is non-zero and creates and/or sets a macro-local
	variable if forceGlobal is 0. To get consistent behavior regardless of how you are
	called, pass 1 for forceGlobal.	

	Returns 0 or error code.
	
	Thread Safety: SetIgorIntVar is thread-safe with Igor Pro 6.20 or later.
*/
int
SetIgorIntVar(const char* numVarName, int value, int forceGlobal)
{
	int varType;
	double d1, d2;
	int result = 0;

	d1 = value;
	d2 = 0.0;
	if (StoreNumVar(numVarName, &d1, &d2) == -1) {			// Does it exist ?
		varType = NT_FP64;									// HR, 10/22/95: Was NT_FP32.
		if (forceGlobal)
			varType |= VAR_GLOBAL;
		if (result = Variable(numVarName, varType))			// Create new variable.
			return result;
		StoreNumVar(numVarName, &d1, &d2);
	}
	
	return result;
}

/*	SetIgorFloatingVar(numVarName, valuePtr, forceGlobal)

	Creates an Igor numeric variable if it does not already exist. Sets the variable to
	a floating-point value.

	When called from the command line or from a user-defined function, SetIgorFloatingVar
	always creates and/or sets a global variable in the current data folder.

	When called from a macro, SetIgorFloatingVar creates and/or sets a global variable in the
	current data folder if forceGlobal is non-zero and creates and/or sets a macro-local
	variable if forceGlobal is 0. To get consistent behavior regardless of how you are
	called, pass 1 for forceGlobal.	

	Returns 0 or error code.
	
	Thread Safety: SetIgorFloatingVar is thread-safe with Igor Pro 6.20 or later.
*/
int
SetIgorFloatingVar(const char* numVarName, const double* valuePtr, int forceGlobal)
{
	int varType;
	double d1, d2;
	int result = 0;

	d1 = *valuePtr;
	d2 = 0.0;
	if (StoreNumVar(numVarName, &d1, &d2) == -1) {			// Does it exist ?
		varType = NT_FP64;
		if (forceGlobal)
			varType |= VAR_GLOBAL;
		if (result = Variable(numVarName, varType))			// Create new variable.
			return result;
		StoreNumVar(numVarName, &d1, &d2);
	}
	
	return result;
}

/*	SetIgorComplexVar(numVarName, realValuePtr, imagValuePtr, int forceGlobal)

	Creates an Igor complex numeric variable if it does not already exist.
	Sets the variable to a complex floating-point value.

	When called from the command line or from a user-defined function, SetIgorComplexVar
	always creates and/or sets a global variable in the current data folder.

	When called from a macro, SetIgorComplexVar creates and/or sets a global variable in the
	current data folder if forceGlobal is non-zero and creates and/or sets a macro-local
	variable if forceGlobal is 0. To get consistent behavior regardless of how you are
	called, pass 1 for forceGlobal.	
	
	Returns 0 or error code.
	
	Thread Safety: SetIgorComplexVar is thread-safe with Igor Pro 6.20 or later.
*/
int
SetIgorComplexVar(const char* numVarName, const double* realValuePtr, const double* imagValuePtr, int forceGlobal)
{
	int varType;
	int result = 0;

	if (StoreNumVar(numVarName, realValuePtr, imagValuePtr) == -1) {	// Does it exist ?
		varType = NT_FP64 | NT_CMPLX;
		if (forceGlobal)
			varType |= VAR_GLOBAL;
		if (result = Variable(numVarName, varType))						// Create new variable.
			return result;
		StoreNumVar(numVarName, realValuePtr, imagValuePtr);
	}
	
	return result;
}

/*	SetIgorStringVar(stringVarName, stringVarValue, forceGlobal)

	Creates an Igor string variable if it does not already exist. Stores the specified string in the variable.

	When called from the command line or from a user-defined function, SetIgorStringVar always creates
	and/or sets a global variable in the current data folder.

	When called from a macro, SetIgorStringVar creates and/or sets a global variable in the
	current data folder if forceGlobal is non-zero and creates and/or sets a macro-local
	variable if forceGlobal is 0. To get consistent behavior regardless of how you are
	called, pass 1 for forceGlobal.

	stringVarValue is a null-terminated C string of any length.	

	Returns 0 or error code.
	
	Thread Safety: SetIgorStringVar is thread-safe with Igor Pro 6.20 or later.
*/
int
SetIgorStringVar(const char* stringVarName, const char* stringVarValue, int forceGlobal)
{
	int varType;
	int result;

	if (result = StoreStrVar(stringVarName, stringVarValue)) {	// Error if string does not exist.
		varType = 0;											// 0 means string variable.
		if (forceGlobal)
			varType |= VAR_GLOBAL;
		if (result = Variable(stringVarName, varType))
			return result;
		result = StoreStrVar(stringVarName, stringVarValue);
	}
	return result;
}


// *** Command Parsing Routines ***
/*	These routines are obsolete and were removed from XOP Toolkit 6.
	See "XOP Toolkit 6 Upgrade Notes" in the XOP Toolkit 6 manual for details.

	AtEndOfCommand
	CheckTerm				FileLoaderGetOperationFlags2
	GetAString				GetAStringInHandle
	GetDataFolder			GetDataFolderAndName
	GetFlag					GetFlagNum
	GetFormat				GetKeyword
	GetLong					GetName
	GetNum					GetNum2
	GetNumVarName			GetStrVarName
	GetSymb					GetTrueOrFalseFlag
	GetWave					GetWaveList
	GetWaveName				GetWaveRange
	IsStringExpression		Keyword
	NextSymb
*/

// *** Name Utility Routines ***

/*	UniqueName(baseName, finalName)

	Given a base name (like "wave") UniqueName returns a name (like "wave3")
	via finalName that does not conflict with any existing names.
	
	Returns the number used to make the name unique.
	
	Thread Safety: UniqueName is thread-safe with Igor Pro 6.20 or later.
*/
int
UniqueName(const char *baseName, char *finalName)
{
	return (int)CallBack2(UNIQUENAME, (void*)baseName, finalName);
}

/*	UniqueName2(nameSpaceCode, baseName, finalName, suffixNumPtr)

	Given a base name (like "wave") UniqueName2 returns a name (like "wave3")
	via finalName that does not conflict with any existing names. The number
	appended to make the name unique will be *suffixNumPtr or greater.
	Igor sets *suffixNumPtr to the number Igor used to make the name unique.
	
	nameSpaceCode is:
		MAIN_NAME_SPACE			for Igor's main name space (waves, variables, windows)
		DATAFOLDER_NAME_SPACE
		See IgorXOP.h for other less frequently-used name space codes.
	
	Returns 0 if OK, -1 for a bad nameSpaceCode or some other non-zero error code.
	
	You should use suffixNumPtr as follows:
		int suffixNum = 0;
		for(i = 0; i < numObjectsToBeCreated; i++) {
			if (err = UniqueName2(nameSpaceCode, baseName, finalName, &suffixNum))
				break;
			MakeObject(finalName);
		}
	
	Thread Safety: UniqueName2 is thread-safe with Igor Pro 6.20 or later.
*/
int
UniqueName2(int nameSpaceCode, const char *baseName, char *finalName, int* suffixNumPtr)
{
	return (int)CallBack4(UNIQUENAME2, (void*)nameSpaceCode, (void*)baseName, finalName, suffixNumPtr);
}

/*	SanitizeWaveName(waveName, column)

	NOTE: If your XOP requires IGOR Pro 3 or later, you should use CleanupName
		  instead of the older SanitizeWaveName.

	Given a pointer to a C string containing a proposed wave name,
	SanitizeWaveName() changes it to make it a valid wave name if necessary.
	
	Returns 1 if it had to make a change, 0 if name was OK to begin with.
	
	First, it chops the string off if it is too long.
	Then, it makes sure that the first character is alphabetic.
	Then it replaces any subsequent characters that are not alphanumeric with underscore.
	
	Thread Safety: SanitizeWaveName is thread-safe with Igor Pro 6.20 or later.
*/
int
SanitizeWaveName(char *waveName, int column)
{
	int len, i, ch;
	int result = 0;
	
	len = (int)strlen(waveName);
	if (len==0) {
		sprintf(waveName, "NoName%d", column);
		return 1;
	}

	ch = waveName[0];
	if (!(ch>0 && isalpha(ch))) {					// First char must be alphabetic.
		memmove(waveName+1, waveName, len+1);		// Shift chars over.
		waveName[0] = 'X';
		len += 1;
		result = 1;
	}
	
	if (len > MAX_WAVE_NAME) {
		waveName[MAX_WAVE_NAME] = '\0';				// Truncate name.
		len = MAX_WAVE_NAME;
		result = 1;
	}
	
	for (i = 1; i < len; i++) {						// Subsequent characters must be.
		ch = waveName[i];							// Alphanumeric or underscore.
		if (!(ch>0 && isalnum(ch))) {
			waveName[i] = '_';
			result = 1;
		}
	}
	
	return result;
}

/*	CheckName(dataFolderH, objectType, name)

	Checks the name for legality and uniqueness.
	
	If dataFolderH is NULL, it looks for conflicts with objects in the current data folder.
	If it is not NULL, it looks for conflicts in the folder specified by dataFolderH.
	
	objectType is one of the following which are defined in IgorXOP.h:
		WAVE_OBJECT
		VAR_OBJECT					(numeric variable)
		STR_OBJECT					(string variable)
		GRAPH_OBJECT
		TABLE_OBJECT
		LAYOUT_OBJECT
		PANEL_OBJECT
		NOTEBOOK_OBJECT
		DATAFOLDER_OBJECT
		PATH_OBJECT
		PICT_OBJECT

	Returns 0 if the name is legal and is not in conflict with an existing object.
	Returns an Igor error code otherwise.
	
	Thread Safety: CheckName is thread-safe with Igor Pro 6.20 or later.
*/
int
CheckName(DataFolderHandle dataFolderH, int objectType, const char* name)
{
	return (int)CallBack3(CHECKNAME, dataFolderH, (void*)objectType, (void*)name);
}

/*	PossiblyQuoteName(name)
	
	name contains an Igor object name.
	PossiblyQuoteName puts single quotes around the name if they would be
	needed to use the name in Igor's command line.

	Igor Pro allows wave and data folder names to contain characters,
	such as space and dot, that were previously illegal in names. We call this
	"liberal" name rules.
	
	If an object has such a name, you must single-quote the name to use it in Igor's
	command line. This includes using it in the Execute operation or in the XOPCommand,
	XOPSilentCommand, or FinishDialogCmd XOPSupport routines. Thus, if you are going
	to use a wave or data folder name for this purpose, you should call PossiblyQuoteName
	to add the quotes if needed.

	Returns true if the name was quoted, false otherwise.
	
	NOTE: name must be able to hold two additional characters. Thus, you
		  should declare name: char name[MAX_OBJ_NAME+2+1];
		  
	NOTE: Liberal rules are still not allowed for string and numeric variable names.
	
	Thread Safety: PossiblyQuoteName is thread-safe with Igor Pro 6.20 or later.
*/
int
PossiblyQuoteName(char *name)
{
	return (int)CallBack1(POSSIBLY_QUOTE_NAME, name);
}

/*	CatPossiblyQuotedName(char* str, char* name)

	Adds the specified Igor object name to the end of the string.
	If necessary, puts single quotes around the name so that it can be
	used in the Igor command line.
	
	Use this to concatenate a wave name to the end of a command string
	when the wave name may be a liberal name that needs to be quoted to
	be used in the command line. 
	
	Example:
		char waveName[MAX_OBJ_NAME+1];				// This contains a wave name.
		char cmd[256];
		strcpy(cmd, "Redimension/N=1000 ");
		CatPossiblyQuotedName(cmd, waveName);
		XOPSilentCommand(cmd);
	
	Thread Safety: CatPossiblyQuotedName is thread-safe with Igor Pro 6.20 or later.
*/
void
CatPossiblyQuotedName(char* str, const char* name)
{
	int len;
	
	len = (int)strlen(str);
	strcpy(str + len, name);
	PossiblyQuoteName(str + len);
}

/*	CleanupName(beLiberal, name, maxNameChars)
	
	name contains an Igor object name.
	CleanupName changes it, if necessary, to make it a legal name.
	
	For most uses, pass MAX_OBJ_NAME for the maxNameChars parameter.

	Igor Pro allows wave and data folder names to contain characters,
	such as space and dot, that were previously illegal in names. We call this
	"liberal" name rules.
	
	If beLiberal is non-zero, CleanupName uses liberal name rules. Liberal rules are still not
	allowed for string and numeric variable names so pass zero for beLiberal for these objects.
	
	If you are going to use the name in Igor's command line (via an Execute operation
	or via the XOPCommand or XOPSilentCommand callbacks), and if the name uses liberal
	rules, the name may need to be single-quoted. In these cases, you should call
	PossiblyQuoteName after calling CleanupName.
	
	Thread Safety: CleanupName is thread-safe with Igor Pro 6.20 or later.
*/
int
CleanupName(int beLiberal, char *name, int maxNameChars)
{
	return (int)CallBack3(CLEANUP_NAME, (void*)beLiberal, name, (void*)maxNameChars);
}

/*	CreateValidDataObjectName(...)
	
	This routine is designed to do all of the nasty work needed to get
	a legal name for a given object. It cleans up illegal names and resolves
	name conflicts.

	It returns in outName a name that can be safely used to create a data object
	of a given type (wave, string, variable, numeric variable).
	
	inName is the proposed name for the object.

	outName is the name after possible cleanup and uniquification.
	
	suffixNumPtr is used to speed up the search for unique names.
	See UniqueName2 for details.
	
	dataFolderH is a handle to a data folder or NULL to use the current data folder.
	
	objectType is one of the following:
		WAVE_OBJECT, VAR_OBJECT, STR_OBJECT, DATAFOLDER_OBJECT.

	beLiberal is 1 if you want to allow liberal names or 0 if not.

	allowOverwrite is 1 if it is OK to for outName to be the name of an
	existing object of the same type.
	
	inNameIsBaseName is 1 if inName is a base name (e.g., "wave") to which a
	suffix (e.g., "0") must always be added to produce the actual name (e.g., "wave0").
	If inNameIsBaseName is 0, then no suffix will be added to inName unless it
	is needed to make the name unique.
	
	printMessage is 1 if you want CreateValidDataObjectName to print a message
	in Igor's history area if an unexpected name conflict occurs. A message is
	printed if you are not using a base name and not allowing overwriting and
	there is a name conflict. A message is also printed if a conflict with
	an object of a different type prevents the normal name from being used.
	
	CreateValidDataObjectName sets *nameChangedPtr to the truth that outName
	is different from inName.
	
	It sets *doOverwritePtr to 1 if outName is the name of an existing object
	of the same type and allowOverwrite is 1.
	
	inName and outName can point to the same array if you don't want to
	preserve the original name. Both must be big enough to hold MAX_OBJ_NAME+1 bytes.
	
	If the object type is VAR_OBJECT or STR_OBJECT, the name will not be
	liberal, even if beLiberal is 1. Igor allows only wave and data folder
	names to be liberal.
	
	Returns 0 if OK or an error code.
	
	Thread Safety: CreateValidDataObjectName is thread-safe with Igor Pro 6.20 or later.
*/
int
CreateValidDataObjectName(
	DataFolderHandle dataFolderH,
	const char* inName, char* outName, int* suffixNumPtr, int objectType,
	int beLiberal, int allowOverwrite, int inNameIsBaseName, int printMessage,
	int* nameChangedPtr, int* doOverwritePtr)
{
	DataFolderHandle oldCurrentDataFolderH;
	char name[MAX_OBJ_NAME+1];
	char inName2[MAX_OBJ_NAME+1];
	int namespaceCode;
	int needToUniquify;
	int objectType2;
	char temp[256];
	int needMessage;
	int err;
	
	err = 0;
	*nameChangedPtr = *doOverwritePtr = 0;
	needMessage = 0;

	// HR, 090629: Was using wrong namespace for data folders.
	namespaceCode = objectType==DATAFOLDER_OBJECT ? DATAFOLDERS_NAME_SPACE : MAIN_NAME_SPACE;
	
	/*	HR, 030701: Set current data folder based on dataFolderH because some of the calls
		below (e.g., UniqueName2), are not data-folder-aware.
	*/
	oldCurrentDataFolderH = NULL;
	if (dataFolderH != NULL) {
		if (GetCurrentDataFolder(&oldCurrentDataFolderH) != 0)
			oldCurrentDataFolderH = NULL;					// Should never happen.
		if (err = SetCurrentDataFolder(dataFolderH))
			return err;										// dataFolderH is invalid?
	}

	// HR, 090629: Was not cleared correctly if inNameIsBaseName was true and inName was long.
	MemClear(inName2, sizeof(inName2));						// Needed to guarantee null terminator if inName >= MAX_OBJ_NAME characters. 

	// HR, 050110, XOP Toolkit 5.04: Avoid crash if inName is too long.
	if (inNameIsBaseName)
		strncpy(inName2, inName, MAX_OBJ_NAME-5);
	else
		strncpy(inName2, inName, MAX_OBJ_NAME);

	strcpy(name, inName2);	// Assume that inName2 is the desired name, not a base name to which we must add a suffix.
	
	// If inName2 is a base name to which we must add a suffix, add the suffix here.
	if (inNameIsBaseName) {
		if (allowOverwrite) {
			sprintf(name, "%s%d", inName2, *suffixNumPtr);				// Use next suffix even if it is already used.
			*suffixNumPtr += 1;
		}
		else {
			UniqueName2(namespaceCode, inName2, name, suffixNumPtr);	// Choose a suffix that is not already used.
		}
	}

	if (objectType!=WAVE_OBJECT && objectType!=DATAFOLDER_OBJECT)
		beLiberal = 0;										// Liberal names not supported for variables.
	if (err = CleanupName(beLiberal, name, MAX_OBJ_NAME))	// Remove illegal characters, etc.
		goto done;
	
	needToUniquify = 0;
	do {
		if (GetDataFolderObject(dataFolderH, name, &objectType2, NULL) == 0) {
			// If here, an object exists with the specified name.
			
			if (allowOverwrite) {
				if (objectType2 == objectType) {
					*doOverwritePtr = 1;
					break;							// OK to overwrite object.
				}
			}
			
			/*	If here, we must choose another name because the name is in use
				for a different type of object or it is in use for a the same type but
				we are not allowed to overwrite it.
			*/
			needToUniquify = 1;
			needMessage = printMessage;
		}
		else {
			/*	Name is not a name of an existing object. Make sure it is not in conflict
				with a function, operation or some other reserved name.
			*/
			if (CheckName(dataFolderH, objectType, name)) {
				// There is a conflict with a different type of object or function or operation.
				if (allowOverwrite) {
					/*	HR, 090218, XOP Toolkit 5.09
						There is a name conflict with an object of a different type (possibly a built-in
						operation or function). Since we are allowed to overwrite, we don't want to create
						a unique name but just a legal name that does not conflict with an object of
						a different type.
					*/
					char suffixStr[16];
					int tempLen, suffixStrLen;

					strcpy(temp, name);
					tempLen = (int)strlen(temp);
					while(1) {
						sprintf(suffixStr, "%d", *suffixNumPtr);		// Use next suffix even if it is already used.
						*suffixNumPtr += 1;
						suffixStrLen = (int)strlen(suffixStr);
						if (tempLen + suffixStrLen > MAX_OBJ_NAME) {	// e.g., tempLen=31, suffixStrLen=1
							temp[MAX_OBJ_NAME-suffixStrLen] = 0;		// e.g., temp[30]=0 leaving temp[30] available for the suffix.
							tempLen = (int)strlen(temp);
						}
						sprintf(name, "%s%s", temp, suffixStr);
						if (GetDataFolderObject(dataFolderH, name, &objectType2, NULL) == 0) {
							if (objectType2 == objectType) {
								*doOverwritePtr = 1;
								break;									// An object of this type exists. Since allowOverwrite is true, we will use this name.
							}
						}
						if (CheckName(dataFolderH, objectType, name) == 0)
							break;
					}
				}
				else {
					needToUniquify = 1;					// There is some kind of conflict and overwrite is not allowed.
					needMessage = printMessage;
				}
			}
		}

		if (!needToUniquify)
			break;

		strcpy(temp, name);
		if (err = UniqueName2(namespaceCode, temp, name, suffixNumPtr))
			goto done;
	} while(0);
	
	*nameChangedPtr = strcmp(inName, name)!=0;
	if (needMessage) {
		char objectTypeStr[64];
		switch(objectType) {
			case WAVE_OBJECT:
				strcpy(objectTypeStr, "Wave");
				break;
			case VAR_OBJECT:
				strcpy(objectTypeStr, "Variable");
				break;
			case STR_OBJECT:
				strcpy(objectTypeStr, "String");
				break;
			case DATAFOLDER_OBJECT:
				strcpy(objectTypeStr, "Data folder");
				break;
			default:
				sprintf(objectTypeStr, "BUG: CreateValidDataObjectName, objectType=%d", objectType);
				break;
		}
		sprintf(temp, "%s \'%s\' changed to \'%s\' because of a conflict.\015", objectTypeStr, inName2, name);
		XOPNotice(temp);
	}
	
	strcpy(outName, name);

done:
	if (oldCurrentDataFolderH != NULL)
		SetCurrentDataFolder(oldCurrentDataFolderH);
	
	return err;
}


// *** IGOR Color Table Routines ***

/*	GetIndexedIgorColorTableName(index, name)

	Returns via name the name of a color table indicated by the index
	or "" if the index is invalid. Valid indices start from zero. You can find
	the maximum valid index by calling this routine with increasing indices
	until it returns an error.

	The function result is 0 if OK or a non-zero error code.
	
	Thread Safety: GetIndexedIgorColorTableName is thread-safe with Igor Pro 6.20 or later.
*/
int
GetIndexedIgorColorTableName(int index, char name[MAX_OBJ_NAME+1])
{
	*name = 0;
	
	return (int)CallBack2(GET_INDEXED_IGORCOLORTABLE_NAME, (void*)index, name);
}

/*	GetNamedIgorColorTableHandle(name, ictHPtr)
	
	Returns via *ictHPtr a handle to an IGOR color table or NULL in case of error.
	
	The returned handle belongs to Igor. Do not modify or delete it.
	
	The IgorColorTableHandle type is defined in IgorXOP.h.
	
	The name parameter is case insensitive.
	
	As of this writing, the following IGOR color tables exist:
		Rainbow, Grays, YellowHot, BlueHot
		BlueRedGreen, RedWhiteBlue, PlanetEarth, Terrain
	You can find the names of all color tables using GetIndexedIgorColorTableName.
	
	The function result is 0 if OK or a non-zero error code.
	
	Thread Safety: GetNamedIgorColorTableHandle is thread-safe with Igor Pro 6.20 or later.
*/
int
GetNamedIgorColorTableHandle(const char *name, IgorColorTableHandle* ictHPtr)
{
	*ictHPtr = NULL;

	return (int)CallBack2(GET_NAMED_IGORCOLORTABLE_HANDLE, (void*)name, ictHPtr);
}

/*	GetIgorColorTableInfo(ictH, name, numColorsPtr)

	Provides access to the name and the number of colors in the IGOR
	color table specified by ictH.

	The IgorColorTableHandle type is defined in IgorXOP.h.
	
	If you don't want to know the name, pass NULL for name.
	If you don't want to know the number of colors, pass NULL for numColorsPtr.
	
	The function result is 0 if OK or a non-zero error code.
	
	Thread Safety: GetIgorColorTableInfo is thread-safe with Igor Pro 6.20 or later.
*/
int
GetIgorColorTableInfo(IgorColorTableHandle ictH, char name[MAX_OBJ_NAME+1], int* numColorsPtr)
{
	// In case IGOR_OBSOLETE.
	if (name != NULL)
		*name = 0;
	if (numColorsPtr != NULL)
		*numColorsPtr = 0;
	
	return (int)CallBack3(GET_IGORCOLORTABLE_INFO, ictH, name, numColorsPtr);
}

/*	GetIgorColorTableValues(ictH, startColorIndex, endColorIndex, updatePixelValues, csPtr)
	
	Returns via csPtr a description of the colors associated with the IGOR color
	table specified by ictH.

	The IgorColorTableHandle type is defined in IgorXOP.h.
	
	startColorIndex and endColorIndex specify the indices of the colors for
	which you want to get a description. startColorIndex must be between 0 and
	the number of colors in the table minus one. endColorIndex must be between
	startColorIndex and the number of colors in the table minus one. You can
	find the number of colors in the table using GetIgorColorTableInfo.
	
	The IgorColorSpec structure contains an RGBColor field which identifies the
	RGB color for a color table entry with a particular index.
	
	The value field of the IgorColorSpec structure tells you the pixel value that
	would need to be written to video RAM to get the associated color to appear on
	the screen when the monitor is in 16 or 256 color mode. It is typically used by
	advanced programmers who are writing directly to offscreen bitmap memory. 
	
	However, when a monitor is running in 16 or 256 color mode, this value is
	invalidated whenever the system changes the hardware color lookup table, which
	can happen at any time. If you pass non-zero for the updateVals parameter, then
	Igor will update the value field for each color before returning it to you and
	it will be accurate until the next time the system changes the hardware color
	lookup table. If you pass zero for the updateVals parameter, then Igor will not
	update the value field and it is likely to be stale. 
	
	Updating the value fields takes time so you should pass non-zero for the
	updateVals parameter only if you really need accurate pixel values. For
	example, if you just want to know what RGB colors appear in a particular color
	table then you don't need the pixel values and should pass 0 for the updateVals
	parameter. On the other hand, if you are writing into an offscreen bitmap in
	preparation for blasting it to the screen, then you need accurate pixel values
	and you should pass 1 for updateVals. 
	
	The function result is 0 if OK or a non-zero error code.
	
	Thread Safety: GetIgorColorTableValues is thread-safe with Igor Pro 6.20 or later.
*/
int
GetIgorColorTableValues(IgorColorTableHandle ictH, int startColorIndex, int endColorIndex, int updatePixelValues, IgorColorSpec* csPtr)
{
	return (int)CallBack5(GET_IGORCOLORTABLE_VALUES, ictH, (void*)startColorIndex, (void*)endColorIndex,  (void*)updatePixelValues, csPtr);
}


// *** Cross-Platform Utilities ***

/*	WinRectToMacRect(wr, mr)

	Thread Safety: WinRectToMacRect is thread-safe with Igor Pro 6.20 or later.
*/
void
WinRectToMacRect(const RECT* wr, Rect* mr)
{
	// In principle, the Windows coordinates could exceed 16 bits. In practice, this should not happen.

	mr->left = (short)wr->left;
	mr->right = (short)wr->right;
	mr->top = (short)wr->top;
	mr->bottom = (short)wr->bottom;
}

/*	MacRectToWinRect(mr, wr)

	Thread Safety: MacRectToWinRect is thread-safe with Igor Pro 6.20 or later.
*/
void
MacRectToWinRect(const Rect *mr, RECT *wr)
{
	wr->left = mr->left;
	wr->top = mr->top;
	wr->right = mr->right;
	wr->bottom = mr->bottom;
}


// *** Miscellaneous Routines ***

/*	XOPBeep()

	Emits a short beep.
	
	Thread Safety: XOPBeep is not thread-safe.
*/
void
XOPBeep(void)
{
	if (!CheckRunningInMainThread("XOPBeep"))
		return;

	#ifdef MACIGOR
		SysBeep(5);
	#endif
	#ifdef WINIGOR
		MessageBeep(MB_ICONEXCLAMATION);
	#endif
}

/*	GetXOPIndString(text, strID, item)

	Tries to get string from a STR# resource in the XOP's resource fork.
	Does not search any other resource forks and does not change current
	resource fork on Macintosh.
	
	text is returned as a C string.
	
	Thread Safety: GetXOPIndString is not thread-safe.
*/
void
GetXOPIndString(
	char* text,									// Should hold up to 256 bytes.
	int strID,									// Resource ID of STR# resource.
	int index)									// String number starting from 1.
{
	*text = 0;					// HR, 981022: For XOP Toolkit 3.1.

	if (!CheckRunningInMainThread("GetXOPIndString"))
		return;
	
	#ifdef MACIGOR				// [
	{	int curResFile;
		
		if (GetXOPResource('STR#', strID) == NULL) 	// No such STR# resource?
			return;
	
		curResFile = CurResFile();
		UseResFile(XOPRefNum());					// Search XOP's resource fork.
		GetIndString((unsigned char*)text, strID, index);
		CopyPascalStringToC((ConstStr255Param)text, text);
		UseResFile(curResFile);
	}
	#endif						// ]

	#ifdef WINIGOR				// [
		GetWinIndString(XOPModule(), text, strID, index);
	#endif						// ]
}

/*	ArrowCursor()
	
	Thread Safety: ArrowCursor is not thread-safe.
*/
void
ArrowCursor(void)
{
	if (!CheckRunningInMainThread("ArrowCursor"))
		return;
	CallBack1(SETCURSOR, (void *)ARROWCURSOR);
}

/*	IBeamCursor()
	
	Thread Safety: IBeamCursor is not thread-safe.
*/
void
IBeamCursor(void)
{
	if (!CheckRunningInMainThread("IBeamCursor"))
		return;
	CallBack1(SETCURSOR, (void *)IBEAMCURSOR);
}

/*	WatchCursor()
	
	Thread Safety: WatchCursor is not thread-safe.
*/
void
WatchCursor(void)
{
	if (!CheckRunningInMainThread("WatchCursor"))
		return;
	CallBack1(SETCURSOR, (void *)WATCHCURSOR);
}

/*	HandCursor()
	
	Thread Safety: HandCursor is not thread-safe.
*/
void
HandCursor(void)
{
	if (!CheckRunningInMainThread("HandCursor"))
		return;
	CallBack1(SETCURSOR, (void *)HANDCURSOR);
}

/*	SpinCursor()
	
	Thread Safety: SpinCursor is not thread-safe.
*/
void
SpinCursor(void)
{
	static TickCountInt lastTime=0;			// Prevents wasting too much time spinning.

	if (!CheckRunningInMainThread("SpinCursor"))
		return;
	
	if (TickCount() > lastTime+6) {
		CallBack1(SETCURSOR, (void *)SPINNINGCURSOR);
		lastTime = TickCount();
	}
}

/*	SpinProcess()

	SpinProcess spins the beach ball cursor AND gives MultiFinder (if active) a chance
	to do background processing. Also, Igor may be moved from the background to the
	foreground or vice versa when SpinProcess is called.
	
	It returns non-zero if the user has pressed cmd-dot recently or zero otherwise.
	
	Thread Safety: SpinProcess is thread-safe with Igor Pro 6.20 or later.
*/
int
SpinProcess(void)
{
	if (igorVersion<620 && !CheckRunningInMainThread("SpinProcess"))	// SpinProcess is threadsafe only
		return -1;														// in Igor Pro 6.20 or later.

	return (int)CallBack0(SPINPROCESS);
}

/*	DoUpdate()

	Causes Igor to do an immediate update.
	An update consists of
		updating any windows that have been uncovered
		updating any graphs, tables or layouts that need it
		recalculating any dependent objects that need it
	
	The DOUPDATE message was added in Igor 2.0.
	If the version of Igor is earlier than 2.0, DoUpdate triggers
	an update by using XOPCommand. The update is a side effect of
	XOPCommand.
	
	Thread Safety: DoUpdate is not thread-safe.
*/
int
DoUpdate(void)
{
	if (!CheckRunningInMainThread("DoUpdate"))
		return NOT_IN_THREADSAFE;

	CallBack0(DOUPDATE);
	return 0;
}

/*	PauseUpdate(savePtr)

	Used to temporarily suppress updating of Igor windows during macro execution.
	savePtr is a pointer to an int to save the previous state of PauseUpdate.
	
	This MUST be balanced by a ResumeUpdate(savePtr).
	
	Thread Safety: PauseUpdate is not thread-safe.
*/
void
PauseUpdate(int *savePtr)
{
	if (!CheckRunningInMainThread("PauseUpdate"))
		return;

	CallBack1(PAUSEUPDATE, savePtr);
}

/*	ResumeUpdate(savePtr)

	Used to undo temporary suppression of updating of Igor windows during macro execution.
	
	This MUST called to balance PauseUpdate(savePtr) calls.
	
	Thread Safety: ResumeUpdate is not thread-safe.
*/
void
ResumeUpdate(int *savePtr)
{
	if (!CheckRunningInMainThread("ResumeUpdate"))
		return;

	CallBack1(RESUMEUPDATE, savePtr);
}

/*	WaveList(listHandle, match, sep, options)
	
	Thread Safety: WaveList is not thread-safe.
*/
int
WaveList(Handle listHandle, const char *match, const char *sep, const char *options)
{
	if (!CheckRunningInMainThread("WaveList"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(WAVELIST, listHandle, (void*)match, (void*)sep, (void*)options);
}

/*	WinList(listHandle, match, sep, options)
	
	Thread Safety: WinList is not thread-safe.
*/
int
WinList(Handle listHandle, const char *match, const char *sep, const char *options)
{
	if (!CheckRunningInMainThread("WinList"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(WINLIST, listHandle, (void*)match, (void*)sep, (void*)options);
}

/*	PathList(listHandle, match, sep, options)
	
	Thread Safety: PathList is not thread-safe.
*/
int
PathList(Handle listHandle, const char *match, const char *sep, const char *options)
{
	if (!CheckRunningInMainThread("PathList"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(PATHLIST, listHandle, (void*)match, (void*)sep, (void*)options);
}

// GetPathInfo is no longer supported. Use GetPathInfo2.

/*	GetPathInfo2(pathName, fullDirPath)

	pathName is the name of an Igor symbolic path.
	
	Returns via fullDirPath the full native path to the directory referenced by pathName.
	The returned path includes a trailing colon on Macintosh and a trailing backslash on Windows.
	
	Returns 0 if OK or an error code if the pathName is not the name of an existing
	Igor symbolic path.
	
	Thread Safety: GetPathInfo2 is not thread-safe.
*/
int
GetPathInfo2(const char* pathName, char fullDirPath[MAX_PATH_LEN+1])
{
	if (!CheckRunningInMainThread("GetPathInfo2"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack2(PATHINFO2, (void*)pathName, (void*)fullDirPath);
}

/*	GetNamedFIFO

	Returns NamedFIFO handle or NULL if none of that name exists
	NamedFIFO data structure and meaning are documented in a separate file.
	
	Thread Safety: GetNamedFIFO is not thread-safe.
*/
struct NamedFIFO **
GetNamedFIFO(const char *name)
{
	if (!CheckRunningInMainThread("GetNamedFIFO"))
		return NULL;
	
	return (NamedFIFO**)CallBack1(GETNAMEDFIFO, (void*)name);
}

/*	MarkFIFOUpdated

	Call this after putting data in a named fifo so chart gadgets will refresh.
	
	Thread Safety: MarkFIFOUpdated is not thread-safe.
*/
void
MarkFIFOUpdated(struct NamedFIFO **fifo)
{
	if (!CheckRunningInMainThread("MarkFIFOUpdated"))
		return;

	CallBack1(MARKFIFOUPDATED, fifo);
}

/*	SaveXOPPrefsHandle(prefsHandle)

	Saves the handle in IGOR's preferences file. You can retrieve the handle
	using GetXOPPrefsHandle.
	
	IGOR makes a copy of the data in the handle, so the handle is still yours
	after you call this. Keep or dispose of it as you wish.
	
	If you pass NULL for the prefsHandle parameter, IGOR removes any existing
	XOP preferences from the IGOR preferences file.
	
	IGOR uses the name of your XOP's file to distinguish your preferences from
	the preferences of other XOPs.
	
	Each time you call this routine, the Igor preferences file is opened and closed.
	Therefore, it is best to call each of it only once. One way to do this is to call
	GetXOPPrefsHandle when your XOPs starts and SaveXOPPrefsHandle when you receive
	the CLEANUP message.
	
	Returns 0 if OK or a non-zero error code.
	
	Thread Safety: SaveXOPPrefsHandle is not thread-safe.
*/
int
SaveXOPPrefsHandle(Handle prefsHandle)
{
	if (!CheckRunningInMainThread("SaveXOPPrefsHandle"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack1(SAVE_XOP_PREFS, prefsHandle);
}

/*	GetXOPPrefsHandle(prefsHandlePtr)

	Retrieve your XOP's preference handle from the IGOR preferences file, if you
	have previously stored it there using SaveXOPPrefsHandle. In this case,
	on return, *prefsHandlePtr will be your preferences handle. This handle
	is allocated by IGOR but belongs to you to keep or dispose as you wish.
	
	If the IGOR preferences file does not contain your preferences, on return,
	*prefsHandlePtr will be NULL and the function result will be 0.
	
	IGOR uses the name of your XOP's file to distinguish your preferences from
	the preferences of other XOPs.
	
	Each time you call this routine, the Igor preferences file is opened and closed.
	Therefore, it is best to call each of it only once. One way to do this is to call
	GetXOPPrefsHandle when your XOPs starts and SaveXOPPrefsHandle when you receive
	the CLEANUP message.
	
	Returns 0 if OK or a non-zero error code.
	
	Thread Safety: GetXOPPrefsHandle is not thread-safe.
*/
int
GetXOPPrefsHandle(Handle* prefsHandlePtr)
{
	*prefsHandlePtr = NULL;

	if (!CheckRunningInMainThread("GetXOPPrefsHandle"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack1(GET_XOP_PREFS, prefsHandlePtr);
}

/*	GetPrefsState(prefsStatePtr)

	Returns via bit 0 of prefsStatePtr the truth that preferences are on.
	See the IGOR Pro manual for information about the preferences on/off state.
	Other bits are reserved for future use.

	Function result is 0 if OK or an error code.
	
	Thread Safety: GetPrefsState is thread-safe with Igor Pro 6.20 or later.
*/
int
GetPrefsState(int* prefsStatePtr)
{
	*prefsStatePtr = 1;				// Default if we are running with old IGOR.
	return (int)CallBack1(GET_PREFS_STATE, prefsStatePtr);
}

/*	XOPDisplayHelpTopic(title, topicStr, flags)

	Displays help for the specified topic.
	
	topicStr is a help topic string that matches a help topic or subtopic in a
	native Igor help file. Igor first searches open help files for the topic.
	If it is not found, Igor then searches all Igor native help files in the folder
	containing the XOP file and subfolders. If it is still not found Igor then
	searches all Igor native help files in the Igor Pro folder and subfolders.
	
	The help file must be compiled in order for Igor to find the topic. Each time
	you open a file as a help file, Igor checks to see if it is compiled and if
	not asks if you want to compile it.
	
	topicStr may have one of the following formats:
		Format							Example
		<topic name>					"GBLoadWave XOP"
		<subtopic name>					"The Load General Binary Dialog"
		<topic name>[<subtopic name>]	"GBLoadWave XOP[The Load General Binary Dialog]"
		
	If the topic that you want to display is a subtopic, you should use the last
	form since it minimizes the chance that Igor will find another help file with
	the same subtopic name. Also, you must choose descriptive topic and
	subtopic names to minimize the chance of a conflict between two help files.
	
	Note that once you reference a help topic or subtopic from your executable code,
	you must be careful to avoid changing the name of the topic or subtopic.
	
	The title parameter is used only for modal help and supplies the title for
	a modal dialog containing the help.
	
	The flags parameter is interpreted bitwise as follows:
		Bit 0			If cleared, Igor displays non-modal help.
						If set, Igor displays modal help.						
						
		Bit 1			If cleared, during modal help Igor displays the entire
						help file (if it is not too big) in the modal help dialog,
						with the specified topic initially in view. If set, during
						modal help Igor displays just the specified topic.
						
		Bit 2			If cleared, if the topic can not be found Igor displays an error
						dialog. If set, if the topic can not be found Igor does not display
						an error dialog.
						
		All other bits	Reserved - set to zero.
						
	You MUST set bit 0 if you call XOPDisplayHelpTopic from a modal dialog. This causes
	Igor do display a dialog containing help on top of your dialog. If you fail to set
	bit zero when calling XOPDisplayHelpTopic from a modal dialog, Igor may behave erratically.
	Unfortunately, links in help files don't work during modal help.
	
	If you are calling XOPDisplayHelpTopic in a non-modal situation, it is appropriate to
	clear bit zero, but not required. If you clear bit zero, Igor displays a normal
	Igor help file. If you set bit zero, Igor displays a modal help dialog.
	
	You must set all other bits to zero.

	Function result is 0 if OK or IGOR_OBSOLETE if you are running with an older
	version of IGOR or some other non-zero code if the topic can not be found.
	
	Thread Safety: XOPDisplayHelpTopic is not thread-safe.
*/
int
XOPDisplayHelpTopic(const char* title, const char* topicStr, int flags)
{
	if (!CheckRunningInMainThread("XOPDisplayHelpTopic"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(DISPLAY_HELP_TOPIC, (void*)title, (void*)topicStr, (void*)flags);
}

/*	XOPSetContextualHelpMessage(theWindow, message, r)

	Displays a message in the Igor Tips help window on Macintosh or in the
	status bar on Windows. Call this when your window is active and the user
	moves the cursor over an icon or other area of the window about which you
	have something to say.
	
	theWindow is your WindowRef on Macintosh or your HWND on Windows.
	This refers to the window containing the control or icon for which
	you are providing help.
	
	message is a C string containing the message to display.
	
	r is a pointer to a Macintosh rectangle, even on Windows, that indicates the
	area of the window that the icon occupies. When the user moves the cursor
	out of this rectangle, Igor will remove the message. On Macintosh, this
	rectangle is in the local coordinates of the window containing the control
	or icon. On Windows, it is in client coordinates of the window containing
	the control or icon. On Windows, use WinRectToMacRect to translate the Windows
	RECT into a Macintosh Rect.

	Function result is 0 if OK or IGOR_OBSOLETE.
	
	The WindowXOP1 sample XOP illustrates the use of this function.
	
	Thread Safety: XOPSetContextualHelpMessage is not thread-safe.
*/
int
XOPSetContextualHelpMessage(XOP_WINDOW_REF theWindow, const char* message, const Rect* r)
{
	if (!CheckRunningInMainThread("XOPSetContextualHelpMessage"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(SET_CONTEXTUAL_HELP_MESSAGE, (void*)theWindow, (void*)message, (void*)r);
}

/*	DoWindowRecreationDialog(procedureName)
	
	This routine is of use only to very advanced XOPs that support window recreation
	macros.

	On input, procedureName is the proposed name for the recreation macro.
	
	DoWindowRecreationDialog displays the standard Igor Close Window dialog that
	allows the user to enter a name for a recreation macro and choose to save the
	macro, replace an existing macro, skip saving the macro, or cancel.
	
	If the user clicks Save or Replace, on output, procedureName is the
	name chosen by the user for the recreation macro.
	
	DoWindowRecreationDialog returns the following codes:
		kCloseWinCancel			The user clicked the cancel button.
								The XOP must cancel the close of the window.
								
		kCloseWinSave			The user clicked the Save button.
								The XOP must save the macro in the Igor procedure
								window and close the window.
								
		kCloseWinReplace		The user clicked the Replace button.
								The XOP must save the macro in the Igor procedure
								window and close the window.
								
		kCloseWinNoSave			The user clicked the No Save button.
								The XOP must close the window without saving any macro.
	
	Thread Safety: DoWindowRecreationDialog is not thread-safe.
*/
enum CloseWinAction
DoWindowRecreationDialog(char* procedureName)
{
	if (!CheckRunningInMainThread("DoWindowRecreationDialog"))
		return kCloseWinCancel;

	return (enum CloseWinAction)CallBack1(WINDOW_RECREATION_DIALOG, procedureName);
}

/*	GetIgorProcedureList(hPtr, flags)

	This routine is intended only for very advanced XOPs, such as XOPs that
	generate their own window recreation macros. An example is the WaveMetrics
	Surface Plotter XOP. It requires Igor Pro 3.13B03 or later
	
	The main use for this routine is to check if a particular macro or function
	exists in the Igor procedure windows.
	
	GetIgorProcedureList returns via *hPtr a handle to a semicolon-separated
	list of procedure names. Depending on the flags parameter, the list may contain
	names of macros, names of user functions, or names of both macros and user functions.
	
	Note that the handle is not null terminated. Use GetHandleSize to determine how
	many characters are in the handle. This handle belongs to you, so call
	DisposeHandle to dispose it when you no longer need it.
	
	If Igor can not build the list, it returns a non-zero error code and sets
	*hPtr to NULL.
	
	The flags parameter is defined as follows:
		Bit 0:	If set, GetIgorProcedureList will list all macros.
				If cleared, it will ignore all macros.
		
		Bit 1:	If set, GetIgorProcedure will list all user functions.
				If cleared, it will ignore all user functions.
		
		All other bits are reserved for future use and must be set to 0.
	
	Igor will be unable to build the list if a syntactical error in the procedure
	files prevents Igor from successfully scanning them. In this case, GetIgorProcedureList
	will return NEED_COMPILE.
	
	GetIgorProcedureList can also return NOMEM if it runs out of memory. This is
	unlikely.
	
	Future versions of GetIgorProcedureList may return other error codes so your
	XOP should not crash or otherwise grossly misbehave if it receives some
	other error code.
	
	Thread Safety: GetIgorProcedureList is not thread-safe.
*/
int
GetIgorProcedureList(Handle* hPtr, int flags)
{
	if (!CheckRunningInMainThread("GetIgorProcedureList"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack2(GET_IGOR_PROCEDURE_LIST, hPtr, (void*)flags);
}

/*	GetIgorProcedure(procedureName, hPtr, flags)

	This routine is intended only for very advanced XOPs, such as XOPs that
	generate their own window recreation macros. An example is the WaveMetrics
	Surface Plotter XOP. It requires Igor Pro 3.13B03 or later.
	
	The main use for this routine is to check if a particular macro or function
	exists in the Igor procedure windows.

	If Igor can find the procedure (macro or function) specified by procedureName,
	GetIgorProcedure returns via *hPtr a handle to the text for the procedure
	and returns a result of 0. The handle will contain the text for the specified
	procedure with a carriage return at the end of each line.
	
	Note that the handle is not null terminated. Use GetHandleSize to determine how
	many characters are in the handle. This handle belongs to you, so call
	DisposeHandle to dispose it when you no longer need it.
	
	If Igor can not find the procedure, it returns a non-zero error code and sets
	*hPtr to NULL.
	
	The flags parameter is defined as follows:
		Bit 0:	If set, GetIgorProcedure will look for macros with the
				specified name. If cleared, it will ignore all macros.
		
		Bit 1:	If set, GetIgorProcedure will look for user functions with the
				specified name. If cleared, it will ignore all user functions.
		
		All other bits are reserved for future use and must be set to 0.
	
	Igor will be unable to find the procedure if there is no such procedure. In
	this case, GetIgorProcedure will return NO_MACRO_OR_FUNCTION.
	
	Igor will be unable to find the procedure if a syntactical error in the procedure
	files prevents Igor from successfully scanning them. In this case, GetIgorProcedure
	will return NEED_COMPILE.
	
	GetIgorProcedure can also return NOMEM if it runs out of memory. This is
	unlikely.
	
	Future versions of GetIgorProcedure may return other error codes so your
	XOP should not crash or otherwise grossly misbehave if it receives some
	other error code.
	
	Thread Safety: GetIgorProcedure is not thread-safe.
*/
int
GetIgorProcedure(const char* procedureName, Handle* hPtr, int flags)
{
	if (!CheckRunningInMainThread("GetIgorProcedure"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(GET_IGOR_PROCEDURE, (void*)procedureName, hPtr, (void*)flags);
}

/*	SetIgorProcedure(procedureName, h, flags)

	This routine is intended only for very advanced XOPs, such as XOPs that
	generate their own window recreation macros. An example is the WaveMetrics
	Surface Plotter XOP. It requires Igor Pro 3.13B03 or later
	
	SetIgorProcedure will fail if Igor procedures are running when it is called.
	It is intended to be called in response to a user action, such as when a user
	clicks the close box of an XOP window and you want to create a recreation macro
	for that window.
	
	SetIgorProcedure stores the procedure text, which is contained in the handle
	h, in a procedure window.
	
	You must make sure that procedureName agrees with the name of the procedure
	as found in the handle. Otherwise, Igor may be left in an unstable state.
	For example, if the text in h is this:
		Macro HelloWorld()
			Print "Hello World"
		End
	then procedureName must be "HelloWorld".
	
	The handle h must contain the lines of the procedure with a carriage return
	after each line. Do not include linefeeds. There must be no extraneous lines
	before or after the procedure text. The text in the handle should not be
	null terminated. There must be no garbage characters in the handle.
	
	The handle h belongs to Igor. Once you pass it to SetIgorProcedure, you
	must not modify it, access it, or dispose of it.

	The flags parameter is defined as follows:
		Bit 0:	Set if the text in the handle is a macro (defined by the keywords
				Macro, Proc or Window).
				
		Bit 1:	Set if the text in the handle is a function (defined by the keyword
				Function).
				
		Bit 2:	If cleared and if there is already a procedure with the same name,
				SetIgorProcedure will return an error. If set and if there is
				already a procedure with the same name, SetIgorProcedure will
				replace the existing procedure and return 0.
		
		All other bits are reserved for future use and must be set to 0.

	Thus, the flags parameter must be one of the following values:
		1		Text is a macro. Return error if conflict.
		2		Text is a function. Return error if conflict.
		5		Text is a macro. Overwrite existing procedure if conflict.
		6		Text is a function. Overwrite existing procedure if conflict.
		
	The flags parameter must be consistent with the text in the handle. For example,
	you must not pass 1 for flags if the text in the handle is a function.
	
	SetIgorProcedure normally stores the procedure in the built-in Procedure window.
	However, if there already exists a procedure with the specified name and if
	that procedure is in another procedure window, SetIgorProcedure stores the
	procedure in that procedure window.
	
	After storing the procedure in the procedure window, SetIgorProcedure causes
	the Igor procedure windows to be compiled if auto-compile is on. If auto-compile
	is off, then the user must cause procedures to be compiled by choosing Compile
	from the Macros menu. Because of this, an XOP can not depend on the procedure
	being available for execution as soon as SetIgorProcedure returns. There is
	currently no way for the XOP to know the state of the auto-compile setting.
	
	SetIgorProcedure returns 0 if it successfully stores the text in a procedure
	window. Otherwise, it returns a non-zero error code.
	
	SetIgorProcedure returns a non-zero error code in the following cases:

	1. SetIgorProcedure was called while Igor was compiling procedures.
	   This could happen only  under weird recursion circumstances.
	   SetIgorProcedure returns -1 in this case.
	   
	2. SetIgorProcedure was called while Igor was running a macro or function.
	   Igor can not change procedures while they are running because this could
	   cause a crash. SetIgorProcedure returns -1 in this case.
	   
	3. Igor can not store the text because a serious syntax error in the
	   procedure windows prevents it from successfully scanning the procedures.
	   In this case, SetIgorProcedure returns NEED_COMPILE.
	  
	4. procedureName is not a valid name (too long or contains bad characters).
	
	5. procedureName is in conflict with a built-in or external operation or function
	   or a wave or a variable. The error returned depends on the type of object
	   that is in conflict.
	
	6. procedureName is in conflict with a macro or user function and bit 2 of
	   flags is not set. SetIgorProcedure returns NAME_MAC_CONFLICT or
	   NAME_USER_FUNC_CONFLICT.
	
	7. SetIgorProcedure ran out of memory. In this case, SetIgorProcedure returns NOMEM.
	
	8. You are doing a replace and the procedure window containing the original
	   procedure is open for read-only. In this case, you will receive a non-zero
	   error code. If the procedure window is open for read/write but the write-protect
	   icon is on, SetIgorProcedure will replace the procedure and return 0. The write-
	   protect icon is intended only to prevent inadvertent manual changes.
	
	Future versions of SetIgorProcedure may return other error codes so your
	XOP should not crash or otherwise grossly misbehave if it receives some
	other error code.
	
	Thread Safety: SetIgorProcedure is not thread-safe.
*/
int
SetIgorProcedure(const char* procedureName, Handle h, int flags)
{
	if (!CheckRunningInMainThread("SetIgorProcedure"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack3(SET_IGOR_PROCEDURE, (void*)procedureName, h, (void*)flags);
}

/*	GetFunctionInfo(name, fip)

	Returns information that you need in order to call an Igor user function
	or an external function from an XOP. This is an advanced feature that most
	XOP programmers will not need.
	
	name is the name of an existing user or external function. If there is no such
	function, GetFunctionInfo returns an error. If the function is a user function
	and procedures are not in a compiled state, GetFunctionInfo returns an
	error. If everything is OK, GetFunctionInfo returns zero.
	
	The information returned by GetFunctionInfo should be used and then discarded.
	If the user does anything to cause procedures to be compiled then the values
	returned by GetFunctionInfo are no longer valid.
	
	GetFunctionInfo returns via fip->compilationIndex a value that you will
	pass to CallFunction. This value is used to make sure that procedures
	have not been changed between the time you call GetFunctionInfo and
	the time you call CallFunction.
	
	GetFunctionInfo returns via fip->functionID a value which you will pass to
	CallFunction to specify which user function you want to execute.
	
	GetFunctionInfo returns via fip->subType a code that identifies certain
	special purpose functions. The value returned currently has no use but may
	be used in the future.
	
	GetFunctionInfo returns via fip->isExternalFunction a value that is non-zero
	for external functions and zero for user functions. This field is for your
	information only. Your code will be the same for external and user functions.
	
	GetFunctionInfo returns via fip->returnType one of the following codes:
		NT_FP64:			Return value is a double-precision number
		NT_FP64 | NT_CMPLX:	Return value is a complex double-precision number
		HSTRING_TYPE:		Return value is a string
		WAVE_TYPE:			Return value is a wave reference (Igor Pro 6.10 or later)
		DATAFOLDER_TYPE:	Return value is a DFREF (Igor Pro 6.10 or later)
	
	GetFunctionInfo returns via fip->numOptionalParameters, fip->numRequiredParameters
	and fip->totalNumParameters the number of optional, required and total parameters for
	the function. Currently, an XOP can call a user function that has optional parameters
	but the XOP can not pass optional parameters to the function. In other words, it
	must pass the required parameters only.
	
	GetFunctionInfo returns via fip->parameterTypes an array of parameter types.
	GetFunctionInfo stores a parameter type value in elements 0 through fip->totalNumParameters-1.
	Elements fip->totalNumParameters and higher are undefined so you must not use them.
	
	You must use the CheckFunctionForm XOPSupport routine to make sure that the
	function is of the form you want. You normally don't need to examine the parameter
	type values directly, but in case you are curious, here is how a parameter type value
	is defined:
		if (parameterType == NT_FP64)
			Parameter is a double-precision number
			
		if (parameterType == (NT_FP64 | NT_CMPLX))
			Parameter is a complex double-precision number
		
		if (parameterType == HSTRING_TYPE)
			Parameter is a string

		if (parameterType == FV_FUNC_TYPE)
			Parameter is a function reference
		
		if (WAVE_TYPE bit is set) {
			Parameter is a numeric or text wave
			if (WAVE_Z_TYPE bit is set)
				Wave declared with /Z flag (used only by Igor debugger)
			
			if (parameterType & 0x01) {
				Parameter is a complex numeric wave
			}
			else {
				if (isExternalFunction) {
					Parameter is a numeric or text wave (can't tell which till runtime)
				}
				else {
					if ((parameterType & 0xFE) == 0)
						Parameter is a text wave
					else
						Parameter is a numeric wave
				}
			
			}
		}
		
		if (DATAFOLDER_TYPE bit is set)
			Parameter is a DFREF
	
	In addition, the parameter type values for numeric, complex numeric and string
	parameters may be ORed with FV_REF_TYPE. This indicates that the corresponding
	parameter is "pass-by-reference", meaning that the function can change the
	value of that parameter.
	
	Thread Safety: GetFunctionInfo is thread-safe with Igor Pro 6.20 or later.
	The FunctionInfo structure contains an isThreadSafe field which tells you if
	the function in question is threadsafe.
*/
int
GetFunctionInfo(const char* name, FunctionInfoPtr fip)
{
	return (int)CallBack2(GET_FUNCTION_INFO, (void*)name, fip);
}

/*	GetFunctionInfoFromFuncRef(fref, fip)

	Returns information that you need in order to call an Igor user function
	or an external function from an XOP. This is an advanced feature that most
	XOP programmers will not need.
	
	This function does the same thing as GetFunctionInfo except that the input parameter
	is the value of a FUNCREF field from an Igor structure instead of the name of the
	function.	
	
	Thread Safety: GetFunctionInfoFromFuncRef is thread-safe with Igor Pro 6.20 or later.
	The FunctionInfo structure contains an isThreadSafe field which tells you if
	the function in question is threadsafe.
*/
int
GetFunctionInfoFromFuncRef(FUNCREF fref, FunctionInfoPtr fip)
{
	return (int)CallBack2(GET_FUNCTION_INFO_FROM_FUNCREF, (void*)fref, fip);
}

/*	CheckFunctionForm(fip, requiredNumParameters, requiredParameterTypes, badParameterNumberPtr, returnType)

	Checks the form (number of parameters, types of parameters and return type) of the function
	against the required form.  This is an advanced feature that most XOP programmers will not need.
	
	You must call CheckFunctionForm before calling CallFunction to make sure that the function
	you are calling has the form you expect. Otherwise you may cause a crash.
	
	fip is pointer to a structure set by calling GetFunctionInfo.
	
	requiredNumParameters is the number of parameters you expect the function to have.
	
	requiredParameterTypes is an array in which you have set each value to one of the following:
		NT_FP64							The parameter must be scalar numeric
		NT_FP64 | NT_CMPLX				The parameter must be complex numeric
		HSTRING_TYPE					The parameter must be a string
		WAVE_TYPE						The parameter must be a scalar numeric wave
		WAVE_TYPE | NT_CMPLX			The parameter must be a complex numeric wave
		TEXT_WAVE_TYPE					The parameter must be a text wave
		DATAFOLDER_TYPE					The parameter must be a DFREF
		FV_FUNC_TYPE					The parameter must be a function reference
		FV_STRUCT_TYPE | FV_REF_TYPE	The parameter must be a structure

	DFREF parameters are supported in Igor Pro 6.10 or later.
		
	The number of elements in requiredParameterTypes must be at least equal to requiredNumParameters. 

	If the parameter must be pass-by-reference, use FV_REF_TYPE in addition to the
	values above. For example:
		NT_FP64 | FV_REF_TYPE
		NT_FP64 | NT_CMPLX | FV_REF_TYPE
		HSTRING_TYPE | FV_REF_TYPE
	
	Pass-by-reference is applicable to numeric and string parameters only.
	
	If you do not want CheckFunctionForm to check a particular parameter, pass
	-1 in the corresponding element of the requiredParameterTypes array.
	
	If the function is an external function, if the function takes a wave parameter,
	there is no way to know if the external function expects a numeric wave or a text wave.
	Consequently, CheckFunctionForm does not distinguish between numeric and text
	waves for external functions. The external function itself is supposed to check
	the type of the wave passed in at runtime and return an error if it is the
	wrong type.
	
	returnType is the required return type of the function which must be one of
	the following:
		NT_FP64
		NT_FP64 | NT_CMPLX
		HSTRING_TYPE
		WAVE_TYPE
		DATAFOLDER_TYPE
		
	Wave reference and DFREF return types are supported in Igor Pro 6.10 or later.

	If you do not want CheckFunctionForm to check the return type, pass -1
	as the returnType parameter.
	
	Sets *badParameterNumberPtr to the zero-based index of the first parameter that does
	not match the required type or to -1 if all parameters match the required type.
	
	Returns 0 if the form matches the requirements or an error code if not.
	
	If a function parameter type does not match the required parameter type, the error
	code returned will indicate the type of parameter required but not which parameter
	type was bad. If you want to inform the user more specifically, use the value returned
	via badParameterNumberPtr to select your own more specific error code. If the
	error was a parameter type mismatch, *badParameterNumberPtr will contain the
	zero-based index of the bad parameter. Otherwise, *badParameterNumberPtr will
	contain -1.
	
	Thread Safety: CheckFunctionForm is thread-safe with Igor Pro 6.20 or later.
*/
int
CheckFunctionForm(FunctionInfoPtr fip, int requiredNumParameters, int requiredParameterTypes[], int* badParameterNumberPtr, int returnType)
{
	return (int)CallBack5(CHECK_FUNCTION_FORM, fip, (void*)requiredNumParameters, requiredParameterTypes, badParameterNumberPtr, (void*)returnType);
}

/*	CallFunction(fip, parameters, resultPtr)

	Calls the function identified by fip. fip is a pointer to a FunctionInfo structure
	whose values were set by calling GetFunctionInfo. This is an advanced feature that most
	XOP programmers will not need.

	fip is a pointer to a FunctionInfo structure whose values were set by calling GetFunctionInfo.
	
	fip->compilationIndex is used by CallFunction to make sure that procedures were
	not recompiled after you called GetFunctionInfo. If procedures were recompiled, then
	the information in the structure may be invalid so CallFunction returns BAD_COMPILATION_INDEX.
	
	parameters is a pointer to a structure containing the values that you want
	to pass to the function. These values must agree in type with the
	function's parameter list, as indicated by the parameter information that
	you obtain via GetFunctionInfo. To guarantee this, you must call CheckFunctionForm
	before calling CallFunction.
	
	NOTE: The parameters structure must use standard XOP structure packing, namely,
	two-byte packing. If you don't set the structure packing correctly, a crash
	is likely.
	
	Parameter type values are discussed in detail in the comments for the GetFunctionInfo
	function. Here is the basic correspondence between function parameter types and the
	respective structure field:
		if (parameterType == NT_FP64)
			structure field is double
			
		if (parameterType == (NT_FP64 | NT_CMPLX))
			structure field is double[2]
		
		if (parameterType == HSTRING_TYPE)
			structure field is Handle
		
		if (WAVE_TYPE bit is set)
			structure field is waveHndl
		
		if (DATAFOLDER_TYPE bit is set)			// Igor Pro 6.10 or later
			structure field is DataFolderHandle
		
		if (FV_FUNC_TYPE bit is set)
			structure field is PSInt
	
	NOTE: For pass-by-value strings parameters, ownership of a handle stored in
	the parameter structure is passed to the function when you call CallFunction.
	The function will dispose the handle and CallFunction will set the field
	to NULL. You must not dispose it or otherwise reference it after calling CallFunction.
	
	NOTE: For pass-by-reference string parameters, the handle stored in the parameter
	structure field may be reused or disposed by the called function. When
	CallFunction returns, you own the handle which may be the same handle
	you passed or a different handle. You own this handle and you must dispose
	it when you no longer need it. If the field is NULL, which could occur in the
	event of an error, you must not dispose of it or otherwise access it.
	
	CallFunction stores the function result at the location indicated by
	resultPtr. Here is the correspondence between the function result type and the
	variable pointed to by resultPtr:
		NT_FP64						double
		NT_FP64 | NT_CMPLX			double[2]
		HSTRING_TYPE				Handle
		WAVE_TYPE					waveHndl
		DATAFOLDER_TYPE				DataFolderHandle

	Wave reference and DFREF return types are supported in Igor Pro 6.10 or later.
		
	NOTE: A function that returns a string, wave or DataFolderHandle can return NULL
	instead of a valid handle. You must test the returned value to make sure it is not
	NULL before using it.
	
	NOTE: For string return values, if the returned handle is not NULL, you own it and
	you must dispose it when you no longer need it.
	
	CallFunction returns 0 as the function result if the function executed.
	If the function could not be executed, it returns a non-zero error code.
	
	Thread Safety: CallFunction is thread-safe with Igor Pro 6.20 or later.
	See the XOP manual for further information on its thread-safety.
*/
int
CallFunction(FunctionInfoPtr fip, void* parameters, void* resultPtr)
{
	return (int)CallBack3(CALL_FUNCTION, fip, parameters, resultPtr);
}

/*	Example of Calling a User or External Function

	In this example, we have written our own curve fitting routine, analogous to Igor's
	FuncFit operation, as an external operation. We want to call a user or external function
	from our external operation. The function that we want to call has this form:
		Function FitFunc(w, x)
			Wave w
			Variable x
			
	To simplify the example, we assume that we known the name of the function that
	we want to execute.

	// Define the parameter structure. These are parameters we will pass to user or external function.
	#pragma pack(2)		// All structures passed between Igor and XOP are two-byte aligned.
	struct OurParams {							// Used to pass parameters to the function.
		waveHndl waveH;							// For the first function parameter.
		double x;								// For the second function parameter.
	};
	typedef struct OurParams OurParams;
	typedef struct OurParams* OurParamsPtr;
	#pragma pack()		// Reset structure alignment to default.

	int
	DoOurOperation(waveHndl coefsWaveH)
	{
		FunctionInfo fi;
		OurParams parameters;
		int badParameterNumber;
		int requiredParameterTypes[2];
		double result;
		int i;
		double values[5];
		int err;
		
		// Make sure the function exists and get information about it.
		if (err = GetFunctionInfo("TestFitFunc", &fi))
			return err;
		
		// Make sure the function has the right form.
		requiredParameterTypes[0] = NT_FP64;			// First parameter is numeric
		requiredParameterTypes[1] = WAVE_TYPE;			// Second parameter is a numeric wave.
		if (err = CheckFunctionForm(&fi, 2, requiredParameterTypes, &badParameterNumber, NT_FP64))
			return err;
		
		// We have a valid function. Let's call it.
		
		parameters.x = 0;
		parameters.waveH = coefsWaveH;
		for(i=0; i<5; i+=1) {
			parameters.x = i;
			if (err = CallFunction(&fi, (void*)&parameters, &result))
				return err;
			values[i] = result;
		}
			
		return 0;	
	}
*/

/*	GetIgorCallerInfo(pathOrTitle, linePtr, routineName, callStackHPtr)

	Returns information about the Igor procedure that is currently running. This may be useful
	for debugging purposes in complex Igor projects.
	
	If the currently running procedure is stored in a standalone file then the full path to the
	procedure file is returned in Macintosh HFS format (using colon separators) via pathOrTitle.
	
	If the currently running procedure is in the built-in procedure window or in a packed procedure
	file or in a new procedure file not yet saved to disk then the window's title is returned via
	pathOrTitle preceded by a colon (e.g., ":Procedure").
	
	If no procedure is running then pathOrTitle is set to an empty string ("").
	
	If a procedure is running then *linePtr is set to the zero-based line number within the file.
	When a macro is running this returns the line number of the macro declaration.
	When a user-defined function is running it returns the line number of the statement which
	is currently executing.
	
	If no procedure is running then *linePtr is set to -1.
	
	If a procedure is running then the name of the procedure is returned via routineName. Otherwise
	routineName is set an empty string ("").

	If GetIgorCallerInfo returns zero then the function has succeeded and you own the
	handle referenced by callStackHPtr, and you must dispose it.
	
	If GetIgorCallerInfo returns a non-zero error code then *callStackHPtr is undefined
	and you must not do anything with it.
	
	If a procedure is running, *callStackHPtr is set to a handle containing a stack crawl.
	If no procedure is running, *callStackHPtr is set to NULL. You should always test
	*callStackHPtr  - do not use it if it is NULL. If *callStackHPtr is not NULL then
	this handle belongs to you and you must dispose it when you no longer need it.

	If *callStackHPtr is not NULL it contains a list of Igor procedures separated by semicolons.
	*callStackHPtr is not null terminated. Use GetHandleSize to determine the size of the stack crawl text.
	
	In some cases the currently running procedure can not be determined, for example because
	because Igor procedures need to be compiled. If so then GetIgorCallerInfo acts as if no
	procedure is running.
	
	Thread Safety: GetIgorCallerInfo is not thread-safe.
*/
int
GetIgorCallerInfo(char pathOrTitle[MAX_PATH_LEN+1], int* linePtr, char routineName[256], Handle* callStackHPtr)
{
	if (!CheckRunningInMainThread("GetIgorCallerInfo"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack4(GET_IGOR_CALLER_INFO, pathOrTitle, linePtr, routineName, callStackHPtr);
}

/*	GetIgorRTStackInfo(code, stackInfoHPtr)

	Returns information about the chain of Igor procedure that are currently running.
	This may be useful for debugging purposes in complex Igor projects.
	
	If GetIgorRTStackInfo returns zero then the function has succeeded and you own the handle
	referenced by stackInfoHPtr, and you must dispose it.
	
	If GetIgorRTStackInfo returns a non-zero error code then *stackInfoHPtr is undefined
	and you must not do anything with it.
	
	The contents of the string returned by GetIgorRTStackInfo is the same as the contents
	of the string returned by the built-in Igor function GetRTStackInfo. See the
	documentation for GetRTStackInfo for details.
	
	The string stored in *stackInfoHPtr is not null-terminated. Use GetHandleSize on *stackInfoHPtr
	to determine the size of the text.
	
	Thread Safety: GetIgorRTStackInfo is not thread-safe.
*/
int
GetIgorRTStackInfo(int code, Handle* stackInfoHPtr)
{
	if (!CheckRunningInMainThread("GetIgorRTStackInfo"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack2(GET_IGOR_RT_STACK_INFO, (void*)code, stackInfoHPtr);
}

/*	RegisterOperation(cmdTemplate, runtimeNumVarList, runtimeStrVarList, runtimeParamStructSize, runtimeAddress, options)

	Registers an XOP operation with Igor's Operation Handler.
	
	cmdTemplate specifies the operation name and syntax.
	
	runtimeNumVarList is a semicolon-separated list of numeric variables that the operation sets
	at runtime or NULL if it sets no numeric variables.
	
	runtimeStrVarList is a semicolon-separated list of string variables that the operation sets
	at runtime or NULL if it sets no string variables.
	
	runtimeParamStructSize is the size of the runtime parameter structure for the operation.
	
	runtimeAddress is the address of the ExecuteOperation function for this operation.
	
	options is reserved for future use. Pass zero for this parameter.
	
	Returns 0 or an error code.
	
	Thread Safety: RegisterOperation is not thread-safe.
*/
int
RegisterOperation(const char* cmdTemplate, const char* runtimeNumVarList, const char* runtimeStrVarList, int runtimeParamStructSize, const void* runtimeAddress, int options)
{
	if (!CheckRunningInMainThread("RegisterOperation"))
		return NOT_IN_THREADSAFE;

	return (int)CallBack6(REGISTER_OPERATION, (void*)cmdTemplate, (void*)runtimeNumVarList, (void*)runtimeStrVarList, (void*)runtimeParamStructSize, (void*)runtimeAddress, (void*)options);
}

/*	SetOperationNumVar(varName, dval)

	This is used only to implement an operation using Igor's Operation Handler.
	It is used to store a value in an external operation output numeric variable such as V_flag.
	
	varName must be the name of a numeric variable that you specified via the runtimeNumVarList
	parameter to RegisterOperation.
	
	dval is the value to be stored in the named variable.
	
	If your operation was invoked from the command line, SetOperationNumVar sets
	a global variable in the current data folder, creating it if necessary.
	
	If your operation was invoked from a macro, SetOperationNumVar sets
	a macro local variable.
	
	If your operation was invoked from a user function, SetOperationNumVar sets
	a function local variable.
	
	Returns 0 or an error code.
	
	Thread Safety: SetOperationNumVar is thread-safe with Igor Pro 6.20 or later.
*/
int
SetOperationNumVar(const char* varName, double dval)
{
	return (int)CallBack2(SET_RUNTIME_NUMERIC_VARIABLE, (void*)varName, (void*)&dval);
}

/*	SetOperationStrVar(varName, str)

	This is used only to implement an operation using Igor's Operation Handler.
	It is used to store a value in an external operation output string variable such as S_fileName.
	
	varName must be the name of a string variable that you specified via the runtimeStrVarList
	parameter to RegisterOperation.
	
	str points to the value to be stored in the named variable.
	
	If your operation was invoked from the command line, SetOperationStrVar sets
	a global variable in the current data folder, creating it if necessary.
	
	If your operation was invoked from a macro, SetOperationStrVar sets
	a macro local variable.
	
	If your operation was invoked from a user function, SetOperationStrVar sets
	a function local variable.
	
	Returns 0 or an error code.
	
	Thread Safety: SetOperationStrVar is thread-safe with Igor Pro 6.20 or later.
*/
int
SetOperationStrVar(const char* varName, const char* str)
{
	return (int)CallBack2(SET_RUNTIME_STRING_VARIABLE, (void*)varName, (void*)str);
}

/*	SetOperationStrVar2(varName, data, dataLength)

	This is the same as SetOperationStrVar except for the dataLength parameter which allows
	you to store binary data in an external operation output string variable.
	
	Added for Igor Pro 6.10B02. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
	
	Thread Safety: SetOperationStrVar2 is thread-safe with Igor Pro 6.20 or later.
*/
int
SetOperationStrVar2(const char* varName, const char* data, BCInt dataLength)
{
	return (int)CallBack3(SET_RUNTIME_STRING_VARIABLE_V2, (void*)varName, (void*)data, (void*)dataLength);
}

/*	VarNameToDataType(varName, dataTypePtr)

	This is used only to implement an operation using Igor's Operation Handler.
	It must be called only from your ExecuteOperation function.
	
	This function is used only for operations which take the name of a variable
	as a parameter and then set the value of that variable. An example is Igor's
	Open operation which takes a "refNum" parameter.
	
	After calling VarNameToDataType you would then call StoreNumericDataUsingVarName
	or StoreStringDataUsingVarName.
	
	varName must be a pointer to a varName field in your operation's runtime parameter
	structure when a pointer to this structure has been passed by Igor to your ExecuteOperation
	function. Igor relies on the varName field being set correctly. If your operation was
	called from a user function, the varName field will actually contain binary data
	used by Igor to access function local variables, NVARs and SVARs.

	On output, *dataTypePtr will contain one of the following:
		0					Means varName refers to a string variable or SVAR.
		NT_FP64				Means varName refers to a scalar local variable or NVAR.
		NT_FP64 | NT_CMPLX	Means varName refers to a complex local variable or NVAR.
	
	Returns 0 or an error code.
	
	Thread Safety: VarNameToDataType is thread-safe with Igor Pro 6.20 or later.
*/
int
VarNameToDataType(const char* varName, int* dataTypePtr)
{
	return (int)CallBack2(VAR_NAME_TO_DATA_TYPE, (void*)varName, (void*)dataTypePtr);
}

/*	FetchNumericDataUsingVarName(varName, realPart, imagPart)

	Retrieves data from a numeric variable which may be local or global.

	This is used only to implement an operation using Igor's Operation Handler.
	It must be called only from your ExecuteOperation function.
	
	This function is used only for operations which take the name of a numeric variable
	as a parameter and then get the value of that variable.
	
	varName must be a pointer to a varName field in your operation's runtime parameter
	structure when a pointer to this structure has been passed by Igor to your ExecuteOperation
	function. Igor relies on the varName field being set correctly. If your operation was
	called from a user function, the varName field will actually contain binary data
	used by Igor to access function local variables, NVARs and SVARs.
	
	You should call this routine only after you have determined that varName refers to
	a numeric variable or NVAR. This will be the case if VarNameToDataType returns
	a non-zero number as the data type.
	
	Returns 0 or an error code.
	
	Added for Igor Pro 6.10. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
	
	Thread Safety: FetchNumericDataUsingVarName is thread-safe with Igor Pro 6.20 or later.
*/
int
FetchNumericDataUsingVarName(const char* varName, double* realPartPtr, double* imagPartPtr)
{
	return (int)CallBack3(FETCH_NUMERIC_DATA_USING_VAR_NAME, (void*)varName, (void*)realPartPtr, (void*)imagPartPtr);
}

/*	StoreNumericDataUsingVarName(varName, realPart, imagPart)

	Stores data in a numeric variable which may be local or global.

	This is used only to implement an operation using Igor's Operation Handler.
	It must be called only from your ExecuteOperation function.
	
	This function is used only for operations which take the name of a numeric variable
	as a parameter and then set the value of that variable. An example is Igor's
	Open operation which takes a "refNum" parameter.
	
	varName must be a pointer to a varName field in your operation's runtime parameter
	structure when a pointer to this structure has been passed by Igor to your ExecuteOperation
	function. Igor relies on the varName field being set correctly. If your operation was
	called from a user function, the varName field will actually contain binary data
	used by Igor to access function local variables, NVARs and SVARs.
	
	You should call this routine only after you have determined that varName refers to
	a numeric variable or NVAR. This will be the case if VarNameToDataType returns
	a non-zero number as the data type.
	
	Returns 0 or an error code.
	
	Thread Safety: StoreNumericDataUsingVarName is thread-safe with Igor Pro 6.20 or later.
*/
int
StoreNumericDataUsingVarName(const char* varName, double realPart, double imagPart)
{
	return (int)CallBack3(STORE_NUMERIC_DATA_USING_VAR_NAME, (void*)varName, (void*)&realPart, (void*)&imagPart);
}

/*	FetchStringDataUsingVarName(varName, hPtr)

	Retrieves data from a string variable which may be local or global.
	
	Returns the string data via a new handle and stores the handle in *hPtr.
	On return, if *hPtr is not NULL then you own the new handle and must
	dispose it.

	This is used only to implement an operation using Igor's Operation Handler.
	It must be called only from your ExecuteOperation function.
	
	This function is used only for operations which take the name of a string variable
	as a parameter and then set the value of that variable.
	
	varName must be a pointer to a varName field in your operation's runtime parameter
	structure when a pointer to this structure has been passed by Igor to your ExecuteOperation
	function. Igor relies on the varName field being set correctly. If your operation was
	called from a user function, the varName field will actually contain binary data
	used by Igor to access function local variables, NVARs and SVARs.
	
	You should call this routine only after you have determined that varName refers to
	a string variable or SVAR. This will be the case if VarNameToDataType returns
	zero as the data type.
	
	Returns 0 or an error code.
	
	Added for Igor Pro 6.10. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
	
	Thread Safety: FetchStringDataUsingVarName is thread-safe with Igor Pro 6.20 or later.
*/
int
FetchStringDataUsingVarName(const char* varName, Handle* hPtr)
{
	return (int)CallBack2(FETCH_STRING_DATA_USING_VAR_NAME, (void*)varName, (void*)hPtr);
}

/*	StoreStringDataUsingVarName(varName, buf, len)

	Stores data in a string variable which may be local or global.

	This is used only to implement an operation using Igor's Operation Handler.
	It must be called only from your ExecuteOperation function.
	
	This function is used only for operations which take the name of a string variable
	as a parameter and then set the value of that variable. An example is Igor's
	FReadLine operation which takes a "stringVarName" parameter.
	
	varName must be a pointer to a varName field in your operation's runtime parameter
	structure when a pointer to this structure has been passed by Igor to your ExecuteOperation
	function. Igor relies on the varName field being set correctly. If your operation was
	called from a user function, the varName field will actually contain binary data
	used by Igor to access function local variables, NVARs and SVARs.
	
	You should call this routine only after you have determined that varName refers to
	a string variable or SVAR. This will be the case if VarNameToDataType returns
	zero as the data type.
	
	Returns 0 or an error code.
	
	Thread Safety: StoreStringDataUsingVarName is thread-safe with Igor Pro 6.20 or later.
*/
int
StoreStringDataUsingVarName(const char* varName, const char* buf, BCInt len)
{
	return (int)CallBack3(STORE_STRING_DATA_USING_VAR_NAME, (void*)varName, (void*)buf, (void*)len);
}

/*	GetOperationWaveRef(tp, dfH, name, waveRefIdentifier, destWaveHPtr)
	
	Returns via destWaveHPtr a reference to the destination wave if it was specified using
	a wave reference or NULL otherwise.

	Except in rare cases, you should not call this function. Use GetOperationDestWave instead.

	This routine supports the /DEST flag for an operation.

	If the destination wave was specified by passing a wave reference to /DEST,
	this routine returns the wave handle via destWaveHPtr.
	
	If the destination wave was not specified via a wave reference, this routine returns NULL
	via destWaveHPtr and you must create the destination wave yourself.
	
	The rule is: If a simple name or the name of a structure wave field are passed to /DEST,
	the specified wave is returned via destWaveHPtr. If $<string> or a path are passed to /DEST,
	NULL is returned. 
	
	Here are typical cases where a wave reference specifies the destination wave:
		SampleOp /DEST=jack	
		Wave w = root:jack
		SampleOp /DEST=w
	In these cases, GetOperationWaveRef would return a handle to wave jack via destWaveHPtr.
	In the first command, Operation Handler creates an automatic wave reference in the calling
	user-defined function if you are running Igor Pro 6.20 or later.
	
	Here are typical cases where a wave reference is not used to specify the destination wave:
		String name = "jack"
		SampleOp /DEST=$name		// Dest specified by $<string>
		SampleOp /DEST=root:jack	// Dest specified by path
		SampleOp /DEST=:jack		// Dest specified by path
	In these cases, there is no wave reference for the destination wave and
	GetOperationWaveRef would return NULL via destWaveHPtr.
	
	The user can specify that a wave reference be created by appending /WAVE after
	the /DEST flag, like this:
		String name = "jack"
		SampleOp /DEST=$name/WAVE=wDest
		SampleOp /DEST=root:jack/WAVE=wDest
		SampleOp /DEST=:jack/WAVE=wDest
	"/WAVE=wDest" is an "inline WAVE reference statement".
	An inline WAVE reference does not determine which wave is the destination wave but
	merely provides a place for the operation to return a reference to the destination wave.

	If a wave reference for the destination wave exists, it must be set to refer
	to the destination wave once that wave has been determined or created.
	You set the wave reference by calling SetOperationWaveRef later in the process.
	
	See GetOperationDestWave for an example.

	Added for Igor Pro 6.20. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
	
	Thread Safety: GetOperationWaveRef is thread-safe with Igor Pro 6.20 or later.
*/
int
GetOperationWaveRef(DataFolderHandle dfH, const char* name, int destWaveRefIdentifier, waveHndl* destWaveHPtr)
{
	return (int)CallBack4(GET_OPERATION_WAVE_REF, dfH, (void*)name, (void*)destWaveRefIdentifier, destWaveHPtr);
}

/*	GetOperationDestWave(dfH, name, waveRefIdentifier, options, dimensionSizes, dataType, destWaveHPtr, destWaveCreatedPtr)

	Returns via destWaveHPtr a reference to the destination wave specified by the
	/DEST flag in an operation.

	GetOperationDestWave optionally creates the wave if it does not already exist. If the
	wave does exist, GetOperationDestWave changes it as specified by the options parameter.

	GetOperationDestWave is typically called to support the /DEST flag for an operation.
	It determines which wave is to be the destination wave, optionally creates it if
	it does not already exist, optionally overwrites it or sets its dimensions and
	data type if it does already exist, and returns the wave handle via destWaveHPtr.
	
	In Igor Pro 6.20 or later, GetOperationDestWave respects wave references
	passed as the destination wave when called from a user-defined function.

	If GetOperationDestWave succeeds, it returns 0 as the function result and
	a non-NULL wave handle via destWaveHPtr. If it fails it returns a non-zero error
	code as the function result and NULL via destWaveHPtr. Possible errors include
	NOMEM, NO_TEXTNUMERIC_WAVE_OVERWRITE, and NAME_WAV_CONFLICT (the destination wave
	already exists and overwrite is not specified).
	
	destWaveRefIdentifier is used by GetOperationDestWave to access a user-defined function
	local wave reference for the destination wave, if such a wave reference was specified
	with the /DEST flag. The example below shows how to obtain destWaveRefIdentifier from
	the operation's runtime parameter structure.
	
	The options parameter is a bitwise combination of the following flags:
		kOpDestWaveOverwriteOK:				Permits returning a pre-existing wave
		kOpDestWaveChangeExistingWave		Data type and dimensions of pre-existing wave will be changed
											but other properties (wave scaling, wave note) will be preserved
		kOpDestWaveOverwriteExistingWave	Pre-existing wave is completely overwritten with a new, pristine wave
		kOpDestWaveMakeFreeWave				Dest wave will be a new free wave (requires Igor Pro 6.20 or later)
		kOpDestWaveMustAlreadyExist			Dest wave must already exist
	
	If kOpDestWaveOverwriteOK is set, the destination wave is returned if it
	already exists. If kOpDestWaveOverwriteOK is not set an error is returned
	if the destination wave already exists.
	
	If kOpDestWaveChangeExistingWave is set, if the destination wave already exists
	its dimensions and data type are changed based on the dimensionSizes and dataType
	parameters but other properties (wave scaling, wave note, etc.) are left unchanged.
	
	If kOpDestWaveOverwriteExistingWave is set, if the destination wave already exists
	it is completely overwritten with a new pristine wave with default properties.
	
	The choice between kOpDestWaveChangeExistingWave and kOpDestWaveOverwriteExistingWave
	must be made on a case-by-case basis. If you are converting an existing operation to
	use GetOperationDestWave you should choose between kOpDestWaveChangeExistingWave
	and kOpDestWaveOverwriteExistingWave so that the behavior of your operation will not change.
	For most new applications, kOpDestWaveOverwriteExistingWave is preferrable since it makes
	the output of the operation independent of pre-existing conditions.
	
	If kOpDestWaveChangeExistingWave and kOpDestWaveOverwriteExistingWave are both cleared,
	then GetOperationDestWave will not change the dimensions and data type of a pre-existing
	destination wave. In this case, it is up to you to make sure that the dimensions and data
	type of the returned wave are appropriate and change them if necessary.

	If kOpDestWaveChangeExistingWave and kOpDestWaveOverwriteExistingWave are both set,
	this is a programmer error and an error is returned.
	
	If your operation has a /FREE flag, signifying that the destination wave should be
	a free wave, set the kOpDestWaveMakeFreeWave bit. If kOpDestWaveMakeFreeWave is set,
	GetOperationDestWave makes a new free wave as the destination wave and returns it
	via destWaveHPtr and all other option flag are ignored. This requires Igor Pro
	6.20 or later.
	
	If kOpDestWaveMakeFreeWave is not set, this procedure is followed: If the
	destination wave was specified via a wave reference by the calling user-defined
	function, GetOperationDestWave returns that wave via destWaveHPtr. Otherwise it
	looks for an existing wave specified by dfH and name. If there is an existing wave,
	GetOperationDestWave reuses it. Otherwise it creates a new wave using dfH and name.
	The rules for obtaining the destination wave via a wave reference are explained
	in the documentation for GetOperationWaveRef.
	
	If kOpDestWaveMustAlreadyExist is set, then, if the destination wave does not
	already exist, the NOWAV ("expected wave name") error is returned. You may want to
	test for NOWAV and return a more specific error code. If you set kOpDestWaveMustAlreadyExist
	you must also set kOpDestWaveOverwriteOK because if kOpDestWaveOverwriteOK is
	cleared GetOperationDestWave returns an error if the destination wave exists.
	
	If you want to return an existing wave only without overwriting or changing it,
	pass (kOpDestWaveMustAlreadyExist | kOpDestWaveOverwriteOK) and leave all other
	bits cleared.
	
	dimensionSizes is an array of wave dimension sizes as shown in the example below.
	If kOpDestWaveMustAlreadyExist is set and kOpDestWaveChangeExistingWave,
	kOpDestWaveOverwriteExistingWave and kOpDestWaveMakeFreeWave are cleared then
	dimensionSizes is not used and you can pass NULL.
	
	dataType is normally a standard Igor data type (e.g., NT_FP64, (NT_FP32 | NT_CMPLX),
	TEXT_WAVE_TYPE). However, if you pass NULL for dimensionSizes, you must pass -1 for
	dataType, indicating that you do not want to change the data type.
	
	If destWaveCreatedPtr is not NULL, *destWaveCreatedPtr is set to 1 if GetOperationDestWave
	created a new wave or to 0 if GetOperationDestWave reused an existing wave.
	
	In the following example, the operation is defined by this Operation Handler syntax:
		SampleOp /DEST=DataFolderAndName:{dest,real}
	This results in an operation runtime parameter structure containing:
		int DESTFlagEncountered;
		DataFolderAndName dest;
		int DESTFlagParamsSet[1];
	
	Here is typical usage for an operation with a /DEST flag. In this example, p is a pointer
	to the operation runtime parameter structure:
	
		waveHndl destWaveH;						// Destination wave handle.
		int destWaveRefIdentifier;				// Used to access the destination wave if specified using a wave reference
		char destWaveName[MAX_OBJ_NAME+1];
		DataFolderHandle dfH;
		int dataType;
		CountInt dimensionSizes[MAX_DIMENSIONS+1];
		int destWaveCreated;
		int options;
		int err;
		
		destWaveH = NULL;
		destWaveRefIdentifier = 0;

		strcpy(destWaveName, "W_SampleOp");			// Default dest wave name
		dfH = NULL;									// Default is current data folder
		
		dataType = NT_FP64;
		MemClear(dimensionSizes, sizeof(dimensionSizes));
		dimensionSizes[ROWS] = 100;
		
		if (p->DESTFlagEncountered) {
			strcpy(destWaveName, p->name);
			dfH = p->dfH;
			// If a wave reference was used, p->destParamsSet[0] contains information that
			// GetOperationDestWave uses to find the address of the wave reference.
			destWaveRefIdentifier = p->destParamsSet[0];
		}
		
		options = kOpDestWaveOverwriteOK | kOpDestWaveOverwriteExistingWave;
		
		err = GetOperationDestWave(tp, dfH, destWaveName, destWaveRefIdentifier, options, dimensionSizes, dataType, &destWaveH, &destWaveCreated);
		if (err != 0)
			return err;
		
		<Store output data in dest wave>
		
		// Set wave reference to refer to destination wave.
		if (destWaveRefIdentifier != 0)
			SetOperationWaveRef(tp, destWaveH, destWaveRefIdentifier);

	The call to SetOperationWaveRef sets the automatic destination wave reference created
	by Operation Handler when you use the DataFolderAndName syntax shown above in your
	operation template.

	GetOperationDestWave was added for Igor Pro 6.20. If you call this with an earlier
	version of Igor, it will still return a destination wave but without recognizing
	any wave referenced passed via the /DEST flag. For example:
		Wave w = root:wave0
		SampleOp /O /DEST=w
	In Igor Pro 6.20 and later, this overwrites a wave named wave0 in the root
	data folder. However, prior to Igor Pro 6.20, it creates a wave named w
	in the current data folder. This is how destination waves worked in XOPs
	prior to Igor Pro 6.20.
	
	If kOpDestWaveMakeFreeWave is set in the options parameter and you are running
	with a version of Igor prior to 6.20, GetOperationDestWave returns IGOR_OBSOLETE.
	
	Thread Safety: GetOperationDestWave is thread-safe with Igor Pro 6.20 or later.
*/
int
GetOperationDestWave(DataFolderHandle dfH, const char* name, int destWaveRefIdentifier, int options, CountInt dimensionSizes[], int dataType, waveHndl* destWaveHPtr, int* destWaveCreatedPtr)
{
	if (igorVersion >= 620)
		return (int)CallBack8(GET_OPERATION_DEST_WAVE, dfH, (void*)name, (void*)destWaveRefIdentifier, (void*)options, dimensionSizes, (void*)dataType, destWaveHPtr, destWaveCreatedPtr);

	{
		// This version of Igor does not support GetOperationDestWave. Drop back to the
		// pre-620 technique of just creating the destination wave using dfH and name
		// without regard to any destination wave reference.
		int overwriteOK = (options & kOpDestWaveOverwriteOK) != 0;
		int changeExistingWave = (options & kOpDestWaveChangeExistingWave) != 0;
		int overwriteExistingWave = (options & kOpDestWaveOverwriteExistingWave) != 0;
		int makeFree = (options & kOpDestWaveMakeFreeWave) != 0;
		int destWaveMustAlreadyExist = (options & kOpDestWaveMustAlreadyExist) != 0;
		int destWavePreExisted = 0;
		int callMakeWave = 0;			// If set we call MDMakeWave rather than MDChangeWave
		int result = 0;
		
		*destWaveHPtr = NULL;
		
		if (destWaveCreatedPtr != NULL)
			*destWaveCreatedPtr = 0;
		
		if (makeFree)
			return IGOR_OBSOLETE;					// Requires Igor Pro 6.20 or later.
			
		if (changeExistingWave && overwriteExistingWave)
			return GENERAL_BAD_VIBS;				// Programmer error.
		
		if (dimensionSizes==NULL && dataType!=-1)
			return GENERAL_BAD_VIBS;				// Programmer error
			
		*destWaveHPtr = FetchWaveFromDataFolder(dfH, name);		// Get existing dest wave, if any.
		if (*destWaveHPtr != NULL)
			destWavePreExisted = 1;
		
		if (!destWavePreExisted && destWaveMustAlreadyExist)
			return NOWAV;		
		
		if (destWavePreExisted && !overwriteOK) {
			*destWaveHPtr = NULL;
			return NAME_WAV_CONFLICT;
		}
		
		if (destWavePreExisted && !changeExistingWave && !overwriteExistingWave)
			return 0;	// Caller just wants to get existing wave and we have it.

		callMakeWave = !destWavePreExisted;						// Truth that we want to call MDMakeWave rather than MDChangeWave
		if (destWavePreExisted && overwriteExistingWave) {
			*destWaveHPtr = NULL;
			callMakeWave = 1;									// Force call of MDMakeWave to get pristine wave
			changeExistingWave = 0;
		}

		if (callMakeWave || changeExistingWave) {
			if (dimensionSizes == NULL) {
				*destWaveHPtr = NULL;
				return GENERAL_BAD_VIBS;
			}
		}
		
		if (callMakeWave) {
			// Make dest wave.
			result = MDMakeWave(destWaveHPtr, name, dfH, dimensionSizes, dataType, overwriteOK);
			if (result==0 && !destWavePreExisted && destWaveCreatedPtr!=NULL)
				*destWaveCreatedPtr = 1;
		}
		else {
			// Dest wave pre-existed
			if (changeExistingWave)
				result = MDChangeWave(*destWaveHPtr, dataType, dimensionSizes);
		}

		return result;
	}
}

/*	SetOperationWaveRef(waveH, waveRefIndentifier)

	Sets a wave reference in an Igor user-defined function to refer to a destination wave
	that your operation created.

	This is used only to implement an operation using Igor's Operation Handler.
	It must be called only from your ExecuteOperation function.
	
	This function is used only for operations which use a DataFolderAndName parameter
	and declare it as a destination wave parameter using this kind of syntax in the
	operation template:
	
		SampleOp DataFolderAndName:{dest,<real> or <complex> or <text>}
		
	When you use this syntax and your operation is compiled in a user-defined function,
	the Igor function compiler may automatically create a wave reference in the function
	for the destination wave. However the automatically created wave reference will be
	NULL until you set it by calling SetOperationWaveRef.
	
	The SetOperationWaveRef callback sets the automatically created wave reference to
	refer to a specific wave, namely the wave that you created in response to the
	DataFolderAndName parameter. You must call it after successfully creating the
	destination wave.
	
	Igor will create an automatic wave reference only if the operation is called
	from a user-defined function and only if the destination wave is specified using
	a simple name (e.g., wave0 but not root:wave0 or $destWave). You have no way to know
	whether an automatic wave reference was created so you must call SetOperationWaveRef
	in all cases. SetOperationWaveRef will do nothing if no automatic wave reference
	exists.	
	
	The waveRefIndentifier parameter allows Igor to determine where in memory the wave
	reference is stored. This information is passed to your XOP in the "ParamSet" field
	of your ExecuteOperation structure.
	
	Your code should look something like this:
	
		// In your RuntimeParams structure
		DataFolderAndName dest;
		int destParamsSet[1];
	
		// In your ExecuteOperation function
		destWaveH = NULL;
		err = MDMakeWave(&destWaveH, p->dest.name, p->dest.dfH, dimensionSizes, type, overwrite);
		if (destWaveH != NULL) {
			int waveRefIndentifier = p->destParamsSet[0];
			err = SetOperationWaveRef(destWaveH, waveRefIndentifier);
		}
	
	Returns 0 or an error code.
	
	Thread Safety: SetOperationWaveRef is thread-safe with Igor Pro 6.20 or later.
*/
int
SetOperationWaveRef(waveHndl waveH, int waveRefIndentifier)
{
	return (int)CallBack2(SET_OPERATION_WAVE_REF, waveH, (void*)waveRefIndentifier);
}

/*	CalcWaveRange(wrp)

	This is used with Operation Handler WaveRange parameters.
	
	wrp is a pointer to a WaveRangeRec.
	
	See declaration of WaveRangeRec structure for details.
	
	Returns error code if error calculating range or 0.
	
	Thread Safety: CalcWaveRange is thread-safe with Igor Pro 6.20 or later.
*/
int
CalcWaveRange(WaveRangeRecPtr wrp)
{
	return (int)CallBack1(CALCWAVERANGE, wrp);
}

/*	DateToIgorDateInSeconds(numValues, year, month, dayOfMonth, secs)
	
	Converts dates into Igor date format (seconds since 1/1/1904).
	
	numValues is the number of dates to convert.
	
	year, month and dayOfMonth and secs are arrays allocated by the calling routine.
	The size of each array is specified by numValues.
	
	On input, year, month and dayOfMonth hold the input values.
	On return, secs holds the output values.
	
	The function result is zero or an error code.
	
	This routine requires Igor Pro 5 or later. Earlier versions will return IGOR_OBSOLETE.
	
	Thread Safety: DateToIgorDateInSeconds is thread-safe with Igor Pro 6.20 or later.
*/
int
DateToIgorDateInSeconds(CountInt numValues, short* year, short* month, short* dayOfMonth, double* secs)
{
	return (int)CallBack5(DATE_TO_IGOR_DATE, (void*)numValues, year, month, dayOfMonth, secs);
}

/*	IgorDateInSecondsToDate(numValues, secs, dates)
	
	Converts dates in Igor date format (seconds since 1/1/1904) into date records.
	
	numValues is the number of Igor dates to convert.
	
	secs is an array of dates in Igor date format. Its size is specified by numValues.

	dates is an array of shorts. It must hold 7*numValues shorts. For each input value,
	7 values are written to the dates array, in the following order:
		year, month, dayOfMonth, hour, minute, second, dayOfWeek		
	
	The function result is zero or an error code.
	
	Example:
		double secs[2];
		short dates[2*7];
		int err;
		
		secs[0] = 0;			// Represents January 1, 1904, midnight.
		secs[1] = 24*60*60;		// Represents January 2, 1904, midnight.
		err = IgorDateInSecondsToDate(2, secs, dates);
	
	This routine requires Igor Pro 5 or later. Earlier versions will return IGOR_OBSOLETE.
	
	Thread Safety: IgorDateInSecondsToDate is thread-safe with Igor Pro 6.20 or later.
*/
int
IgorDateInSecondsToDate(CountInt numValues, double* secs, short* dates)
{
	return (int)CallBack3(IGOR_DATE_TO_DATE, (void*)numValues, secs, dates);
}

/*	GetNVAR(nvp, realPartPtr, imagPartPtr, numTypePtr)

	Retrieves the data and type of a global numeric variable referenced by an NVAR
	field in an Igor Pro structure.
	
	nvp is a pointer to an NVAR field in an Igor structure.
	The Igor structure pointer would be passed into the XOP as a parameter.

	If GetNVAR returns 0 then *realPartPtr will be the contents of the real part of
	the global variable and *imagPartPtr will be the contents of the imaginary part of
	the global variable, if it is complex.
	
	realPartPtr and imagPartPtr must each point to storage for a double whether
	the global variable is complex or not.
	
	*numTypePtr is set to the numeric type of the global. This will be either NT_FP64
	or (NT_FP64 | NT_CMPLX).
	
	Returns 0 or an error code.
	
	Thread Safety: GetNVAR is thread-safe with Igor Pro 6.20 or later.
*/
int
GetNVAR(const NVARRec* nvp, double* realPartPtr, double* imagPartPtr, int* numTypePtr)
{
	return (int)CallBack4(GET_NVAR, (void*)nvp, realPartPtr, imagPartPtr, numTypePtr);
}

/*	SetNVAR(nvp, realPartPtr, imagPartPtr)

	Sets the value of a global numeric variable referenced by an NVAR field
	in an Igor Pro structure.

	nvp is a pointer to an NVAR field in an Igor structure.
	The Igor structure pointer would be passed into the XOP as a parameter.

	*realPartPtr is the value to store in the real part of the global variable and
	*imagPartPtr is the value to store in the imaginary part of the global variable,
	if it is complex.
	
	realPartPtr and imagPartPtr must each point to storage for a double whether the global
	variable is complex or not.

	Returns 0 or an error code.
	
	Thread Safety: SetNVAR is thread-safe with Igor Pro 6.20 or later.
*/
int
SetNVAR(const NVARRec* nvp, const double* realPartPtr, const double* imagPartPtr)
{
	return (int)CallBack3(SET_NVAR, (void*)nvp, (void*)realPartPtr, (void*)imagPartPtr);
}

/*	GetSVAR(svp, strHPtr)

	Retrieves the handle for a global string variable referenced by an SVAR field
	in an Igor Pro structure.

	svp is a pointer to an SVAR field in an Igor structure.
	The Igor structure pointer would be passed into the XOP as a parameter.

	If GetSVAR returns 0 then *strHPtr will be the handle for an Igor global string variable.

	NOTE:	*strHPtr can be NULL if the global string variable contains no characters.

	NOTE:	*strHPtr belongs to Igor. Do not dispose it or alter it in any way.
			Use this function only to read the contents of the string.
		
	Returns 0 or an error code.
	
	Thread Safety: GetSVAR is thread-safe with Igor Pro 6.20 or later.
*/
int
GetSVAR(const SVARRec* svp, Handle* strHPtr)
{
	return (int)CallBack2(GET_SVAR, (void*)svp, strHPtr);
}

/*	SetSVAR(svp, strH)

	Sets the value of a global string variable referenced by an SVAR field
	in an Igor Pro structure.

	svp is a pointer to an SVAR field in an Igor structure.
	The Igor structure pointer would be passed into the XOP as a parameter.

	strH is a handle that you own. It can be NULL to set the global string variable
	to empty.

	NOTE:	Igor copies the contents of strH. You retain ownership of it and
			must dispose it.
	
	Returns 0 or an error code.
	
	Thread Safety: SetSVAR is thread-safe with Igor Pro 6.20 or later.
*/
int
SetSVAR(const SVARRec* svp, Handle strH)
{
	return (int)CallBack2(SET_SVAR, (void*)svp, strH);
}
