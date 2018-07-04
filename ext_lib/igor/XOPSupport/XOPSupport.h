/*	XOPSupport.h

	Declares routines, global variables and other items needed to use XOPSupport routines.
*/

#ifdef __cplusplus
extern "C" {						// This allows C++ to call the XOPSupport routines.
#endif

// Global variables used by XOPSupport.c.
extern IORecHandle XOPRecHandle;
extern int XOPModified;
extern int igorVersion;

// NaN represents a missing value.

#ifdef MACIGOR
	// Visual C++ does not accept this syntax.
	#define SINGLE_NAN ((float)NAN)
	#define DOUBLE_NAN ((double)NAN)
#endif

#ifdef WINIGOR
	// Apple's ProjectBuilder (GNU C) does not accept this syntax.
	static unsigned char NAN_BYTES[] = {0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0x7F};
	#define SINGLE_NAN *(float *)NAN_BYTES
	#define DOUBLE_NAN *(double *)NAN_BYTES
#endif

/*	"Private" routines.
	The XOPSupport files use these to call Igor. You should not call them directly.
*/
XOPIORecResult CallBack0(int message);
XOPIORecResult CallBack1(int message, void* item0);
XOPIORecResult CallBack2(int message, void* item0, void* item1);
XOPIORecResult CallBack3(int message, void* item0, void* item1, void* item2);
XOPIORecResult CallBack4(int message, void* item0, void* item1, void* item2, void* item3);
XOPIORecResult CallBack5(int message, void* item0, void* item1, void* item2, void* item3, void* item4);
XOPIORecResult CallBack6(int message, void* item0, void* item1, void* item2, void* item3, void* item4, void* item5);

// Notice Routines (in XOPSupport.c).
void XOPNotice(const char* noticePtr);
void XOPNotice2(const char* noticePtr, UInt32 options);
void XOPNotice3(const char *noticePtr, const char* rulerName, UInt32 options);
void XOPResNotice(int strListID, int index);

// Utility routines (in XOPSupport.c).
// Capitalize was removed from XOP Toolkit 6. Use CmpStr instead of calling Capitalize and strcmp.
int CmpStr(const char* str1, const char* str2);
const char* strchr2(const char* str, int ch);
const char* strrchr2(const char* str, int ch);
int EscapeSpecialCharacters(const char* input, int inputLength, char* output, int outputBufferSize, int* numCharsOutPtr);
int UnEscapeSpecialCharacters(const char* input, int inputLength, char* output, int outputBufferSize, int* numCharsOutPtr);
void MemClear(void* p, BCInt n);
// MoveLockHandle is obsolete and was removed from XOP Toolkit 6. Remove all calls to MoveLockHandle, HLock, HGetState, HSetState and HUnlock from your source code.
int GetCStringFromHandle(Handle h, char* str, int maxChars);
int PutCStringInHandle(const char* str, Handle h);
int CheckAbort(TickCountInt timeOutTicks);
int IsNaN32(float* floatPtr);
int IsNaN64(double* doublePtr);
void SetNaN32(float* fPtr);
void SetNaN64(double* dPtr);
int IsINF32(float* floatPtr);
int IsINF64(double* doublePtr);
// int IgorVersion(void);	// HR, 091118: Made static in XOP Toolkit 6. Use the igorVersion global instead.
void XOPInit(IORecHandle ioRecHandle);
int RunningInMainThread(void);
int CheckRunningInMainThread(const char* routineName);
void SetXOPType(int type);
void SetXOPEntry(void (*entryPoint)(void));
void SetXOPResult(XOPIORecResult result);
XOPIORecResult GetXOPResult(void);
void SetXOPMessage(int message);
int GetXOPMessage(void);
// SetXOPRefCon is obsolete and was removed from XOP Toolkit 6.
// GetXOPRefCon is obsolete and was removed from XOP Toolkit 6.
int GetXOPStatus(void);
XOPIORecParam GetXOPItem(int itemNumber);
void IgorError(const char* title, int errCode);
int GetIgorErrorMessage(int errCode, char errorMessage[256]);
int WinInfo(int index, int typeMask, char* name, XOP_WINDOW_REF* windowRefPtr);

// Numeric conversion utilities (in XOPNumericConversion.c).
#define SIGNED_INT 1
#define UNSIGNED_INT 2
#define IEEE_FLOAT 3

void DoubleToFloat(const double* inPtr, float* outPtr, CountInt numValues);
void DoubleToSInt32(const double* inPtr, SInt32* outPtr, CountInt numValues);
void DoubleToShort(const double* inPtr, short* outPtr, CountInt numValues);
void DoubleToByte(const double* inPtr, char* outPtr, CountInt numValues);
void DoubleToUInt32(const double* inPtr, UInt32* outPtr, CountInt numValues);
void DoubleToUnsignedShort(const double* inPtr, unsigned short* outPtr, CountInt numValues);
void DoubleToUnsignedByte(const double* inPtr, unsigned char* outPtr, CountInt numValues);
int ConvertDouble(const double* src, void* dest, CountInt numValues, int destFormat, int destBytes);

void FloatToDouble(const float* inPtr, double* outPtr, CountInt numValues);
void FloatToSInt32(const float* inPtr, SInt32* outPtr, CountInt numValues);
void FloatToShort(const float* inPtr, short* outPtr, CountInt numValues);
void FloatToByte(const float* inPtr, char* outPtr, CountInt numValues);
void FloatToUInt32(const float* inPtr, UInt32* outPtr, CountInt numValues);
void FloatToUnsignedShort(const float* inPtr, unsigned short* outPtr, CountInt numValues);
void FloatToUnsignedByte(const float* inPtr, unsigned char* outPtr, CountInt numValues);
int ConvertFloat(const float* src, void* dest, CountInt numValues, int destFormat, int destBytes);

void SInt32ToDouble(const SInt32* inPtr, double* outPtr, CountInt numValues);
void SInt32ToFloat(const SInt32* inPtr, float* outPtr, CountInt numValues);
void SInt32ToShort(const SInt32* inPtr, short* outPtr, CountInt numValues);
void SInt32ToByte(const SInt32* inPtr, char* outPtr, CountInt numValues);
void SInt32ToUInt32(const SInt32* inPtr, UInt32* outPtr, CountInt numValues);
void SInt32ToUnsignedShort(const SInt32* inPtr, unsigned short* outPtr, CountInt numValues);
void SInt32ToUnsignedByte(const SInt32* inPtr, unsigned char* outPtr, CountInt numValues);
int ConvertSInt32(const SInt32* src, void* dest, CountInt numValues, int destFormat, int destBytes);

void ShortToDouble(const short* inPtr, double* outPtr, CountInt numValues);
void ShortToFloat(const short* inPtr, float* outPtr, CountInt numValues);
void ShortToSInt32(const short* inPtr, SInt32* outPtr, CountInt numValues);
void ShortToByte(const short* inPtr, char* outPtr, CountInt numValues);
void ShortToUInt32(const short* inPtr, UInt32* outPtr, CountInt numValues);
void ShortToUnsignedShort(const short* inPtr, unsigned short* outPtr, CountInt numValues);
void ShortToUnsignedByte(const short* inPtr, unsigned char* outPtr, CountInt numValues);
int ConvertShort(const short* src, void* dest, CountInt numValues, int destFormat, int destBytes);

void ByteToDouble(const char* inPtr, double* outPtr, CountInt numValues);
void ByteToFloat(const char* inPtr, float* outPtr, CountInt numValues);
void ByteToSInt32(const char* inPtr, SInt32* outPtr, CountInt numValues);
void ByteToShort(const char* inPtr, short* outPtr, CountInt numValues);
void ByteToUInt32(const char* inPtr, UInt32* outPtr, CountInt numValues);
void ByteToUnsignedShort(const char* inPtr, unsigned short* outPtr, CountInt numValues);
void ByteToUnsignedByte(const char* inPtr, unsigned char* outPtr, CountInt numValues);
int ConvertByte(const char* src, void* dest, CountInt numValues, int destFormat, int destBytes);

void UInt32ToDouble(const UInt32* inPtr, double* outPtr, CountInt numValues);
void UInt32ToFloat(const UInt32* inPtr, float* outPtr, CountInt numValues);
void UInt32ToSInt32(const UInt32* inPtr, SInt32* outPtr, CountInt numValues);
void UInt32ToShort(const UInt32* inPtr, short* outPtr, CountInt numValues);
void UInt32ToByte(const UInt32* inPtr, char* outPtr, CountInt numValues);
void UInt32ToUnsignedShort(const UInt32* inPtr, unsigned short* outPtr, CountInt numValues);
void UInt32ToUnsignedByte(const UInt32* inPtr, unsigned char* outPtr, CountInt numValues);
int ConvertUInt32(const UInt32* src, void* dest, CountInt numValues, int destFormat, int destBytes);

void UnsignedShortToDouble(const unsigned short* inPtr, double* outPtr, CountInt numValues);
void UnsignedShortToFloat(const unsigned short* inPtr, float* outPtr, CountInt numValues);
void UnsignedShortToSInt32(const unsigned short* inPtr, SInt32* outPtr, CountInt numValues);
void UnsignedShortToShort(const unsigned short* inPtr, short* outPtr, CountInt numValues);
void UnsignedShortToByte(const unsigned short* inPtr, char* outPtr, CountInt numValues);
void UnsignedShortToUInt32(const unsigned short* inPtr, UInt32* outPtr, CountInt numValues);
void UnsignedShortToUnsignedByte(const unsigned short* inPtr, unsigned char* outPtr, CountInt numValues);
int ConvertUnsignedShort(const unsigned short* src, void* dest, CountInt numValues, int destFormat, int destBytes);

void UnsignedByteToDouble(const unsigned char* inPtr, double* outPtr, CountInt numValues);
void UnsignedByteToFloat(const unsigned char* outPtr, float* fPtr, CountInt numValues);
void UnsignedByteToSInt32(const unsigned char* inPtr, SInt32* outPtr, CountInt numValues);
void UnsignedByteToShort(const unsigned char* inPtr, short* outPtr, CountInt numValues);
void UnsignedByteToByte(const unsigned char* inPtr, char* outPtr, CountInt numValues);
void UnsignedByteToUInt32(const unsigned char* inPtr, UInt32* outPtr, CountInt numValues);
void UnsignedByteToUnsignedShort(const unsigned char* inPtr, unsigned short* outPtr, CountInt numValues);
int ConvertUnsignedByte(const unsigned char* src, void* dest, CountInt numValues, int destFormat, int destBytes);

int NumTypeToNumBytesAndFormat(int numType, int* numBytesPerPointPtr, int* dataFormatPtr, int* isComplexPtr);
int NumBytesAndFormatToNumType(int numBytesPerValue, int dataFormat, int* numTypePtr);
void FixByteOrder(void* p, int bytesPerPoint, CountInt count);
int ConvertData(const void* src, void* dest, CountInt numValues, int srcBytes, int srcFormat, int destBytes, int destFormat);
int ConvertData2(const void* src, void* dest, CountInt numValues, int srcDataType, int destDataType);

void ScaleDouble(double* dPtr, double* offset, double* multiplier, CountInt numValues);
void ScaleFloat(float* fPtr, double* offset, double* multiplier, CountInt numValues);
void ScaleSInt32(SInt32* iPtr, double* offset, double* multiplier, CountInt numValues);
void ScaleShort(short* iPtr, double* offset, double* multiplier, CountInt numValues);
void ScaleByte(char* iPtr, double* offset, double* multiplier, CountInt numValues);
void ScaleUInt32(UInt32* iPtr, double* offset, double* multiplier, CountInt numValues);
void ScaleUnsignedShort(unsigned short* iPtr, double* offset, double* multiplier, CountInt numValues);
void ScaleUnsignedByte(unsigned char* iPtr, double* offset, double* multiplier, CountInt numValues);
void ScaleData(int dataType, void* dataPtr, double* offsetPtr, double* multiplierPtr, CountInt numValues);
void ScaleClipAndRoundData(int dataType, void* dataPtr, CountInt numValues, double offset, double multiplier, double dMin, double dMax, int doRound);

// Wave access routines (in XOPWaveAccess.c).
int FetchNumericValue(int type, const char* dataStartPtr, IndexInt index, double value[2]);
int StoreNumericValue(int type, char* dataStartPtr, IndexInt index, double value[2]);
void WaveHandleModified(waveHndl waveHandle);
// WaveHandlesModified is obsolete and was removed from XOP Toolkit 6. Use WaveHandleModified.
void WaveModified(const char* waveName);
waveHndl FetchWave(const char* waveName);
waveHndl FetchWaveFromDataFolder(DataFolderHandle dataFolderH, const char* waveName);
int WaveType(waveHndl waveHandle);
CountInt WavePoints(waveHndl waveHandle);
void WaveName(waveHndl waveHandle, char* namePtr);
void* WaveData(waveHndl waveHandle);
void WaveScaling(waveHndl waveHandle, double* hsAPtr, double* hsBPtr, double* topPtr, double* botPtr);
void SetWaveScaling(waveHndl waveHandle, const double* hsAPtr, const double* hsBPtr, const double* topPtr, const double* botPtr);
void WaveUnits(waveHndl waveHandle, char* xUnits, char* dataUnits);
void SetWaveUnits(waveHndl waveHandle, const char* xUnits, const char* dataUnits);
Handle WaveNote(waveHndl waveHandle);
void SetWaveNote(waveHndl waveHandle, Handle noteHandle);
TickCountInt WaveModDate(waveHndl waveH);
int WaveLock(waveHndl waveH);
int SetWaveLock(waveHndl waveH, int lockState);
int WaveModState(waveHndl waveH);
int WaveModCount(waveHndl waveH);
// GetWavesInfo is obsolete and was removed from XOP Toolkit 6. Use WaveType, WavePoints and WaveData as needed.
// long GetWavesInfo(waveHndl waves[], int numWaves, int waveTypes[], long wavePoints[], int waveStates[], void* wavePtrs[]);
// SetWavesStates is obsolete and was removed from XOP Toolkit 6. It is obsolete and no longer needed.
// void SetWavesStates(waveHndl waves[], int numWaves, int waveStates[]);
int MakeWave(waveHndl* waveHandlePtr, const char* waveName, CountInt numPoints, int type, int overwrite);
int ChangeWave(waveHndl waveHandle, CountInt numPoints, int type);
int KillWave(waveHndl waveHandle);

// Data folder access routines (in XOPDataFolderAccess.c).
int GetDataFolderNameOrPath(DataFolderHandle dataFolderH, int flags, char dataFolderPathOrName[MAXCMDLEN+1]);
int GetDataFolderIDNumber(DataFolderHandle dataFolderH, int* IDNumberPtr);
int GetDataFolderProperties(DataFolderHandle dataFolderH, int* propertiesPtr);
int SetDataFolderProperties(DataFolderHandle dataFolderH, int properties);
int GetDataFolderListing(DataFolderHandle dataFolderH, int optionsFlag, Handle h);
int GetRootDataFolder(int refNum, DataFolderHandle* rootFolderHPtr);
int GetCurrentDataFolder(DataFolderHandle* currentFolderHPtr);
int SetCurrentDataFolder(DataFolderHandle dataFolderH);
int GetNamedDataFolder(DataFolderHandle startingDataFolderH, const char dataFolderPath[MAXCMDLEN+1], DataFolderHandle* dataFolderHPtr);
int GetDataFolderByIDNumber(int IDNumber, DataFolderHandle* dataFolderHPtr);
int GetParentDataFolder(DataFolderHandle dataFolderH, DataFolderHandle* parentFolderHPtr);
int GetNumChildDataFolders(DataFolderHandle parentDataFolderH, int* numChildDataFolderPtr);
int GetIndexedChildDataFolder(DataFolderHandle parentDataFolderH, int index, DataFolderHandle* childDataFolderHPtr);
int GetWavesDataFolder(waveHndl waveH, DataFolderHandle* dataFolderHPtr);
int NewDataFolder(DataFolderHandle parentFolderH, const char newDataFolderName[MAX_OBJ_NAME+1], DataFolderHandle* newDataFolderHPtr);
int KillDataFolder(DataFolderHandle dataFolderH);
int DuplicateDataFolder(DataFolderHandle sourceDataFolderH, DataFolderHandle parentDataFolderH, const char newDataFolderName[MAX_OBJ_NAME+1]);
int MoveDataFolder(DataFolderHandle sourceDataFolderH, DataFolderHandle newParentDataFolderH);
int RenameDataFolder(DataFolderHandle dataFolderH, const char newName[MAX_OBJ_NAME+1]);
int GetNumDataFolderObjects(DataFolderHandle dataFolderH, int objectType, int* numObjectsPtr);
int GetIndexedDataFolderObject(DataFolderHandle dataFolderH, int objectType, int index, char objectName[MAX_OBJ_NAME+1], DataObjectValuePtr objectValuePtr);
int GetDataFolderObject(DataFolderHandle dataFolderH, const char objectName[MAX_OBJ_NAME+1], int* objectTypePtr, DataObjectValuePtr objectValuePtr);
int SetDataFolderObject(DataFolderHandle dataFolderH, const char objectName[MAX_OBJ_NAME+1], int objectType, DataObjectValuePtr objectValuePtr);
int KillDataFolderObject(DataFolderHandle dataFolderH, int objectType, const char objectName[MAX_OBJ_NAME+1]);
int MoveDataFolderObject(DataFolderHandle sourceDataFolderH, int objectType, const char objectName[MAX_OBJ_NAME+1], DataFolderHandle destDataFolderH);
int RenameDataFolderObject(DataFolderHandle dataFolderH, int objectType, const char objectName[MAX_OBJ_NAME+1], const char newObjectName[MAX_OBJ_NAME+1]);
int DuplicateDataFolderObject(DataFolderHandle dataFolderH, int objectType, const char objectName[MAX_OBJ_NAME+1], DataFolderHandle destFolderH, const char newObjectName[MAX_OBJ_NAME+1], int overwrite);
void ClearDataFolderFlags(void);
int GetDataFolderChangesCount(void);
int GetDataFolderChangeFlags(DataFolderHandle dataFolderH, int* flagsP);
int HoldDataFolder(DataFolderHandle dfH, DataFolderHandle* dfRefPtr);	// Added for Igor Pro 6.20B04.
int ReleaseDataFolder(DataFolderHandle* dfRefPtr);						// Added for Igor Pro 6.20B04.

// Multi-dimension wave access routines (in XOPWaveAccess.c).
int MDMakeWave(waveHndl* waveHPtr, const char* waveName, DataFolderHandle dataFolderH, CountInt dimSizes[MAX_DIMENSIONS+1], int type, int overwrite);
int MDGetWaveDimensions(waveHndl waveH, int* numDimensionsPtr, CountInt dimSizes[MAX_DIMENSIONS+1]);
int MDChangeWave(waveHndl waveH, int dataType, CountInt dimSizes[MAX_DIMENSIONS+1]);
int MDChangeWave2(waveHndl waveH, int dataType, CountInt dimSizes[MAX_DIMENSIONS+1], int mode);
int MDGetWaveScaling(waveHndl waveH, int dimension, double* sfA, double* sfB);
int MDSetWaveScaling(waveHndl waveH, int dimension, const double* sfA, const double* sfB);
int MDGetWaveUnits(waveHndl waveH, int dimension, char units[MAX_UNIT_CHARS+1]);
int MDSetWaveUnits(waveHndl waveH, int dimension, const char units[MAX_UNIT_CHARS+1]);
int MDGetDimensionLabel(waveHndl waveH, int dimension, IndexInt element, char label[MAX_DIM_LABEL_CHARS+1]);
int MDSetDimensionLabel(waveHndl waveH, int dimension, IndexInt element, char label[MAX_DIM_LABEL_CHARS+1]);
int MDAccessNumericWaveData(waveHndl waveH, int accessMode, BCInt* dataOffsetPtr);
int MDGetNumericWavePointValue(waveHndl waveH, IndexInt indices[MAX_DIMENSIONS], double value[2]);
int MDSetNumericWavePointValue(waveHndl waveH, IndexInt indices[MAX_DIMENSIONS], double value[2]);
int MDGetDPDataFromNumericWave(waveHndl waveH, double* dPtr);
int MDStoreDPDataInNumericWave(waveHndl waveH, const double* dPtr);
int MDGetTextWavePointValue(waveHndl waveH, IndexInt indices[MAX_DIMENSIONS], Handle textH);
int MDSetTextWavePointValue(waveHndl waveH, IndexInt indices[MAX_DIMENSIONS], Handle textH);

// Command line routines (in XOPSupport.c).
int XOPCommand(const char* cmdPtr);
int XOPSilentCommand(const char* cmdPtr);
int XOPCommand2(const char *cmdPtr, int silent, int sendToHistory);
int XOPCommand3(const char *cmdPtr, int silent, int sendToHistory, Handle* historyTextHPtr);
void PutCmdLine(const char* cmd, int mode);
void FinishDialogCmd(const char* cmd, int mode);

// Variable access routines (in XOPSupport.c).
int FetchNumVar(const char* varName, double* doublePtr1, double* doublePtr2);
int StoreNumVar(const char* varName, const double* doublePtr1, const double* doublePtr2);
int FetchStrVar(const char* varName, char* stringPtr);
Handle FetchStrHandle(const char* varName);
int StoreStrVar(const char* varName, const char* stringPtr);
int Variable(const char* varName, int varType);
int VariableList(Handle listHandle, const char* match, const char* sep, int varTypeCode);
int StringList(Handle listHandle, const char* match, const char* sep);
int SetIgorIntVar(const char* numVarName, int value, int forceGlobal);
int SetIgorFloatingVar(const char* numVarName, const double* valuePtr, int forceGlobal);
int SetIgorComplexVar(const char* numVarName, const double* realValuePtr, const double* imagValuePtr, int forceGlobal);
int SetIgorStringVar(const char* stringVarName, const char* stringVarValue, int forceGlobal);

// Command parsing routines (in XOPSupport.c).
/*	These routines are obsolete and were removed in XOP Toolkit 6.
	See "XOP Toolkit 6 Upgrade Notes" in the XOP Toolkit 6 manual for details.

	NextSymb, GetSymb, IsStringExpression
	GetFlag, GetFlagNum, GetTrueOrFalseFlag
	GetAString, GetAStringInHandle, GetFormat
	GetName, Keyword, GetKeyword
	GetNum, GetNum2, GetLong
	GetWaveName, GetWave, GetWaveList, GetWaveRange
	GetNumVarName, GetStrVarName
	GetDataFolder, GetDataFolderAndName
	AtEndOfCommand, CheckTerm
*/

// Name utilities (in XOPSupport.c).
int UniqueName(const char* baseName, char* finalName);
int UniqueName2(int nameSpaceCode, const char* baseName, char* finalName, int* suffixNumPtr);
int SanitizeWaveName(char* waveName, int column);
int CheckName(DataFolderHandle dataFolderH, int objectType, const char* name);
int PossiblyQuoteName(char* name);
void CatPossiblyQuotedName(char* str, const char* name);
int CleanupName(int beLiberal, char* name, int maxNameChars);
int CreateValidDataObjectName(DataFolderHandle dataFolderH, const char* inName, char* outName, int* suffixNumPtr, int objectType, int beLiberal, int allowOverwrite, int inNameIsBaseName, int printMessage, int* nameChangedPtr, int* doOverwritePtr);

// As of XOP Toolkit 6, this is no longer supported.
// int DoCHIO(CHIORecPtr CHIOPtr);

// Utilities for XOPs with menu items (in XOPMenus.c).
#ifdef MACIGOR		// The Windows versions of these routines are implemented inside Igor and are declared in XOPWinMacSupport.h.
	void WMDrawMenuBar(void);
	void WMDeleteMenu(short menuID);
	void WMInsertMenu(MenuHandle menuH, short beforeID);
	// void ReleaseMenu(MenuHandle menuH);		// HR, 010427: This name ReleaseMenu was usurped by Apple. It is no longer supplied by the XOP Toolkit.
#endif
MenuHandle WMGetMenu(short resourceID);

int	ResourceToActualMenuID(int resourceMenuID);
MenuHandle ResourceMenuIDToMenuHandle(int resourceMenuID);
int	ActualToResourceMenuID(int menuID);
int ActualToResourceItem(int igorMenuID, int actualItemNumber);
int ResourceToActualItem(int igorMenuID, int resourceItemNumber);
int SetIgorMenuItem(int message, int enable, const char* text, int param);
void FillMenu(MenuHandle theMenu, const char* itemList, int itemListLen, int afterItem);
void FillMenuNoMeta(MenuHandle theMenu, const char* itemList, int itemListLen, int afterItem);
int FillWaveMenu(MenuHandle theMenu, const char* match, const char* options, int afterItem);
int FillPathMenu(MenuHandle theMenu, const char* match, const char* options, int afterItem);
int FillWinMenu(MenuHandle theMenu, const char* match, const char* options, int afterItem);
void WMDeleteMenuItems(MenuHandle theMenu, int afterItem);	// HR, 010427: This was previously called DeleteMenuItems but Apple usurped that name.

// Utilities for XOPs with windows (in XOPWindows.c).
XOP_WINDOW_REF GetActiveWindowRef(void);
int IsXOPWindowActive(XOP_WINDOW_REF w);
void ShowXOPWindow(XOP_WINDOW_REF w);
void HideXOPWindow(XOP_WINDOW_REF w);
void ShowAndActivateXOPWindow(XOP_WINDOW_REF w);
void HideAndDeactivateXOPWindow(XOP_WINDOW_REF w);
void SetXOPWindowTitle(XOP_WINDOW_REF windowRef, const char* title);
void GetXOPWindowPositionAndState(XOP_WINDOW_REF windowRef, Rect* r, int* winStatePtr);
void SetXOPWindowPositionAndState(XOP_WINDOW_REF windowRef, Rect* r, int winState);
void TransformWindowCoordinates(int mode, double coords[4]);
void GetXOPWindowIgorPositionAndState(XOP_WINDOW_REF windowRef, double coords[4], int* winStatePtr);
void SetXOPWindowIgorPositionAndState(XOP_WINDOW_REF windowRef, double coords[4], int winState);
int TellIgorWindowStatus(XOP_WINDOW_REF windowRef, int status, int options);

// Utilities for XOPs with text windows (in XOPWindows.c).
// TUNew is obsolete and was removed from XOP Toolkit 6. Use TUNew2.
// Handle TUNew(WindowPtr winPtr, Rect* borderRectPtr, int font, int size, int crOnly);
int TUNew2(const char* winTitle, const Rect* winRectPtr, Handle* TUPtr, XOP_WINDOW_REF* windowRefPtr);
void TUDispose(TUStuffHandle TU);
void TUDisplaySelection(TUStuffHandle TU);
void TUGrow(TUStuffHandle TU, int size);
void TUDrawWindow(TUStuffHandle TU);
void TUUpdate(TUStuffHandle TU);
void TUFind(TUStuffHandle TU, int code);
void TUReplace(TUStuffHandle TU);
void TUIndentLeft(TUStuffHandle TU);
void TUIndentRight(TUStuffHandle TU);
void TUClick(TUStuffHandle TU, EventRecord* eventPtr);
void TUActivate(TUStuffHandle TU, int flag);
void TUIdle(TUStuffHandle TU);
void TUNull(TUStuffHandle TU, EventRecord* eventPtr);
void TUCopy(TUStuffHandle TU);
void TUCut(TUStuffHandle TU);
void TUPaste(TUStuffHandle TU);
void TUClear(TUStuffHandle TU);
void TUKey(TUStuffHandle TU, EventRecord* eventPtr);
void TUInsert(TUStuffHandle TU, const char* dataPtr, int dataLen);
void TUDelete(TUStuffHandle TU);
// TUSetSelect is obsolete and was removed from XOP Toolkit 6. Use TUSetSelLocs instead.
// void TUSetSelect(TUStuffHandle TU, int start, int end);
void TUSelectAll(TUStuffHandle TU);
void TUUndo(TUStuffHandle TU);
void TUPrint(TUStuffHandle TU);
void TUFixEditMenu(TUStuffHandle TU);
void TUFixFileMenu(TUStuffHandle TU);
// TUGetText is obsolete and was removed from XOP Toolkit 6. Use TUFetchParagraphText instead.
// Handle TUGetText(TUStuffHandle TU);
// TUFetchText is obsolete and was removed from XOP Toolkit 6. Use TUFetchParagraphText instead.
// void TUFetchText(TUStuffHandle TU, int offset, int numChars, char* buffer);
// TULength is obsolete and was removed from XOP Toolkit 6. Use TUGetDocInfo instead.
// int TULength(TUStuffHandle TU);
int TULines(TUStuffHandle TU);
// TUSelStart is obsolete and was removed from XOP Toolkit 6. Use TUGetSelLocs instead.
// int TUSelStart(TUStuffHandle TU);
// TUSelEnd is obsolete and was removed from XOP Toolkit 6. Use TUGetSelLocs instead.
// int TUSelEnd(TUStuffHandle TU);
// TUSelectionLength is obsolete and was removed from XOP Toolkit 6. Use TUGetSelLocs instead.
// int TUSelectionLength(TUStuffHandle TU);
// TUInsertFile is obsolete and was removed from XOP Toolkit 6. There is no direct replacement.
// int TUInsertFile(TUStuffHandle TU, const char* fileName, int wdRefNum);
// TUWriteFile is obsolete and was removed from XOP Toolkit 6. There is no direct replacement.
// int TUWriteFile(TUStuffHandle TU, const char* fileName, int wdRefNum, int allFlag);
int TUSFInsertFile(TUStuffHandle TU, const char* prompt, OSType fileTypes[], int numTypes);
int TUSFWriteFile(TUStuffHandle TU, const char* prompt, OSType fileType, int allFlag);
void TUPageSetupDialog(TUStuffHandle TU);
int TUGetDocInfo(TUStuffHandle TU, TUDocInfoPtr dip);
int TUGetSelLocs(TUStuffHandle TU, TULocPtr startLocPtr, TULocPtr endLocPtr);
int TUSetSelLocs(TUStuffHandle TU, TULocPtr startLocPtr, TULocPtr endLocPtr, int flags);
int TUFetchParagraphText(TUStuffHandle TU, int paragraph,  Ptr* textPtrPtr, int* lengthPtr);
int TUFetchSelectedText(TUStuffHandle TU, Handle* textHandlePtr, void* reservedForFuture, int flags);
int TUFetchText2(TUStuffHandle TU, TULocPtr startLocPtr, TULocPtr endLocPtr, Handle* textHandlePtr, void* reservedForFuture, int flags);
int TUSetStatusArea(TUStuffHandle TU, const char* message, int eraseFlags, int statusAreaWidth);
void TUMoveToPreferredPosition(TUStuffHandle TU);
void TUMoveToFullSizePosition(TUStuffHandle TU);
void TURetrieveWindow(TUStuffHandle TU);

// Utilities for accessing the history area
void HistoryDisplaySelection(void);
void HistoryInsert(const char* dataPtr, int dataLen);
void HistoryDelete(void);
int HistoryLines(void);
int HistoryGetSelLocs(TULocPtr startLocPtr, TULocPtr endLocPtr);
int HistorySetSelLocs(TULocPtr startLocPtr, TULocPtr endLocPtr, int flags);
int HistoryFetchParagraphText(int paragraph,  Ptr* textPtrPtr, int* lengthPtr);
int HistoryFetchText(TULocPtr startLocPtr, TULocPtr endLocPtr, Handle* textHPtr);

// Cross-platform dialog routines (in XOPDialogsMac.c and XOPDialogsWin.c)
void XOPEmergencyAlert(const char* message);
void XOPOKAlert(const char* title, const char* message);
int XOPOKCancelAlert(const char* title, const char* message);
int XOPYesNoAlert(const char* title, const char* message);
int XOPYesNoCancelAlert(const char* title, const char* message);
void SetDialogBalloonHelpID(int balloonHelpID);
void GetDBox(XOP_DIALOG_REF theDialog, int itemID, Rect *box);
void HiliteDControl(XOP_DIALOG_REF theDialog, int itemID, int enable);
void DisableDControl(XOP_DIALOG_REF theDialog, int itemID);
void EnableDControl(XOP_DIALOG_REF theDialog, int itemID);
void SetRadBut(XOP_DIALOG_REF theDialog, int firstID, int lastID, int theID);
int ToggleCheckBox(XOP_DIALOG_REF theDialog, int itemID);
int GetCheckBox(XOP_DIALOG_REF theDialog, int itemID);
int SetCheckBox(XOP_DIALOG_REF theDialog, int itemID, int val);
void DisplayDialogCmd(XOP_DIALOG_REF theDialog, int itemID, const char* cmd);
#ifdef MACIGOR
	CGrafPtr SetDialogPort(XOP_DIALOG_REF theDialog);
#endif
#ifdef WINIGOR
	XOP_DIALOG_REF SetDialogPort(XOP_DIALOG_REF theDialog);		// This is a NOP on Windows.
#endif
void ShowDialogWindow(XOP_DIALOG_REF theDialog);
int GetRadBut(XOP_DIALOG_REF theDialog, int itemID);
int GetDText(XOP_DIALOG_REF theDialog, int theItem, char *theText);
void SetDText(XOP_DIALOG_REF theDialog, int theItem, const char *theText);
int GetDInt(XOP_DIALOG_REF theDialog, int theItem, int *theInt);
void SetDInt(XOP_DIALOG_REF theDialog, int theItem, int theInt);
int GetDInt64(XOP_DIALOG_REF theDialog, int theItem, SInt64* valPtr);		// Replaces GetDLong in XOP Toolkit 6
void SetDInt64(XOP_DIALOG_REF theDialog, int theItem, SInt64 val);			// Replaces SetDLong in XOP Toolkit 6
int GetDDouble(XOP_DIALOG_REF theDialog, int theItem, double *theDouble);
void SetDDouble(XOP_DIALOG_REF theDialog, int theItem, double *theDouble);
void SelEditItem(XOP_DIALOG_REF theDialog, int itemID);
void SelMacEditItem(XOP_DIALOG_REF theDialog, int itemID);
int ItemIsPopMenu(XOP_DIALOG_REF theDialog, int itemID);
void InitPopMenus(XOP_DIALOG_REF theDialog);
int CreatePopMenu(XOP_DIALOG_REF theDialog, int popupItemNumber, int titleItemNumber, const char* itemList, int initialItem);	// Added in Igor Pro 3.1 but works with any version.
#ifdef MACIGOR
	// This routine is not supported on Windows.
	MenuHandle GetPopMenuHandle(XOP_DIALOG_REF theDialog, int itemID);
#endif
void GetPopMenu(XOP_DIALOG_REF theDialog, int itemID, int *selItem, char *selStr);
int SetPopMatch(XOP_DIALOG_REF theDialog, int itemID, const char *selStr);
void SetPopItem(XOP_DIALOG_REF theDialog, int itemID, int theItem);
void KillPopMenus(XOP_DIALOG_REF theDialog);
void AddPopMenuItems(XOP_DIALOG_REF theDialog, int itemID, const char *itemList);
void FillPopMenu(XOP_DIALOG_REF theDialog, int itemID, const char *itemList, int itemListLen, int afterItem);
int FillWavePopMenu(XOP_DIALOG_REF theDialog, int itemID, const char *match, const char *options, int afterItem);
int FillPathPopMenu(XOP_DIALOG_REF theDialog, int itemID, const char *match, const char *options, int afterItem);
int FillWindowPopMenu(XOP_DIALOG_REF theDialog, int itemID, const char *match, const char *options, int afterItem);
void DeletePopMenuItems(XOP_DIALOG_REF theDialog, int itemID, int afterItem);

// Cross-platform file handling routines (in XOPFiles.c).
#ifdef MACIGOR
	int HFSToPosixPath(const char* hfsPath, char posixPath[MAX_PATH_LEN+1], int isDirectory);
#endif
int XOPCreateFile(const char* fullFilePath, int overwrite, int macCreator, int macFileType);
int XOPDeleteFile(const char* fullFilePath);
int XOPOpenFile(const char* fullFilePath, int readOrWrite, XOP_FILE_REF* fileRefPtr);
int XOPCloseFile(XOP_FILE_REF fileRef);
int XOPReadFile(XOP_FILE_REF fileRef, UInt32 count, void* buffer, UInt32* numBytesReadPtr);
int XOPReadFile2(XOP_FILE_REF fileRef, UInt32 count, void* buffer, UInt32* numBytesReadPtr);
int XOPReadFile64(XOP_FILE_REF fileRef, SInt64 count, void* buffer, SInt64* numBytesReadPtr);
int XOPWriteFile(XOP_FILE_REF fileRef, UInt32 count, const void* buffer, UInt32* numBytesWrittenPtr);
int XOPWriteFile64(XOP_FILE_REF fileRef, SInt64 count, const void* buffer, SInt64* numBytesWrittenPtr);
int XOPGetFilePosition(XOP_FILE_REF fileRef, UInt32* filePosPtr);
int XOPSetFilePosition(XOP_FILE_REF fileRef, SInt32 filePos, int mode);
int XOPGetFilePosition2(XOP_FILE_REF fileRef, SInt64* dFilePosPtr);
int XOPSetFilePosition2(XOP_FILE_REF fileRef, SInt64 dFilePos);
int XOPAtEndOfFile(XOP_FILE_REF fileRef);
int XOPNumberOfBytesInFile(XOP_FILE_REF fileRef, UInt32* numBytesPtr);
int XOPNumberOfBytesInFile2(XOP_FILE_REF fileRef, SInt64* dNumBytesPtr);
int XOPReadLine(XOP_FILE_REF fileRef, char* buffer, UInt32 bufferLength, UInt32* numBytesReadPtr);
int FullPathPointsToFile(const char* fullPath);
int FullPathPointsToFolder(const char* fullPath);
int WinToMacPath(char path[MAX_PATH_LEN+1]);
int MacToWinPath(char path[MAX_PATH_LEN+1]);
int GetNativePath(const char* filePathIn, char filePathOut[MAX_PATH_LEN+1]);
int EscapeBackslashesInUNCVolumeName(char macFilePath[MAX_PATH_LEN+1]);
int GetDirectoryAndFileNameFromFullPath(const char* fullFilePath, char dirPath[MAX_PATH_LEN+1], char fileName[MAX_FILENAME_LEN+1]);
int GetLeafName(const char* filePath, char name[MAX_FILENAME_LEN+1]);
int GetFullPathFromSymbolicPathAndFilePath(const char* symbolicPathName, const char filePath[MAX_PATH_LEN+1], char fullFilePath[MAX_PATH_LEN+1]);
int ConcatenatePaths(const char* pathIn1, const char* nameOrPathIn2, char pathOut[MAX_PATH_LEN+1]);
int XOPOpenFileDialog(const char* prompt, const char* fileFilterStr, int* fileIndexPtr, const char* initialDir, char filePath[MAX_PATH_LEN+1]);
int XOPSaveFileDialog(const char* prompt, const char* fileFilterStr, int* fileIndexPtr, const char* initialDir, const char* defaultExtensionStr, char filePath[MAX_PATH_LEN+1]);
int ParseFilePath(int mode, const char* pathIn, const char* separator, int whichEnd, int whichElement, char pathOut[MAX_PATH_LEN+1]);	// Added for Igor Pro 6.20B03
int SpecialDirPath(const char* pathID, int domain, int flags, int createDir, char pathOut[MAX_PATH_LEN+1]);								// Added for Igor Pro 6.20B03

// File loader utilities (in XOPFiles.c).
// FileLoaderGetOperationFlags2 is obsolete and was removed from XOP Toolkit 6.
// int FileLoaderGetOperationFlags2(const char* additionalFlagsChars, int* flagsPtr, char* baseName, char symbolicPathName[MAX_OBJ_NAME+1]);
int FileLoaderMakeWave(int column, char* waveName, CountInt numPoints, int fileLoaderFlags, waveHndl* waveHandlePtr);
int SetFileLoaderOutputVariables(const char* fileNameOrPath, int numWavesLoaded, const char* waveNames);
int SetFileLoaderOperationOutputVariables(int runningInUserFunction, const char* fileNameOrPath, int numWavesLoaded, const char* waveNames);

// Data loading and saving utilities for internal WaveMetrics use only (used by WaveMetrics Browser). (In XOPFiles.c).
struct LoadDataInfo;		// GNU C requires this.
struct LoadFileInfo;		// GNU C requires this.
struct SaveDataInfo;		// GNU C requires this.
int PrepareLoadIgorData(struct LoadDataInfo* ldiPtr, int* refNumPtr, struct LoadFileInfo*** topFIHPtr);
int LoadIgorData(struct LoadDataInfo* ldiPtr, int refNum, struct LoadFileInfo** topFIH, DataFolderHandle destDataFolderH);
int EndLoadIgorData(struct LoadDataInfo* ldiPtr, int refNum, struct LoadFileInfo** topFIH);
int SaveIgorData(struct SaveDataInfo* sdiPtr, DataFolderHandle topDataFolderH);

// IGOR color table routines (in XOPSupport.c).
int GetIndexedIgorColorTableName(int index, char name[MAX_OBJ_NAME+1]);
int GetNamedIgorColorTableHandle(const char* name, IgorColorTableHandle* ictHPtr);
int GetIgorColorTableInfo(IgorColorTableHandle ictH, char name[MAX_OBJ_NAME+1], int* numColorsPtr);
int GetIgorColorTableValues(IgorColorTableHandle ictH, int startColorIndex, int endColorIndex, int updatePixelValues, IgorColorSpec* csPtr);

// Cross-Platform Utilities (in XOPSupport.c).
void WinRectToMacRect(const RECT* wr, Rect* mr);
void MacRectToWinRect(const Rect *mr, RECT *wr);

// Miscellaneous routines (in XOPSupport.c).
void XOPBeep(void);
void GetXOPIndString(char* text, int strID, int index);
void ArrowCursor(void);
void IBeamCursor(void);
void WatchCursor(void);
void HandCursor(void); 
void SpinCursor(void);
int SpinProcess(void);
int DoUpdate(void);
void PauseUpdate(int* savePtr);
void ResumeUpdate(int* savePtr);
int WaveList(Handle listHandle, const char* match, const char* sep, const char* options);
int WinList(Handle listHandle, const char* match, const char* sep, const char* options);
int PathList(Handle listHandle, const char* match, const char* sep, const char* options);
int GetPathInfo2(const char* pathName, char fullDirPath[MAX_PATH_LEN+1]);
struct NamedFIFO** GetNamedFIFO(const char* name);
void MarkFIFOUpdated(struct NamedFIFO** fifo);
int SaveXOPPrefsHandle(Handle prefsHandle);
int GetXOPPrefsHandle(Handle* prefsHandlePtr);
int GetPrefsState(int* prefsStatePtr);
int XOPDisplayHelpTopic(const char* title, const char* topicStr, int flags);
enum CloseWinAction DoWindowRecreationDialog(char* procedureName);
int GetIgorProcedureList(Handle* hPtr, int flags);
int GetIgorProcedure(const char* procedureName, Handle* hPtr, int flags);
int SetIgorProcedure(const char* procedureName, Handle h, int flags);
int XOPSetContextualHelpMessage(XOP_WINDOW_REF theWindow, const char* message, const Rect* r);

int GetFunctionInfo(const char* name, FunctionInfoPtr fip);
int CheckFunctionForm(struct FunctionInfo* fip, int requiredNumParameters, int requiredParameterTypes[], int* badParameterNumberPtr, int returnType);
int CallFunction(struct FunctionInfo* fip, void* parameters, void* resultPtr);
int GetFunctionInfoFromFuncRef(FUNCREF fref, FunctionInfoPtr fip);
int GetIgorCallerInfo(char pathOrTitle[MAX_PATH_LEN+1], int* linePtr, char routineName[256], Handle* callStackHPtr);
int GetIgorRTStackInfo(int code, Handle* stackInfoHPtr);	// Added for Igor Pro 6.02B01.

int RegisterOperation(const char* cmdTemplate, const char* runtimeNumVarList, const char* runtimeStrVarList, int runtimeParamStructSize, const void* runtimeAddress, int options);
int SetOperationNumVar(const char* varName, double dval);
int SetOperationStrVar(const char* varName, const char* str);
int SetOperationStrVar2(const char* varName, const char* data, BCInt dataLength);				// Added for Igor Pro 6.10B02.
int VarNameToDataType(const char* varName, int* dataTypePtr);
int FetchNumericDataUsingVarName(const char* varName, double* realPartPtr, double* imagPartPtr);// Added for Igor Pro 6.10.
int StoreNumericDataUsingVarName(const char* varName, double realPart, double imagPart);
int FetchStringDataUsingVarName(const char* varName, Handle* hPtr);								// Added for Igor Pro 6.10.
int StoreStringDataUsingVarName(const char* varName, const char* buf, BCInt len);
int GetOperationWaveRef(DataFolderHandle dfH, const char* name, int destWaveRefIdentifier, waveHndl* destWaveHPtr);	// Added in Igor Pro 6.20
int GetOperationDestWave(DataFolderHandle dfH, const char* name, int destWaveRefIdentifier, int options, CountInt dimensionSizes[], int dataType, waveHndl* destWaveHPtr, int* destWaveCreatedPtr);	// Added in Igor Pro 6.20
int SetOperationWaveRef(waveHndl waveH, int waveRefIndentifier);
int CalcWaveRange(WaveRangeRecPtr wrp);
int HoldWave(waveHndl waveH, waveHndl* waveRefPtr);												// Added for Igor Pro 6.20.
int ReleaseWave(waveHndl* waveRefPtr);															// Added for Igor Pro 6.20.

int DateToIgorDateInSeconds(CountInt numValues, short* year, short* month, short* dayOfMonth, double* secs);
int IgorDateInSecondsToDate(CountInt numValues, double* secs, short* dates);

int GetNVAR(const NVARRec* nvp, double* realPartPtr, double* imagPartPtr, int* numTypePtr);
int SetNVAR(const NVARRec* nvp, const double* realPartPtr, const double* imagPartPtr);
int GetSVAR(const SVARRec* nvp, Handle* strHPtr);
int SetSVAR(const SVARRec* nvp, Handle strH);

int GetTextWaveData(waveHndl waveH, int mode, Handle* textDataHPtr);
int SetTextWaveData(waveHndl waveH, int mode, Handle textDataH);
int GetWaveDimensionLabels(waveHndl waveH, Handle dimLabelsHArray[MAX_DIMENSIONS]);
int SetWaveDimensionLabels(waveHndl waveH, Handle dimLabelsHArray[MAX_DIMENSIONS]);

#ifdef __cplusplus
}
#endif
