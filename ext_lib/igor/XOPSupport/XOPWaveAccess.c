/*	XOPWaveAccess.c
	
	Routines for Igor XOPs that provide access to multi-dimensional waves.
*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h

/*	ABOUT WAVE ACCESS
	
	Waves are stored in handles that contain a wave header structure followed
	by the wave data. There have been 3 main wave structures. The current wave
	structure has been in use since Igor Pro 3.
	
	Routines in this file provide access to wave properties and wave data.
	Some of these routines call back to Igor and others access the wave structure
	directly. These routines provide compatibility with future versions of Igor.	
	
	IGOR64 Versus IGOR32

	As of Igor Pro 6.2, all waves in both IGOR32 and IGOR64 use the version 3
	wave structure. However, the structures are different when running with
	IGOR64 versus IGOR32 because pointer fields are 8 bytes instead of 4.
	The routines in this file take care of this distinction.
*/

#define MAXDIMS1 4					// Maximum dimensions in wave struct version 3.
#define kWaveStruct3Version 1		// The value in the version field of the version 3 wave structure is 1.

#ifdef IGOR64
	#define WAVE_VERSION(waveH) (*(short*)((char*)*waveH+34))
#else
	#define WAVE_VERSION(waveH) (*(short*)((char*)*waveH+26))
#endif

/*	*** Wave Access Routines *** */

/*	WaveHandleModified(waveH)

	WaveHandleModified does a callBack to Igor to notify it that the specified
	wave has been modified.
	
	Thread Safety: WaveHandleModified is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
void
WaveHandleModified(waveHndl waveH)
{
	CallBack1(WAVEMODIFIED, waveH);
}

// WaveHandlesModified has been removed. Use WaveHandleModified.

/*	WaveModified(waveName)

	WaveModified() does a callBack to Igor to notify it that the specified wave
	has been modified. It ASSUMES that waveName is the name of a valid wave.
	
	Thread Safety: WaveModified is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
void
WaveModified(const char *waveName)
{
	CallBack1(WAVEMODIFIED, FetchWave(waveName));
}

/*	FetchWave(waveName)

	FetchWave returns a handle to the data structure for the specified wave
	or NULL if the wave does not exist.
	
	Thread Safety: FetchWave is thread-safe with Igor Pro 6.20 or later.
*/
waveHndl
FetchWave(const char *waveName)
{
	waveHndl waveH;
	
	waveH = (waveHndl)CallBack1(FETCHWAVE, (void*)waveName);
	return(waveH);
}

/*	FetchWaveFromDataFolder(dataFolderH, waveName)

	FetchWaveFromDataFolder returns a handle to the data structure for the
	specified wave in the specified data folder or NULL if the wave does not exist.

	If dataFolderH is NULL, it uses the current data folder.
	
	Thread Safety: FetchWaveFromDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
waveHndl
FetchWaveFromDataFolder(DataFolderHandle dataFolderH, const char* waveName)
{
	waveHndl waveH;
	
	waveH = (waveHndl)CallBack2(FETCHWAVE_FROM_DATAFOLDER, dataFolderH, (void*)waveName);
	return(waveH);
}

/*	WaveType(waveH)

	Returns wave type which is:
		NT_FP32, NT_FP64 for single or double precision floating point
		NT_I8, NT_I16, NT_I32 for 8, 16 or 32 bit signed integer
	This is ORed with NT_COMPLEX if the wave is complex
	and ORed with NT_UNSIGNED if the wave is unsigned integer.
	
	The wave type can also be one of the following:
		TEXT_WAVE_TYPE			text wave
		WAVE_TYPE				wave wave - holds wave references
		DATAFOLDER_TYPE			DFREF wave - holds data folder references
	These types of waves can not be complex.
	
	Note: Future versions may support other wave types.
		  Always check the wave type to make sure it is a type you can handle.
	
	Thread Safety: WaveType is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
WaveType(waveHndl waveH)
{
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:		// type field
			#ifdef IGOR64
				return (*(short*)((char*)*waveH+24));
			#else
				return (*(short*)((char*)*waveH+16));
			#endif
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (int)CallBack1(WAVETYPE, waveH);
			break;
	}
}

/*	WavePoints(waveH)

	Returns number of points in wave.
	
	For multi-dimensional waves, WavePoints returns the total number of
	points in all dimensions.
	
	Thread Safety: WavePoints is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
CountInt
WavePoints(waveHndl waveH)
{
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:		// npnts field
			#ifdef IGOR64
				return (*(CountInt*)((char*)*waveH+16));
			#else
				return (*(CountInt*)((char*)*waveH+12));
			#endif
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (CountInt)CallBack1(WAVEPOINTS, waveH);
			break;
	}
}

/*	WaveName(waveH, namePtr)

	Returns name of the wave via namePtr.
	
	Thread Safety: WaveName is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
void
WaveName(waveHndl waveH, char *namePtr)
{
	CallBack2(WAVENAME, waveH, namePtr);
}

/*	WaveData(waveH)

	Returns pointer to start of wave's data.
	
	WARNING: Do not use WaveData to access a text wave.
	
	Thread Safety: WaveData is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
void*
WaveData(waveHndl waveH)
{
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:		// Address of wData field
			#ifdef IGOR64
				return (void*)((char*)*waveH+400);
			#else
				return (void*)((char*)*waveH+320);
			#endif
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (void*)CallBack1(WAVEDATA, waveH);
			break;
	}
}

/*	WaveScaling(waveH, hsAPtr, hsBPtr, topPtr, botPtr)

	Returns the wave's X and Y scaling information.
	
	hsA and hsB define the transformation from point number to X value where
	  X value = hsA*Point# + hsB.
	  
	top and bottom are the values that the user entered for the wave's Y Full Scale.
	If both are zero, there is no Y Full Scale for this wave.
	
	The wave can be multi-dimensional. See MDGetWaveScaling.
	
	Thread Safety: WaveScaling is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
void
WaveScaling(waveHndl waveH, double *hsAPtr, double *hsBPtr, double *topPtr, double *botPtr)
{
	CallBack5(WAVESCALING, waveH, hsAPtr, hsBPtr, topPtr, botPtr);
}

/*	SetWaveScaling(waveH, hsAPtr, hsBPtr, topPtr, botPtr)

	Sets the wave's X and Y scaling information.
	If hsAPtr and/or hsBPtr are NULL, does not set X scaling.
	If topPtr and/or botPtr are NULL, does not set Y scaling.
	
	For multidimensional waves use MDSetWaveScaling.
	
	Thread Safety: SetWaveScaling is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
void
SetWaveScaling(waveHndl waveH, const double *hsAPtr, const double *hsBPtr, const double *topPtr, const double *botPtr)
{
	CallBack5(SETWAVESCALING, waveH, (void*)hsAPtr, (void*)hsBPtr, (void*)topPtr, (void*)botPtr);
}

/*	WaveUnits(waveH, xUnits, dataUnits)

	Returns the wave's X and data units string.
	
	In Igor Pro 3.0, the number of characters allowed was increased from
	3 to 49 (MAX_UNIT_CHARS). For backward compatibility, WaveUnits will return
	no more than 3 characters (plus the null terminator). To get the full units and
	to access units for higher dimensions, new XOPs should use MDGetWaveUnits instead
	of WaveUnits.
	
	Thread Safety: WaveUnits is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
void
WaveUnits(waveHndl waveH, char *xUnits, char *dataUnits)
{
	CallBack3(WAVEUNITS, waveH, xUnits, dataUnits);
}

/*	SetWaveUnits(waveH, xUnits, dataUnits)

	Sets the wave's X and data units string.

	If xUnits is NULL, does not set X units.
	If dataUnits is NULL, does not set data units.

	When running with Igor Pro 3.0 or later, units can be up to 49 (MAX_UNIT_CHARS)
	characters long. In earlier versions of Igor, units are limited to 3 characters.
	You can pass a longer units string but only the first 3 characters will be used.
	
	To access units for higher dimensions, use MDSetWaveUnits instead
	of SetWaveUnits.
	
	Thread Safety: SetWaveUnits is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
void
SetWaveUnits(waveHndl waveH, const char *xUnits, const char *dataUnits)
{
	CallBack3(SETWAVEUNITS, waveH, (void*)xUnits, (void*)dataUnits);
}

/*	WaveNote(waveH)

	Returns handle to wave's note text or NULL if wave has no note.
	
	NOTE: You should not modify or dispose this handle. Use SetWaveNote() to change note.
	
	Thread Safety: WaveNote is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
Handle
WaveNote(waveHndl waveH)
{
	return (Handle)CallBack1(WAVENOTE, waveH);
}

/*	SetWaveNote(waveH, noteHandle)

	Sets wave's note text.
	noteHandle is handle to text or NULL if you want to kill wave's note.
	
	NOTE: once you pass the noteHandle to Igor, don't modify or dispose of it.
	
	Thread Safety: SetWaveNote is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
void
SetWaveNote(waveHndl waveH, Handle noteHandle)
{
	CallBack2(SETWAVENOTE, waveH, noteHandle);
}

/*	WaveModDate(waveH)

	Returns wave modification date. This is an unsigned long in Macintosh
	date/time format, namely the number of seconds since midnight, January 1, 1904.
	
	The main use for this is to allow an XOP to check if a particular wave has been
	changed. For example, an XOP that displays a wave can check the wave's mod date
	in the XOP's idle routine. If the date has changed since the last time the
	XOP updated its display, it can update the display again.
	
	Modification date tracking was added to Igor in Igor 1.2. If a wave
	is loaded from a file created by an older version of Igor, the mod date field
	will be zero and this routine will return zero.
	
	Thread Safety: WaveModDate is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
TickCountInt
WaveModDate(waveHndl waveH)
{
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:		// modDate field
			#ifdef IGOR64
				return (*(UInt32*)((char*)*waveH+12));
			#else
				return (*(UInt32*)((char*)*waveH+8));
			#endif
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (TickCountInt)CallBack1(WAVEMODDATE, waveH);
			break;
	}
}

/*	WaveLock(waveH)

	Returns the lock state of the wave.
	
	A return value of 0 signifies that the wave is not locked.
	
	A return value of 1 signifies that the wave is locked. In that case, you should
	not kill the wave or modify it in any way.
	
	Thread Safety: WaveLock is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
WaveLock(waveHndl waveH)
{
	return (int)CallBack1(WAVELOCK, waveH);
}

/*	SetWaveLock(waveH, lockState)

	Sets wave's lock state.
	
	If lockState is 0, the wave will be unlocked. If it is 1, the wave will be locked.
	
	All other bits are reserved.
	
	Returns the previous state of the wave lock setting.
	
	Thread Safety: SetWaveLock is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
SetWaveLock(waveHndl waveH, int lockState)
{
	return (int)CallBack2(SETWAVELOCK, waveH, (void*)lockState);
}

/*	WaveModState(waveH)

	Returns the truth that the wave has been modified since the last save to disk.

	This routine works with all versions of Igor.
	
	Thread Safety: WaveModState is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
WaveModState(waveHndl waveH)
{
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:		// swModified field
			#ifdef IGOR64
				return (*(short*)((char*)*waveH+364));
			#else
				return (*(short*)((char*)*waveH+296));
			#endif
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (int)CallBack1(WAVEMODSTATE, waveH);
			break;
	}
}

/*	WaveModCount(waveH)

	Returns a value that can be used to tell if a wave has been changed between
	one call to WaveModCount and another. This function was created so that the
	Igor Data Browser could tell when a wave was changed. Previously, the Data
	Browser used the WaveModDate function, but that function can not identify
	changes that happen closer than 1 second apart.
	
	The exact value returned by WaveModCount has no significance. The only valid
	use for it is to compare the values returned by two calls to WaveModCount. If
	they are the different, the wave was changed in the interim.
	
	Example:
		waveModCount1 = WaveModCount(waveH);
		. . .
		waveModCount2 = WaveModCount(waveH);
		if (waveModCount2 != waveModCount1)
			// Wave has changed.
	
	Thread Safety: WaveModCount is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
WaveModCount(waveHndl waveH)
{
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:		// modCount
			#ifdef IGOR64
				return (*(short*)((char*)*waveH+380));
			#else
				return (*(short*)((char*)*waveH+308));
			#endif
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (int)CallBack1(WAVEMODCOUNT, waveH);
			break;
	}
}

// GetWavesInfo has been removed. Use WaveType, WavePoints and WaveData as needed.

// SetWavesStates has been removed. It is obsolete and no longer needed.

/*	MakeWave(waveHPtr,waveName,numPoints,type,overwrite)

	Tries to make wave with specified name, number of points, numeric type.
	
	NOTE: For historical reasons from ancient times, prior to XOP Toolkit 6,
	if numPoints was zero, MakeWave created a 128 point wave. As of XOP Toolkit 6,
	passing 0 for numPoints creates a zero-point wave.

	Returns error code or 0 if wave was made.
	
	Thread Safety: MakeWave is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MakeWave(waveHndl *waveHPtr, const char *waveName, CountInt numPoints, int type, int overwrite)
{
	*waveHPtr = NULL;
	
	if (numPoints == 0) {	// See note above.
		CountInt dimSizes[MAX_DIMENSIONS+1];
		MemClear(dimSizes, sizeof(dimSizes));
		return (int)MDMakeWave(waveHPtr,waveName,NULL,dimSizes,type,overwrite);	
	}
	
	return (int)CallBack5(MAKEWAVE, waveHPtr, (void*)waveName, (void *)numPoints, (void *)type, (void *)overwrite);
}

/*	ChangeWave(waveH, numPoints, type)

	Changes wave length to specified number of points and numeric type.
	Returns error code or 0 if OK.
	
	Thread Safety: ChangeWave is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
ChangeWave(waveHndl waveH, CountInt numPoints, int type)
{
	return (int)CallBack3(CHANGEWAVE,waveH,(void *)numPoints,(void *)type);
}

/*	KillWave(waveH)

	Kills wave.
	Returns error code or 0 if OK.
	
	NOTE: if wave is in use, returns error code and wave is not killed.
	
	Thread Safety: KillWave is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
KillWave(waveHndl waveH)
{
	return (int)CallBack1(KILLWAVE, waveH);
}

/*	*** Multi-dimension Wave Access Routines *** */

/*	MDMakeWave(waveHPtr,waveName,dataFolderH,dimSizes,type,overwrite)

	Tries to make wave with specified name and type in the specified data folder.
	
	If dataFolderH is NULL, it uses the current data folder.
	
	For each dimension, dimSizes[i] specifies the number of points
	in that dimension. For a wave of dimension n, i goes from 0 to n-1.
	
	NOTE: dimSizes[n] must be zero. This is how Igor determines
		  how many dimensions the wave is to have.
	
	If you are running with Igor Pro 6.20 or later, you can pass -1 for
	dataFolderH to create a free wave. This would be appropriate if, for example,
	your external function were called from a user-defined function with a flag
	parameter indicating that the user wants to create a free wave.
	You must be sure that you are running with Igor Pro 6.20 or a later.
	If you are running with an earlier version, this will cause a crash.
	For example:
		if (igorVersion < 620)
			return IGOR_OBSOLETE;	// Can't create free wave in old Igor
		result = MDMakeWave(&waveH, "freejack", (DataFolderHandle)-1, dimensionSizes, type, overwrite);
	When making a free wave, the overwrite parameter is irrelevant and is ignored.
	
	You can also create a free wave using GetOperationDestWave.

	Returns error code or 0 if wave was made.
	
	Thread Safety: MDMakeWave is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDMakeWave(waveHndl *waveHPtr, const char *waveName, DataFolderHandle dataFolderH, CountInt dimSizes[MAX_DIMENSIONS+1], int type, int overwrite)
{
	*waveHPtr = NULL;
	return (int)CallBack6(MD_MAKEWAVE,waveHPtr, (void*)waveName,dataFolderH, dimSizes, (void *)type, (void *)overwrite);
}

/*	MDGetWaveDimensions(waveH, numDimensionsPtr, dimSizes)

	Returns number of used dimensions in wave via numDimensionsPtr
	and the number of points in each used dimension via dimSizes.

	If you only want o know the number of dimensions, you can pass NULL for dimSizes.
	
	NOTE: dimSizes (if not NULL) should have room for MAX_DIMENSIONS+1 values.

	For an n dimensional wave, MDGetWaveDimensions sets dimSizes[0..n-1] 
	to the number of elements in the corresponding dimension and sets
	dimSizes[n..MAX_DIMENSIONS] to zero, indicating that they are unused
	dimensions. This guarantees that there will always be an element containing
	zero in the dimSizes array.
	
	The function result is 0 or an Igor error code.
	
	Thread Safety: MDGetWaveDimensions is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDGetWaveDimensions(waveHndl waveH, int* numDimensionsPtr, CountInt dimSizes[MAX_DIMENSIONS+1])
{
	CountInt dimSize;
	CountInt* nDim;
	int i;
	
	if (dimSizes != NULL)
		MemClear(dimSizes, (MAX_DIMENSIONS+1)*sizeof(CountInt));		/* All unused dimensions are zeroed. */
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:		// nDim field
			#ifdef IGOR64
				nDim = (CountInt*)((char*)*waveH + 80);
			#else
				nDim = (CountInt*)((char*)*waveH + 68);
			#endif
			*numDimensionsPtr = 0;
			for(i=0; i<MAXDIMS1; i++) {
				dimSize = nDim[i];
				if (dimSize==0)
					break;
				*numDimensionsPtr += 1;
				if (dimSizes!=NULL)
					dimSizes[i] = dimSize;
			}
			return 0;
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (int)CallBack3(MD_GETWAVEDIMENSIONS, waveH, numDimensionsPtr, dimSizes);
			break;
	}
}

/*	MDChangeWave(waveH, type, dimSizes)

	Changes one or more of the following:
		the wave's data type
		the number of dimensions in the wave
		the number of points in one or more dimensions
	
	type is one of the following:
		-1 for no change in data type.
		
		NT_FP32, NT_FP64 for single or double precision floating point
		NT_I8, NT_I16, NT_I32 for 8, 16 or 32 bit signed integer.
		These may be ORed with NT_COMPLEX to make the wave complex
		and ORed with NT_UNSIGNED to make the wave unsigned integer.
	
		TEXT_WAVE_TYPE.
	However converting a text wave to numeric or vice versa is currently not
	supported and will result in an error.
	
	dimSizes[i] contains the desired number of points for dimension i.
	For n dimensions, dimSizes[n] must be zero. Then the size of each
	dimension is set by dimSizes[0..n-1].
	
	If dimSizes[i] == -1, then the size of dimension i will be unchanged.
	
	Returns 0 or an error code.
	
	Thread Safety: MDChangeWave is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDChangeWave(waveHndl waveH, int type, CountInt dimSizes[MAX_DIMENSIONS+1])
{
	return (int)CallBack3(MD_CHANGEWAVE, waveH, (void *)type, dimSizes);
}

/*	MDChangeWave2(waveH, type, dimSizes, mode)

	This is the same as MDChangeWave except for the added mode parameter.
	
	mode = 0:	Does a normal redimension.
	
	mode = 1:	Changes the wave's dimensions without changing the wave data.
				This is useful, for example, when you have a 2D wave
				consisting of 5 rows and 3 columns which you want to treat
				as a 2D wave consisting of 3 rows and 5 columns.
	
	mode = 2:	Changes the wave data from big-endian to little-endian
				or vice versa. This is useful when you have loaded data
				from a file that uses a byte ordering different from that
				of the platform on which you are running.
	
	Returns 0 or an error code.
	
	Thread Safety: MDChangeWave2 is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDChangeWave2(waveHndl waveH, int type, CountInt dimSizes[MAX_DIMENSIONS+1], int mode)
{
	if (mode == 0)
		return (int)CallBack3(MD_CHANGEWAVE, waveH, (void *)type, dimSizes);
	
	return (int)CallBack4(MD_CHANGEWAVE2, waveH, (void *)type, dimSizes, (void *)mode);
}

/*	MDGetWaveScaling(waveH, dimension, sfAPtr, sfBPtr)

	Returns the dimension scaling values or the full scale values for the wave.
	If dimension is -1, it returns the full scale values. Otherwise, it returns the dimension
	scaling.
	
	For dimension d, the scaled index of point p is:
	  scaled index = sfA[d]*p + sfB[d].
	
	If dimension is -1, this gets the wave's data full scale setting instead of
	dimension scaling. *sfAPtr points to the top full scale value and *sfBPtr
	points to the bottom full scale value.
	
	The function result is 0 or an Igor error code.
	
	Thread Safety: MDGetWaveScaling is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDGetWaveScaling(waveHndl waveH, int dimension, double* sfAPtr, double* sfBPtr)
{
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:
		default:		/* Future version of Igor with a different wave structure. */
			return (int)CallBack4(MD_GETWAVESCALING, waveH, (void*)dimension, sfAPtr, sfBPtr);	// HR, 091201: This was CallBack3 which was wrong.
			break;
	}
}

/*	MDSetWaveScaling(waveH, dimension, sfAPtr, sfBPtr)

	Sets the dimension scaling values or the full scale values for the wave.
	If dimension is -1, it sets the full scale values. *sfAPtr is the top full
	scale value and *sfBPtr is the bottom full scale value.
	
	If dimension is 0 or greater, it sets the dimension scaling.
	For dimension d, the scaled index of point p is:
	  scaled index = sfA[d]*p + sfB[d].
	The sfA value can never be set to zero. If you pass 0.0 for sfA, this routine
	will use 1.0 instead.
	
	If dimension is -1, this sets the wave's data full scale setting instead of
	dimension scaling. *sfAPtr points to the top full scale value and *sfBPtr
	points to the bottom full scale value.
	
	The function result is 0 or an Igor error code.
	
	Thread Safety: MDSetWaveScaling is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDSetWaveScaling(waveHndl waveH, int dimension, const double* sfAPtr, const double* sfBPtr)
{
	/*	HR, 031112, XOP Toolkit 5.00:
		This routine previously set wave fields directly instead of calling
		back to Igor. This was not sufficient because Igor internally sets
		other flags in addition to the wave fields. So now it always does
		a callback to Igor.
	*/
	return (int)CallBack4(MD_SETWAVESCALING, waveH, (void*)dimension, (void*)sfAPtr, (void*)sfBPtr);
}

/*	MDGetWaveUnits(waveH, dimension, units)

	Returns the units string for the specified dimension in the wave via units.
	
	To get the data units (as opposed to dimension units), pass -1 for dimension.
	
	In Igor Pro 3.0 or later, units can be up to 49 (MAX_UNIT_CHARS) characters. In
	earlier versions, units were limited to 3 characters.
		
	The function result is 0 or an Igor error code.
	
	Thread Safety: MDGetWaveUnits is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDGetWaveUnits(waveHndl waveH, int dimension, char units[MAX_UNIT_CHARS+1])
{
	*units = 0;
	
	switch (WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:
		default:
			return (int)CallBack3(MD_GETWAVEUNITS, waveH, (void*)dimension, units);
			break;
	}
}

/*	MDSetWaveUnits(waveH, int dimension, units[MAX_UNIT_CHARS+1])

	Sets the units string for the specified dimension.
	
	To set the data units (as opposed to dimension units), pass -1 for dimension.
	
	In Igor Pro 3.0 or later, units can be up to 49 (MAX_UNIT_CHARS) characters. In
	earlier versions, units were limited to 3 characters. In either case, if the
	string you pass is too long, Igor will store a truncated version of it. 
	
	The function result is 0 or an Igor error code.
	
	Thread Safety: MDSetWaveUnits is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDSetWaveUnits(waveHndl waveH, int dimension, const char units[MAX_UNIT_CHARS+1])
{
	switch (WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:
		default:
			return (int)CallBack3(MD_SETWAVEUNITS, waveH, (void*)dimension, (void*)units);
			break;
	}
}

/*	MDGetDimensionLabel(waveH, dimension, element, label)

	Returns the label for the specified element of the specified dimension via label.
	
	If element is -1, this specifies a label for the entire dimension.
	If element is between 0 to n-1, where n is the size of the dimension,
	then element specifies a label for that element of the dimension only.
	
	label should have room for MAX_DIM_LABEL_CHARS+1 bytes.
	At present, Igor Pro truncates dimension labels at 31 characters.
	A future version may increase this.
		
	The function result is 0 or an Igor error code.
	
	Thread Safety: MDGetDimensionLabel is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDGetDimensionLabel(waveHndl waveH, int dimension, IndexInt element, char label[MAX_DIM_LABEL_CHARS+1])
{
	switch (WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:
		default:
			return (int)CallBack4(MD_GETDIMLABELS, waveH, (void*)dimension, (void*)element, label);
			break;
	}
}

/*	MDSetDimensionLabel(waveH, dimension, element, label)

	Sets the label for the specified element of the specified dimension via label.
	
	If element is -1, this specifies a label for the entire dimension.
	If element is between 0 to n-1, where n is the size of the dimension,
	then element specifies a label for that element of the dimension only.
	
	At present, Igor Pro truncates dimension labels at 31 characters.
	A future version may increase this.
		
	The function result is 0 or an Igor error code.
	
	Thread Safety: MDSetDimensionLabel is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDSetDimensionLabel(waveHndl waveH, int dimension, IndexInt element, char label[MAX_DIM_LABEL_CHARS+1])
{
	switch (WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:
		default:
			return (int)CallBack4(MD_SETDIMLABELS, waveH, (void*)dimension, (void*)element, label);
			break;
	}
}

/*	MDAccessNumericWaveData(waveH, accessMode, dataOffsetPtr)
	
	MDAccessNumericWaveData provides access to the data for numeric waves.
	
	waveH is the wave handle containing the data you want to access.
	
	accessMode is a code that tells Igor how you plan to access the wave data
	and is used for a future compatibility check. At present, there is only one
	accessMode. You should use the symbol kMDWaveAccessMode0 for the accessMode parameter.
	
	On output, if there is no error, *dataOffsetPtr contains the offset in bytes
	from the start of the wave handle to the data.
	
	MDAccessNumericWaveData returns 0 or an error code.
	
	If it returns a non-zero error code, you should not attempt to access
	the wave data but merely return the error code to Igor.
	
	At present, there is only one case in which MDAccessNumericWaveData will return an
	error code. This is if the wave is a text wave.
	
	It is possible that a future version of Igor Pro will store wave data in
	a different way, such that the current method of accessing wave data will no
	longer work. If your XOP ever runs with such a future Igor, MDAccessNumericWaveData
	will return an error code indicating the incompatibility. Your XOP will refrain
	from attempting to access the wave data and return the error code to Igor.
	This will prevent a crash and indicate the nature of the problem to the user.
	
	Numeric wave data is stored contiguously in the wave handle in one of the
	supported data types (NT_I8, NT_I16, NT_I32, NT_FP32, NT_FP64). These types
	will be ORed with NT_CMPLX if the wave is complex and ORed with NT_UNSIGNED if
	the wave is unsigned integer.
	
	Although they are not truly numeric waves, MDAccessNumericWaveData also allows
	you to access wave reference waves (WAVE_TYPE) and data folder reference waves
	(DATAFOLDER_TYPE) which store waveHndls and DataFolderHandles respectively.
	
	To access the a particular point, you need to know the number of data points in each
	dimension. To find this, you must call MDGetWaveDimensions. This returns the number
	of used dimensions in the wave and an array of dimension lengths. The dimension lengths
	are interpreted as follows:
		dimSizes[0]		number of rows in a column
		dimSizes[1]		number of columns in a layer
		dimSizes[2]		number of layers in a chunk
		dimSizes[3]		number of chunks in the wave
	
	The data is stored in row/column/layer/chunk order. This means that,
	as you step linearly through memory one point at a time, you first pass the
	value for each row in the first column. At the end of the first column,
	you reach the start of the second column. After you have passed the data for
	each column in the first layer, you reach the data for the first column
	in the second layer. After you have passed the data for each layer, you
	reach the data for the first layer of the second chunk. And so on.
	
	Thread Safety: MDAccessNumericWaveData is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDAccessNumericWaveData(waveHndl waveH, int accessMode, BCInt* dataOffsetPtr)
{
	if (accessMode!=kMDWaveAccessMode0)		/* Call back to Igor if this is an unknown (future) access mode. */
		return (int)CallBack3(MD_ACCESSNUMERICWAVEDATA, waveH, (void*)accessMode, dataOffsetPtr);
	
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:		// offset of wData field
			if (WaveType(waveH) == TEXT_WAVE_TYPE)
				return NUMERIC_ACCESS_ON_TEXT_WAVE;
			#ifdef IGOR64
				*dataOffsetPtr = 400;
			#else
				*dataOffsetPtr = 320;
			#endif
			return 0;
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (int)CallBack3(MD_ACCESSNUMERICWAVEDATA, waveH, (void*)accessMode, dataOffsetPtr);
			break;
	}
}

/*	MDPointIndexV3(waveH, indices[MAX_DIMENSIONS], indexPtr)
	
	For version 3 waves, returns a linear index number to access the point
	indicated by the indices.
	
	waveH must be a version 3 wave.
	
	NOTE: This routine ignores indices for dimensions that do not exist in the wave.
	
	Thread Safety: MDPointIndexV3 is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
static int
MDPointIndexV3(waveHndl waveH, CountInt indices[MAX_DIMENSIONS], IndexInt* indexPtr)
{
	CountInt* nDim;
	IndexInt p;
	CountInt dimtot= 1;
	CountInt dimSize;
	int i;
	
	#ifdef IGOR64
		nDim = (CountInt*)((char*)*waveH + 80);
	#else
		nDim = (CountInt*)((char*)*waveH + 68);
	#endif
	
	*indexPtr = 0;
	for(i=0; i<MAXDIMS1; i++) {
		dimSize = nDim[i];
		if (dimSize == 0)
			break;							/* We've done all used dimensions. */
		p = indices[i];
		if (p<0 || p>=dimSize)
			return MD_WAVE_BAD_INDEX;
		if (dimSize != 0) {
			*indexPtr += p*dimtot;
			dimtot *= dimSize;
		}
	}

	return 0;		
}

/*	FetchNumericValue(type, dataStartPtr, index, value)

	type is an Igor numeric type.
	
	dataStartPtr points to the start of the data of the specified type.
	
	index is the "point number" of the point of interest, considering
	the data as one long vector.
	
	Thread Safety: FetchNumericValue is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
FetchNumericValue(int type, const char* dataStartPtr, IndexInt index, double value[2])
{
	int isComplex;
	double* dp;
	float* fp;
	SInt32* lp;
	UInt32* ulp;
	short* sp;
	unsigned short* usp;
	signed char* cp;
	unsigned char* ucp;

	isComplex = type & NT_CMPLX;
	if (isComplex)
		index *= 2;
	type &= ~NT_CMPLX;
	
	switch (type) {
		case NT_FP64:
			dp = (double*)dataStartPtr + index;
			value[0] = *dp++;
			if (isComplex)
				value[1] = *dp;
			break;
	
		case NT_FP32:
			fp = (float*)dataStartPtr + index;
			value[0] = *fp++;
			if (isComplex)
				value[1] = *fp;
			break;
	
		case NT_I32:
			lp = (SInt32*)dataStartPtr + index;
			value[0] = *lp++;
			if (isComplex)
				value[1] = *lp;
			break;
	
		case NT_I32 | NT_UNSIGNED:
			ulp = (UInt32*)dataStartPtr + index;
			value[0] = *ulp++;
			if (isComplex)
				value[1] = *ulp;
			break;
	
		case NT_I16:
			sp = (short*)dataStartPtr + index;
			value[0] = *sp++;
			if (isComplex)
				value[1] = *sp;
			break;
	
		case NT_I16 | NT_UNSIGNED:
			usp = (unsigned short*)dataStartPtr + index;
			value[0] = *usp++;
			if (isComplex)
				value[1] = *usp;
			break;
	
		case NT_I8:
			cp = (signed char*)dataStartPtr + index;
			value[0] = *cp++;
			if (isComplex)
				value[1] = *cp;
			break;
	
		case NT_I8 | NT_UNSIGNED:
			ucp = (unsigned char*)dataStartPtr + index;
			value[0] = *ucp++;
			if (isComplex)
				value[1] = *ucp;
			break;
		
		case TEXT_WAVE_TYPE:
			return NUMERIC_ACCESS_ON_TEXT_WAVE;

		case WAVE_TYPE:
		case DATAFOLDER_TYPE:
			return NUMERIC_ACCESS_ON_NON_NUMERIC_WAVE;
			
		default:
			return WAVE_TYPE_INCONSISTENT;		/* Should never happen. */
	}
	
	return 0;
}

/*	MDGetNumericWavePointValue(waveH, indices, value)
	
	Returns via value the value of a particular point in the specified wave.
	The value returned is always double precision floating point, regardless
	of the precision of the wave.
	
	indices is an array of dimension indices. For example, for a 3D wave,
		indices[0] should contain the row number
		indices[1] should contain the column number
		indices[2] should contain the layer number
	
	NOTE: This routine ignores indices for dimensions that do not exist in the wave.
	
	The real part of the value of specified point is returned in value[0].
	If the wave is complex, the imaginary part of the value of specified point
	is returned in value[1]. If the wave is not complex, value[1] is undefined.
	
	The function result is 0 or an error code.

	Currently the only error code returned is MD_WAVE_BAD_INDEX, indicating
	that you have passed invalid indices. Future versions of Igor may return
	other error codes. If you receive an error, just return it to Igor so
	that it will be reported to the user.
	
	Thread Safety: MDGetNumericWavePointValue is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDGetNumericWavePointValue(waveHndl waveH, IndexInt indices[MAX_DIMENSIONS], double value[2])
{
	int waveType;
	IndexInt index;
	int result;
	
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:
			waveType = WaveType(waveH);
			switch(waveType) {
				case TEXT_WAVE_TYPE:
					return NUMERIC_ACCESS_ON_TEXT_WAVE;
				case WAVE_TYPE:
				case DATAFOLDER_TYPE:
					return NUMERIC_ACCESS_ON_NON_NUMERIC_WAVE;
					break;
			}
			if (result = MDPointIndexV3(waveH, indices, &index))
				return result;
			return FetchNumericValue(waveType, (char*)WaveData(waveH), index, value);
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (int)CallBack3(MD_GETWAVEPOINTVALUE, waveH, indices, value);
			break;
	}
}

/*	StoreNumericValue(type, dataStartPtr, index, value)

	type is an Igor numeric type.
	
	dataStartPtr points to the start of the data of the specified type.
	
	index is the "point number" of the point of interest, considering
	the data as one long vector.

	HR, 091001: Added casts to prevent VC warnings.
	
	Thread Safety: StoreNumericValue is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
StoreNumericValue(int type, char* dataStartPtr, CountInt index, double value[2])
{
	int isComplex;
	double* dp;
	float* fp;
	SInt32* lp;
	UInt32* ulp;
	short* sp;
	unsigned short* usp;
	signed char* cp;
	unsigned char* ucp;

	isComplex = type & NT_CMPLX;
	if (isComplex)
		index *= 2;
	type &= ~NT_CMPLX;
	
	switch (type) {
		case NT_FP64:
			dp = (double*)dataStartPtr + index;
			*dp++ = value[0];
			if (isComplex)
				*dp = value[1];
			break;
	
		case NT_FP32:
			fp = (float*)dataStartPtr + index;
			*fp++ = (float)value[0];
			if (isComplex)
				*fp = (float)value[1];
			break;
	
		case NT_I32:
			lp = (SInt32*)dataStartPtr + index;
			*lp++ = (SInt32)value[0];
			if (isComplex)
				*lp = (SInt32)value[1];
			break;
	
		case NT_I32 | NT_UNSIGNED:
			ulp = (UInt32*)dataStartPtr + index;
			*ulp++ = (UInt32)value[0];
			if (isComplex)
				*ulp = (UInt32)value[1];
			break;
	
		case NT_I16:
			sp = (short*)dataStartPtr + index;
			*sp++ = (short)value[0];
			if (isComplex)
				*sp = (short)value[1];
			break;
	
		case NT_I16 | NT_UNSIGNED:
			usp = (unsigned short*)dataStartPtr + index;
			*usp++ = (unsigned short)value[0];
			if (isComplex)
				*usp = (unsigned short)value[1];
			break;
	
		case NT_I8:
			cp = (signed char*)dataStartPtr + index;
			*cp++ = (signed char)value[0];
			if (isComplex)
				*cp = (signed char)value[1];
			break;
	
		case NT_I8 | NT_UNSIGNED:
			ucp = (unsigned char*)dataStartPtr + index;
			*ucp++ = (unsigned char)value[0];
			if (isComplex)
				*ucp = (unsigned char)value[1];
			break;
		
		case TEXT_WAVE_TYPE:
			return NUMERIC_ACCESS_ON_TEXT_WAVE;

		case WAVE_TYPE:
		case DATAFOLDER_TYPE:
			return NUMERIC_ACCESS_ON_NON_NUMERIC_WAVE;
		
		default:
			return WAVE_TYPE_INCONSISTENT;		/* Should never happen. */
	}
	
	return 0;
}

/*	MDSetNumericWavePointValue(waveH, indices, value)
	
	Sets the value of a particular point in the specified wave.
	The value that you supply is always double precision floating point,
	regardless of the precision of the wave.
	
	indices is an array of dimension indices. For example, for a 3D wave,
		indices[0] should contain the row number
		indices[1] should contain the column number
		indices[2] should contain the layer number
	
	NOTE: This routine ignores indices for dimensions that do not exist in the wave.
	
	You should pass in value[0] the real part of the value.
	If the wave is complex, you should pass the complex part in value[1].
	If the wave is not complex, Igor ignores value[1].
	
	The function result is 0 or an error code.

	Currently the only error code returned is MD_WAVE_BAD_INDEX, indicating
	that you have passed invalid indices. Future versions of Igor may return
	other error codes. If you receive an error, just return it to Igor so
	that it will be reported to the user.
	
	Thread Safety: MDSetNumericWavePointValue is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDSetNumericWavePointValue(waveHndl waveH, IndexInt indices[MAX_DIMENSIONS], double value[2])
{
	int waveType;
	IndexInt index;
	int result;

	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:
			waveType = WaveType(waveH);
			switch(waveType) {
				case TEXT_WAVE_TYPE:
					return NUMERIC_ACCESS_ON_TEXT_WAVE;
				case WAVE_TYPE:
				case DATAFOLDER_TYPE:
					return NUMERIC_ACCESS_ON_NON_NUMERIC_WAVE;
					break;
			}
			if (result = MDPointIndexV3(waveH, indices, &index))
				return result;
			return StoreNumericValue(waveType, (char*)WaveData(waveH), index, value);
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (int)CallBack3(MD_SETWAVEPOINTVALUE, waveH, indices, value);
			break;
	}
}

/*	MDGetDPDataFromNumericWave(waveH, dPtr)
	
	MDGetDPDataFromNumericWave stores a double-precision representation of
	the specified wave's data in the memory pointed to by dPtr. dPtr must
	point to a block of memory that you have allocated and which must be
	at least (WavePoints(waveH)*sizeof(double)) bytes.
	
	This routine is intended for use with MDStoreDPDataInNumericWave.
	
	The function result is zero or an error code.
	
	Thread Safety: MDGetDPDataFromNumericWave is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDGetDPDataFromNumericWave(waveHndl waveH, double* dPtr)
{
	CountInt numNumbers;
	int bytesPerPoint;
	int waveType, type2;
	int srcFormat;
	
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:
			waveType = WaveType(waveH);
			switch(waveType) {
				case TEXT_WAVE_TYPE:
					return NUMERIC_ACCESS_ON_TEXT_WAVE;
				case WAVE_TYPE:
				case DATAFOLDER_TYPE:
					return NUMERIC_ACCESS_ON_NON_NUMERIC_WAVE;
					break;
			}
			type2 = waveType & ~NT_CMPLX;
			switch(type2) {
				case NT_I8:
					bytesPerPoint = sizeof(char);
					srcFormat = SIGNED_INT;
					break;
				case NT_I8 | NT_UNSIGNED:
					bytesPerPoint = sizeof(char);
					srcFormat = UNSIGNED_INT;
					break;
				case NT_I16:
					bytesPerPoint = sizeof(short);
					srcFormat = SIGNED_INT;
					break;
				case NT_I16 | NT_UNSIGNED:
					bytesPerPoint = sizeof(short);
					srcFormat = UNSIGNED_INT;
					break;
				case NT_I32:
					bytesPerPoint = sizeof(SInt32);
					srcFormat = SIGNED_INT;
					break;
				case NT_I32 | NT_UNSIGNED:
					bytesPerPoint = sizeof(UInt32);
					srcFormat = UNSIGNED_INT;
					break;
				case NT_FP32:
					bytesPerPoint = sizeof(float);
					srcFormat = IEEE_FLOAT;
					break;
				case NT_FP64:
					bytesPerPoint = sizeof(double);
					srcFormat = IEEE_FLOAT;
					break;
				default:
					return WAVE_TYPE_INCONSISTENT;		/* Corrupted wave. */
			}
			numNumbers = WavePoints(waveH);
			if (waveType & NT_CMPLX)
				numNumbers *= 2;
			ConvertData(WaveData(waveH), dPtr, numNumbers, bytesPerPoint, srcFormat, sizeof(double), IEEE_FLOAT);
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (int)CallBack2(MD_GETDPDATAFROMNUMERICWAVE, waveH, dPtr);
			break;
	}
	return 0;
}

/*	MDStoreDPDataInNumericWave(waveH, dPtr)
	
	MDStoreDPDataInNumericWave stores the data pointed to by dPtr in the specified wave.
	During the transfer, it converts the data from double precision to the numeric type
	of the wave. The conversion is done on-the-fly and the data pointed to by dPtr is not
	changed.
	
	This routine is intended for use with MDGetDPDataInNumericWave.
	
	The function result is zero or an error code.
	
	Thread Safety: MDStoreDPDataInNumericWave is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDStoreDPDataInNumericWave(waveHndl waveH, const double* dPtr)
{
	CountInt numNumbers;
	int bytesPerPoint;
	int waveType, type2;
	int destFormat;
	
	switch(WAVE_VERSION(waveH)) {
		case kWaveStruct3Version:
			waveType = WaveType(waveH);
			switch(waveType) {
				case TEXT_WAVE_TYPE:
					return NUMERIC_ACCESS_ON_TEXT_WAVE;
				case WAVE_TYPE:
				case DATAFOLDER_TYPE:
					return NUMERIC_ACCESS_ON_NON_NUMERIC_WAVE;
					break;
			}
			type2 = waveType & ~NT_CMPLX;
			switch(type2) {
				case NT_I8:
					bytesPerPoint = sizeof(char);
					destFormat = SIGNED_INT;
					break;
				case NT_I8 | NT_UNSIGNED:
					bytesPerPoint = sizeof(char);
					destFormat = UNSIGNED_INT;
					break;
				case NT_I16:
					bytesPerPoint = sizeof(short);
					destFormat = SIGNED_INT;
					break;
				case NT_I16 | NT_UNSIGNED:
					bytesPerPoint = sizeof(short);
					destFormat = UNSIGNED_INT;
					break;
				case NT_I32:
					bytesPerPoint = sizeof(SInt32);
					destFormat = SIGNED_INT;
					break;
				case NT_I32 | NT_UNSIGNED:
					bytesPerPoint = sizeof(UInt32);
					destFormat = UNSIGNED_INT;
					break;
				case NT_FP32:
					bytesPerPoint = sizeof(float);
					destFormat = IEEE_FLOAT;
					break;
				case NT_FP64:
					bytesPerPoint = sizeof(double);
					destFormat = IEEE_FLOAT;
					break;
				default:
					return WAVE_TYPE_INCONSISTENT;		/* Corrupted wave. */
			}
			numNumbers = WavePoints(waveH);
			if (waveType & NT_CMPLX)
				numNumbers *= 2;
			ConvertData(dPtr, WaveData(waveH), numNumbers, sizeof(double), IEEE_FLOAT, bytesPerPoint, destFormat);
			break;
		
		default:		/* Future version of Igor with a different wave structure. */
			return (int)CallBack2(MD_STOREDPDATAINNUMERICWAVE, waveH, (void*)dPtr);
			break;
	}
	return 0;
}

/*	MDGetTextWavePointValue(waveH, indices, textH)
	
	Returns via textH the value of a particular point in the specified wave.
	Any previous contents of textH are overwritten.
	
	If the wave is not a text wave, returns an error code and does not
	alter textH.
	
	indices is an array of dimension indices. For example, for a 3D wave,
		indices[0] should contain the row number
		indices[1] should contain the column number
		indices[2] should contain the layer number
	
	NOTE: This routine ignores indices for dimensions that do not exist in the wave.
	
	You must create textH before calling MDGetTextWavePointValue.
	For example:
		textH = NewHandle(0L);
	
	On output, if there is no error, textH contains a copy of the characters
	in the specified wave point. A point in an Igor text wave can contain
	any number of characters, including zero. Therefore, the handle can
	contain any number of characters. Igor text waves can contain any characters,
	including control characters. No characters codes are considered illegal.
	
	The characters in the handle are not null terminated. If you want
	to treat them as a C string, you should add a null character at the end.
	Make sure to remove any null termination if you pass the handle back to Igor.
		
	The function result is 0 or an error code.
	
	Thread Safety: MDGetTextWavePointValue is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDGetTextWavePointValue(waveHndl waveH, IndexInt indices[MAX_DIMENSIONS], Handle textH)
{
	return (int)CallBack3(MD_GETTEXTWAVEPOINTVALUE, waveH, indices, textH);
}

/*	MDSetTextWavePointValue(waveH, indices, textH)
	
	Transfers the characters in textH to the specified point in the specified
	wave. The contents of textH is not altered.
	
	If the wave is not a text wave, returns an error code.
	
	indices is an array of dimension indices. For example, for a 3D wave,
		indices[0] should contain the row number
		indices[1] should contain the column number
		indices[2] should contain the layer number
	
	NOTE: This routine ignores indices for dimensions that do not exist in the wave.
	
	A point in an Igor text wave can contain any number of characters, including zero.
	Therefore, the handle can contain any number of characters. Igor text waves can
	contain any characters, including control characters. No characters codes are
	considered illegal.
	
	The characters in the handle should not null terminated. If you have
	put a null terminator in the handle, remove it before calling MDSetTextWavePointValue.
	
	After calling MDSetTextWavePointValue, the handle is still yours so
	you should dispose it when you no longer need it.
		
	The function result is 0 or an error code.
	
	Thread Safety: MDSetTextWavePointValue is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
MDSetTextWavePointValue(waveHndl waveH, IndexInt indices[MAX_DIMENSIONS], Handle textH)
{
	return (int)CallBack3(MD_SETTEXTWAVEPOINTVALUE, waveH, indices, textH);
}

/*	GetTextWaveData(waveH, mode, textDataHPtr)

	Returns all of the text for the specified text wave via textDataHPtr.
	
	NOTE: This routine is for advanced programmers who are comfortable with
	pointer arithmetic and handles. Less experienced programmers should use
	MDGetTextWavePointValue to get the wave text values one at a time.
	
	If the function result is 0 then *textDataHPtr is a handle that you own.
	When you are finished, dispose of it using DisposeHandle.
	
	In the event of an error, the function result will be non-zero and
	*textDataHPtr will be NULL.
	
	The returned handle will contain the text for all of the wave's elements
	in one of several formats explained below. The format depends on the mode
	parameter.
	
	Modes 0 and 1 use a null byte to mark the end of a string and thus
	will not work if 0 is considered to be a legal character value.

	mode = 0
		The returned handle contains one C string (null-terminated) for each element
		of the wave.

		Example:
			"Zero"<null>
			"One"<null>
			"Two"<null>
	
	mode = 1
		The returned handle contains a list of 32-bit (in IGOR32) or 64-bit (in IGOR64)
		offsets to strings followed by the string data. There is one extra offset which
		is the offset to where the string would be for the next element if the wave had
		one more element.
		
		The text for each element in the wave is represented by a C string (null-terminated).
		
		Example:
			<Offset to "Zero">
			<Offset to "One">
			<Offset to "Two">
			<Extra offset>
			"Zero"<null>
			"One"<null>
			"Two"<null>
	
	mode = 2
		The returned handle contains a list of 32-bit (in IGOR32) or 64-bit (in IGOR64)
		offsets to strings followed by the string data.	There is one extra offset which
		is the offset to where the string would be for the next element if the wave had
		one more element.
		
		The text for each element in the wave is not null-terminated.
		
		Example:
			<Offset to "Zero">
			<Offset to "One">
			<Offset to "Two">
			<Extra offset>
			"Zero"
			"One"
			"Two"
			
	Using modes 1 and 2, you can determine the length of element i by subtracting
	offset i from offset i+1.
	
	You can convert the offsets into pointers to strings by adding
	**textDataHPtr to each of the offsets. However, since the handle in
	theory can be relocated in memory, you should lock the handle before
	converting to pointers and unlock it when you are done with it.
	
	For the purposes of GetTextWaveData, the wave is treated as a 1D wave
	regardless of its true dimensionality. If waveH a 2D text wave, the
	data returned via textDataHPtr is in column-major order. This means that
	the data for each row of the first column appears first in memory, followed
	by the data for the each row of the next column, and so on.
		
	Returns 0 or an error code.
	
	For an example using this routine, see TestGetAndSetTextWaveData below.
	
	Thread Safety: GetTextWaveData is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
GetTextWaveData(waveHndl waveH, int mode, Handle* textDataHPtr)
{
	return (int)CallBack3(GET_TEXT_WAVE_DATA, waveH, (void*)mode, textDataHPtr);
}

/*	SetTextWaveData(waveH, mode, textDataH)

	Sets all of the text for the specified text wave according to textDataH.
	
	NOTE: This routine is for advanced programmers who are comfortable with
	pointer arithmetic and handles. Less experienced programmers should use
	MDSetTextWavePointValue to set the wave text values one at a time.
	
	WARNING: If you pass inconsistent data in textDataH you will cause Igor to crash.
	
	SetTextWaveData can not change the number of points in the text wave.
	Therefore the data in textDataH must be consistent with the number of
	points in text wave. Otherwise a crash will occur.

	Also, when using modes 1 or 2, the offsets must be correct. Otherwise a
	crash will occur.
	
	Crashes caused by inconsistent data may occur at unpredictable times making
	it hard to trace it to the problem. So triple-check your code.
	
	You own the textDataH handle. When you are finished with it, dispose of it
	using DisposeHandle.
	
	The format of textDataH depends on the mode parameter. See the documentation
	for GetTextWaveData for a description of these formats.
	
	Modes 0 and 1 use a null byte to mark the end of a string and thus
	will not work if 0 is considered to be a legal character value.
		
	Returns 0 or an error code.
	
	For an example using this routine, see TestGetAndSetTextWaveData below.
	
	Thread Safety: SetTextWaveData is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
SetTextWaveData(waveHndl waveH, int mode, Handle textDataH)
{
	return (int)CallBack3(SET_TEXT_WAVE_DATA, waveH, (void*)mode, textDataH);
}

/*	TestGetAndSetTextWaveData(sourceWaveH, destWaveH, mode, echo)

	This routine is here just to give you and example of using GetTextWaveData
	and SetTextWaveData.
	
	If echo is true the contents of the source wave are printed in the history area.

	Then the data from the source wave is copied to the dest wave.
	
	Thread Safety: TestGetAndSetTextWaveData is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
#if 0	// This is just for testing
static int
TestGetAndSetTextWaveData(waveHndl sourceWaveH, waveHndl destWaveH, int mode, int echo)
{
	Handle textDataH;
	CountInt npnts;
	IndexInt* pTableOffset;
	char* pTextData;
	char message[256];
	int dataLen, prefixLen, availableBytes;
	int i;
	int err;
	
	if (WaveType(sourceWaveH) != TEXT_WAVE_TYPE)
		return TEXT_ACCESS_ON_NUMERIC_WAVE;
	
	if (WaveType(destWaveH) != TEXT_WAVE_TYPE)
		return TEXT_ACCESS_ON_NUMERIC_WAVE;
	
	npnts = WavePoints(sourceWaveH);
	
	if (err = ChangeWave(destWaveH, npnts, TEXT_WAVE_TYPE))
		return err;
	
	if (err = GetTextWaveData(sourceWaveH, mode, &textDataH))
		return err;
		
	pTableOffset = (PSInt*)*textDataH;					// Pointer to table of offsets if mode==1 or mode==2
	pTextData = *textDataH;
	if (mode > 0)
		pTextData += (npnts+1) * sizeof(PSInt);
		
	for(i=0; i<npnts; i+=1) {
		switch(mode) {
			case 0:
				dataLen = strlen(pTextData);
				break;
			
			case 1:
				dataLen = pTableOffset[i+1] - pTableOffset[i];
				dataLen -= 1;							// Exclude null terminator.
				break;

			case 2:
				dataLen = pTableOffset[i+1] - pTableOffset[i];
				break;	
		}

		if (echo) {
			sprintf(message, "Element %d: ", i);
			prefixLen = strlen(message);
			availableBytes = sizeof(message) - prefixLen - 1 - 1;		// Allow 1 for CR and 1 for null terminator.
			if (dataLen > availableBytes)
				dataLen = availableBytes;
			memcpy(message+prefixLen, pTextData, dataLen);
			message[prefixLen + dataLen] = 0x0D;
			message[prefixLen + dataLen+1] = 0;
			XOPNotice(message);
		}
		
		switch(mode) {
			case 0:
				pTextData += dataLen + 1;
				break;
			
			case 1:
				pTextData += dataLen + 1;
				break;

			case 2:
				pTextData += dataLen;
				break;	
		}		
	}
	
	err = SetTextWaveData(destWaveH, mode, textDataH);

	DisposeHandle(textDataH);
	return err;
}
#endif

/*	GetWaveDimensionLabels(waveH, dimLabelsHArray)

	dimLabelsHArray points to an array of MAX_DIMENSIONS handles. GetWaveDimensionLabels
	sets each element of this array to a handle containing dimension labels or to NULL.
	
	On output, if the function result is 0 (no error), dimLabelsHArray[i] will be a
	handle containing dimension labels for dimension i or NULL if dimension i has no
	dimension labels.
	
	If the function result is non-zero then all handles in dimLabelsHArray will be NULL.
	
	Any non-NULL output handles belong to you. Dispose of them with DisposeHandle
	when you are finished with them.
	
	For each dimension, the corresponding dimension label handle consists of
	an array of N+1 C strings, each in a field of (MAX_DIM_LABEL_CHARS+1) bytes.
	
	The first label is the overall dimension label for that dimension.
	
	Label i+1 is the dimension label for element i of the dimension.
	
	N is the smallest number such that the last non-empty dimension label
	for a given dimension and all dimension labels before it, whether empty
	or not, can be stored in the handle.
	
	For example, if a 5 point 1D wave has dimension labels for rows 0 and 2
	with all other dimension labels being empty then dimLabelsHArray[0] will
	contain four dimension labels, one for the overall dimension and three
	for rows 0 through 2. dimLabelsHArray[0] will not contain any storage
	for any point after row 2 because the remaining dimension labels for
	that dimension are empty.
		
	Returns 0 or an error code.
	
	For an example using this routine, see TestGetAndSetWaveDimensionLabels below.
	
	Thread Safety: GetWaveDimensionLabels is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
GetWaveDimensionLabels(waveHndl waveH, Handle dimLabelsHArray[MAX_DIMENSIONS])
{
	return (int)CallBack2(GET_WAVE_DIMENSION_LABELS, waveH, dimLabelsHArray);
}

/*	SetWaveDimensionLabels(waveH, dimLabelsHArray)

	dimLabelsHArray points to an array of MAX_DIMENSIONS handles. SetWaveDimensionLabels
	sets the dimension labels for each existing dimension of waveH based on the
	corresponding handle in dimLabelsHArray.
	
	The handles in dimLabelsHArray belong to you. Dispose of them with DisposeHandle
	when you are finished with them.
	
	See the documentation for GetWaveDimensionLabels for a discussion of how
	the dimension labels are stored in the handles.
		
	Returns 0 or an error code.
	
	For an example using this routine, see TestGetAndSetWaveDimensionLabels below.
	
	Thread Safety: SetWaveDimensionLabels is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
SetWaveDimensionLabels(waveHndl waveH, Handle dimLabelsHArray[MAX_DIMENSIONS])
{
	return (int)CallBack2(SET_WAVE_DIMENSION_LABELS, waveH, dimLabelsHArray);
}

/*	TestGetAndSetWaveDimensionLabels(sourceWaveH, destWaveH, echo)

	This routine is here just to give you and example of using GetWaveDimensionLabels
	and SetWaveDimensionLabels.
	
	If echo is true the source wave's dimension labels are printed in the history area.
	
	Then the dimension labels are copied from the source to the dest wave.
	
	Thread Safety: TestGetAndSetWaveDimensionLabels is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
#if 0	// This is just for testing
static int
TestGetAndSetWaveDimensionLabels(waveHndl sourceWaveH, waveHndl destWaveH, int echo)
{
	Handle dimLabelsH;
	CountInt dimSizes[MAX_DIMENSIONS+1];
	CountInt numLabels;
	IndexInt element;
	int dim, numDimensions;
	Handle dimLabelsHArray[MAX_DIMENSIONS];
	char label[MAX_DIM_LABEL_CHARS+1];
	char message[256];
	int err;
	
	if (err = MDGetWaveDimensions(sourceWaveH, &numDimensions, dimSizes))
		return err;
	
	if (err = MDChangeWave(destWaveH, WaveType(destWaveH), dimSizes))
		return err;
	
	if (err = GetWaveDimensionLabels(sourceWaveH, dimLabelsHArray))
		return err;
		
	if (echo) {
		for(dim=0; dim<numDimensions; dim+=1) {
			dimLabelsH = dimLabelsHArray[dim];
			if (dimLabelsH == NULL) {
				sprintf(message, "Dimension %d has no labels"CR_STR, dim);
				XOPNotice(message);
			}
			else {
				numLabels = GetHandleSize(dimLabelsH) / (MAX_DIM_LABEL_CHARS+1);
				for(element=-1; element<(numLabels-1); element++) {
					strcpy(label, *dimLabelsH + (element+1)*(MAX_DIM_LABEL_CHARS+1));
					sprintf(message, "Dimension %d, element %ld = '%s'"CR_STR, dim, element, label);
					XOPNotice(message);				
				}
			}
		}
	}
		
	err = SetWaveDimensionLabels(destWaveH, dimLabelsHArray);

	// We own the label handles and thus must dispose of them.
	for(dim=0; dim<numDimensions; dim+=1) {
		dimLabelsH = dimLabelsHArray[dim];
		if (dimLabelsH != NULL)
			DisposeHandle(dimLabelsH);
	}
	
	return err;
}
#endif

/*	HoldWave(waveH, waveRefPtr)

	Tells Igor that you are holding a wave and that it should not be killed.
	
	waveH is a wave handle that you have obtained from Igor.
	
	waveRefPtr contains the address of a waveHndl variable that you will use
	to refer to the wave going forward.
	
	HoldWave calls Igor to notify it that you are holding the wave and stores
	waveH in the waveHndl pointed to by waveRefPtr. Use the waveHndl variable
	pointed to by waveRefPtr to refer to the wave going forward.
	
	An XOP may "hold" a wave handle. "Hold" means that you are storing a waveHndl
	over a period of time during which Igor could possibly kill it. HoldWave
	allows you to tell Igor that the wave should not be killed until further notice.
	
	For example, a data acquisition XOP that stores data in a wave during IDLE
	messages would typically store the wave handle in a global variable and write
	to that wave each time it receives an IDLE message. It is important that Igor
	not kill a wave while the XOP is holding it.

	Before killing a wave, Igor sends an OBJINUSE message to each XOP. If the XOP
	reports that the wave is in use (i.e., the XOP is holding the wave handle), Igor
	refuses to kill it. However, the OBJINUSE message is not sent to XOPs when
	a preemptive thread is executing because the XOP message protocol is not threadsafe.
	
	Igor Pro 6 and later implement a more efficient method for preventing the untimely
	killing of a wave: wave reference counting. Igor keeps a reference count for
	each wave. This count is used to determine if a wave is in use, for example,
	used in a graph or table.

	As of Igor Pro 6.20, XOPs can participate in wave reference counting.
	If your XOP obtains a wave handle and holds it after the XOP returns to Igor,
	you should call HoldWave to increment Igor's internal wave reference for that
	wave. When you no longer need to access the wave, you must call ReleaseWave.
	ReleaseWave decrements the wave reference count and, if the count reaches zero,
	kills the wave.

	For example:

	static waveHndl gWaveH = NULL;	// Global to hold data acq wave

	int StartAcq()		// Called from your StartAcq operation or function
	{
		waveHndle waveH;
		
		waveH = FetchWave("jack");
		if (waveH == NULL)
		 return NOWAV;
		
		// Tell Igor we are holding this wave handle
		HoldWave(waveH, &gWaveH);	// gWaveH now holds wave reference 
		
		// Global gWaveH now refers to the wave
		waveH = NULL;	// Make it clear that waveH is no longer to be used.
		. . .
	}

	int DoAcq()			// Called from IDLE message
	{
		if (gWaveH == NULL)
			return NOWAV;		// gWaveH should not be NULL during data acq
		. . .
	}
 
	int EndAcq()		// Called from your EndAcq operation or function
	{
		if (gWaveH == NULL)
			return DATA_ACQ_NOT_RUNNING;
		
		// Tell Igor we are no longer holding this wave handle
		ReleaseWave(&gWaveH);	// NOTE: Sets gWaveH to NULL
		. . .
	}

	Note that ReleaseWave sets your wave handle (gWaveH in this example)
	to NULL and it is no longer valid.
	
	You also should call HoldWave and ReleaseWave if you are doing a callback
	to Igor and the callback could possibly kill the wave - for example, if
	you are doing an XOPCommand callback to run a user-defined function that
	could kill the wave. This constitutes "holding the wave" because you are
	holding a reference to the wave over a period of time in which it
	could be killed. Therefore you should use HoldWave to tell Igor you are
	using it and ReleaseWave when you are no longer using it.

	If you are just using the wave handle temporarily during the execution of your
	external function or operation and you make no calls that could potentially
	kill the wave then you do not need to and should not call HoldWave and ReleaseWave.

	Once you have called HoldWave on a wave, Igor will not allow the user to kill
	it until ReleaseWave has been called on it. Normally you will call ReleaseWave
	yourself but here is a case where Igor will call ReleaseWave on it. If you
	define an external function with a structure parameter and the structure contains
	a WAVE field and your XOP returns a wave handle to the calling user-defined function
	by setting the WAVE field, you must call HoldWave passing a pointer to the WAVE field
	as the waveRefPtr parameter. This signifies that the WAVE field is holding the wave.
	When the calling user-defined function returns and the structure goes out of scope,
	Igor will automatically call ReleaseWave on the structure WAVE field. In this example,
	by calling HoldWave with a pointer to the structure WAVE field, you are passing
	ownership of the wave to the calling user-defined function which is then responsible
	for calling ReleaseWave.
		
	When an experiment is closed, either by doing New Experiment, opening another
	experiment or quitting Igor, all waves are killed whether you are holding them
	or not. Consequently, you must do something like this when your XOPEntry
	routine receives the CLEANUP message:
	
	case CLEANUP:
		// Wave about to be killed. Make sure global reference is NULL.
		if (gWaveH != NULL)
			ReleaseWave(&gWaveH);

	If your XOP will run only with Igor Pro 6.20 or later then you can implement
	just reference counting using HoldWave and ReleaseWave and you do not need to
	respond to the OBJINUSE message. If your XOP will run with earlier versions
	of Igor, you should implement reference counting and also respond to the
	OBJINUSE message.

	Returns 0 or IGOR_OBSOLETE.
	
	Added for Igor Pro 6.20. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
	
	Thread Safety: HoldWave is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
HoldWave(waveHndl waveH, waveHndl* waveRefPtr)
{
	return (int)CallBack2(HOLD_WAVE, waveH, waveRefPtr);
}

/*	ReleaseWave(waveRefPtr)

	Tells Igor that you are no longer holding a wave.
	
	waveRefPtr contains the address of your waveHndl variable that refers to a wave.
	ReleaseWave sets *waveRefPtr to NULL so your waveHndl variable is not valid
	after you call ReleaseWave.
	
	See HoldWave for a detailed discussion.

	Returns 0 or IGOR_OBSOLETE.
	
	Added for Igor Pro 6.20. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
	
	Thread Safety: ReleaseWave is thread-safe with Igor Pro 6.20 or later except for
	waves passed to threads as parameters.
*/
int
ReleaseWave(waveHndl* waveRefPtr)
{
	return (int)CallBack1(RELEASE_WAVE, waveRefPtr);
}

