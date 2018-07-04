/*	XOPDataFolderAccess.c
	
	Routines for Igor XOPs that provide access to Igor Pro data folders.
*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h

// *** Data Folder Access Routines ***

/*	GetDataFolderNameOrPath(dataFolderH, flags, dataFolderPathOrName)

	Given a handle to a data folder, returns
		the name of the folder if bit 0 of flags is zero
		a full path to the folder if bit 0 of flags is set
	
	If bit 1 of flags is set, Igor returns the dataFolderPathOrName with
	single quotes if they would be needed to use the name or path in Igor's
	command line.
	
	If bit 1 if flags is zero, dataFolderPathOrName will have no quotes.
	
	Set bit 1 of flags if you are going to use the path or name in a command
	that you submit to Igor via the XOPCommand or XOPSilentCommand callbacks.
	Clear bit 1 of flags for other purpose, for example, if you are getting the
	name or path just to display to the user.
	
	All other bits of the flags parameter are reserved; you must set them to zero.
	
	If dataFolderH is NULL, it uses the current data folder.
	
	A data folder name can be up to MAX_OBJ_NAME characters while a full path
	can be up to MAXCMDLEN characters. To be safe, allocate MAXCMDLEN+1 characters
	for dataFolderPathOrName.
	
	Returns 0 or error code.
	
	Thread Safety: GetDataFolderNameOrPath is thread-safe with Igor Pro 6.20 or later.
*/
int
GetDataFolderNameOrPath(DataFolderHandle dataFolderH, int flags, char dataFolderPathOrName[MAXCMDLEN+1])
{
	return (int)CallBack3(GET_DATAFOLDER_NAMEORPATH, dataFolderH, (void*)flags, dataFolderPathOrName);
}

/*	GetDataFolderIDNumber(dataFolderH, IDNumberPtr)

	Returns the unique ID number for the data folder via *IDNumberPtr.
	
	If dataFolderH is NULL, it uses the current data folder.
	
	Each data folder has a unique ID number that stays the same as long as the data
	folder exists, even if it is renamed or moved. If you need to reference a data
	folder over a period of time during which it could be killed, then you should
	store the data folder's ID number.
	
	Given the ID number, you can call GetDataFolderByIDNumber to check if the data
	folder still exists and to get a handle it.
	
	The ID number is valid until the user creates a new Igor experiment or quits Igor.
	ID numbers are not remembered from one running of Igor to the next.
	
	Thread Safety: GetDataFolderIDNumber is thread-safe with Igor Pro 6.20 or later.
*/
int
GetDataFolderIDNumber(DataFolderHandle dataFolderH, int* IDNumberPtr)
{
	return (int)CallBack2(GET_DATAFOLDER_IDNUMBER, dataFolderH, IDNumberPtr);
}

/*	GetDataFolderProperties(dataFolderH, propertiesPtr)

	Returns the bit-flag properties of the specified data folder.
	
	If dataFolderH is NULL, it uses the current data folder.
	
	At present, Igor does not support any properties and this routine will always return 0
	in *propertiesPtr. In the future, it might support properties such as "locked".
	
	Returns 0 or error code.
	
	Thread Safety: GetDataFolderProperties is thread-safe with Igor Pro 6.20 or later.
*/
int
GetDataFolderProperties(DataFolderHandle dataFolderH, int* propertiesPtr)
{
	return (int)CallBack2(GET_DATAFOLDER_PROPERTIES, dataFolderH, propertiesPtr);
}

/*	SetDataFolderProperties(dataFolderH, properties)

	Sets the bit-flag properties of the specified data folder.
	
	If dataFolderH is NULL, it uses the current data folder.
	
	At present, Igor does not support any properties and there is no reason
	to call this routine. It will return a -1 error code for any value
	of properties other than 0.
	
	In the future, it might support properties such as "locked".
	
	Returns 0 or error code.
	
	Thread Safety: SetDataFolderProperties is thread-safe with Igor Pro 6.20 or later.
*/
int
SetDataFolderProperties(DataFolderHandle dataFolderH, int properties)
{
	return (int)CallBack2(SET_DATAFOLDER_PROPERTIES, dataFolderH, (void*)properties);
}

/*	GetDataFolderListing(dataFolderH, optionsFlag, h)

	Returns via the handle h a listing of the contents of the specified data folder.
	You must create h and dispose it when you no longer need it. Its contents will
	be replaced by the listing.
	
	The listing does not include a null terminator character. Use GetHandleSize to
	find the length of the text and add a null terminator if you want to treat it
	as a C string.
	
	If dataFolderH is NULL, it uses the current data folder.
	
	optionsFlag determines what is in the listing.
		If bit 0 of optionsFlag is set:
			includes "FOLDERS:<subFolder0>,<subFolder1>...,<subFolderN>;<CR>"
		If bit 1 of optionsFlag is set:
			includes "WAVES:<waveName0>,<waveName1>...,<waveNameN>;<CR>"
		If bit 2 of optionsFlag is set:
			includes "VARIABLES:<variableName0>,<variableName1>...,<variableNameN>;<CR>"
		If bit 3 of optionsFlag is set:
			includes "STRINGS:<stringVariableName0>,<stringVariableName1>...,<stringVariableNameN>;<CR>"
		All other bits are reserved and should be set to zero.

	Names in the listing of waves, variables and strings are quoted with single
	quotes if this is necessary to make them suitable for use in the Igor command line.

	Returns 0 or error code.

	Thread Safety: GetDataFolderListing is thread-safe with Igor Pro 6.20 or later.
*/
int
GetDataFolderListing(DataFolderHandle dataFolderH, int optionsFlag, Handle h)
{
	return (int)CallBack3(GET_DATAFOLDER_LISTING, dataFolderH, (void*)optionsFlag, h);
}

/*	int GetRootDataFolder(refNum, rootDataFolderHPtr);

	Returns a handle to the root data folder.
	
	Data folder handles belong to Igor so you should not modify or dispose them.
	
	NOTE: Always pass 0 for refNum. It is for future use.

	Returns 0 or error code.

	Thread Safety: GetRootDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
GetRootDataFolder(int refNum, DataFolderHandle* rootDataFolderHPtr)
{
	return (int)CallBack2(GETROOT_DATAFOLDER, (void*)refNum, rootDataFolderHPtr);
}

/*	GetCurrentDataFolder(currentDataFolderHPtr)

	Returns a handle to the current data folder in *currentDataFolderHPtr.
	
	Data folder handles belong to Igor so you should not modify or dispose them.
	
	The only use for this handle is to pass it to other data-folder-related
	XOPSupport routines.
	
	Returns 0 or error code.

	Thread Safety: GetCurrentDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
GetCurrentDataFolder(DataFolderHandle* currentDataFolderHPtr)
{
	return (int)CallBack1(GETCURRENT_DATAFOLDER, currentDataFolderHPtr);
}

/*	SetCurrentDataFolder(dataFolderH)

	Sets the current data folder to the data folder referenced by dataFolderH.

	Returns 0 or error code.

	Thread Safety: SetCurrentDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
SetCurrentDataFolder(DataFolderHandle dataFolderH)
{
	return (int)CallBack1(SETCURRENT_DATAFOLDER, dataFolderH);
}

/*	GetNamedDataFolder(startingDataFolderH, dataFolderPath, dataFolderHPtr)

	Returns in *dataFolderHPtr the data folder specified by startingDataFolderH and dataFolderPath.

	Data folder handles belong to Igor so you should not modify or dispose them.
	
	dataFolderPath can be an absolute path (e.g., "root:FolderA:FolderB:"), a relative
	path (e.g., ":FolderA:FolderB:") or a data folder name (e.g., "FolderA").
	
	If dataFolderPath is an absolute path then startingDataFolderH is immaterial - you
	can pass any data folder handle or NULL. An absolute path must always start with "root:".
	It should include a trailing colon but GetNamedDataFolder tolerates an absolute path
	without the trailing colon. Note that "root" is an not an absolute path whereas "root:" is.

	If dataFolderPath is a relative path or a data folder name then dataFolderPath is relative
	to startingDataFolderH. However, if startingDataFolderH is NULL then dataFolderPath is relative
	to the current folder.
	
	Passing "root" as dataFolderPath will not find the root data folder because "root" is a
	data folder name. Igor will try to find a data folder named "root" relative to the current
	data folder. The actual root data folder is never relative to any data folder so it can
	not be found this way. Use "root:" instead.

	Returns 0 or error code.

	Thread Safety: GetNamedDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
GetNamedDataFolder(DataFolderHandle startingDataFolderH, const char dataFolderPath[MAXCMDLEN+1], DataFolderHandle* dataFolderHPtr)
{
	return (int)CallBack3(GETNAMED_DATAFOLDER, startingDataFolderH, (void*)dataFolderPath, dataFolderHPtr);
}

/*	GetDataFolderByIDNumber(IDNumber, dataFolderHPtr)

	Returns via *dataFolderHPtr the data folder handle associated with the
	specified ID number.
	
	Data folder handles belong to Igor so you should not modify or dispose them.
	
	Returns 0 if OK or a non-zero error code if the data folder doesn't exist,
	which would be the case if the data folder were killed since you got its
	ID number.
	
	Each data folder has a unique ID number that stays the same as long as the data
	folder exists. You can get the ID number for a data folder using GetDataFolderIDNumber().
	
	If you need to reference a data folder over a period of time during which it
	could be killed, then you should store the data folder's ID number. Given the ID
	number, GetDataFolderByIDNumber tells you if the data folder still exists and,
	if it does, gives you the data folder handle.
	
	The ID number is valid until the user creates a new Igor experiment or quits Igor.
	ID numbers are not remembered from one running of Igor to the next.

	Thread Safety: GetDataFolderByIDNumber is thread-safe with Igor Pro 6.20 or later.
*/
int
GetDataFolderByIDNumber(int IDNumber, DataFolderHandle* dataFolderHPtr)
{
	return (int)CallBack2(GET_DATAFOLDER_BYIDNUMBER, (void*)IDNumber, dataFolderHPtr);
}

/*	GetParentDataFolder(dataFolderH, parentFolderHPtr)

	Returns the parent of the specified data folder via *parentFolderHPtr.
	
	Data folder handles belong to Igor so you should not modify or dispose them.
	
	If dataFolderH is NULL, it uses the current data folder.
	
	Passing the root data folder as dataFolderH is an error. In this
	case GetParentDataFolder returns NO_PARENT_DATAFOLDER.
	
	Returns 0 or error code.

	Thread Safety: GetParentDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
GetParentDataFolder(DataFolderHandle dataFolderH, DataFolderHandle* parentFolderHPtr)
{
	return (int)CallBack2(GETPARENT_DATAFOLDER, dataFolderH, parentFolderHPtr);
}

/*	int GetNumChildDataFolders(parentDataFolderH, numChildDataFolderPtr)

	Returns the number of child data folders in the specified parent data folder.

	If parentDataFolderH is NULL, it uses the current data folder.
	
	Returns 0 or error code.

	Thread Safety: GetNumChildDataFolders is thread-safe with Igor Pro 6.20 or later.
*/
int
GetNumChildDataFolders(DataFolderHandle parentDataFolderH, int* numChildDataFolderPtr)
{
	return (int)CallBack2(GETNUMCHILD_DATAFOLDERS, parentDataFolderH, (void*)numChildDataFolderPtr);
}

/*	int GetIndexedChildDataFolder(parentDataFolderH, index, childDataFolderHPtr)

	Returns a handle to the child data folder specified by the index.
	index starts from 0.
	
	Data folder handles belong to Igor so you should not modify or dispose them.
	
	If parentDataFolderH is NULL, it uses the current data folder.
	
	Returns 0 or error code.

	Thread Safety: GetIndexedChildDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
GetIndexedChildDataFolder(DataFolderHandle parentDataFolderH, int index, DataFolderHandle* childDataFolderHPtr)
{
	return (int)CallBack3(GETINDEXEDCHILD_DATAFOLDER, parentDataFolderH, (void*)index, childDataFolderHPtr);
}

/*	GetWavesDataFolder(waveH, dataFolderHPtr)

	Returns the handle to the data folder containing the specified wave.
	
	Data folder handles belong to Igor so you should not modify or dispose them.
	
	Returns 0 or error code.

	Thread Safety: GetWavesDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
GetWavesDataFolder(waveHndl waveH, DataFolderHandle* dataFolderHPtr)
{
	return (int)CallBack2(GETWAVES_DATAFOLDER, waveH, dataFolderHPtr);
}

/*	NewDataFolder(parentFolderH, newDataFolderName, newDataFolderHPtr)

	Creates a new data folder in the data folder specified by parentFolderH.
	
	parentFolderH can be
		a handle to an Igor data folder
		NULL to use the current data folder
	
	On output, *newDataFolderHPtr will contain a handle to the new data folder
	or NULL if an error occurred.
	
	Data folder handles belong to Igor so you should not modify or dispose them.
	
	NewDataFolder does not change the current data folder. If you want to make the
	new folder the current folder, call SetCurrentDataFolder after NewDataFolder.
	
	Returns 0 or error code.

	Thread Safety: NewDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
NewDataFolder(DataFolderHandle parentFolderH, const char newDataFolderName[MAX_OBJ_NAME+1], DataFolderHandle* newDataFolderHPtr)
{
	return (int)CallBack3(NEW_DATAFOLDER, parentFolderH, (void*)newDataFolderName, newDataFolderHPtr);
}

/*	KillDataFolder(dataFolderH)

	Kills an existing data folder, removing it and its contents, including
	any child data folders, from memory.
	
	dataFolderH is a handle to an existing Igor data folder or NULL to use
	the current data folder.

	You will receive an error and the data folder will not be killed if it
	contains waves or variables that are in use (e.g. displayed in tables or graphs).
	
	If you kill the current data folder or a data folder that contains the current
	data folder, Igor will set the current data folder to the parent of the killed
	data folder.
	
	If you kill the root data folder, its contents will be killed but not
	the root data folder itself.
	
	Returns 0 or error code.
	
	NOTE: Once a data folder is successfully killed, dataFolderH is no longer
		  valid. You should not reference it for any purpose.

	Thread Safety: KillDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
KillDataFolder(DataFolderHandle dataFolderH)
{
	return (int)CallBack1(KILL_DATAFOLDER, dataFolderH);
}

/*	DuplicateDataFolder(sourceDataFolderH, parentDataFolderH, newDataFolderName)

	Creates a clone of the source data folder. The contents of the destination
	will be clones of the contents of the source.
	
	sourceDataFolderH is a handle to the data folder to be duplicated
	or NULL to use the current data folder.

	parentDataFolderH is a handle to the data folder in which the new data folder
	is to be created or NULL to use the current data folder.
	
	newDataFolderName is the name to be given to the new data folder.

	Returns 0 or error code.

	Thread Safety: DuplicateDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
DuplicateDataFolder(DataFolderHandle sourceDataFolderH, DataFolderHandle parentDataFolderH, const char newDataFolderName[MAX_OBJ_NAME+1])
{
	return (int)CallBack3(DUPLICATE_DATAFOLDER, sourceDataFolderH, parentDataFolderH, (void*)newDataFolderName);
}

/*	MoveDataFolder(sourceDataFolderH, newParentDataFolderH)
	
	Moves the source data folder into a new location in the hierarchy.
	It is an error to attempt to move a parent folder into itself or
	one of its children.
	
	sourceDataFolderH is a handle to the data folder to be moved
	or NULL to use the current data folder.

	newParentDataFolderH is a handle to the data folder in which the source data folder
	is to be moved or NULL to use the current data folder.

	Returns 0 or error code.

	Thread Safety: MoveDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
MoveDataFolder(DataFolderHandle sourceDataFolderH, DataFolderHandle newParentDataFolderH)
{
	return (int)CallBack2(MOVE_DATAFOLDER, sourceDataFolderH, newParentDataFolderH);
}

/*	RenameDataFolder(dataFolderH, newName)
	
	Renames the data folder.
	
	dataFolderH is a handle to the data folder to be renamed or NULL to use
	the current data folder.

	Returns 0 or error code.

	Thread Safety: RenameDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
RenameDataFolder(DataFolderHandle dataFolderH, const char newName[MAX_OBJ_NAME+1])
{
	return (int)CallBack2(RENAME_DATAFOLDER, dataFolderH, (void*)newName);
}

/*	GetNumDataFolderObjects(dataFolderH, int objectType, numObjectsPtr)

	Returns via numObjectsPtr the number of objects of the specified type in the
	specified data folder.

	If dataFolderH is NULL, it uses the current data folder.
	
	objectType is one of the following:
		WAVE_OBJECT		for waves
		VAR_OBJECT		for numeric variables
		STR_OBJECT		for string variables
	
	Returns 0 or error code.

	Thread Safety: GetNumDataFolderObjects is thread-safe with Igor Pro 6.20 or later.
*/
int
GetNumDataFolderObjects(DataFolderHandle dataFolderH, int objectType, int* numObjectsPtr)
{
	return (int)CallBack3(GETNUM_DATAFOLDER_OBJECTS, dataFolderH, (void*)objectType, (void*)numObjectsPtr);
}

/*	GetIndexedDataFolderObject(dataFolderH, objectType, index, objectName, objectValuePtr)

	Returns information that allows you to access an object of the specified type in the
	specified data folder.

	index starts from 0.

	If dataFolderH is NULL, it uses the current data folder.
	
	objectType is one of the following:
		WAVE_OBJECT		for waves
		VAR_OBJECT		for numeric variables
		STR_OBJECT		for string variables
	
	You can pass NULL for objectName if you don't need to know the name of the object.

	If you do not want to get the "value" of the object, pass NULL for objectValuePtr.
	If objectValuePtr is not NULL, GetIndexedDataFolderObject sets fields depending
	on the object's type:
		WAVE_OBJECT		sets objectValuePtr->waveH field to wave's handle
		VAR_OBJECT		stores numeric variable's value in objectValuePtr->nv field
		STR_OBJECT		sets objectValuePtr->strH field to strings's handle
	
	The handles returned via the waveH and strH fields belong to Igor. Do not
	modify or dispose them.
	
	Returns 0 or error code.

	Thread Safety: GetIndexedDataFolderObject is thread-safe with Igor Pro 6.20 or later.
*/
int
GetIndexedDataFolderObject(DataFolderHandle dataFolderH, int objectType, int index, char objectName[MAX_OBJ_NAME+1], DataObjectValuePtr objectValuePtr)
{
	return (int)CallBack5(GETINDEXED_DATAFOLDER_OBJECT, dataFolderH, (void*)objectType, (void*)index, objectName, objectValuePtr);
}

/*	GetDataFolderObject(dataFolderH, objectName, objectTypePtr, objectValuePtr)
	
	Returns 0 if the object exists or an error code if it does not exist.
	
	Returns the object type (WAVE_OBJECT, VAR_OBJECT, STR_OBJECT or DATAFOLDER_OBJECT
	via objectTypePtr.

	If dataFolderH is NULL, it uses the current data folder.
	
	If objectValuePtr is not NULL, it returns:
		For WAVE_OBJECT:		A handle to the wave via the objectValuePtr->waveH field.
		For VAR_OBJECT:			The value of the variable via the objectValuePtr->nv field.
		For STR_OBJECT:			A handle to the strings contents via the objectValuePtr->strH field.
		For DATAFOLDER_OBJECT:	A handle to the data folder via the objectValuePtr->dfH field.
	
	Note that the wave, string and data folder handles returned belong to Igor.
	The caller should not modify or dispose them.
	
	Remember also that strings in handles do not contain a null terminator (they are not C strings).
	To find the number of characters, call GetHandleSize on the handle.

	Thread Safety: GetDataFolderObject is thread-safe with Igor Pro 6.20 or later.
*/
int
GetDataFolderObject(DataFolderHandle dataFolderH, const char objectName[MAX_OBJ_NAME+1], int* objectTypePtr, DataObjectValuePtr objectValuePtr)
{
	int err;

	err= (int)CallBack4(GET_DATAFOLDER_OBJECT, dataFolderH, (void*)objectName, objectTypePtr, objectValuePtr);
	return err;
}

/*	SetDataFolderObject(dataFolderH, objectName, objectType, objectValuePtr)
	
	Returns 0 if the object exists or an error code if it does not exist
	or is of the wrong type.

	If dataFolderH is NULL, it uses the current data folder.
	
	If objectValuePtr is not NULL, it does the following:
		For WAVE_OBJECT:		Does nothing.
		
		For VAR_OBJECT:			Sets the value of the numeric variable based on objectValuePtr->nv.
								Note that you can't change between real and complex by changing
								the nv numType field.
								
		For STR_OBJECT:			Sets the value of the string variable based on objectValuePtr->strH.
								Igor just copies the data from the handle. The handle is yours to dispose
								(if you created it).
								Remember also that strings in handles do not contain a null terminator
								(they are not C strings). To find the number of characters, call
								GetHandleSize on the handle.
								
		For DATAFOLDER_OBJECT:	Does nothing.

	Thread Safety: SetDataFolderObject is thread-safe with Igor Pro 6.20 or later.
*/
int
SetDataFolderObject(DataFolderHandle dataFolderH, const char objectName[MAX_OBJ_NAME+1], int objectType, DataObjectValuePtr objectValuePtr)
{
	return (int)CallBack4(SET_DATAFOLDER_OBJECT, dataFolderH, (void*)objectName, (void*)objectType, objectValuePtr);
}

/*	KillDataFolderObject(dataFolderH, objectType, objectName)

	Kills the named object of the specified type in the specified data folder.

	If dataFolderH is NULL, it uses the current data folder.
	
	objectType is one of the following:
		WAVE_OBJECT		for waves
		VAR_OBJECT		for numeric variables
		STR_OBJECT		for string variables
	
	NOTE: If you attempt to kill a wave that is in use (e.g. in a graph, table or user function)
		  the wave will not be killed and you will receive a non-zero error code.
		  Igor does not check if numeric and string variables are in use. You can
		  kill a numeric or string variable at any time without receiving an error.
	
	Returns 0 or error code.

	Thread Safety: KillDataFolderObject is thread-safe with Igor Pro 6.20 or later.
*/
int
KillDataFolderObject(DataFolderHandle dataFolderH, int objectType, const char objectName[MAX_OBJ_NAME+1])
{
	return (int)CallBack3(KILL_DATAFOLDER_OBJECT, dataFolderH, (void*)objectType, (void*)objectName);
}

/*	MoveDataFolderObject(sourceDataFolderH, objectType, objectName, destDataFolderH)

	Moves the named object of the specified type from the source data folder to the
	destination data folder.

	If sourceDataFolderH is NULL, it uses the current data folder.

	If destDataFolderH is NULL, it uses the current data folder.
	
	objectType is one of the following:
		WAVE_OBJECT		for waves
		VAR_OBJECT		for numeric variables
		STR_OBJECT		for string variables
	
	NOTE: If an object with the same name exists in the destination data folder,
		  the object will not be moved and you will receive a non-zero error code.
	
	Returns 0 or error code.

	Thread Safety: MoveDataFolderObject is thread-safe with Igor Pro 6.20 or later.
*/
int
MoveDataFolderObject(DataFolderHandle sourceDataFolderH, int objectType, const char objectName[MAX_OBJ_NAME+1], DataFolderHandle destDataFolderH)
{
	return (int)CallBack4(MOVE_DATAFOLDER_OBJECT, sourceDataFolderH, (void*)objectType, (void*)objectName, destDataFolderH);
}

/*	RenameDataFolderObject(dataFolderH, objectType, objectName, newObjectName)

	Renames the named object of the specified type in the specified data folder.

	If dataFolderH is NULL, it uses the current data folder.
	
	objectType is one of the following:
		WAVE_OBJECT		for waves
		VAR_OBJECT		for numeric variables
		STR_OBJECT		for string variables
	
	NOTE: If the new name is illegal or in use the object will not be renamed and you will
		  receive a non-zero error code.
	
	Returns 0 or error code.

	Thread Safety: RenameDataFolderObject is thread-safe with Igor Pro 6.20 or later.
*/
int
RenameDataFolderObject(DataFolderHandle dataFolderH, int objectType, const char objectName[MAX_OBJ_NAME+1], const char newObjectName[MAX_OBJ_NAME+1])
{
	return (int)CallBack4(RENAME_DATAFOLDER_OBJECT, dataFolderH, (void*)objectType, (void*)objectName, (void*)newObjectName);
}

/*	DuplicateDataFolderObject(dataFolderH, objectType, objectName, destFolderH, newObjectName, overwrite)

	Duplicates the named object of the specified type.
	
	If dataFolderH and/or destFolderH is NULL, it uses the current data folder.
	
	objectType is one of the following:
		WAVE_OBJECT		for waves
		VAR_OBJECT		for numeric variables
		STR_OBJECT		for string variables
	
	If the new name is illegal you will receive a non-zero error code.
	If the new name is in use and overwrite is false, you will receive a non-zero error code.
	If the new name is in use for a different kind of object, you will receive a non-zero error code.
	To avoid these errors, you can check and if necessary fix the new name using the CheckName,
	CleanupName and UniqueName2 routines.

	Returns 0 or error code.

	Thread Safety: DuplicateDataFolderObject is thread-safe with Igor Pro 6.20 or later.
*/
int
DuplicateDataFolderObject(
	DataFolderHandle dataFolderH, int objectType, const char objectName[MAX_OBJ_NAME+1],
	DataFolderHandle destFolderH, const char newObjectName[MAX_OBJ_NAME+1], int overwrite)
{
	return (int)CallBack6(DUPLICATE_DATAFOLDER_OBJECT, dataFolderH, (void*)objectType, (void*)objectName, destFolderH, (void*)newObjectName, (void*)overwrite);
}

/*	ClearDataFolderFlags()

	This routine is for use by the WaveMetrics Data Browser only.

	Thread Safety: ClearDataFolderFlags is not thread-safe.
*/
void
ClearDataFolderFlags(void)
{
	if (!CheckRunningInMainThread("ClearDataFolderFlags"))
		return;
	
	CallBack0(CLEAR_DATAFOLDER_FLAGS);
}

/*	GetDataFolderChangesCount()

	This routine is for use by the WaveMetrics Data Browser only.

	Thread Safety: GetDataFolderChangesCount is not thread-safe.
*/
int
GetDataFolderChangesCount(void)
{
	if (!CheckRunningInMainThread("GetDataFolderChangesCount"))
		return 0;
	
	return (int)CallBack0(GET_DATAFOLDER_CHANGESCOUNT);
}

/*	GetDataFolderChangeFlags(dataFolderH, flagsP)

	This routine is for use by the WaveMetrics Data Browser only.

	Thread Safety: GetDataFolderChangeFlags is not thread-safe.
*/
int
GetDataFolderChangeFlags(DataFolderHandle dataFolderH, int *flagsP)
{
	if (!CheckRunningInMainThread("GetDataFolderChangeFlags"))
		return 0;
	
	return (int)CallBack2(GET_DATAFOLDER_CHANGEFLAGS, dataFolderH, flagsP);
}

/*	HoldDataFolder(dfH, dfRefPtr)

	Tells Igor that you are holding a data folder handle and that the corresponding
	data folder should not be killed.
	
	dfH is a data folder handle that you have obtained from Igor.
	
	dfRefPtr contains the address of a DataFolderHandle variable that you will use
	to refer to the data folder going forward.
	
	HoldDataFolder calls Igor to notify it that you are holding the data folder
	and stores dfH in the DataFolderHandle pointed to by dfRefPtr. Use the
	DataFolderHandle variable pointed to by dfRefPtr to refer to the data
	folder going forward.
	
	An XOP may "hold" a data folder handle. "Hold" means that you are storing a
	DataFolderHandle over a period of time during which Igor could possibly kill it.
	HoldDataFolder allows you to tell Igor that the data folder should not be killed
	until further notice.
	
	For example, a data acquisition XOP that stores data in a data folder during IDLE
	messages would typically store the data folder handle in a global variable and write
	to that data folder each time it receives an IDLE message. It is important that Igor
	not kill a data folder while the XOP is holding it.
	
	Igor Pro 6 and later implement data folder reference counting. Igor keeps a
	reference count for each data folder. This count is used to determine if a 
	data folder is in use.

	As of Igor Pro 6.20, XOPs can participate in data folder reference counting.
	If your XOP obtains a data folder handle and holds it after the XOP returns
	to Igor, you should call HoldDataFolder to increment Igor's internal reference
	count for that data folder. When you no longer need to access the data folder,
	you must call ReleaseDataFolder. ReleaseDataFolder decrements the data folder
	reference count and, if the count reaches zero, kills the data folder.
	
	You also should call HoldDataFolder and ReleaseDataFolder if you are doing a
	callback to Igor and the callback could possibly kill the data folder - for
	example, if you are doing an XOPCommand callback to run a user-defined function
	that could kill the data folder. This constitutes "holding the data folder"
	because you are holding a reference to the data folder over a period of time
	in which it could be killed. Therefore you should use HoldDataFolder to tell
	Igor you are using it and ReleaseDataFolder when you are no longer using it.

	If you are just using the data folder handle temporarily during the execution
	of your external function or operation and you make no calls that could potentially
	kill the data folder then you do not need to and should not call HoldDataFolder
	and ReleaseDataFolder.
	
	If you are running with an Igor version prior to 6.20, HoldDataFolder does nothing
	and returns IGOR_OBSOLETE. If you must run with a version earlier than 6.20, you have
	no way to prevent the killing of a data folder. Therefore you should use a different
	strategy. See GetDataFolderIDNumber and GetDataFolderByIDNumber for further information.

	Once you have called HoldDataFolder on a data folder, Igor will not allow the user
	to kill it until ReleaseDataFolder has been called on it. Normally you will call
	ReleaseDataFolder yourself but here is a case where Igor will call ReleaseDataFolder
	on it. If you define an external function with a structure parameter and the structure
	contains a DFREF field and your XOP returns a data folder handle to the calling
	user-defined function by setting the DFREF field, you must call HoldDataFolder
	passing a pointer to the DFREF field as the dfRefPtr parameter. This signifies
	that the DFREF field is holding the data folder. When the calling user-defined
	function returns and the structure goes out of scope, Igor will automatically call
	ReleaseDataFolder on the structure DFREF field. In this example, by calling HoldDataFolder
	with a pointer to the structure DFREF field, you are passing ownership of the data
	folder to the calling user-defined function which is then responsible for calling
	ReleaseDataFolder.
		
	When an experiment is closed, either by doing New Experiment, opening another
	experiment or quitting Igor, all data folders are killed whether you are holding them
	or not. Consequently, you must release any data folder that you are holding by
	calling ReleaseDataFolder.	

	HoldDataFolder returns 0 or IGOR_OBSOLETE.
	
	HoldDataFolder was added for Igor Pro 6.20. If you call this with an earlier
	version of Igor, it will do nothing and return IGOR_OBSOLETE.

	Thread Safety: HoldDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
HoldDataFolder(DataFolderHandle dfH, DataFolderHandle* dfRefPtr)
{
	return (int)CallBack2(HOLD_DATAFOLDER, dfH, dfRefPtr);
}

/*	ReleaseDataFolder(dfRefPtr)

	Tells Igor that you are no longer holding a data folder.
	
	dfRefPtr contains the address of your DataFolderHandle variable that refers to
	a data folder. ReleaseDataFolder sets *dfRefPtr to NULL so your DataFolderHandle
	variable is not valid after you call ReleaseDataFolder.
	
	See HoldDataFolder for a detailed discussion.

	Returns 0 or IGOR_OBSOLETE.
	
	ReleaseDataFolder was added for Igor Pro 6.20. If you call this with an earlier
	version of Igor, it will do nothing and return IGOR_OBSOLETE.

	Thread Safety: ReleaseDataFolder is thread-safe with Igor Pro 6.20 or later.
*/
int
ReleaseDataFolder(DataFolderHandle* dfRefPtr)
{
	return (int)CallBack1(RELEASE_DATAFOLDER, dfRefPtr);
}

