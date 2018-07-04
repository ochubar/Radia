/*-------------------------------------------------------------------------
*
* File name:      mathlink_wrap.cpp
*
* Project:        RADIA
*
* Description:    Wraps of some MathLink functions ensuring backwards-compatibility with previous MathLink versions (based of fixes suggested by D. Hidas on 070715)
*                 Use these functions in the code(s) instead of the corresponding MathLink functions.
*
* Author(s):      Dean Hidas and Oleg Chubar
*
* First release:  2015
* 
* Copyright (C):  2015 by European Synchrotron Radiation Facility (France) and Brookhaven National Laboratory (USA)
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __MATHLINK_WRAP_H
#define __MATHLINK_WRAP_H

#ifndef MLINTERFACE
#define MLINTERFACE 3 //4 //Set this to 3 to compile with ML3 for Math9; or to 4 to compile with ML3 for Math10
#endif

extern "C" {
#include <mathlink.h>
}

//-------------------------------------------------------------------------
#if(MLVERSION >= 3)
static void MLWrapDeleteSymbol(MLINK link, const char *s) 
#else
static void MLWrapDeleteSymbol(MLINK link, char *s) 
#endif
{//DH070715
#if MLINTERFACE >= 4
	MLReleaseSymbol(link, s);
#else
	MLDisownSymbol(link, s);
#endif
}

//-------------------------------------------------------------------------
#if(MLVERSION >= 3)
static void MLWrapDeleteString(MLINK link, const char *s)
#else
static void MLWrapDeleteString(MLINK link, char *s)
#endif
{//DH070715
#if MLINTERFACE >= 4
	MLReleaseString(link, s);
#else
	MLDisownString(link, s);
#endif
}

//-------------------------------------------------------------------------
static void MLWrapDeleteDoubleArray(MLINK link, double *data, long *dims, char **heads, long depth)
{//DH070715
#if MLINTERFACE >= 4
	MLReleaseDoubleArray(link, data, dims, heads, depth);
	//MLReleaseRealArray(link, data, dims, heads, depth);
#else
	MLDisownDoubleArray(link, data, dims, heads, depth);
#endif
}

//-------------------------------------------------------------------------
static void MLWrapDeleteRealArray(MLINK link, double *data, long *dims, char **heads, long depth)
{
	MLWrapDeleteDoubleArray(link, data, dims, heads, depth);
}

//-------------------------------------------------------------------------
static void MLWrapDeleteIntegerList(MLINK link, int *data, long n)
{//DH070715
#if MLINTERFACE >= 4
	MLReleaseIntegerList(link, data, n);
#else
	MLDisownIntegerList(link, data, n);
#endif
}

//-------------------------------------------------------------------------

#endif
