/*-------------------------------------------------------------------------
*
* File name:      radyield.h
*
* Project:        RADIA
*
* Description:    Yield functionality for platforms with no preemtive multi-tasking
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADYIELD_H
#define __RADYIELD_H

#include <time.h>

#ifdef __MATHEMATICA__
extern "C" {
#include <mathlink.h>
}
#endif

#ifndef __RADSEND_H
#include "radsend.h"
#endif

#if defined ALPHA__DLL__ || defined ALPHA__LIB__
extern int (*pgRadYieldExternFunc)();
#endif

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTYield {
private:
	clock_t oldtime;
	clock_t delta;
public:
	radTYield() { delta=0;}

	inline void YieldInit(double t);
	inline int Check();
};

//-------------------------------------------------------------------------

inline int radTYield::Check() 
{
// takes 1.6, 2.3 µs or 0.2 s of CPU time on Performa 5200
	if(delta <= 0) return 1;

#ifdef __MATHEMATICA__
	if(MLAbort) 
	{ 
		MLPutFunction(stdlink, "Abort", 0); 
		return 0; // aborting (more  smart action needed?)
	}  
	if(clock() > oldtime) 
	{
		MLCallYieldFunction(MLYieldFunction(stdlink), stdlink, (MLYieldParameters)0);
		oldtime = clock()+delta;
	}
#endif
#if defined ALPHA__DLL__ || defined ALPHA__LIB__

	if(clock() > oldtime) 
	{
		if(pgRadYieldExternFunc != 0) 
		{
			if((*pgRadYieldExternFunc)())
			{
				radTSend::ErrorMessage("Radia::Error998"); 
				throw 0;
				return 0;
			}
		}
		oldtime = clock() + delta;
	}

#endif

	return 1; // normal
}

//-------------------------------------------------------------------------

inline void radTYield::YieldInit(double t)
{
	if(t<=0)
	{
		delta=0;
		return;
	}
	delta = (clock_t)(CLOCKS_PER_SEC*t);
	oldtime = clock()+delta;
	return;
}

//-------------------------------------------------------------------------

#endif
