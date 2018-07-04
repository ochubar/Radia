/*-------------------------------------------------------------------------
*
* File name:      radg3da1.h
*
* Project:        RADIA
*
* Description:    Some inline functions of the base class for 3D objects
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADG3DA1_H
#define __RADG3DA1_H

#include "radg3d.h"
//#include "radtrans.h"
#include "gmtrans.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

inline void radTg3d::FindResTransfWithMultOne(radTrans& ResTrans, short& SomethingFound)
{// This multiplies all the transformations with m = 1
	radTrans BufTrans;
	short NoOneFoundYet = 1;
	for(radTlphg::iterator iter = g3dListOfTransform.begin(); iter != g3dListOfTransform.end(); ++iter)
	{
		int Mult = (*iter).m;
		if(Mult==1)
		{
			if(NoOneFoundYet) 
			{ 
				BufTrans = *((radTrans*)((*iter).Handler_g.rep));
				NoOneFoundYet = 0;
			}
			else BufTrans = Product(BufTrans, *((radTrans*)((*iter).Handler_g.rep)));
		}
	}
	if(NoOneFoundYet) { SomethingFound = 0; return;}
	else { ResTrans = BufTrans; SomethingFound = 1; return;}
}

//-------------------------------------------------------------------------

inline void radTg3d::FindInnerTransfWithMultOne(radTrans& ResTrans, short& SomethingFound)
{// This multiplies all the transformations with m = 1
	short LocSomethingFound = 0;
	for(radTlphg::reverse_iterator iter = g3dListOfTransform.rbegin(); iter != g3dListOfTransform.rend(); ++iter)
	{
		int Mult = (*iter).m;
		if(Mult==1)
		{
			ResTrans = *((radTrans*)((*iter).Handler_g.rep));
			LocSomethingFound = 1;
			break;
		}
	}
	SomethingFound = LocSomethingFound;
}

//-------------------------------------------------------------------------

#endif
