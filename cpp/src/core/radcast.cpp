/*-------------------------------------------------------------------------
*
* File name:      radcast.cpp
*
* Project:        RADIA
*
* Description:    "Dynamic cast" for RADIA classes
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radcast.h"

#ifndef __RADINTRC_H
#include "radintrc.h"
#endif

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTInteraction* radTCast::InteractCast(radTg* gPtr)
{
	radTInteraction Interact;
	if(gPtr->Type_g()==Interact.Type_g()) return (radTInteraction*)gPtr;
	else return 0;
	// Do not move it to the .h file: compilation problems may appear!
}
