/*-------------------------------------------------------------------------
*
* File name:      radapl1.h
*
* Project:        RADIA
*
* Description:    Wrapping RADIA application function calls
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADAPL1_H
#define __RADAPL1_H

#ifndef __GMVECT_H
#include "gmvect.h"
#endif

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radTIntPtrAndInt {
	int* pInt;
	int AnInt;
	radTIntPtrAndInt(int* In_pInt =0, int InAnInt =0) { pInt=In_pInt; AnInt=InAnInt;}
};

//-------------------------------------------------------------------------

#ifdef __GCC__
typedef vector<radTIntPtrAndInt> radTVectIntPtrAndInt;
#else
typedef vector<radTIntPtrAndInt, allocator<radTIntPtrAndInt> > radTVectIntPtrAndInt;
#endif

//-------------------------------------------------------------------------

struct radTPtrsToPgnAndVect2d {
	radTPolygon* pPgn;
	TVector2d* pVect2d;
	int AmOfPoints;
	radTPtrsToPgnAndVect2d() { pPgn = 0; pVect2d = 0; AmOfPoints = 0;}
};

//-------------------------------------------------------------------------

struct radTVertexPointLiaison {
#ifdef __GCC__
	vector<int> FirstIndVect;
	vector<int> SecondIndVect;
#else
	vector<int, allocator<int> > FirstIndVect;
	vector<int, allocator<int> > SecondIndVect;
#endif
	char AdjSegmentUsed;

	radTVertexPointLiaison() { AdjSegmentUsed = 0;}
};

//-------------------------------------------------------------------------

#endif
