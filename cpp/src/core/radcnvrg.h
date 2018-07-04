/*-------------------------------------------------------------------------
*
* File name:      radcnvrg.h
*
* Project:        RADIA
*
* Description:    Randomization
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADCNVRG_H
#define __RADCNVRG_H

#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "gmvect.h"

//-------------------------------------------------------------------------

class radTConvergRepair {
public:

	double AbsRand;
	double RelRand;
	double ZeroRand;
	
	short ActOnDoubles;

	radTConvergRepair() 
	{ 
		ActOnDoubles = 1;

		AbsRand = 1.E-09;
		RelRand = 1.E-11;
		ZeroRand = AbsRand;

		//srand((unsigned)time(NULL)); // To seed random generator by current time
	}

	void SwitchActOnDoubles(short InActOnDoubles, double InAbsRand =0., double InRelRand =0., double InZeroRand =0.)
	{
		ActOnDoubles = InActOnDoubles;

		AbsRand = InAbsRand;
		RelRand = InRelRand;
		ZeroRand = InZeroRand;

		srand((unsigned)time(NULL));
	}

	double Double(double A)
	{
		if(!ActOnDoubles) return A;
		A += AbsRand*(double(rand())/double(RAND_MAX) - 0.5);
		A *= (1. + RelRand*(double(rand())/double(RAND_MAX) - 0.5));
		if(A==0.) A = ZeroRand*(double(rand())/double(RAND_MAX) - 0.5);
		return A;
	}
	double DoublePlus(double A)
	{
		if(!ActOnDoubles) return A;
		A += AbsRand*(double(rand())/double(RAND_MAX));
		A *= (1. + RelRand*(double(rand())/double(RAND_MAX)));
		if(A == 0.) A = ZeroRand*(double(rand())/double(RAND_MAX));
		return A;
	}
	double DoubleMinus(double A)
	{
		if(!ActOnDoubles) return A;
		A -= AbsRand*(double(rand())/double(RAND_MAX));
		A *= (1. - RelRand*(double(rand())/double(RAND_MAX)));
		if(A == 0.) A = -ZeroRand*(double(rand())/double(RAND_MAX));
		return A;
	}

	double AbsRandMagnitude(double A)
	{
		if(!ActOnDoubles) return 0.;
		double AbsFromRel = RelRand*A;
		double AbsMax = (AbsFromRel<AbsRand)? AbsRand : AbsFromRel;
		return (A!=0.)? AbsMax : ((AbsMax<ZeroRand)? ZeroRand : AbsMax);
	}
};

//-------------------------------------------------------------------------

#endif

