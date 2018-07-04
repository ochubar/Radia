/*-------------------------------------------------------------------------
*
* File name:      radflm.h
*
* Project:        RADIA
*
* Description:    Magnetic field source: filament conductor
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADFLM_H
#define __RADFLM_H

#include "radg3d.h"
#include "radtrans.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTFlmLinCur : public radTg3d {
	radTrans NativeRotation;
	double Length;

public:
	double I;
	TVector3d StartPoint, EndPoint;

	radTFlmLinCur(const TVector3d&, const TVector3d&, double);

	radTFlmLinCur(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{//Instantiates from string according to DumpBin
		DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);
		NativeRotation.DumpBinParse_Trans(inStr);

		//double Length;
		inStr >> Length;

		//double I;
		inStr >> I;

		//TVector3d StartPoint, EndPoint;
		inStr >> StartPoint;
		inStr >> EndPoint;
	}

	radTFlmLinCur() {}

	int Type_g3d() { return 4;}

	void SetNativeRotation(const TVector3d&, double);

	void B_comp(radTField*);
	void B_intComp(radTField*);

	void SimpleEnergyComp(radTField* FieldPtr)
	{
		const double ConstForJ = -1.E-06;
		radTFieldKey LocFieldKey; LocFieldKey.A_ = 1;
		radTField LocField(LocFieldKey, FieldPtr->CompCriterium);
		LocField.P = 0.5*(StartPoint + EndPoint);
		((radTg3d*)(FieldPtr->HandleEnergyForceTorqueCompData.rep->hSource.rep))->B_genComp(&LocField);
		FieldPtr->Energy += (ConstForJ*I)*((EndPoint - StartPoint)*LocField.A);
	}

	void Dump(std::ostream&, int ShortSign =0); // Porting
	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey);
	radTg3dGraphPresent* CreateGraphPresent();

	int SubdivideItself(double*, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideItselfByOneSetOfParPlanes(TVector3d&, TVector3d*, int, radThg&, radTApplication*, radTSubdivOptions* pSubdivOptions, radTvhg*) 
	{
		radTSend Send;
		if(!pSubdivOptions->SubdivideCoils) return 1;
		else { Send.ErrorMessage("Radia::Error109"); return 0;}
	}
	int SubdivideItselfByParPlanes(double*, int, radThg&, radTApplication*, radTSubdivOptions* pSubdivOptions) 
	{
		radTSend Send;
		if(!pSubdivOptions->SubdivideCoils) return 1;
		else { Send.ErrorMessage("Radia::Error109"); return 0;}
	}
	int CutItself(TVector3d* InCuttingPlane, radThg& In_hg, radTPair_int_hg&, radTPair_int_hg&, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{
		radTSend Send;
		if(!pSubdivOptions->SubdivideCoils) return 1;
		else { Send.ErrorMessage("Radia::Error109"); return 0;}
	}
	int SubdivideItselfByEllipticCylinder(double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions* pSubdivOptions)
	{ 
		radTSend Send;
		if(!pSubdivOptions->SubdivideCoils) return 1;
		else { Send.ErrorMessage("Radia::Error109"); return 0;}
	}

	int FindLowestAndUppestVertices(TVector3d&, radTSubdivOptions* pSubdivOptions, TVector3d&, TVector3d&, radTrans&, char&, char& Ignore)
	{
		Ignore = 1;
		if(!pSubdivOptions->SubdivideCoils) return 1;
		else 
		{ 
			radTSend Send;
			Send.ErrorMessage("Radia::Error109"); return 0;
		}
	}

	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{ 
		return FinishDuplication(new radTFlmLinCur(*this), hg);
	}

	int NumberOfDegOfFreedom() { return 0;}
	int SizeOfThis() { return sizeof(radTFlmLinCur);}

	int ScaleCurrent(double scaleCoef) //virtual in g3d
	{
		I *= scaleCoef; 
		return 1;
	}
};

//-------------------------------------------------------------------------

#endif
