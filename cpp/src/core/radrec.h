/*-------------------------------------------------------------------------
*
* File name:      radrec.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 rectangular parallelepiped with constant magnetization 
*                 or currect density
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADREC_H
#define __RADREC_H

#include "radg3d.h"

#include "radvlpgn.h"
#include "radcnvrg.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

extern radTConvergRepair& radCR;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radTParallelepSurfIntData {
	TVector3d PointOnSurface;
	int SurfBoundInd;
// Surface bound indicator: 1 - lower, 2 - upper, 3 - left, 4 - right, 5 - back, 6 - front, 0 - No
	int IntegrandLen;
	void (radTg3d::*IntegrandFunPtr)(radTField*);
	radTField Field;

	double* InnerAbsPrecAndLimitsArray;
	short* InnerElemCompNotFinished;
	TVector3d** InnerIntegVal;
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTRecMag : public radTg3dRelax {
	radTParallelepSurfIntData* SurfIntDataPtr;
public:
	TVector3d Dimensions;
	TVector3d J;
	short J_IsNotZero;
	short InternalFacesAfterCut;

	radTRecMag(const TVector3d& InCPoiVect, const TVector3d& InDimsVect, 
			   const TVector3d& InMagnVect, const TVector3d& InJ_vect, const radThg& InMaterHandle, short InJ_IsNotZero =0) 
			   : radTg3dRelax(InCPoiVect, InMagnVect, InMaterHandle)
	{
		Dimensions=InDimsVect; J=InJ_vect;
		if(InMaterHandle.rep != 0) J_IsNotZero = 0;
		InternalFacesAfterCut = 0;

		J_IsNotZero = InJ_IsNotZero;
	}
	radTRecMag(const TVector3d& InCPoiVect, const TVector3d& InDimsVect, 
			   const TVector3d& InMagnVect, 
			   const TVector3d& InJ_vect, short InJ_IsNotZero)
			   : radTg3dRelax(InMagnVect)
	{
		CentrPoint=InCPoiVect; Dimensions=InDimsVect; //Magn=InMagnVect; 
		J=InJ_vect; J_IsNotZero = InJ_IsNotZero;
		InternalFacesAfterCut = 0;
	}
	radTRecMag(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{//Instantiates from string according to DumpBin
		DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);
		DumpBinParse_g3dRelax(inStr, mKeysOldNew, gMapOfHandlers);
		DumpBinParse_RecMag(inStr);
	}
	radTRecMag() : radTg3dRelax()
	{ 
		InternalFacesAfterCut = 0;
	}

	int Type_g3dRelax() { return 1;}
	virtual int Type_RecMag() { return 0;}

	void B_comp(radTField*);
	void B_compMultipole(radTField*, double*);
	void B_intComp(radTField*);
	void B_intUtilSpecCaseZeroVxVy(const TVector3d&, const TVector3d&, short, TMatrix3d&, TVector3d&);

	void IntOverShape(radTField* FieldPtr) 
	{
		if(FieldPtr->ShapeIntDataPtr->IntOverSurf_) IntOverSurf(FieldPtr);
		else if(FieldPtr->ShapeIntDataPtr->IntOverVol_) IntOverVol(FieldPtr);
	}
	void IntOverSurf(radTField*);
	void FunForOuterIntAtSurfInt(double, TVector3d*);
	inline void FunForInnerIntAtSurfInt(double, TVector3d*);
	void IntOverVol(radTField*) {}

	double Volume() { return Dimensions.x*Dimensions.y*Dimensions.z;}
	void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique);

	void SimpleEnergyComp(radTField* FieldPtr)
	{
		const double PI = 3.14159265358979;
		const double ConstForM = -1./(4.*PI*100.);
		const double ConstForJ = -1.E-06;
		const double MagnCurDensTol = 1.E-09;
		char LocJ_IsNotZero = Abs(J.x)>MagnCurDensTol || Abs(J.y)>MagnCurDensTol || Abs(J.z)>MagnCurDensTol || J_IsNotZero;
		char LocM_IsNotZero = Abs(Magn.x)>MagnCurDensTol || Abs(Magn.y)>MagnCurDensTol || Abs(Magn.z)>MagnCurDensTol;
		radTFieldKey LocFieldKey;
		if(LocM_IsNotZero) LocFieldKey.B_ = 1;
		if(LocJ_IsNotZero) LocFieldKey.A_ = 1;
		radTField LocField(LocFieldKey, FieldPtr->CompCriterium);
		LocField.P = CentrPoint;
		((radTg3d*)(FieldPtr->HandleEnergyForceTorqueCompData.rep->hSource.rep))->B_genComp(&LocField);
		if(LocM_IsNotZero) FieldPtr->Energy += (ConstForM*Volume())*(Magn*LocField.B);
		if(LocJ_IsNotZero) FieldPtr->Energy += (ConstForJ*Volume())*(J*LocField.A);
	}
	void UniformlyDistrPoints(double* q, TVector3d& P)
	{// This is not used
		P.x = CentrPoint.x + Dimensions.x*((*(q++))-0.5); 
		P.y = CentrPoint.y + Dimensions.y*((*(q++))-0.5); 
		P.z = CentrPoint.z + Dimensions.z*((*q)-0.5);
	}

	void Push_backCenterPointAndField(radTFieldKey*, radTVectPairOfVect3d*, radTrans*, radTg3d*, radTApplication*);
	
	void Dump(std::ostream&, int ShortSign =0);
	void DumpPureObjInfo(std::ostream&, int);
	//void DumpBin(CAuxBinStrVect& oStr, radTmhg& mEl, radThg& hg); 
	//void DumpBin(CAuxBinStrVect& oStr, map<int, radTHandle<radTg>, less<int> >& mEl, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey);
	void DumpBin_RecMag(CAuxBinStrVect& oStr);
	void DumpBinParse_RecMag(CAuxBinStrVect& inStr);
	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey);

	radTg3dGraphPresent* CreateGraphPresent();

	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{
		return FinishDuplication(new radTRecMag(*this), hg);
	}

	int SubdivideItself(double*, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideItselfByOneSetOfParPlanes(TVector3d&, TVector3d*, int, radThg&, radTApplication*, radTSubdivOptions*, radTvhg*);
	int SubdivideItselfByPlanesParToFace(short, TVector3d*, int, radThg&, radTApplication*, char, char);
	int SubdivideItselfByParPlanes(double* SubdivArray, int AmOfDir, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{
		char ThereIsCurrent = (J.x!=0. || J.y!=0. || J.z!=0. || J_IsNotZero);
		if(ThereIsCurrent)
		{
			radTSend Send;
			if(pSubdivOptions->SubdivideCoils) { Send.ErrorMessage("Radia::Error109"); return 0;}
			else return 1;
		}
		radThg& NewHandle = In_hg;
		radThg OldHandle = In_hg;
		if(!((radTg3d*)(OldHandle.rep))->ConvertToPolyhedron(NewHandle, radPtr, pSubdivOptions->PutNewStuffIntoGenCont)) return 0;
		OldHandle = NewHandle;
		if(!((radTg3d*)(OldHandle.rep))->SubdivideItselfByParPlanes(SubdivArray, AmOfDir, NewHandle, radPtr, pSubdivOptions)) return 0;
		return 1;
	}
	int SubdivideItselfByEllipticCylinder(double* SubdivArray, radTCylindricSubdivSpec* pSubdivSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{ 
		char ThereIsCurrent = (J.x!=0. || J.y!=0. || J.z!=0. || J_IsNotZero);
		if(ThereIsCurrent)
		{
			radTSend Send;
			if(pSubdivOptions->SubdivideCoils) { Send.ErrorMessage("Radia::Error109"); return 0;}
			else return 1;
		}
		radThg& NewHandle = In_hg;
		radThg OldHandle = In_hg;
		if(!((radTg3d*)(OldHandle.rep))->ConvertToPolyhedron(NewHandle, radPtr, pSubdivOptions->PutNewStuffIntoGenCont)) return 0;
		OldHandle = NewHandle;
		if(!((radTg3d*)(OldHandle.rep))->SubdivideItselfByEllipticCylinder(SubdivArray, pSubdivSpec, NewHandle, radPtr, pSubdivOptions)) return 0;
		return 1;
	}

	int CutItself(TVector3d*, radThg&, radTPair_int_hg&, radTPair_int_hg&, radTApplication*, radTSubdivOptions*);
	int FindLowestAndUppestVertices(TVector3d&, radTSubdivOptions*, TVector3d&, TVector3d&, radTrans&, char&, char&);

	int SetMaterial(radThg& InMatHandle, radTApplication* ApPtr) 
	{ 
		if(!J_IsNotZero) return radTg3dRelax::SetMaterial(InMatHandle, ApPtr);
		else return 1;
	}

	int ScaleCurrent(double scaleCoef) //virtual in g3d
	{//note: if(scaleCoef == 0) this still doesn't change J_IsNotZero
		if(J_IsNotZero) 
		{
			J *= scaleCoef; return 1;
		}
		else return 0;
	}

	int NumberOfDegOfFreedom() { return (MaterHandle.rep == 0)? 0 : 3;}
	int SizeOfThis() { return sizeof(radTRecMag);}

	int ConvertToPolyhedron(radThg&, radTApplication*, char);
	void CheckVertexPtsPositionsWithRespectToPlane(TVector3d*, char&);

	void DefineRelAndAbsTol(double* RelAbsTol)
	{
		double RelZeroToler = 1.E-09;
		RelZeroToler = 500.*((RelZeroToler>radCR.RelRand)? RelZeroToler : radCR.RelRand);

		TVector3d VectToCenter = 0.5*Dimensions;
		RelAbsTol[1] = RelZeroToler*NormAbs(VectToCenter);
		RelAbsTol[0] = RelZeroToler;
	}

	int ItemIsNotFullyInternalAfterCut()
	{
		return (InternalFacesAfterCut == 63)? 0 : 1;
	}
	void MapFaceAsInternalAfterCut(short FaceNo)
	{
		short ExtraFaceCode;
		switch(FaceNo)
		{
			case 1:
				ExtraFaceCode = 1;
				break;
			case 2:
				ExtraFaceCode = 2;
				break;
			case 3:
				ExtraFaceCode = 4;
				break;
			case 4:
				ExtraFaceCode = 8;
				break;
			case 5:
				ExtraFaceCode = 16;
				break;
			case 6:
				ExtraFaceCode = 32;
				break;
		}
		InternalFacesAfterCut |= ExtraFaceCode;
	}
	void MapFaceAsExternal(short FaceNo)
	{
		short ExtraFaceCode;
		switch(FaceNo)
		{
			case 1:
				ExtraFaceCode = 1;
				break;
			case 2:
				ExtraFaceCode = 2;
				break;
			case 3:
				ExtraFaceCode = 4;
				break;
			case 4:
				ExtraFaceCode = 8;
				break;
			case 5:
				ExtraFaceCode = 16;
				break;
			case 6:
				ExtraFaceCode = 32;
				break;
		}
		InternalFacesAfterCut &= (!ExtraFaceCode);
	}
	void ListFacesInternalAfterCut(short* FacesState)
	{
		short BufNum = InternalFacesAfterCut;
		for(int k=0; k<6; k++) { *(FacesState++) = BufNum & 1; BufNum >>= 1;}
	}
	void SetFacesInternalAfterCut(short* FacesState)
	{
		for(int k=0; k<6; k++) 
		{
			if(*(FacesState++)) MapFaceAsInternalAfterCut(k+1);
		}
	}
};

//-------------------------------------------------------------------------

inline void radTRecMag::FunForInnerIntAtSurfInt(double Arg, TVector3d* VectArray)
{
	if(SurfIntDataPtr->SurfBoundInd==1 || SurfIntDataPtr->SurfBoundInd==2 || 
	   SurfIntDataPtr->SurfBoundInd==3 || SurfIntDataPtr->SurfBoundInd==4)
		SurfIntDataPtr->PointOnSurface.x = Arg;
	else if(SurfIntDataPtr->SurfBoundInd==5 || SurfIntDataPtr->SurfBoundInd==6)
		SurfIntDataPtr->PointOnSurface.y = Arg;

	SurfIntDataPtr->Field.P = SurfIntDataPtr->PointOnSurface;
	(((radTg3d*)this)->*(SurfIntDataPtr->IntegrandFunPtr))(&(SurfIntDataPtr->Field));

	for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++) 
		VectArray[i] = (SurfIntDataPtr->Field.ShapeIntDataPtr->VectArray)[i];
}

//-------------------------------------------------------------------------

#endif