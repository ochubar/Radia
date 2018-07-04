/*-------------------------------------------------------------------------
*
* File name:      radsbdac.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 subdivided rectangular cross-section arc with azimuthal current
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADSBDAC_H
#define __RADSBDAC_H

#include "radgroup.h"
#include "radarccu.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTSubdividedArcCur : public radTGroup, public radTArcCur {
public:
	int AmOfSubElem;

	radTSubdividedArcCur(const radTArcCur* ArcCurPtr)
	{
		radTArcCur::CentrPoint = radTGroup::CentrPoint = ArcCurPtr->CentrPoint; //OC101008
		CircleCentrPoint = ArcCurPtr->CircleCentrPoint;
		R_min = ArcCurPtr->R_min; R_max = ArcCurPtr->R_max;
		Phi_min = ArcCurPtr->Phi_min; Phi_max = ArcCurPtr->Phi_max;
		Height = ArcCurPtr->Height;
		J_azim = ArcCurPtr->J_azim;
		NumberOfSectors = ArcCurPtr->NumberOfSectors;
		BasedOnPrecLevel = ArcCurPtr->BasedOnPrecLevel;
		InternalFacesAfterCut = ArcCurPtr->InternalFacesAfterCut;
		J_IsNotZero = ArcCurPtr->J_IsNotZero;

		radTGroup::IsGroupMember = ArcCurPtr->IsGroupMember;
		radTGroup::g3dListOfTransform = ArcCurPtr->g3dListOfTransform;
		radTGroup::HandleAuxCompData = ArcCurPtr->HandleAuxCompData;
		radTGroup::ConsiderOnlyWithTrans = ArcCurPtr->ConsiderOnlyWithTrans;

		radTGroup::MessageChar = ArcCurPtr->MessageChar;
	}
	radTSubdividedArcCur() {}

	int Type_Group() { return 4;}

	void B_comp(radTField* FieldPtr) { radTGroup::B_comp(FieldPtr);}
	void B_intComp(radTField* FieldPtr) { radTGroup::B_intComp(FieldPtr);}

	radTg3dGraphPresent* CreateGraphPresent();

	int SubdivideItself(double* SubdivArray, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions) 
	{
		return radTGroup::SubdivideItself(SubdivArray, In_hg, radPtr, pSubdivOptions);
	}
	int SubdivideItselfByOneSetOfParPlanes(TVector3d& InPlanesNormal, TVector3d* InPointsOnCuttingPlanes, int AmOfPieces_mi_1, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions, radTvhg* pVectOfHgChanged)
	{
		return radTGroup::SubdivideItselfByOneSetOfParPlanes(InPlanesNormal, InPointsOnCuttingPlanes, AmOfPieces_mi_1, In_hg, radPtr, pSubdivOptions, pVectOfHgChanged);
	}
	int SubdivideItselfByParPlanes(double* SubdivArray, int AmOfDir, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{
		return radTGroup::SubdivideItselfByParPlanes(SubdivArray, AmOfDir, In_hg, radPtr, pSubdivOptions);
	}
	int CutItself(TVector3d* InCuttingPlane, radThg& In_hg, radTPair_int_hg& LowerNewPair_int_hg, radTPair_int_hg& UpperNewPair_int_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{
		return radTGroup::CutItself(InCuttingPlane, In_hg, LowerNewPair_int_hg, UpperNewPair_int_hg, radPtr, pSubdivOptions);
	}
	int SubdivideItselfByEllipticCylinder(double* SubdivArray, radTCylindricSubdivSpec* pSubdivSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{ 
		return radTGroup::SubdivideItselfByEllipticCylinder(SubdivArray, pSubdivSpec, In_hg, radPtr, pSubdivOptions);
	}
	int CreateFromSym(radThg& In_hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
	{
		return radTGroup::CreateFromSym(In_hg, radPtr, PutNewStuffIntoGenCont);
	}
	int FindLowestAndUppestVertices(TVector3d& PlanesNormal, radTSubdivOptions* pSubdivOptions, TVector3d& LowestVertexPoint, TVector3d& UppestVertexPoint, radTrans& Trans, char& TransWasSet, char& Ignore)
	{
		return radTGroup::FindLowestAndUppestVertices(PlanesNormal, pSubdivOptions, LowestVertexPoint, UppestVertexPoint, Trans, TransWasSet, Ignore);
	}

	void Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl)
	{
		radTGroup::Push_backCenterPointAndField(pFieldKey, pVectPairOfVect3d, pBaseTrans, g3dSrcPtr, pAppl);
	}

	int DuplicateItself(radThg& hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
	{
		radTGroup* PtrToNewGroup = new radTSubdividedArcCur(*this);
		return DuplicateGroupStuff(PtrToNewGroup, hg, radPtr, PutNewStuffIntoGenCont);
	}
	int DuplicateWithoutDuplicatingGroupStuff(radThg& hgGroup)
	{
		radTSend Send;
		radTGroup* pNewGroup = new radTSubdividedArcCur(*this);
		if(pNewGroup==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgLoc(pNewGroup);
		hgGroup = hgLoc;
		return 1;
	}

	int NumberOfDegOfFreedom() { return 0;}
	int ItemIsNotFullyInternalAfterCut() { return radTArcCur::ItemIsNotFullyInternalAfterCut();}
	int SizeOfThis() 
	{
		int GenSize = sizeof(radTArcCur);
		return GenSize + radTGroup::SizeOfThis();
	}

	void Dump(std::ostream&, int ShortSign =0);

	double Volume() { return radTArcCur::Volume();}
	void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique) { radTArcCur::VerticesInLocFrame(OutVect, EnsureUnique);}

	void SimpleEnergyComp(radTField* FieldPtr) { radTGroup::SimpleEnergyComp(FieldPtr);}
	void ActualEnergyForceTorqueCompWithAdd(radTField* FieldPtr) { radTGroup::ActualEnergyForceTorqueCompWithAdd(FieldPtr);}
	void MarkFurtherSubdNeed(char SubdNeedX, char SubdNeedY, char SubdNeedZ) { radTGroup::MarkFurtherSubdNeed(SubdNeedX, SubdNeedY, SubdNeedZ);}
	void MarkFurtherSubdNeed1D(char SubdNeed, char XorYorZ) { radTGroup::MarkFurtherSubdNeed1D(SubdNeed, XorYorZ);}
	int NextStepEnergyForceTorqueComp(double* TotSubdArr, radThg& HandleOfThis, radTField* FieldPtr, char& MoreSubdNeeded) 
	{
		return radTGroup::NextStepEnergyForceTorqueComp(TotSubdArr, HandleOfThis, FieldPtr, MoreSubdNeeded);
	}
	int ProceedNextStepEnergyForceTorqueComp(double* SubdArr, radThg& HandleOfThis, radTField* LocFieldPtr, radTField* FieldPtr, char& OutSubdNeed, char XorYorZ)
	{
		return radTGroup::ProceedNextStepEnergyForceTorqueComp(SubdArr, HandleOfThis, LocFieldPtr, FieldPtr, OutSubdNeed, XorYorZ);
	}
	void SetupFurtherSubdInd(char InSubdInd) { radTGroup::SetupFurtherSubdInd(InSubdInd);}
	void SetMessageChar(char InMessageChar) { radTGroup::SetMessageChar(InMessageChar);}
};

//-------------------------------------------------------------------------

#endif
