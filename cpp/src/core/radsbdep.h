/*-------------------------------------------------------------------------
*
* File name:      radsbdep.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 subdivided extruded polygon (prism)
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADSBDEP_H
#define __RADSBDEP_H

#include "radgroup.h"
#include "radexpgn.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTSubdividedExtrPolygon : public radTGroup, public radTExtrPolygon {
public:
	int kx, ky, kz;
	double qx, qy, qz;
	int AmOfSubElem;
	bool AlgsBasedOnKsQsMayNotWork;

	radTSubdividedExtrPolygon(const radTExtrPolygon* ExtrPolygPtr, const double* InSubdivArray)
	{
		FirstPoint = ExtrPolygPtr->FirstPoint; 
		AxOrnt = ExtrPolygPtr->AxOrnt; 
		Thickness = ExtrPolygPtr->Thickness;
		BasePolygonHandle = ExtrPolygPtr->BasePolygonHandle;

		radTExtrPolygon::CentrPoint = radTGroup::CentrPoint = ExtrPolygPtr->CentrPoint; // ??

		Magn = ExtrPolygPtr->Magn;
		MaterHandle = ExtrPolygPtr->MaterHandle;

		radTGroup::IsGroupMember = ExtrPolygPtr->IsGroupMember;
		radTGroup::g3dListOfTransform = ExtrPolygPtr->g3dListOfTransform;
		radTGroup::HandleAuxCompData = ExtrPolygPtr->HandleAuxCompData;
		radTGroup::ConsiderOnlyWithTrans = ExtrPolygPtr->ConsiderOnlyWithTrans;

		radTGroup::MessageChar = ExtrPolygPtr->MessageChar;
		
		const double Small = 1.E-09;

		kx = int(InSubdivArray[0]+Small); ky = int(InSubdivArray[2]+Small); kz = int(InSubdivArray[4]+Small);
		qx = InSubdivArray[1]; qy = InSubdivArray[3]; qz = InSubdivArray[5];

		AlgsBasedOnKsQsMayNotWork = false;
	}
	radTSubdividedExtrPolygon(const radTExtrPolygon* ExtrPolygPtr)
	{
		FirstPoint = ExtrPolygPtr->FirstPoint; 
		AxOrnt = ExtrPolygPtr->AxOrnt; 
		Thickness = ExtrPolygPtr->Thickness;
		BasePolygonHandle = ExtrPolygPtr->BasePolygonHandle;

		radTExtrPolygon::CentrPoint = radTGroup::CentrPoint = ExtrPolygPtr->CentrPoint; // ??

		Magn = ExtrPolygPtr->Magn;
		MaterHandle = ExtrPolygPtr->MaterHandle;

		radTGroup::IsGroupMember = ExtrPolygPtr->IsGroupMember;
		radTGroup::g3dListOfTransform = ExtrPolygPtr->g3dListOfTransform;
		radTGroup::HandleAuxCompData = ExtrPolygPtr->HandleAuxCompData;
		radTGroup::ConsiderOnlyWithTrans = ExtrPolygPtr->ConsiderOnlyWithTrans;

		radTGroup::MessageChar = ExtrPolygPtr->MessageChar;

		AlgsBasedOnKsQsMayNotWork = true;
	}
	radTSubdividedExtrPolygon(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers);
	radTSubdividedExtrPolygon() { AlgsBasedOnKsQsMayNotWork = true;};

	int Type_Group() { return 2;}
	int Type_g3dRelax() { return 4;}

	radTg3dGraphPresent* CreateGraphPresent();

	int SubdivideItself(double* SubdivArray, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions) 
	{
		return radTGroup::SubdivideItself(SubdivArray, In_hg, radPtr, pSubdivOptions);
	}
	int SubdivideItselfByOneSetOfParPlanes(TVector3d& InPlanesNormal, TVector3d* InPointsOnCuttingPlanes, int AmOfPieces_mi_1, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions, radTvhg* pVectOfHgChanged)
	{
		AlgsBasedOnKsQsMayNotWork = true;
		return radTGroup::SubdivideItselfByOneSetOfParPlanes(InPlanesNormal, InPointsOnCuttingPlanes, AmOfPieces_mi_1, In_hg, radPtr, pSubdivOptions, pVectOfHgChanged);
	}
	int SubdivideItselfByParPlanes(double* SubdivArray, int AmOfDir, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{
		AlgsBasedOnKsQsMayNotWork = true;
		return radTGroup::SubdivideItselfByParPlanes(SubdivArray, AmOfDir, In_hg, radPtr, pSubdivOptions);
	}
	int CutItself(TVector3d* InCuttingPlane, radThg& In_hg, radTPair_int_hg& LowerNewPair_int_hg, radTPair_int_hg& UpperNewPair_int_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{
		AlgsBasedOnKsQsMayNotWork = true;
		return radTGroup::CutItself(InCuttingPlane, In_hg, LowerNewPair_int_hg, UpperNewPair_int_hg, radPtr, pSubdivOptions);
	}
	int SubdivideItselfByEllipticCylinder(double* SubdivArray, radTCylindricSubdivSpec* pSubdivSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{ 
		AlgsBasedOnKsQsMayNotWork = true;
		return radTGroup::SubdivideItselfByEllipticCylinder(SubdivArray, pSubdivSpec, In_hg, radPtr, pSubdivOptions);
	}
	int CreateFromSym(radThg& In_hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
	{
		AlgsBasedOnKsQsMayNotWork = true;
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
		radTGroup* PtrToNewGroup = new radTSubdividedExtrPolygon(*this);
		return DuplicateGroupStuff(PtrToNewGroup, hg, radPtr, PutNewStuffIntoGenCont);
	}
	int DuplicateWithoutDuplicatingGroupStuff(radThg& hgGroup)
	{
		radTSend Send;
		radTGroup* pNewGroup = new radTSubdividedExtrPolygon(*this);
		if(pNewGroup==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgLoc(pNewGroup);
		hgGroup = hgLoc;
		return 1;
	}

	int NumberOfDegOfFreedom() { return radTGroup::NumberOfDegOfFreedom();}
	int SizeOfThis() { return (radTGroup::SizeOfThis() + radTExtrPolygon::SizeOfThis());}

	void Dump(std::ostream&, int ShortSign =0);
	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey);

	int SetMaterial(radThg& InMatHandle, radTApplication* ApPtr) 
	{ 
		//return radTg3dRelax::SetMaterial(InMatHandle, ApPtr) && radTGroup::SetMaterial(InMatHandle, ApPtr);
		//OC081008 bug fix: The above sequence of statements caused buggy behaviour in those cases when 
		//the orientation of M vector in group members happened to be different from that of the entire subdivided polyhedron.
		//The reason was in eventual change of the EasyAxisDefined switch after calling SetMaterial(...)

		int res1 = radTGroup::SetMaterial(InMatHandle, ApPtr); //OC081008 first go through all elements of the group
		int res2 = radTg3dRelax::SetMaterial(InMatHandle, ApPtr); //OC081008 then set it for the entire subdivided polyhedron
		return res1 && res2;
	}
	void SetM(TVector3d& M) 
	{ 
		radTg3dRelax::SetM(M);
		radTGroup::SetM(M);
	}

	double Volume() { return radTExtrPolygon::Volume();}
	void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique) { radTExtrPolygon::VerticesInLocFrame(OutVect, EnsureUnique);}

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
