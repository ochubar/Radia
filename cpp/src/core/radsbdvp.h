/*-------------------------------------------------------------------------
*
* File name:      radsbdvp.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 subdivided polyhedron with constant magnetization
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADSBDVP_H
#define __RADSBDVP_H

#include "radgroup.h"
#include "radvlpgn.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTSubdividedPolyhedron : public radTGroup, public radTPolyhedron {
public:
	int AmOfSubElem;

	radTSubdividedPolyhedron(const radTPolyhedron* PolyhedronPtr)
	{
		AmOfFaces = PolyhedronPtr->AmOfFaces;
		int VectSize = (int)(PolyhedronPtr->VectHandlePgnAndTrans.size());
		if(AmOfFaces != VectSize) 
		{
			radTSend::ErrorMessage("Radia::Error117"); 
			throw 0;
		}

		for(int i=0; i<AmOfFaces; i++) 
		{
			VectHandlePgnAndTrans.push_back(PolyhedronPtr->VectHandlePgnAndTrans[i]);
		}

		Magn = PolyhedronPtr->Magn;
		MaterHandle = PolyhedronPtr->MaterHandle;

		radTPolyhedron::CentrPoint = radTGroup::CentrPoint = PolyhedronPtr->CentrPoint; //OC101008

		radTGroup::IsGroupMember = PolyhedronPtr->IsGroupMember;
		radTGroup::g3dListOfTransform = PolyhedronPtr->g3dListOfTransform;
		radTGroup::HandleAuxCompData = PolyhedronPtr->HandleAuxCompData;
		radTGroup::ConsiderOnlyWithTrans = PolyhedronPtr->ConsiderOnlyWithTrans;

		radTGroup::MessageChar = PolyhedronPtr->MessageChar;
	}
	radTSubdividedPolyhedron(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{
		//Members of radTg3d 
		radTGroup::DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers); //all G3D members should be accessed via radTGroup

		//Members of radTGroup 
		DumpBinParse_Group(inStr, mKeysOldNew, gMapOfHandlers);

		//Members of radTg3dRelax 
		DumpBinParse_g3dRelax(inStr, mKeysOldNew, gMapOfHandlers);

		//Members of radTPolyhedron 
		DumpBinParse_Polyhedron(inStr);

		//Members of radTSubdividedPolyhedron
		//int AmOfSubElem;
		inStr >> AmOfSubElem;
	}

	radTSubdividedPolyhedron() {};

	int Type_Group() { return 3;}
	int Type_g3dRelax() { return 6;}

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
		radTGroup* PtrToNewGroup = new radTSubdividedPolyhedron(*this);
		return DuplicateGroupStuff(PtrToNewGroup, hg, radPtr, PutNewStuffIntoGenCont);
	}
	int DuplicateWithoutDuplicatingGroupStuff(radThg& hgGroup)
	{
		radTSend Send;
		radTGroup* pNewGroup = new radTSubdividedPolyhedron(*this);
		if(pNewGroup==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgLoc(pNewGroup);
		hgGroup = hgLoc;
		return 1;
	}

	int NumberOfDegOfFreedom() { return radTGroup::NumberOfDegOfFreedom();}
	int ItemIsNotFullyInternalAfterCut() { return radTPolyhedron::ItemIsNotFullyInternalAfterCut();}
	int SizeOfThis() 
	{
		int GenSize = sizeof(radTPolyhedron);
		int BufSize = sizeof(radTrans);
		GenSize += AmOfFaces*BufSize;
		for(int i=0; i<AmOfFaces; i++) GenSize += (VectHandlePgnAndTrans[i].PgnHndl.rep)->SizeOfThis();
		return GenSize + radTGroup::SizeOfThis();
	}

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

	void Dump(std::ostream&, int ShortSign =0);
	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey);

	double Volume() { return radTPolyhedron::Volume();}
	void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique) { radTPolyhedron::VerticesInLocFrame(OutVect, EnsureUnique);}

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
