/*-------------------------------------------------------------------------
*
* File name:      radgroup.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 group (/container) of magnetic field sources
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radappl.h"
#include "radgroup.h"
#include "radg3dgr.h"
#include "radg3da1.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTGroup::Dump(std::ostream& o, int ShortSign) // Porting
{
	radTg3d::Dump(o);
	DumpPureObjInfo(o, ShortSign);
	if(ShortSign==1) return;

	DumpTransApplied(o);

	o << endl;
	o << "   Memory occupied (incl. the content): " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

void radTGroup::DumpPureObjInfo(std::ostream& o, int ShortSign)
{
	o << "Container";

	if(ShortSign==1) return;

	o << endl;
	o << "   Content:";

	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
	{
		o << endl;
		o << "      Index " << (*iter).first << ": ";
		(((*iter).second).rep)->Dump(o, 1);
	}
}

//-------------------------------------------------------------------------

void radTGroup::DumpBin_Group_TreatMembers(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, vector<int>& vGroupMemKeys)
{
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		int existKey = 0;
		const radThg &cur_hg = iter->second;

		for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
		{
			if(mit->second == cur_hg)
			{
				existKey = mit->first; break;
			}
		}

		if(existKey == 0)
		{
			existKey = gUniqueMapKey;
			gMapOfHandlers[gUniqueMapKey++] = cur_hg;
		}

		int indExist = CAuxParse::FindElemInd(existKey, vElemKeysOut);
		if(indExist < 0)
		{
			cur_hg.rep->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, existKey); //adding element to the map (mEl) should happen here
		}
		vGroupMemKeys.push_back(existKey);
	}
}

//-------------------------------------------------------------------------

void radTGroup::DumpBin_Group_OutMemKeys(CAuxBinStrVect& oStr, vector<int>& vGroupMemKeys)
{
	int nGroupMem = (int)vGroupMemKeys.size();
	oStr << nGroupMem;
	for(int i=0; i<nGroupMem; i++)
	{
		oStr << vGroupMemKeys[i];
	}
}

//-------------------------------------------------------------------------

//void radTGroup::DumpBin(CAuxBinStrVect& oStr, radTmhg& mEl, radThg& hg)
void radTGroup::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	//Dumping objects that may be used by this object
	vector<pair<int, int> > vTrfKeys;
	DumpBin_g3d_TreatTrfs(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vTrfKeys);

	vector<int> vGroupMemKeys;
	DumpBin_Group_TreatMembers(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vGroupMemKeys);

	//Start dumping this object
	//elemCount++;
	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)Type_g();
	oStr << (char)Type_g3d();
	oStr << (char)Type_Group();
	oStr << (char)0;
	oStr << (char)0;

	//Members of radTg3d
	DumpBin_g3d(oStr, vTrfKeys);

	DumpBin_Group_OutMemKeys(oStr, vGroupMemKeys);
}

//-------------------------------------------------------------------------

void radTGroup::DumpBinParse_Group(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
{
	int nGroupMem = 0;
	inStr >> nGroupMem;
	for(int i=0; i<nGroupMem; i++)
	{
		int oldKey = 0;
		inStr >> oldKey;

		map<int, int>::const_iterator itKey = mKeysOldNew.find(oldKey);
		int newKey = 0;
		if(itKey == mKeysOldNew.end()) throw 0;

		newKey = itKey->second;
		if(newKey > 0)
		{
			radTmhg::const_iterator iter = gMapOfHandlers.find(newKey);
			if(iter == gMapOfHandlers.end()) throw 0;
			radThg hg = (*iter).second;
			
			if(radTCast::g3dCast(hg.rep)==0) throw 0;
			AddElement(newKey, hg);
		}
	}
}

//-------------------------------------------------------------------------

radTGroup::radTGroup(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
{//Instantiates from string according to DumpBin
	DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);
	DumpBinParse_Group(inStr, mKeysOldNew, gMapOfHandlers);
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTGroup::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTGroupGraphPresent(this);
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------

int radTGroup::DuplicateGroupStuff(radTGroup* NewGroupPtr, radThg& hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
{
	radTSend Send;
	if(NewGroupPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
	NewGroupPtr->IsGroupMember = 0;
	NewGroupPtr->GroupMapOfHandlers.erase(NewGroupPtr->GroupMapOfHandlers.begin(), NewGroupPtr->GroupMapOfHandlers.end());
	radThg hgLoc(NewGroupPtr);

	int NewElemKey, NewStuffCounter = 0;
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg hgSubLoc;
		if(!((radTg3d*)(((*iter).second).rep))->DuplicateItself(hgSubLoc, radPtr, PutNewStuffIntoGenCont)) return 0;

		if(PutNewStuffIntoGenCont)
		{
			NewElemKey = radPtr->AddElementToContainer(hgSubLoc);
			radPtr->CopyDrawAttr((*iter).first, NewElemKey);
			NewGroupPtr->AddElement(NewElemKey, hgSubLoc);
		}
		else NewGroupPtr->AddElement(++NewStuffCounter, hgSubLoc);
	}
	hg = hgLoc; return 1;
	// This really creates new copies of all the GroupMembers
}

//-------------------------------------------------------------------------

int radTGroup::SetMaterial(radThg& InMatHandle, radTApplication* ApPtr)
{
	char PutNewStuffIntoGenCont = 1; // For Material: Maybe not necessary?
	radTMaterial* pMat = (radTMaterial*)(InMatHandle.rep);
	char EasyAxisDefinedInMat = pMat->EasyAxisDefined;

	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg hgMat = InMatHandle;
		if(!EasyAxisDefinedInMat) 
		{
			if(!pMat->DuplicateItself(hgMat, ApPtr, PutNewStuffIntoGenCont)) return 0;
			if(PutNewStuffIntoGenCont) ApPtr->AddElementToContainer(hgMat); // Maybe not necessary
		}
		if(!((radTg3d*)(((*iter).second).rep))->SetMaterial(hgMat, ApPtr)) return 0;
	}
	return 1;
}

//-------------------------------------------------------------------------

void radTGroup::SetM(TVector3d& M)
{
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		((radTg3d*)(((*iter).second).rep))->SetM(M);
	}
}

//-------------------------------------------------------------------------

int radTGroup::SubdivideItself(double* SubdivArray, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	if(pSubdivOptions->SubdivisionFrame == 1) // AsWhole differs even with no transformations (distinguish AsWhole and Lab. frame ?)
		return SubdivideItselfAsWholeInLabFrame(SubdivArray, In_hg, radPtr, pSubdivOptions);

	if((pSubdivOptions->SubdivisionFrame == 2) && (!g3dListOfTransform.empty())) // In Lab Frame, but each object separately
	{
		double LocSubdivArray[] = {1.,0.,0., SubdivArray[0], SubdivArray[1], 0.,1.,0., SubdivArray[2], SubdivArray[3], 0.,0.,1., SubdivArray[4], SubdivArray[5]};
		return SubdivideItselfByParPlanes(LocSubdivArray, 3, In_hg, radPtr, pSubdivOptions);
	}

	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg& NewHandle = (*iter).second;
		radThg OldHandle = NewHandle;
		int SubdOK = ((radTg3d*)(OldHandle.rep))->SubdivideItself(SubdivArray, NewHandle, radPtr, pSubdivOptions);
		if(!SubdOK) return 0;
		if(pSubdivOptions->PutNewStuffIntoGenCont)
		{
			radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
			if(((radTg3d*)(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::SubdivideItselfAsWholeInLabFrame(double* SubdivArray, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	TVector3d Directions[3];
	*Directions = TVector3d(1.,0.,0.); Directions[1] = TVector3d(0.,1.,0.); Directions[2] = TVector3d(0.,0.,1.);

	radTSubdivOptions LocSubdivOptions = *pSubdivOptions;
	LocSubdivOptions.SubdivisionFrame = 1; // In Laboratory fr.
	LocSubdivOptions.PutNewStuffIntoGenCont = 1;
	
	TVector3d LowestVertexPoint, UppestVertexPoint;
	radTrans Trans;
	char TransWasSet = 0, Ignore = 0;
	const double RelZeroTol = 5.E-13;

	radTvhg AuxVectOfHgChanged;

	radTSend Send;
	//radTCast Cast;

	SetMessageChar(0);

	for(int k=0; k<3; k++)
	{
		int k2 = k << 1;
		double SubdivisionParam[] = { *(SubdivArray + k2), *(SubdivArray + k2 + 1)};

		FindLowestAndUppestVertices(Directions[k], &LocSubdivOptions, LowestVertexPoint, UppestVertexPoint, Trans, TransWasSet, Ignore);
		TVector3d V = UppestVertexPoint - LowestVertexPoint;
		double Size = fabs(V*Directions[k]);

		double& kk = *SubdivisionParam;
		double& qq = SubdivisionParam[1];
		if(pSubdivOptions->SubdivisionParamCode == 1)
		{
			kk = (kk < Size)? Round(Size/kk) : 1.;
		}
		int AmOfPieces = int(*SubdivisionParam + 1.E-10);

		if(AmOfPieces > 1)
		{
			int AmOfPieces_mi_1 = AmOfPieces - 1;
			TVector3d* PointsOnCuttingPlanes = new TVector3d[AmOfPieces_mi_1];
			if(PointsOnCuttingPlanes == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

			double q0 = (fabs(kk-1.)>RelZeroTol)? pow(qq, 1./(kk-1.)) : qq;
			double Buf = qq*q0 - 1.;
			double dTau = (fabs(Buf) > RelZeroTol)? (q0 - 1.)/Buf : 1./kk;
			double Tau = dTau;
			for(int j=0; j<AmOfPieces_mi_1; j++)
			{
				PointsOnCuttingPlanes[j] = LowestVertexPoint + Tau*V;
				dTau *= q0; Tau += dTau;
			}

			for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
			{
				radThg& NewHandle = (*iter).second;
				radThg OldHandle = NewHandle;
				radTg* gPtrOld = OldHandle.rep;

				int SubdOK = ((radTg3d*)(OldHandle.rep))->SubdivideItselfByOneSetOfParPlanes(Directions[k], PointsOnCuttingPlanes, AmOfPieces_mi_1, NewHandle, radPtr, &LocSubdivOptions, &AuxVectOfHgChanged);
				if(!SubdOK) return 0;
				
				if(gPtrOld != NewHandle.rep) 
				{
					if(pSubdivOptions->PutNewStuffIntoGenCont) 
					{
						radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
						if(((radTg3d*)(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
					}
				}
			}
			delete[] PointsOnCuttingPlanes;
		}
	}

	for(radTvhg::iterator It = AuxVectOfHgChanged.begin(); It != AuxVectOfHgChanged.end(); ++It)
	{
		radTGroup* pGroup = radTCast::GroupCast((radTg3d*)((*It).rep));
		if(pGroup != 0) 
		{
			char RespectKeys = 0;
			pGroup->FlattenNestedStructure(0, RespectKeys); // 0 is essential
			radTmhg AuxMapOfElem;
			for(radTmhg::iterator iter = pGroup->GroupMapOfHandlers.begin(); iter != pGroup->GroupMapOfHandlers.end(); ++iter)
			{
				radThg& hg = (*iter).second;
				AuxMapOfElem[radPtr->AddElementToContainer(hg)] = hg;
				((radTg3d*)(hg.rep))->IsGroupMember = 1;
			}
			pGroup->GroupMapOfHandlers = AuxMapOfElem;
		}
	}
	SetMessageChar(0);
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::SetUpCuttingPlanes(TVector3d& PlanesNormal, double* SubdivisionParam, radTSubdivOptions* pSubdivOptions, TVector3d* PointsOnCuttingPlanes)
{// This is not used
	TVector3d LowestVertexPoint, UppestVertexPoint;
	radTrans Trans;
	char TransWasSet = 0, Ignore = 0;
	FindLowestAndUppestVertices(PlanesNormal, pSubdivOptions, LowestVertexPoint, UppestVertexPoint, Trans, TransWasSet, Ignore);

	TVector3d V = UppestVertexPoint - LowestVertexPoint;
	double Size = sqrt(V.x*V.x + V.y*V.y + V.z*V.z);

	const double RelZeroTol = 5.E-13;

	double& kk = *SubdivisionParam;
	double& qq = SubdivisionParam[1];
	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		kk = (kk < Size)? Round(Size/kk) : 1.;
	}

	double q0 = (fabs(kk-1.)>RelZeroTol)? pow(qq, 1./(kk-1.)) : qq;
	double Buf = qq*q0 - 1.;
	double dTau = (fabs(Buf) > RelZeroTol)? (q0 - 1.)/Buf : 1./kk;

	int AmOfPieces_mi_1 = int(kk + 1.E-10 - 1.);
	double Tau = dTau;
	for(int j=0; j<AmOfPieces_mi_1; j++)
	{
		PointsOnCuttingPlanes[j] = LowestVertexPoint + Tau*V;
		dTau *= q0; Tau += dTau;
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::FindLowestAndUppestVertices(TVector3d& PlanesNormal, radTSubdivOptions* pSubdivOptions, TVector3d& LowestVertexPoint, TVector3d& UppestVertexPoint, radTrans& Trans, char& TransWasSet, char& Ignore)
{// This may change PlanesNormal !
	short SubdInLocFrame = (pSubdivOptions->SubdivisionFrame == 0)? 1 : 0;

	TransWasSet = 0;
	TVector3d& ActualPlanesNormal = PlanesNormal;

	radTrans ResTransf;
	short SomethingFound = 0;
	if(!SubdInLocFrame)
	{
		FindResTransfWithMultOne(ResTransf, SomethingFound);
	}
	else if(ConsiderOnlyWithTrans)
	{
		FindInnerTransfWithMultOne(ResTransf, SomethingFound);
	}
	if(SomethingFound) 
	{
		ActualPlanesNormal = ResTransf.TrBiPoint_inv(ActualPlanesNormal);
		Trans = ResTransf; 
		TransWasSet = 1;
	}

	TVector3d LowestPo = (1.E+23)*PlanesNormal, UppestPo = (-1.E+23)*PlanesNormal;
	TVector3d LocLowestPo, LocUppestPo;

	Ignore = 1;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		TVector3d LocNormal = PlanesNormal; // To suppress changing the PlanesNormal
		radTrans LocTrans;
		char LocTransWasSet = 0, LocIgnore = 0;
		if(!((radTg3d*)(((*iter).second).rep))->FindLowestAndUppestVertices(LocNormal, pSubdivOptions, LocLowestPo, LocUppestPo, LocTrans, LocTransWasSet, LocIgnore)) return 0;
		if(!LocIgnore)
		{
			Ignore = 0;
			if(LocTransWasSet) 
			{
				LocLowestPo = LocTrans.TrPoint(LocLowestPo);
				LocUppestPo = LocTrans.TrPoint(LocUppestPo);
			}
			TVector3d TestLoV = LocLowestPo - LowestPo;
			TVector3d TestUpV = LocUppestPo - UppestPo;
			if(TestLoV*PlanesNormal < 0.) LowestPo = LocLowestPo;
			if(TestUpV*PlanesNormal > 0.) UppestPo = LocUppestPo;
		}
	}
	LowestVertexPoint = LowestPo;
	UppestVertexPoint = UppestPo;
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::SubdivideItselfByOneSetOfParPlanes(
	TVector3d& InPlanesNormal, TVector3d* InPointsOnCuttingPlanes, int AmOfPieces_mi_1, 
	radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions, radTvhg* VectOfHgChanged)
{
	TVector3d PlanesNormal, *PointsOnCuttingPlanes;
	if(!TransferSubdivisionStructToLocalFrame(InPlanesNormal, InPointsOnCuttingPlanes, AmOfPieces_mi_1, pSubdivOptions, PlanesNormal, PointsOnCuttingPlanes)) return 0;

	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg& NewHandle = (*iter).second;
		radThg OldHandle = NewHandle;
		int SubdOK = ((radTg3d*)(OldHandle.rep))->SubdivideItselfByOneSetOfParPlanes(PlanesNormal, PointsOnCuttingPlanes, AmOfPieces_mi_1, NewHandle, radPtr, pSubdivOptions, VectOfHgChanged);
		if(!SubdOK) return 0;
		if(pSubdivOptions->PutNewStuffIntoGenCont) 
		{
			radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
			if(((radTg3d*)(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
		}
	}

	if(PointsOnCuttingPlanes != InPointsOnCuttingPlanes) delete[] PointsOnCuttingPlanes;
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::CutItself(TVector3d* CuttingPlane, radThg& In_hg, radTPair_int_hg& LowerNewPair_int_hg, radTPair_int_hg& UpperNewPair_int_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	short SubdInLocFrame = (pSubdivOptions->SubdivisionFrame == 0)? 1 : 0;
	char AddNewElemsToGenCont = pSubdivOptions->PutNewStuffIntoGenCont;
	char ReplaceOldStuff = pSubdivOptions->ReplaceOldStuff;
	char SeparatePiecesAtCutting = pSubdivOptions->SeparatePiecesAtCutting;

	TVector3d PointOnCutPlane = *CuttingPlane;
	TVector3d CutPlaneNormal = CuttingPlane[1];
	double SqLen = CutPlaneNormal.x*CutPlaneNormal.x + CutPlaneNormal.y*CutPlaneNormal.y + CutPlaneNormal.z*CutPlaneNormal.z;
	double Norm = 1./sqrt(SqLen);
	CutPlaneNormal = Norm*CutPlaneNormal;
	
	radTrans ResTransf;
	short SomethingFound = 0;
	if(!SubdInLocFrame)
	{
		FindResTransfWithMultOne(ResTransf, SomethingFound);
	}
	else if(ConsiderOnlyWithTrans)
	{
		FindInnerTransfWithMultOne(ResTransf, SomethingFound);
	}
	if(SomethingFound) 
	{
		CutPlaneNormal = ResTransf.TrBiPoint_inv(CutPlaneNormal);
		PointOnCutPlane = ResTransf.TrPoint_inv(PointOnCutPlane);
	}

	TVector3d LocCuttingPlane[] = { PointOnCutPlane, CutPlaneNormal };

	radTGroup *pNewLowerGroup = 0, *pNewUpperGroup = 0;
	int LowerElemKeyGenerator = 0, UpperElemKeyGenerator = 0;

	radTSend Send;
	//radTCast Cast;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radTPair_int_hg LocLowerNewPair_int_hg, LocUpperNewPair_int_hg;

		radThg& NewHandle = (*iter).second;
		radThg OldHandle = NewHandle;
		int SubdOK = ((radTg3d*)(OldHandle.rep))->CutItself(LocCuttingPlane, NewHandle, LocLowerNewPair_int_hg, LocUpperNewPair_int_hg, radPtr, pSubdivOptions); // Change this !!!
		if(!SubdOK) return 0;
		if(ReplaceOldStuff) 
		{
			radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
			if(((radTg3d*)(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
		}
		if(SeparatePiecesAtCutting)
		{
			int NewElemKey;
			if(LocLowerNewPair_int_hg.Handler_g.rep != 0)
			{
				if(pNewLowerGroup == 0)
				{
					radTSubdividedPolyhedron* pSubdividedPolyhedron = radTCast::SubdPolyhedronCastFromGroup(this);
					if(pSubdividedPolyhedron != 0) pNewLowerGroup = new radTSubdividedPolyhedron(*pSubdividedPolyhedron);
					else
					{
						radTSubdividedRecMag* pSubdividedRecMag = radTCast::SubdividedRecMagCast(this);
						if(pSubdividedRecMag != 0) pNewLowerGroup = new radTSubdividedRecMag(*pSubdividedRecMag);
						else
						{
							pNewLowerGroup = new radTGroup(*this);
						}
					}
					if(pNewLowerGroup == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
					radTmhg* pNewLowerGroupMapOfHandlers = &(pNewLowerGroup->GroupMapOfHandlers);
					pNewLowerGroupMapOfHandlers->erase(pNewLowerGroupMapOfHandlers->begin(), pNewLowerGroupMapOfHandlers->end());
				}
				if(AddNewElemsToGenCont)
				{
					NewElemKey = LocLowerNewPair_int_hg.m;
					if(NewElemKey == 0) NewElemKey = radPtr->AddElementToContainer(LocLowerNewPair_int_hg.Handler_g);
				}
				else NewElemKey = ++LowerElemKeyGenerator;
				pNewLowerGroup->AddElement(NewElemKey, LocLowerNewPair_int_hg.Handler_g);
			}
			if(LocUpperNewPair_int_hg.Handler_g.rep != 0)
			{
				if(pNewUpperGroup == 0)
				{
					radTSubdividedPolyhedron* pSubdividedPolyhedron = radTCast::SubdPolyhedronCastFromGroup(this);
					if(pSubdividedPolyhedron != 0) pNewUpperGroup = new radTSubdividedPolyhedron(*pSubdividedPolyhedron);
					else
					{
						radTSubdividedRecMag* pSubdividedRecMag = radTCast::SubdividedRecMagCast(this);
						if(pSubdividedRecMag != 0) pNewUpperGroup = new radTSubdividedRecMag(*pSubdividedRecMag);
						else
						{
							pNewUpperGroup = new radTGroup(*this);
						}
					}
					if(pNewUpperGroup == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
					radTmhg* pNewUpperGroupMapOfHandlers = &(pNewUpperGroup->GroupMapOfHandlers);
					pNewUpperGroupMapOfHandlers->erase(pNewUpperGroupMapOfHandlers->begin(), pNewUpperGroupMapOfHandlers->end());
				}
				if(AddNewElemsToGenCont)
				{
					NewElemKey = LocUpperNewPair_int_hg.m;
					if(NewElemKey == 0) NewElemKey = radPtr->AddElementToContainer(LocUpperNewPair_int_hg.Handler_g);
				}
				else NewElemKey = ++UpperElemKeyGenerator;
				pNewUpperGroup->AddElement(NewElemKey, LocUpperNewPair_int_hg.Handler_g);
			}
		}
	}
	if(SeparatePiecesAtCutting)
	{
		if(pNewLowerGroup != 0)
		{
			radThg hgNewLowerGroup(pNewLowerGroup);
			LowerNewPair_int_hg.m = AddNewElemsToGenCont? radPtr->AddElementToContainer(hgNewLowerGroup) : 1;
			LowerNewPair_int_hg.Handler_g = hgNewLowerGroup;
		}
		if(pNewUpperGroup != 0)
		{
			radThg hgNewUpperGroup(pNewUpperGroup);
			UpperNewPair_int_hg.m = AddNewElemsToGenCont? radPtr->AddElementToContainer(hgNewUpperGroup) : 1;
			UpperNewPair_int_hg.Handler_g = hgNewUpperGroup;
		}

		if(AddNewElemsToGenCont)
		{
			int GroupElemKey = radPtr->RetrieveElemKey(this);
			if(pNewLowerGroup != 0) radPtr->CopyDrawAttr(GroupElemKey, LowerNewPair_int_hg.m);
			if(pNewUpperGroup != 0) radPtr->CopyDrawAttr(GroupElemKey, UpperNewPair_int_hg.m);
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::SubdivideItselfByParPlanes(double* SubdivArray, int AmOfDir, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	if(pSubdivOptions->SubdivisionFrame == 1) // AsWhole differs even with no transformations (distinguish AsWhole and Lab. frame ?)
		return SubdivideItselfByParPlanesAsWholeInLabFrame(SubdivArray, AmOfDir, In_hg, radPtr, pSubdivOptions);

	short SubdInLocFrame = (pSubdivOptions->SubdivisionFrame == 0)? 1 : 0;
	char SubdivideCoils = pSubdivOptions->SubdivideCoils;
	char PutNewStuffIntoGenCont = pSubdivOptions->PutNewStuffIntoGenCont;

	double LocSubdivArray[15];
	for(int i=0; i<AmOfDir; i++)
	{
		int i_mu_5 = i*5;
		int i_mu_5_p_1 = i_mu_5+1, i_mu_5_p_2 = i_mu_5+2, i_mu_5_p_3 = i_mu_5+3, i_mu_5_p_4 = i_mu_5+4;
		TVector3d PlanesNormal(SubdivArray[i_mu_5], SubdivArray[i_mu_5_p_1], SubdivArray[i_mu_5_p_2]);
	
		radTrans ResTransf;
		short SomethingFound = 0;
		if(!SubdInLocFrame)
		{
			FindResTransfWithMultOne(ResTransf, SomethingFound);
		}
		else if(ConsiderOnlyWithTrans)
		{
			FindInnerTransfWithMultOne(ResTransf, SomethingFound);
		}
		if(SomethingFound) PlanesNormal = ResTransf.TrBiPoint_inv(PlanesNormal);

		LocSubdivArray[i_mu_5] = PlanesNormal.x;
		LocSubdivArray[i_mu_5_p_1] = PlanesNormal.y;
		LocSubdivArray[i_mu_5_p_2] = PlanesNormal.z;
		LocSubdivArray[i_mu_5_p_3] = SubdivArray[i_mu_5_p_3];
		LocSubdivArray[i_mu_5_p_4] = SubdivArray[i_mu_5_p_4];
	}
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg& NewHandle = (*iter).second;
		radThg OldHandle = NewHandle;
		int SubdOK = ((radTg3d*)(OldHandle.rep))->SubdivideItselfByParPlanes(LocSubdivArray, AmOfDir, NewHandle, radPtr, pSubdivOptions);
		if(!SubdOK) return 0;

		if(PutNewStuffIntoGenCont) 
		{
			radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
			if(((radTg3d*)(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::SubdivideItselfByParPlanesAsWholeInLabFrame(double* SubdivArray, int AmOfDir, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	radTSend Send;
	TVector3d* Directions = new TVector3d[AmOfDir];
	if(Directions == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

	TVector3d* tDirections = Directions;
	double* tSubdivArray = SubdivArray;
	for(int i=0; i<AmOfDir; i++)
	{
		tDirections->x = *(tSubdivArray++);
		tDirections->y = *(tSubdivArray++);
		tDirections->z = *(tSubdivArray++);
		tDirections++;
		tSubdivArray += 2;
	}

	radTSubdivOptions LocSubdivOptions = *pSubdivOptions;
	LocSubdivOptions.SubdivisionFrame = 1; // In Laboratory fr.
	LocSubdivOptions.PutNewStuffIntoGenCont = 1;
	
	TVector3d LowestVertexPoint, UppestVertexPoint;
	radTrans Trans;
	char TransWasSet = 0, Ignore = 0;
	const double RelZeroTol = 5.E-13;

	SetMessageChar(0);
	radTvhg AuxVectOfHgChanged;

	double* StartKs = SubdivArray + 3;

	for(int k=0; k<AmOfDir; k++)
	{
		int k5 = k*5;
		double SubdivisionParam[] = {*(StartKs + k5), *(StartKs + k5 + 1)};

		FindLowestAndUppestVertices(Directions[k], &LocSubdivOptions, LowestVertexPoint, UppestVertexPoint, Trans, TransWasSet, Ignore);
		TVector3d V = UppestVertexPoint - LowestVertexPoint;
		double Size = fabs(V*Directions[k]);

		double& kk = *SubdivisionParam;
		double& qq = SubdivisionParam[1];
		if(pSubdivOptions->SubdivisionParamCode == 1)
		{
			kk = (kk < Size)? Round(Size/kk) : 1.;
		}
		int AmOfPieces = int(*SubdivisionParam + 1.E-10);
		if(AmOfPieces > 1)
		{
			int AmOfPieces_mi_1 = AmOfPieces - 1;
			TVector3d* PointsOnCuttingPlanes = new TVector3d[AmOfPieces_mi_1];
			if(PointsOnCuttingPlanes == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

			double q0 = (fabs(kk-1.)>RelZeroTol)? pow(qq, 1./(kk-1.)) : qq;
			double Buf = qq*q0 - 1.;
			double dTau = (fabs(Buf) > RelZeroTol)? (q0 - 1.)/Buf : 1./kk;
			double Tau = dTau;
			for(int j=0; j<AmOfPieces_mi_1; j++)
			{
				PointsOnCuttingPlanes[j] = LowestVertexPoint + Tau*V;
				dTau *= q0; Tau += dTau;
			}

			for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
			{
				radThg& NewHandle = (*iter).second;
				radThg OldHandle = NewHandle;
				radTg* gPtrOld = OldHandle.rep;

				int SubdOK = ((radTg3d*)(OldHandle.rep))->SubdivideItselfByOneSetOfParPlanes(Directions[k], PointsOnCuttingPlanes, AmOfPieces_mi_1, NewHandle, radPtr, &LocSubdivOptions, &AuxVectOfHgChanged);
				if(!SubdOK) return 0;
				
				if(gPtrOld != NewHandle.rep) 
				{
					if(pSubdivOptions->PutNewStuffIntoGenCont) 
					{
						radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
						if(((radTg3d*)(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
					}
				}
			}
			delete[] PointsOnCuttingPlanes;
		}
	}
	//radTCast Cast;
	for(radTvhg::iterator It = AuxVectOfHgChanged.begin(); It != AuxVectOfHgChanged.end(); ++It)
	{
		radTGroup* pGroup = radTCast::GroupCast((radTg3d*)((*It).rep));
		if(pGroup != 0) 
		{
			char RespectKeys = 0;
			pGroup->FlattenNestedStructure(0, RespectKeys); // 0 is essential
			radTmhg AuxMapOfElem;
			for(radTmhg::iterator iter = pGroup->GroupMapOfHandlers.begin(); iter != pGroup->GroupMapOfHandlers.end(); ++iter)
			{
				radThg& hg = (*iter).second;
				AuxMapOfElem[radPtr->AddElementToContainer(hg)] = hg;
				((radTg3d*)(hg.rep))->IsGroupMember = 1;
			}
			pGroup->GroupMapOfHandlers = AuxMapOfElem;
		}
	}

	SetMessageChar(0);
	if(Directions != 0) delete[] Directions;
	return 1;
}

//-------------------------------------------------------------------------

void radTGroup::FlattenNestedStructure(radTApplication* radPtr, char RespectKeys)
{// If radPtr != 0, operations on the Global Container occur !
	radTmhg MapOfNonGroupElem;

	radTlphg Aux_g3dListOfTransform;
	char GroupHasTransforms = (!g3dListOfTransform.empty());
	if(GroupHasTransforms) 
	{
		Aux_g3dListOfTransform = g3dListOfTransform;
		g3dListOfTransform.erase(g3dListOfTransform.begin(), g3dListOfTransform.end());
	}

	int MemberCount = 0;
	CollectNonGroupElements(&MapOfNonGroupElem, MemberCount, radPtr, RespectKeys);
	GroupMapOfHandlers = MapOfNonGroupElem;

	if(GroupHasTransforms) g3dListOfTransform = Aux_g3dListOfTransform;
}

//-------------------------------------------------------------------------

void radTGroup::CollectNonGroupElements(radTmhg* pMapOfNonGroupElem, int& MemberCount, radTApplication* radPtr, char RespectKeys)
{
	char GroupHasTransforms = (!g3dListOfTransform.empty());

	//radTCast Cast;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg hg = (*iter).second;
		int ElemKey = (*iter).first;
		radTg3d* g3dPtr = (radTg3d*)(hg.rep);

		if(GroupHasTransforms)
		{
			for(radTlphg::reverse_iterator riter = g3dListOfTransform.rbegin(); riter != g3dListOfTransform.rend(); ++riter)
			{
				g3dPtr->g3dListOfTransform.push_front(*riter);
			}
		}

		radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
		if(pGroup == 0)
		{
			if(RespectKeys) (*pMapOfNonGroupElem)[ElemKey] = hg;
			else (*pMapOfNonGroupElem)[MemberCount] = hg;
			MemberCount++;
		}
		else
		{
			pGroup->CollectNonGroupElements(pMapOfNonGroupElem, MemberCount, radPtr, RespectKeys);

			if(radPtr != 0)
			{
				short OldSendingIsRequired = radPtr->SendingIsRequired;
				radPtr->SendingIsRequired = 0;
				radPtr->DeleteElement(ElemKey);
				radPtr->SendingIsRequired = OldSendingIsRequired;
			}
		}
	}
}

//-------------------------------------------------------------------------

int radTGroup::SubdivideItselfByEllipticCylinder(double* SubdivArray, radTCylindricSubdivSpec* pCylSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	radTCylindricSubdivSpec LocCylSpec = *pCylSpec;

	radTrans ResTransf;
	short SomethingFound = 0;
	if((pSubdivOptions->SubdivisionFrame != 0) && (!g3dListOfTransform.empty()))
	{
		FindResTransfWithMultOne(ResTransf, SomethingFound);
	}
	else if(ConsiderOnlyWithTrans && (!g3dListOfTransform.empty()))
	{
		FindInnerTransfWithMultOne(ResTransf, SomethingFound);
	}
	if(SomethingFound)
	{
		LocCylSpec.PointOnCylAx = ResTransf.TrPoint_inv(pCylSpec->PointOnCylAx);
		LocCylSpec.CylAxVect = ResTransf.TrAxialVect_inv(pCylSpec->CylAxVect);
		LocCylSpec.PointOnEllAx = ResTransf.TrPoint_inv(pCylSpec->PointOnEllAx);
	}

	if(pSubdivOptions->SubdivisionFrame == 1) // AsWhole differs even with no transformations (distinguish AsWhole and Lab. frame ?)
		return SubdivideItselfByEllipticCylinderAsWholeInLabFrame(SubdivArray, &LocCylSpec, In_hg, radPtr, pSubdivOptions);

	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg& NewHandle = (*iter).second;
		radThg OldHandle = NewHandle;
		int SubdOK = ((radTg3d*)(OldHandle.rep))->SubdivideItselfByEllipticCylinder(SubdivArray, &LocCylSpec, NewHandle, radPtr, pSubdivOptions);
		if(!SubdOK) return 0;
		if(pSubdivOptions->PutNewStuffIntoGenCont)
		{
			radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
			if(((radTg3d*)(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::SubdivideItselfByEllipticCylinderAsWholeInLabFrame(double* SubdivArray, radTCylindricSubdivSpec* pSubdivSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{
	short OldRecognizeRecMagsInPolyhedrons = radPtr->RecognizeRecMagsInPolyhedrons;
	radPtr->RecognizeRecMagsInPolyhedrons = 0;
	
	char ReplaceOldStuff = pSubdivOptions->ReplaceOldStuff;
	if(!ConvertToPolyhedron(In_hg, radPtr, ReplaceOldStuff)) return 0;

	double kPhi = SubdivArray[2], kz = SubdivArray[4];
	double qPhi = SubdivArray[3], qz = SubdivArray[5];

	const double PI = 3.1415926535898;

	TVector3d EdgePointsOverPhiAndAxForCylSubd[6];
	double Limits[6];

	SetMessageChar(0);

	int ThereArePolyhedrons = FindEdgePointsOverPhiAndAxForCylSubd(pSubdivSpec, EdgePointsOverPhiAndAxForCylSubd, Limits);
	if(!ThereArePolyhedrons) return 1;
	char AxisCrossesVolume = MessageChar;

	radTSend Send;
	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		double SizePhi;
		if(!AxisCrossesVolume)
		{
			double aCenPt, PhiCenPt;
			TVector3d CenterPoint;
			int ThereAreRelaxables = EstimateCenterPointOverRelaxables(CenterPoint);
			if(!ThereAreRelaxables) return 1;

			FindEllipticCoordOfPoint(pSubdivSpec, CenterPoint, aCenPt, PhiCenPt);
			SizePhi = EstimateLengthAlongEllipse(aCenPt, pSubdivSpec->EllAxRatio, Limits[0], Limits[1]);
		}
		else
		{
			TVector3d MaxOffsetVertPt = EdgePointsOverPhiAndAxForCylSubd[5];
			double APt, PhiPt;
			FindEllipticCoordOfPoint(pSubdivSpec, MaxOffsetVertPt, APt, PhiPt);
			double a = 0.5*APt;
			SizePhi = PI*a*(1.5*(pSubdivSpec->EllAxRatio + 1.) - sqrt(pSubdivSpec->EllAxRatio));
		}
		double SizeAx = Limits[3] - Limits[2];

		kPhi = (kPhi < SizePhi)? Round(SizePhi/kPhi) : 1.;
		kz = (kz < SizeAx)? Round(SizeAx/kz) : 1.;
	}
	if(AxisCrossesVolume && (int(kPhi)==1)) { Send.ErrorMessage("Radia::Error068"); return 0;}

	double kPhi_qPhi[2];
	kPhi_qPhi[0] = kPhi; kPhi_qPhi[1] = qPhi;

	radTSubdivOptions LocSubdivOptions = *pSubdivOptions;
	LocSubdivOptions.SubdivisionFrame = 0;
	LocSubdivOptions.SubdivisionParamCode = 0;
	LocSubdivOptions.SubdivideCoils = 0;
	LocSubdivOptions.PutNewStuffIntoGenCont = 0;
	LocSubdivOptions.ReplaceOldStuff = pSubdivOptions->ReplaceOldStuff;
	LocSubdivOptions.SeparatePiecesAtCutting = 1;
	LocSubdivOptions.MapInternalFacesAfterCut = 1;

	if(!SubdivideItselfOverAzimuth(kPhi_qPhi, Limits, pSubdivSpec, In_hg, radPtr, &LocSubdivOptions)) return 0;

	TVector3d EdgePointsOverEllipseSet[2];
	double EllipticCoordOfEdgePoints[2];

	int EdgePointsOverEllipseSetFound = 0;
	if(LocSubdivOptions.MethForRadialSegmAtEllCylSubd == 0)
		EdgePointsOverEllipseSetFound = FindEdgePointsOverEllipseSet0(SubdivArray, pSubdivSpec, In_hg, EdgePointsOverEllipseSet, EllipticCoordOfEdgePoints, &LocSubdivOptions);
	else if(LocSubdivOptions.MethForRadialSegmAtEllCylSubd == 1)
		EdgePointsOverEllipseSetFound = FindEdgePointsOverEllipseSet(SubdivArray, pSubdivSpec, In_hg, EdgePointsOverEllipseSet, EllipticCoordOfEdgePoints, &LocSubdivOptions);

	if(!EdgePointsOverEllipseSetFound) { Send.ErrorMessage("Radia::Error999"); return 0;}

	double ka = SubdivArray[0];
	double qa = SubdivArray[1];

	if((int(ka)==1) && (int(kPhi)==1) && (int(kz)==1)) return 1;

	double ka_qa[2];
	ka_qa[0] = ka; ka_qa[1] = qa;
	if(!SubdivideByEllipses(ka_qa, EllipticCoordOfEdgePoints, pSubdivSpec, In_hg, radPtr, &LocSubdivOptions)) return 0;

	double kz_qz[] = {kz, qz};
	double LimitsAx[] = {Limits[2], Limits[3]};
	TVector3d EdgePointsAx[] = { EdgePointsOverPhiAndAxForCylSubd[2], EdgePointsOverPhiAndAxForCylSubd[3]};

	LocSubdivOptions.PutNewStuffIntoGenCont = pSubdivOptions->PutNewStuffIntoGenCont;
	if(!SubdivideByPlanesPerpToCylAx(kz_qz, LimitsAx, EdgePointsAx, pSubdivSpec, In_hg, radPtr, &LocSubdivOptions)) return 0;

	SetMessageChar(0);
	radPtr->RecognizeRecMagsInPolyhedrons = OldRecognizeRecMagsInPolyhedrons;

	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::FindEdgePointsOverPhiAndAxForCylSubd(radTCylindricSubdivSpec* pSubdivSpec, TVector3d* EdgePointsForCylSubd, double* PhiAndAxLimits)
{// This is only used with Frame->LabTot
	double MinRe2=1.E+23, MaxRe2=0., MinZ=1.E+23, MaxZ=-1.E+23;
	const double TwoPI = 2.*3.1415926535898;

	radTVectPairOfDouble CleanPhiSectors, CleanPhiSectorsDupl;

	//radTCast Cast;
	int PolyhdrFound = 0;
	char AxisCrossesVolume = 0;
	MessageChar = 0;

	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		TVector3d LocEdgePoints[6];
		double LocPhiAndAxLimits[6];
		char LocAxisCrossesVolume = 0;

		char SomethingWasExecuted = 0;
		radThg hg = (*iter).second;
		radTg3d* g3dPtr = (radTg3d*)(hg.rep);

		radTCylindricSubdivSpec LocCylSpec = *pSubdivSpec;
		char BackTransformIsNeeded = 0;
		radTrans ResTransf;
		if(!g3dPtr->g3dListOfTransform.empty())
		{
			short SomethingFound = 0;
			g3dPtr->FindResTransfWithMultOne(ResTransf, SomethingFound);
			if(SomethingFound)
			{
				LocCylSpec.PointOnCylAx = ResTransf.TrPoint_inv(pSubdivSpec->PointOnCylAx);
				LocCylSpec.CylAxVect = ResTransf.TrAxialVect_inv(pSubdivSpec->CylAxVect);
				LocCylSpec.PointOnEllAx = ResTransf.TrPoint_inv(pSubdivSpec->PointOnEllAx);
				BackTransformIsNeeded = 1;
			}
		}

		radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
		if(pGroup != 0)
		{
			pGroup->FindEdgePointsOverPhiAndAxForCylSubd(&LocCylSpec, LocEdgePoints, LocPhiAndAxLimits);
			LocAxisCrossesVolume = pGroup->MessageChar;
			SomethingWasExecuted = 1;
		}
		else
		{
			radTg3dRelax* g3dRelaxPtr = radTCast::g3dRelaxCast(g3dPtr);
			if(g3dRelaxPtr != 0)
			{
				radTPolyhedron* pPolyhdr = radTCast::PolyhedronCast(g3dRelaxPtr);
				if(pPolyhdr != 0)
				{
					pPolyhdr->FindEdgePointsOverPhiAndAxForCylSubd(&LocCylSpec, LocEdgePoints, LocPhiAndAxLimits);
					LocAxisCrossesVolume = pPolyhdr->MessageChar;
					SomethingWasExecuted = 1;
				}
			}
		}
		if(BackTransformIsNeeded)
		{
			LocEdgePoints[2] = ResTransf.TrPoint(LocEdgePoints[2]);
			LocEdgePoints[3] = ResTransf.TrPoint(LocEdgePoints[3]);
			LocEdgePoints[4] = ResTransf.TrPoint(LocEdgePoints[4]);
			LocEdgePoints[5] = ResTransf.TrPoint(LocEdgePoints[5]);
		}

		if(SomethingWasExecuted)
		{
			if(MinZ > LocPhiAndAxLimits[2])
			{
				MinZ = LocPhiAndAxLimits[2]; EdgePointsForCylSubd[2] = LocEdgePoints[2];
			}
			if(MaxZ < LocPhiAndAxLimits[3])
			{
				MaxZ = LocPhiAndAxLimits[3]; EdgePointsForCylSubd[3] = LocEdgePoints[3];
			}
			if(MinRe2 > LocPhiAndAxLimits[4])
			{
				MinRe2 = LocPhiAndAxLimits[4]; EdgePointsForCylSubd[4] = LocEdgePoints[4];
			}
			if(MaxRe2 < LocPhiAndAxLimits[5])
			{
				MaxRe2 = LocPhiAndAxLimits[5]; EdgePointsForCylSubd[5] = LocEdgePoints[5];
			}

			if(LocAxisCrossesVolume && (!AxisCrossesVolume))
			{
				AxisCrossesVolume = 1;
				MessageChar = 1; // To map AxisCrossesVolume; Do not use SetMessageChar here: it will override this value in group members

				PhiAndAxLimits[0] = LocPhiAndAxLimits[0];
				PhiAndAxLimits[1] = LocPhiAndAxLimits[1];
			}

			if(!AxisCrossesVolume)
			{
				if(LocPhiAndAxLimits[0] >= TwoPI) LocPhiAndAxLimits[0] -= TwoPI;
				if(LocPhiAndAxLimits[1] >= TwoPI) LocPhiAndAxLimits[1] -= TwoPI;

				double CurPhiMin = LocPhiAndAxLimits[0];
				double CurPhiMax = LocPhiAndAxLimits[1];

				if(!PolyhdrFound)
				{
					radTPairOfDouble NewCleanSector(CurPhiMax, CurPhiMin);
					CleanPhiSectors.push_back(NewCleanSector);

					PolyhdrFound = 1;
				}
				else
				{
					CleanPhiSectorsDupl.erase(CleanPhiSectorsDupl.begin(), CleanPhiSectorsDupl.end());
					for(int j=0; j<(int)(CleanPhiSectors.size()); j++)
					{
						radTPairOfDouble& CurSect = CleanPhiSectors[j];
						double& CurSectSt = CurSect.First;
						double& CurSectFi = CurSect.Second;

						char PhiMinIsBetween = AngleIsBetween(CurPhiMin, CurSectSt, CurSectFi);
						char PhiMaxIsBetween = AngleIsBetween(CurPhiMax, CurSectSt, CurSectFi);

						if(PhiMinIsBetween && PhiMaxIsBetween)
						{
							char OldSectorIsInsideNew = (AngleIsBetween(CurSectSt, CurPhiMax, CurPhiMin) && AngleIsBetween(CurSectFi, CurPhiMax, CurPhiMin));
							if(OldSectorIsInsideNew)
							{
								radTPairOfDouble NewCleanSector1(CurSectSt, CurPhiMin), NewCleanSector2(CurPhiMax, CurSectFi);
								CleanPhiSectorsDupl.push_back(NewCleanSector1);
								CleanPhiSectorsDupl.push_back(NewCleanSector2);
							}
							else
							{
								radTPairOfDouble NewCleanSector(CurPhiMax, CurPhiMin);
								CleanPhiSectorsDupl.push_back(NewCleanSector);
							}
						}
						else if(PhiMinIsBetween && (!PhiMaxIsBetween))
						{
							radTPairOfDouble NewCleanSector(CurSectSt, CurPhiMin);
							CleanPhiSectorsDupl.push_back(NewCleanSector);
						}
						else if((!PhiMinIsBetween) && PhiMaxIsBetween)
						{
							radTPairOfDouble NewCleanSector(CurPhiMax, CurSectFi);
							CleanPhiSectorsDupl.push_back(NewCleanSector);
						}
						else if((!PhiMinIsBetween) && (!PhiMaxIsBetween))
						{
							CleanPhiSectorsDupl.push_back(CurSect);
						}
					}
					CleanPhiSectors.erase(CleanPhiSectors.begin(), CleanPhiSectors.end());
					CleanPhiSectors = CleanPhiSectorsDupl;
				}
			}
			else
			{
				PolyhdrFound = 1;
			}
		}
	}
	if(!PolyhdrFound) return 0;

	PhiAndAxLimits[2] = MinZ; PhiAndAxLimits[3] = MaxZ;
	PhiAndAxLimits[4] = MinRe2; PhiAndAxLimits[5] = MinRe2;
	if(AxisCrossesVolume) return 1;

	double MaxAngDiff = -1.E+23;
	int MaxAngSectNo = -1;
	for(int i=0; i<(int)(CleanPhiSectors.size()); i++)
	{
		radTPairOfDouble& CurSector = CleanPhiSectors[i];
		double AngDiff = AngularDifference(CurSector.Second, CurSector.First);
		if(MaxAngDiff < AngDiff)
		{
			MaxAngDiff = AngDiff; MaxAngSectNo = i;
		}
	}
	radTPairOfDouble& MaxSector = CleanPhiSectors[MaxAngSectNo];
	PhiAndAxLimits[0] = MaxSector.Second; PhiAndAxLimits[1] = MaxSector.First;

	CleanPhiSectors.erase(CleanPhiSectors.begin(), CleanPhiSectors.end());
	CleanPhiSectorsDupl.erase(CleanPhiSectorsDupl.begin(), CleanPhiSectorsDupl.end());

	if(PhiAndAxLimits[1] < PhiAndAxLimits[0]) PhiAndAxLimits[1] += TwoPI;

	return 1;
// Attention: This leaves EdgePointsForCylSubd [0] and [1] (over Phi) undefined. They seem not to be used anyway.
}

//-------------------------------------------------------------------------

int radTGroup::ConvertToPolyhedron(radThg& In_hg, radTApplication* radPtr, char ReplaceOldStuff)
{
	//radTCast Cast;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg& NewHandle = (*iter).second;
		radTGroup* pGroup = radTCast::GroupCast((radTg3d*)(NewHandle.rep));
		if(pGroup != 0)
		{
			radThg OldHandle = NewHandle;
			if(!pGroup->ConvertToPolyhedron(NewHandle, radPtr, ReplaceOldStuff)) return 0;
		}
		else
		{
			radTg3dRelax* g3dRelaxPtr = radTCast::g3dRelaxCast((radTg3d*)(NewHandle.rep));
			if(g3dRelaxPtr != 0)
			{
				radThg OldHandle = NewHandle;
				if(!g3dRelaxPtr->ConvertToPolyhedron(NewHandle, radPtr, ReplaceOldStuff)) return 0;

				if(ReplaceOldStuff)
				{
					radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
					if(((radTg3d*)(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
				}
			}
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::EstimateCenterPointOverRelaxables(TVector3d& CenterPoint)
{// This is only used with Frame->LabTot
	TVector3d Sum(0.,0.,0.);
	long RelaxCount = 0;
	//radTCast Cast;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg& NewHandle = (*iter).second;
		radTg3d* g3dPtr = (radTg3d*)(NewHandle.rep);

		char SomethingWasExecuted = 0;
		TVector3d LocCP;

		radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
		if(pGroup != 0)
		{
			if(EstimateCenterPointOverRelaxables(LocCP) != 0)
			{
				SomethingWasExecuted = 1;
				RelaxCount++;
			}
		}
		else
		{
			radTg3dRelax* g3dRelaxPtr = radTCast::g3dRelaxCast(g3dPtr);
			if(g3dRelaxPtr != 0)
			{
				char ThisIsNotRecCur = 1;
				radTRecMag* pRecMag = radTCast::RecMagCast(g3dRelaxPtr);
				if(pRecMag != 0) ThisIsNotRecCur = !(pRecMag->J_IsNotZero);

				if(ThisIsNotRecCur)
				{
					SomethingWasExecuted = 1;
					RelaxCount++;
					LocCP = g3dRelaxPtr->CentrPoint;
				}
			}
		}
		if(SomethingWasExecuted) 
		{
			radTrans ResTransf;
			if(!g3dPtr->g3dListOfTransform.empty())
			{
				short SomethingFound = 0;
				g3dPtr->FindResTransfWithMultOne(ResTransf, SomethingFound);
				if(SomethingFound)
				{
					LocCP = ResTransf.TrPoint(LocCP);
				}
			}
			Sum += LocCP;
		}
	}
	if(RelaxCount == 0) return 0;
	CenterPoint = (1./double(RelaxCount))*Sum;
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::SubdivideItselfOverAzimuth(double* kPhi_qPhi, double* Limits, radTCylindricSubdivSpec* pCylSubdivSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{// This is only used with Frame->LabTot
	char AxisCrossesVolume = MessageChar;

	radTSend Send;
	int kPh = int(*kPhi_qPhi);
	if(AxisCrossesVolume && (kPh==1)) { Send.ErrorMessage("Radia::Error068"); return 0;}

	const double PI = 3.1415926535898;

	if((Limits[1] < Limits[0]) && (!AxisCrossesVolume)) 
	{
		Limits[1] += 2.*PI;
		if(Limits[1] < Limits[0]) Limits[0] -= 2.*PI;
	}
	if(AxisCrossesVolume) 
	{
		Limits[0] = 0.;
		Limits[1] = 2.*PI;
	}
	if(kPh==1)
	{
		return 1;
	}

	//radTCast Cast;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg& NewHandle = (*iter).second;
		radTg3d* g3dPtr = (radTg3d*)(NewHandle.rep);

		radTCylindricSubdivSpec LocCylSpec = *pCylSubdivSpec;
		if(!g3dPtr->g3dListOfTransform.empty())
		{
			radTrans ResTransf;
			short SomethingFound = 0;
			g3dPtr->FindResTransfWithMultOne(ResTransf, SomethingFound);
			if(SomethingFound)
			{
				LocCylSpec.PointOnCylAx = ResTransf.TrPoint_inv(pCylSubdivSpec->PointOnCylAx);
				LocCylSpec.CylAxVect = ResTransf.TrAxialVect_inv(pCylSubdivSpec->CylAxVect);
				LocCylSpec.PointOnEllAx = ResTransf.TrPoint_inv(pCylSubdivSpec->PointOnEllAx);
			}
		}

		radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
		if(pGroup != 0)
		{
			radThg OldHandle = NewHandle;
			int SubdOK = pGroup->SubdivideItselfOverAzimuth(kPhi_qPhi, Limits, &LocCylSpec, NewHandle, radPtr, pSubdivOptions);
			if(!SubdOK) return 0;
		}
		else
		{
			radTg3dRelax* g3dRelaxPtr = radTCast::g3dRelaxCast(g3dPtr);
			if(g3dRelaxPtr != 0)
			{
				radTPolyhedron* pPolyhdr = radTCast::PolyhedronCast(g3dRelaxPtr);
				if(pPolyhdr != 0)
				{
					radThg OldHandle = NewHandle;
					int SubdOK = pPolyhdr->SubdivideItselfOverAzimuth(kPhi_qPhi, Limits, &LocCylSpec, NewHandle, radPtr, pSubdivOptions);
					if(!SubdOK) return 0;

					if(pSubdivOptions->ReplaceOldStuff && (OldHandle.rep != NewHandle.rep))
					{
						radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
						if(((radTg3d*)(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
					}
				}
			}
		}
	}
	MessageChar = 0; // Essential for further processing: Releasing this member variable; But do not use SetMessageChar here !!!
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::FindEdgePointsOverEllipseSet0(double* SubdivArray, radTCylindricSubdivSpec* pSubdivSpec, radThg In_hg, TVector3d* EdgePoints, double* Limits, radTSubdivOptions* pSubdivOptions)
{
	double& ka = SubdivArray[0];
	double& qa = SubdivArray[1];
	if((pSubdivOptions->SubdivisionParamCode == 0) && (int(ka) == 1)) return 1;

	TVector3d& CylAxVect = pSubdivSpec->CylAxVect;
	TVector3d& PointOnCylAx = pSubdivSpec->PointOnCylAx;
	TVector3d& PointOnEllAx = pSubdivSpec->PointOnEllAx;

	EdgePoints[0] = PointOnCylAx;
	EdgePoints[1] = PointOnEllAx;

	TVector3d PointOnEllAx_mi_PointOnCylAx = PointOnEllAx - PointOnCylAx;
	TVector3d Ex = PointOnEllAx_mi_PointOnCylAx - (PointOnEllAx_mi_PointOnCylAx*CylAxVect)*CylAxVect;

	double DelA = sqrt(Ex.x*Ex.x + Ex.y*Ex.y + Ex.z*Ex.z);

	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		ka = (ka < DelA)? Round(DelA/ka) : 1.;
	}

	Limits[0] = 0.;
	Limits[1] = DelA;

	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::FindEdgePointsOverEllipseSet(double* SubdivArray, radTCylindricSubdivSpec* pSubdivSpec, radThg In_hg, TVector3d* EdgePoints, double* Limits, radTSubdivOptions* pSubdivOptions)
{// This is only used with Frame->LabTot
	double& ka = SubdivArray[0];
	double& qa = SubdivArray[1];
	if((pSubdivOptions->SubdivisionParamCode == 0) && (int(ka) == 1)) return 1;
	
	double aMin = 1.E+23, aMax = 0.;
	//radTCast Cast;
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg hgLoc = (*iter).second;
		radTg3d* g3dPtr = (radTg3d*)(hgLoc.rep);

		radTCylindricSubdivSpec LocCylSpec = *pSubdivSpec;
		if(!g3dPtr->g3dListOfTransform.empty())
		{
			radTrans ResTransf;
			short SomethingFound = 0;
			g3dPtr->FindResTransfWithMultOne(ResTransf, SomethingFound);
			if(SomethingFound)
			{
				LocCylSpec.PointOnCylAx = ResTransf.TrPoint_inv(pSubdivSpec->PointOnCylAx);
				LocCylSpec.CylAxVect = ResTransf.TrAxialVect_inv(pSubdivSpec->CylAxVect);
				LocCylSpec.PointOnEllAx = ResTransf.TrPoint_inv(pSubdivSpec->PointOnEllAx);
			}
		}

		TVector3d LocEdgePoints[2];
		double LocLimits[] = { 1.E+23, 0.};

		radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
		if(pGroup != 0)
		{
			if(!pGroup->FindEdgePointsOverEllipseSet(SubdivArray, &LocCylSpec, hgLoc, LocEdgePoints, LocLimits, pSubdivOptions)) return 0;
		}
		else
		{
			radTg3dRelax* g3dRelaxPtr = radTCast::g3dRelaxCast(g3dPtr);
			if(g3dRelaxPtr != 0)
			{
				radTPolyhedron* pPolyhdr = radTCast::PolyhedronCast(g3dRelaxPtr);
				if(pPolyhdr != 0)
				{
					if(!pPolyhdr->FindEdgePointsOverEllipseSet(SubdivArray, &LocCylSpec, hgLoc, LocEdgePoints, LocLimits, pSubdivOptions)) return 0;
				}
			}
		}
		if(*LocLimits < aMin)
		{
			aMin = *LocLimits; *EdgePoints = *LocEdgePoints;
		}
		if(LocLimits[1] > aMax)
		{
			aMax = LocLimits[1]; EdgePoints[1] = LocEdgePoints[1];
		}
	}
	*Limits = aMin; Limits[1] = aMax;

	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		double SizeA = Limits[1] - Limits[0];
		ka = (ka < SizeA)? Round(SizeA/ka) : 1.;
	}

	return (Limits[1] >= *Limits)? 1 : 0;
}

//-------------------------------------------------------------------------

int radTGroup::SubdivideByEllipses(double* ka_qa, double* aLimits, radTCylindricSubdivSpec* pCylSubdivSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{// This is only used with Frame->LabTot
	if(int(*ka_qa)==1) return 1;

	//radTCast Cast;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg& NewHandle = (*iter).second;
		radTg3d* g3dPtr = (radTg3d*)(NewHandle.rep);

		radTCylindricSubdivSpec LocCylSpec = *pCylSubdivSpec;
		if(!g3dPtr->g3dListOfTransform.empty())
		{
			radTrans ResTransf;
			short SomethingFound = 0;
			g3dPtr->FindResTransfWithMultOne(ResTransf, SomethingFound);
			if(SomethingFound)
			{
				LocCylSpec.PointOnCylAx = ResTransf.TrPoint_inv(pCylSubdivSpec->PointOnCylAx);
				LocCylSpec.CylAxVect = ResTransf.TrAxialVect_inv(pCylSubdivSpec->CylAxVect);
				LocCylSpec.PointOnEllAx = ResTransf.TrPoint_inv(pCylSubdivSpec->PointOnEllAx);
			}
		}
		
		radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
		if(pGroup != 0)
		{
			radThg OldHandle = NewHandle;
			if(!pGroup->SubdivideByEllipses(ka_qa, aLimits, &LocCylSpec, NewHandle, radPtr, pSubdivOptions)) return 0;

			if(pGroup->MessageChar == 1)
			{
				char RespectKeys = 0;
				radTApplication* Loc_radPtr = 0; // No operations on Global Container !
				pGroup->FlattenNestedStructure(Loc_radPtr, RespectKeys);
			}
		}
		else
		{
			radTg3dRelax* g3dRelaxPtr = radTCast::g3dRelaxCast(g3dPtr);
			if(g3dRelaxPtr != 0)
			{
				radTPolyhedron* pPolyhdr = radTCast::PolyhedronCast(g3dRelaxPtr);
				if(pPolyhdr != 0)
				{
					radThg OldHandle = NewHandle;
					if(!pPolyhdr->SubdivideByEllipses(ka_qa, aLimits, &LocCylSpec, NewHandle, radPtr, pSubdivOptions)) return 0;

					if(pSubdivOptions->ReplaceOldStuff && (OldHandle.rep != NewHandle.rep))
					{
						radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
						if(((radTg3d*)(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
					}
				}
			}
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::SubdivideByPlanesPerpToCylAx(double* kz_qz, double* Limits, TVector3d* EdgePoints, radTCylindricSubdivSpec* pSubdivSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
{// This is only used with Frame->LabTot !!!
	double kz = *kz_qz, qz = kz_qz[1];
	if(int(kz) == 1) 
	{
		char RespectKeys = 0;
		FlattenNestedStructureIfMessageCharIsSet(radPtr, RespectKeys, pSubdivOptions);
		return 1; 
	}

	radTSend Send;
	int AmOfPieces_mi_1 = int(kz) - 1;
	TVector3d* PointsOnCuttingPlanes = new TVector3d[AmOfPieces_mi_1];
	if(PointsOnCuttingPlanes==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

	const double AbsZeroTol = 5.E-12;
	double DelZ = Limits[1] - Limits[0];
	double q0z = (fabs(kz-1.) > AbsZeroTol)? pow(qz, 1./(kz-1.)) : qz;
	double BufZ = qz*q0z - 1.;
	double NewDelZ = (fabs(BufZ) > AbsZeroTol)? DelZ*(q0z - 1.)/BufZ : DelZ/kz;

	TVector3d Pz = EdgePoints[0];
	TVector3d* tPointsOnCuttingPlanes = PointsOnCuttingPlanes;

	for(int k=0; k<AmOfPieces_mi_1; k++)
	{
		Pz += NewDelZ*pSubdivSpec->CylAxVect;
		NewDelZ *= q0z;
		*(tPointsOnCuttingPlanes++) = Pz;
	}

	radTSubdivOptions SubdivOptionsForCutByParPlanes = *pSubdivOptions;
	SubdivOptionsForCutByParPlanes.SubdivisionFrame = 1;
	SubdivOptionsForCutByParPlanes.PutNewStuffIntoGenCont = 0;

	radTSubdivOptions SubdivOptionsForFlatten = *pSubdivOptions;

	radTvhg DummyVectOfHgChanged;
	//radTCast Cast;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg& NewHandle = (*iter).second;
		radThg OldHandle = NewHandle;
		radTg3d* g3dPtr = (radTg3d*)(OldHandle.rep);

		int SubdOK = g3dPtr->SubdivideItselfByOneSetOfParPlanes(pSubdivSpec->CylAxVect, PointsOnCuttingPlanes, AmOfPieces_mi_1, NewHandle, radPtr, &SubdivOptionsForCutByParPlanes, &DummyVectOfHgChanged);
		if(!SubdOK) return 0;

		radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
		if(pGroup != 0)
		{
			char RespectKeys = 0;
			pGroup->FlattenNestedStructureIfMessageCharIsSet(radPtr, RespectKeys, &SubdivOptionsForFlatten);
		}

		if(pSubdivOptions->ReplaceOldStuff && (OldHandle.rep != NewHandle.rep)) 
		{
			radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
			if(((radTg3d*)(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
		}
	}

	if(PointsOnCuttingPlanes!=0) delete[] PointsOnCuttingPlanes;
	return 1;
}
		
//-------------------------------------------------------------------------

void radTGroup::FlattenNestedStructureIfMessageCharIsSet(radTApplication* radPtr, char RespectKeys, radTSubdivOptions* pSubdivOptions)
{
	if(MessageChar == 1)
	{
		char LocRespectKeys = 0;
		radTApplication* Loc_radPtr = 0; // No operations on Global Container here!
		FlattenNestedStructure(Loc_radPtr, LocRespectKeys);

		if(pSubdivOptions->PutNewStuffIntoGenCont)
		{
			radTmhg AuxGroupMapOfHandlers;
			for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
			{
				radThg hg = (*iter).second;

				AuxGroupMapOfHandlers[radPtr->AddElementToContainer(hg)] = hg;
				((radTg3d*)(hg.rep))->IsGroupMember = 1;
			}
			GroupMapOfHandlers.erase(GroupMapOfHandlers.begin(), GroupMapOfHandlers.end());
			GroupMapOfHandlers = AuxGroupMapOfHandlers;
		}
	}
	else
	{
		//radTCast Cast;
		for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
		{
			radThg OldHandle = (*iter).second;
			radTg3d* g3dPtr = (radTg3d*)(OldHandle.rep);

			radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
			if(pGroup != 0)
			{
				pGroup->FlattenNestedStructureIfMessageCharIsSet(radPtr, RespectKeys, pSubdivOptions);
			}
		}
	}
}

//-------------------------------------------------------------------------

void radTGroup::JustTraverse()
{
	//radTCast Cast;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg hg = (*iter).second;
		radTg3d* g3dPtr = (radTg3d*)(hg.rep);

		radTGroup* pGroup = radTCast::GroupCast(g3dPtr);
		if(pGroup != 0)
		{
			pGroup->JustTraverse();
		}
	}
}

//-------------------------------------------------------------------------

int radTGroup::CreateFromSym(radThg& In_hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
{
	radThg hgDpl;
	if(!DuplicateWithoutDuplicatingGroupStuff(hgDpl)) return 0;
	radTGroup* pGroupDpl = (radTGroup*)((radTg3d*)(hgDpl.rep));

	pGroupDpl->GroupMapOfHandlers.erase(pGroupDpl->GroupMapOfHandlers.begin(), pGroupDpl->GroupMapOfHandlers.end());

	int NewStuffCounter = 0;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg hgNew = (*iter).second;
		radThg hgOld = hgNew;
		((radTg3d*)(hgOld.rep))->CreateFromSym(hgNew, radPtr, PutNewStuffIntoGenCont);

		if(PutNewStuffIntoGenCont)
		{
			int NewElemKey = radPtr->AddElementToContainer(hgNew);
			radPtr->CopyDrawAttr((*iter).first, NewElemKey);
			pGroupDpl->AddElement(NewElemKey, hgNew);
		}
		else pGroupDpl->AddElement(++NewStuffCounter, hgNew);
	}

	int MaxMult = 0;
	for(radTlphg::const_iterator iterTr = g3dListOfTransform.begin(); iterTr != g3dListOfTransform.end(); ++iterTr)
		if((*iterTr).m > MaxMult) MaxMult = (*iterTr).m;

	In_hg = hgDpl;
	if(MaxMult > 1) 
	{
		radTSend Send;
		radTGroup* pGroup = new radTGroup();
		if(pGroup==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgGroup(pGroup);

		radTrans BaseTrans;
		BaseTrans.SetupIdent();
		if(!pGroupDpl->NestedFor_CreateFromSym(pGroup, radPtr, PutNewStuffIntoGenCont, &BaseTrans, pGroupDpl->g3dListOfTransform.begin())) return 0;
		if(!pGroup->GroupMapOfHandlers.empty()) 
		{
			pGroup->IsGroupMember = pGroupDpl->IsGroupMember;
			pGroup->ConsiderOnlyWithTrans = pGroupDpl->ConsiderOnlyWithTrans;
			pGroup->MessageChar = pGroupDpl->MessageChar;

			In_hg = hgGroup;
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

void radTGroup::Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl)
{// Attention: this assumes no more than one transformation with mult. no more than 1 !!!
	radTrans* pTrans = (g3dListOfTransform.empty())? 0 : (radTrans*)((*(g3dListOfTransform.begin())).Handler_g.rep);
	radTrans TotTrans;
	if(pTrans != 0)
	{
		if(pBaseTrans != 0) 
		{
			TrProduct(pBaseTrans, pTrans, TotTrans);
			pTrans = &TotTrans;
		}
	}
	else
	{
		if(pBaseTrans != 0) pTrans = pBaseTrans;
	}
	
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg hg = (*iter).second;
		radTg3d* g3dPtr = (radTg3d*)(hg.rep);
		//g3dPtr->Push_backCenterPointAndField(pFieldKey, pVectPairOfVect3d, pTrans, g3dSrcPtr);

		radThg hgDplWithoutSym; //OC061007_BNL
		char PutNewStuffIntoGenCont = 0;
		if(!g3dPtr->CreateFromSym(hgDplWithoutSym, pAppl, PutNewStuffIntoGenCont)) return;

		radTg3d* g3dDplWithoutSymPtr = (radTg3d*)(hgDplWithoutSym.rep);

		radTvhg vhFlatTransforms; //OC061007_BNL
		g3dDplWithoutSymPtr->FlattenSpaceTransforms(vhFlatTransforms);
		if(vhFlatTransforms.size() > 0)
		{
			g3dDplWithoutSymPtr->EraseAllTransformations();
			g3dDplWithoutSymPtr->AddTransform(1, vhFlatTransforms[0]);
		}

		g3dDplWithoutSymPtr->Push_backCenterPointAndField(pFieldKey, pVectPairOfVect3d, pTrans, g3dSrcPtr, pAppl);
	}
}

//-------------------------------------------------------------------------

int radTGroup::NextStepEnergyForceTorqueComp(double* TotSubdArr, radThg& HandleOfThis, radTField* FieldPtr, char& OutMoreSubdNeeded)
{
	char SubdNeedX, SubdNeedY, SubdNeedZ, MoreSubdNeeded;
	if(HandleAuxCompData.rep != 0)
	{
		HandleAuxCompData.rep->ShowSubdNeed(SubdNeedX, SubdNeedY, SubdNeedZ);
		MoreSubdNeeded = SubdNeedX || SubdNeedY || SubdNeedZ;
		if(!MoreSubdNeeded) return 1;
	}

	char ThisListOfTransformNotEmpty = !g3dListOfTransform.empty(); // Modif. 01.04.99

	MoreSubdNeeded = 0;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		char LocMoreSubdNeeded = 0;
		radThg &hg = (*iter).second;
		radThg hgOld = hg;

		radTg3d* OldPtr = (radTg3d*)(hgOld.rep);
		radTlphg Oldg3dListOfTransform = OldPtr->g3dListOfTransform;
		if(ThisListOfTransformNotEmpty) // Modif. 01.04.99
		{ // Modif. 01.04.99
			radTlphg CopyOfThisListOfTransform = g3dListOfTransform; // Modif. 01.04.99
			(OldPtr->g3dListOfTransform).splice((OldPtr->g3dListOfTransform).begin(), CopyOfThisListOfTransform); // Modif. 01.04.99
		} // Modif. 01.04.99

		if(!OldPtr->NextStepEnergyForceTorqueComp(TotSubdArr, hg, FieldPtr, LocMoreSubdNeeded)) return 0;

		((radTg3d*)(hg.rep))->g3dListOfTransform = Oldg3dListOfTransform;
		MoreSubdNeeded |= LocMoreSubdNeeded;
	}
	OutMoreSubdNeeded = MoreSubdNeeded;
	return 1;
}

//-------------------------------------------------------------------------

int radTGroup::ProceedNextStepEnergyForceTorqueComp(double* SubdArr, radThg& HandleOfThis, radTField* LocFieldPtr, radTField* FieldPtr, char& OutSubdNeed, char XorYorZ)
{
	char SubdNeed;
	if(HandleAuxCompData.rep != 0)
	{
		HandleAuxCompData.rep->ShowSubdNeed1D(SubdNeed, XorYorZ);
		if(!SubdNeed) return 1;
	}

	char ThisListOfTransformNotEmpty = !g3dListOfTransform.empty();
	SubdNeed = 0;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		char LocSubdNeed = 0;
		radThg &hg = (*iter).second;
		radThg hgOld = hg;

		radTg3d* OldPtr = (radTg3d*)(hgOld.rep);
		radTlphg Oldg3dListOfTransform = OldPtr->g3dListOfTransform;

		if(ThisListOfTransformNotEmpty)
		{
			radTlphg CopyOfThisListOfTransform = g3dListOfTransform;
			(OldPtr->g3dListOfTransform).splice((OldPtr->g3dListOfTransform).begin(), CopyOfThisListOfTransform);
		}

		if(!OldPtr->ProceedNextStepEnergyForceTorqueComp(SubdArr, hg, LocFieldPtr, FieldPtr, LocSubdNeed, XorYorZ)) return 0;

		((radTg3d*)(hg.rep))->g3dListOfTransform = Oldg3dListOfTransform;
		SubdNeed |= LocSubdNeed;
	}

	radTg3d::MarkFurtherSubdNeed1D(SubdNeed, XorYorZ);
	OutSubdNeed = SubdNeed; 
	return 1;
}

//-------------------------------------------------------------------------

void radTGroup::ActualEnergyForceTorqueCompWithAdd(radTField* FieldPtr)
{
	char ThisListOfTransformNotEmpty = !g3dListOfTransform.empty();

	for(radTmhg::iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg &hg = (*iter).second;
		radTg3d* Ptr = (radTg3d*)(hg.rep);
		radTlphg Oldg3dListOfTransform = Ptr->g3dListOfTransform;

		if(ThisListOfTransformNotEmpty)
		{
			radTlphg CopyOfThisListOfTransform = g3dListOfTransform;
			(Ptr->g3dListOfTransform).splice((Ptr->g3dListOfTransform).begin(), CopyOfThisListOfTransform);
		}

		Ptr->ActualEnergyForceTorqueCompWithAdd(FieldPtr);
		Ptr->g3dListOfTransform = Oldg3dListOfTransform;
	}
}

//-------------------------------------------------------------------------

radTGroup* radTGroup::CreateGroupIncludingAllMembersExceptIt(const radTmhg::const_iterator& it)
{
	radTGroup* OutGroupPtr = new radTGroup();

	int LocInd = -1;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		LocInd++;
		if(iter == it) continue;

		radThg cur_hg = (*iter).second;
        OutGroupPtr->AddElement(LocInd, cur_hg); 
	}
	return OutGroupPtr;
}

//-------------------------------------------------------------------------
