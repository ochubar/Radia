/*-------------------------------------------------------------------------
*
* File name:      radgroup.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 group (/container) of magnetic field sources
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADGROUP_H
#define __RADGROUP_H

#include "radg3d.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

#ifdef __GCC__
typedef vector<int> radTVectOfInt;
#else
typedef vector<int, allocator<int> > radTVectOfInt;
#endif

//-------------------------------------------------------------------------

class radTGroup : public radTg3d {
public:
	radTmhg GroupMapOfHandlers;

	radTGroup(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers);
	radTGroup() {}
	~radTGroup() {}

	int Type_g3d() { return 2;}
	virtual int Type_Group() { return 0;}

	inline void AddElement(int, const radThg&);
	
	inline void B_comp(radTField*);  // Modified by P. Elleaume 8 Nov 96
	void B_intComp(radTField* FieldPtr) { B_comp(FieldPtr);} // This is not an Error!!!

	void Dump(std::ostream&, int =0);
	void DumpPureObjInfo(std::ostream&, int);
	void DumpBin_Group_TreatMembers(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, vector<int>& vGroupMemKeys);
	void DumpBin_Group_OutMemKeys(CAuxBinStrVect& oStr, vector<int>& vGroupMemKeys);
	void DumpBinParse_Group(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers);
	//void DumpBin(CAuxBinStrVect& oStr, radTmhg& mEl, radThg& hg);
	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey);

	radTg3dGraphPresent* CreateGraphPresent();

	int DuplicateItself(radThg& hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
	{
		return DuplicateGroupStuff(new radTGroup(*this), hg, radPtr, PutNewStuffIntoGenCont);
	}
	int DuplicateGroupStuff(radTGroup*, radThg&, radTApplication*, char);
	virtual int DuplicateWithoutDuplicatingGroupStuff(radThg& hgGroup)
	{
		radTSend Send;
		radTGroup* pNewGroup = new radTGroup(*this);
		if(pNewGroup==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgLoc(pNewGroup);
		hgGroup = hgLoc;
		return 1;
	}

	int CreateFromSym(radThg&, radTApplication*, char);

	int SubdivideItself(double*, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideItselfAsWholeInLabFrame(double*, radThg&, radTApplication*, radTSubdivOptions*);
	int SetUpCuttingPlanes(TVector3d&, double*, radTSubdivOptions*, TVector3d*);
	int FindLowestAndUppestVertices(TVector3d&, radTSubdivOptions*, TVector3d&, TVector3d&, radTrans&, char&, char&);

	int SubdivideItselfByEllipticCylinder(double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideItselfByEllipticCylinderAsWholeInLabFrame(double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	int FindEdgePointsOverPhiAndAxForCylSubd(radTCylindricSubdivSpec*, TVector3d*, double*);
	int ConvertToPolyhedron(radThg&, radTApplication*, char);
	int EstimateCenterPointOverRelaxables(TVector3d&);
	int SubdivideItselfOverAzimuth(double*, double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	int FindEdgePointsOverEllipseSet0(double*, radTCylindricSubdivSpec*, radThg, TVector3d*, double*, radTSubdivOptions*);
	int FindEdgePointsOverEllipseSet(double*, radTCylindricSubdivSpec*, radThg, TVector3d*, double*, radTSubdivOptions*);
	int SubdivideByEllipses(double*, double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideByPlanesPerpToCylAx(double*, double*, TVector3d*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	void FlattenNestedStructureIfMessageCharIsSet(radTApplication*, char, radTSubdivOptions*);
	void JustTraverse();

	int CutItself(TVector3d*, radThg&, radTPair_int_hg&, radTPair_int_hg&, radTApplication*, radTSubdivOptions*);

	int SubdivideItselfByParPlanes(double*, int, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideItselfByParPlanesAsWholeInLabFrame(double*, int, radThg&, radTApplication*, radTSubdivOptions*);

	int SubdivideItselfByOneSetOfParPlanes(TVector3d&, TVector3d*, int, radThg&, radTApplication*, radTSubdivOptions*, radTvhg*);

	int SetMaterial(radThg&, radTApplication*);
	void SetM(TVector3d& M); //virtual
	inline int ScaleCurrent(double); //virtual in radTg3d

	void FlattenNestedStructure(radTApplication*, char);
	void CollectNonGroupElements(radTmhg*, int&, radTApplication*, char);

	inline int ItemIsNotFullyInternalAfterCut();

	inline int NumberOfDegOfFreedom();
	inline int SizeOfThis();

	void Push_backCenterPointAndField(radTFieldKey*, radTVectPairOfVect3d*, radTrans*, radTg3d*, radTApplication*);

	inline double Volume();
	inline void LimitsAtTransform(radTrans* pExtTr, double* LimArr);
	inline void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique);

	inline void SimpleEnergyComp(radTField*);
	void ActualEnergyForceTorqueCompWithAdd(radTField*);
	inline void MarkFurtherSubdNeed(char, char, char);
	inline void MarkFurtherSubdNeed1D(char, char);
	int NextStepEnergyForceTorqueComp(double*, radThg&, radTField*, char&);
	int ProceedNextStepEnergyForceTorqueComp(double*, radThg&, radTField*, radTField*, char&, char);
	inline void SetupFurtherSubdInd(char);

	inline void SetMessageChar(char);

	radTGroup* CreateGroupIncludingAllMembersExceptIt(const radTmhg::const_iterator&);
};

//-------------------------------------------------------------------------

inline void radTGroup::AddElement(int ElemKey, const radThg& hg) 
{ 
	GroupMapOfHandlers[ElemKey] = hg;
	((radTg3d*)(hg.rep))->IsGroupMember = 1;
}

//-------------------------------------------------------------------------

inline void radTGroup::B_comp(radTField* FieldPtr)
{
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		((radTg3d*)(((*iter).second).rep))->B_genComp(FieldPtr);  // To check carefully!!!
}

//-------------------------------------------------------------------------

inline int radTGroup::NumberOfDegOfFreedom() 
{
	int DegFrCount = 0;
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		DegFrCount += ((radTg3d*)(((*iter).second).rep))->NumberOfDegOfFreedom();  // To check carefully!!!
	return DegFrCount;
}

//-------------------------------------------------------------------------

inline int radTGroup::SizeOfThis()
{
	int GenSize = sizeof(*this);
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		GenSize += (((*iter).second).rep)->SizeOfThis();
	return GenSize;
}

//-------------------------------------------------------------------------

inline int radTGroup::ItemIsNotFullyInternalAfterCut()
{
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		if(((radTg3d*)(((*iter).second).rep))->ItemIsNotFullyInternalAfterCut()) return 1;
	return 0;
}

//-------------------------------------------------------------------------

inline double radTGroup::Volume()
{
	double SumVol = 0.;
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
		//SumVol += ((radTg3d*)(((*iter).second).rep))->Volume();
		SumVol += ((radTg3d*)(((*iter).second).rep))->VolumeWithSym();

	return SumVol;
}

//-------------------------------------------------------------------------

inline void radTGroup::LimitsAtTransform(radTrans* pExtTr, double* LimArr)
{
	double &xMin = LimArr[0], &xMax = LimArr[1];
	double &yMin = LimArr[2], &yMax = LimArr[3];
	double &zMin = LimArr[4], &zMax = LimArr[5];
	xMin = 1e+23; xMax = -1e+23;
	yMin = 1e+23; yMax = -1e+23;
	zMin = 1e+23; zMax = -1e+23;

	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		double LocLimArr[6];
		((radTg3d*)(((*iter).second).rep))->Limits(pExtTr, LocLimArr);
		if(xMin > LocLimArr[0]) xMin = LocLimArr[0];
		if(xMax < LocLimArr[1]) xMax = LocLimArr[1];
		if(yMin > LocLimArr[2]) yMin = LocLimArr[2];
		if(yMax < LocLimArr[3]) yMax = LocLimArr[3];
		if(zMin > LocLimArr[4]) zMin = LocLimArr[4];
		if(zMax < LocLimArr[5]) zMax = LocLimArr[5];
	}
}

//-------------------------------------------------------------------------

inline void radTGroup::VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique)
{
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		((radTg3d*)(((*iter).second).rep))->VerticesInLocFrame(OutVect, EnsureUnique);
}

//-------------------------------------------------------------------------

inline void radTGroup::SimpleEnergyComp(radTField* FieldPtr)
{
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
	{
		radTg3d* g3dPtr = (radTg3d*)(((*iter).second).rep);
		if(g3dPtr->g3dListOfTransform.empty()) g3dPtr->SimpleEnergyComp(FieldPtr);
		else
		{
			g3dPtr->NestedFor_Energy(FieldPtr, g3dPtr->g3dListOfTransform.begin());
		}
	}
}

//-------------------------------------------------------------------------

inline void radTGroup::MarkFurtherSubdNeed(char SubdNeedX, char SubdNeedY, char SubdNeedZ)
{
	radTg3d::MarkFurtherSubdNeed(SubdNeedX, SubdNeedY, SubdNeedZ);

	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		((radTg3d*)(((*iter).second).rep))->MarkFurtherSubdNeed(SubdNeedX, SubdNeedY, SubdNeedZ);
}

//-------------------------------------------------------------------------

inline void radTGroup::MarkFurtherSubdNeed1D(char SubdNeed, char XorYorZ)
{
	radTg3d::MarkFurtherSubdNeed1D(SubdNeed, XorYorZ);

	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		((radTg3d*)(((*iter).second).rep))->MarkFurtherSubdNeed1D(SubdNeed, XorYorZ);
}

//-------------------------------------------------------------------------

inline void radTGroup::SetupFurtherSubdInd(char InSubdInd)
{
	radTg3d::SetupFurtherSubdInd(InSubdInd);

	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		((radTg3d*)(((*iter).second).rep))->SetupFurtherSubdInd(InSubdInd);
}

//-------------------------------------------------------------------------

inline void radTGroup::SetMessageChar(char InMessageChar)
{
	MessageChar = InMessageChar;
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		((radTg3d*)(((*iter).second).rep))->SetMessageChar(InMessageChar);
}

//-------------------------------------------------------------------------

inline int radTGroup::ScaleCurrent(double scaleCoef)
{
	int scalingWasApplied = 0;
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		if(((radTg3d*)(((*iter).second).rep))->ScaleCurrent(scaleCoef)) scalingWasApplied = 1;
	}
	return scalingWasApplied;
}

//-------------------------------------------------------------------------

#endif
