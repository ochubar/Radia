/*-------------------------------------------------------------------------
*
* File name:      radg3d.h
*
* Project:        RADIA
*
* Description:    Base class for 3D objects - magnetic field sources;
*                 auxiliary classes/structures for field computation
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

//-------------------------------------------------------------------------
//	Definition of class radTg3d - a class of objects capable
//	to generate magnetic field.
//	radTg3d is a parent for radTg3dRelax, radTArcCur, ...
//-------------------------------------------------------------------------

#ifndef __RADG3D_H
#define __RADG3D_H

#include "radmater.h"
#include "radg.h"
#include "gmvect.h"
#include "radmamet.h"

#include <list>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

#ifdef __GCC__
typedef list <radTPair_int_hg> radTlphg;
#else
typedef list <radTPair_int_hg, allocator<radTPair_int_hg> > radTlphg; // Porting
#endif

//-------------------------------------------------------------------------

class radTField;

//-------------------------------------------------------------------------

struct radTAuxCompDataG3D {
	TVector3d Force, Torque;
	double Energy;
	char SubdNeedInd;

	radTAuxCompDataG3D() { SubdNeedInd = 0;}
	~radTAuxCompDataG3D() {}

	void MarkSubdNeed(char SubdNeedX, char SubdNeedY, char SubdNeedZ)
	{// This accepts only 0 or 1
		char SubdNeedX_Code = SubdNeedX << 2;
		char SubdNeedY_Code = SubdNeedY << 1;
		SubdNeedInd = SubdNeedX_Code | SubdNeedY_Code | SubdNeedZ;
	}
	void MarkSubdNeed1D(char SubdNeed, char XorYorZ)
	{// This accepts only 0 or 1
		char SubdNeedCode, OtherSubdNeed;
		if(XorYorZ == 'x')
		{
			SubdNeedCode = SubdNeed << 2;
			OtherSubdNeed = SubdNeedInd & 3;
		}
		else if(XorYorZ == 'y')
		{
			SubdNeedCode = SubdNeed << 1;
			OtherSubdNeed = SubdNeedInd & 5;
		}
		else if(XorYorZ == 'z')
		{
			SubdNeedCode = SubdNeed;
			OtherSubdNeed = SubdNeedInd & 6;
		}
		SubdNeedInd = SubdNeedCode | OtherSubdNeed;
	}
	void ShowSubdNeed(char& SubdNeedX, char& SubdNeedY, char& SubdNeedZ)
	{// This gives 0 or 1
		SubdNeedZ = SubdNeedInd & 1;
		SubdNeedY = (SubdNeedInd >> 1) & 1;
		SubdNeedX = SubdNeedInd >> 2;
	}
	void ShowSubdNeed1D(char& SubdNeed, char XorYorZ)
	{
		if(XorYorZ == 'x') SubdNeed = SubdNeedInd >> 2;
		else if(XorYorZ == 'y') SubdNeed = (SubdNeedInd >> 1) & 1;
		else if(XorYorZ == 'z') SubdNeed = SubdNeedInd & 1;
	}
	inline void StoreDataFromField(radTField*);
	inline void PutStoredDataToField(radTField*);
};

//-------------------------------------------------------------------------

struct radTSubdivOptions {
	char SubdivisionFrame; // 0- Local; 1- Laboratory, Group as whole; 2- Laboratory, Group Members individually;
	char SubdivisionParamCode; // 0- kx,ky,kz are subdiv. numbers; 1- kx,ky,kz are average sizes of pieces;
	char SubdivideCoils; // 0- do not subdivide coils; 1- subdivide coils;
	char PutNewStuffIntoGenCont; // 0- do not put; 1- put;
	char ReplaceOldStuff;
	char SeparatePiecesAtCutting; // 0- do not separate; 1- separate;
	char MapInternalFacesAfterCut;
	char MethForRadialSegmAtEllCylSubd; // 0- Radial segmentation based on position of point defining the Ellipse axis; 1- Radial segmentation based on dimensions of the object with respect to the Elliptic cylinder of subdivision;

	radTSubdivOptions() { SeparatePiecesAtCutting = 0; MapInternalFacesAfterCut = 1;}
};

//-------------------------------------------------------------------------

struct radTCylindricSubdivSpec {
	TVector3d PointOnCylAx, CylAxVect;
	TVector3d PointOnEllAx;
	double EllAxRatio;
	char EllAxNotDefined;

	radTCylindricSubdivSpec(double* SubdivSpecData, int LenSubdivSpecData) 
	{
		if(LenSubdivSpecData == 0) return;
		PointOnCylAx.x = *(SubdivSpecData++);
		PointOnCylAx.y = *(SubdivSpecData++);
		PointOnCylAx.z = *(SubdivSpecData++);
		CylAxVect.x = *(SubdivSpecData++);
		CylAxVect.y = *(SubdivSpecData++);
		CylAxVect.z = *(SubdivSpecData++);
		if(LenSubdivSpecData > 6)
		{
			PointOnEllAx.x = *(SubdivSpecData++);
			PointOnEllAx.y = *(SubdivSpecData++);
			PointOnEllAx.z = *(SubdivSpecData++);
			EllAxRatio = *SubdivSpecData;
			EllAxNotDefined = 0;
		}
		else
		{
			PointOnEllAx = PointOnCylAx;
			EllAxRatio = 1.;
			EllAxNotDefined = 1;
		}
		double InvLen = 1./sqrt(CylAxVect.x*CylAxVect.x + CylAxVect.y*CylAxVect.y + CylAxVect.z*CylAxVect.z);
		CylAxVect = InvLen*CylAxVect;
	}
	radTCylindricSubdivSpec() {}
};

//-------------------------------------------------------------------------

typedef radTHandle<radTAuxCompDataG3D> radTHandleAuxCompDataG3D;
class radTApplication;
class radTg3dGraphPresent;
class radTField;
struct radTFieldKey;
class radTrans;
class radTGroup;

//-------------------------------------------------------------------------

class radTg3d : public radTg {
public:
	radTlphg g3dListOfTransform; // Don't make it private!!!
	int IsGroupMember;
	char ConsiderOnlyWithTrans;
	TVector3d CentrPoint; //moved from derived classes OC061008

	radTHandleAuxCompDataG3D HandleAuxCompData;
	char MessageChar; 
	//double gCurrentScaleCoef; //required for current-carrying objects?

	radTg3d() 
	{ 
		IsGroupMember = 0; ConsiderOnlyWithTrans = 0; //gCurrentScaleCoef = 1;
	}
	~radTg3d() {}

	int Type_g() { return 1;}
	virtual int Type_g3d() { return 0;}

	inline void AddTransform(int, const radThg&);
	inline void AddTransform_OtherSide(int Multiplicity, const radThg& hg);

	inline void FindResTransfWithMultOne(radTrans&, short&);
	inline void FindInnerTransfWithMultOne(radTrans&, short&);
	inline void EraseOuterTransform();
	inline void EraseInnerTransform();
	inline void EraseAllTransformations();

	void FlattenSpaceTransforms(radTvhg&);
	int TotalMultiplicity();

	inline void B_genComp(radTField*);
	virtual void B_comp(radTField*) {}
	void NestedFor_B(radTField*, const radTlphg::iterator&);
	inline void B_comp_Or_NestedFor(radTField*, const radTlphg::iterator&);

	virtual void B_intComp(radTField*) {}
	void B_intCompFinNum(radTField*);

	void EnergyForceTorqueComp(radTField*);
	void ActualEnergyForceTorqueComp(radTField*);
	void NestedFor_Energy(radTField*, const radTlphg::iterator&);
	void Energy_Or_NestedFor(radTField* FieldPtr, const radTlphg::iterator& Iter)
	{
		if(Iter == g3dListOfTransform.end()) SimpleEnergyComp(FieldPtr);
		else NestedFor_Energy(FieldPtr, Iter);
	}
	void CreateAuxCompData() { HandleAuxCompData = radTHandleAuxCompDataG3D(new radTAuxCompDataG3D());}
	void EnergyForceTorqueCompAutoDestSubd(radTField*);
	char CheckIfMoreEnrFrcTrqCompNeededAndUpdate(radTField*, radTField*);
	virtual void SimpleEnergyComp(radTField*) {}
	virtual int NextStepEnergyForceTorqueComp(double*, radThg&, radTField*, char&);
	virtual void ActualEnergyForceTorqueCompWithAdd(radTField*);
	virtual int ProceedNextStepEnergyForceTorqueComp(double*, radThg&, radTField*, radTField*, char&, char);
	virtual void MarkFurtherSubdNeed(char SubdNeedX, char SubdNeedY, char SubdNeedZ)
	{// This is for everything except for Groups and their childs.
		if(HandleAuxCompData.rep != 0) HandleAuxCompData.rep->MarkSubdNeed(SubdNeedX, SubdNeedY, SubdNeedZ);
	}
	virtual void MarkFurtherSubdNeed1D(char SubdNeed, char XorYorZ)
	{// This is for everything except for Groups and their childs.
		if(HandleAuxCompData.rep != 0) HandleAuxCompData.rep->MarkSubdNeed1D(SubdNeed, XorYorZ);
	}
	virtual void SetupFurtherSubdInd(char InSubdInd)
	{
		if(HandleAuxCompData.rep != 0) HandleAuxCompData.rep->SubdNeedInd = InSubdInd;
	}

	void NormStressTensor(radTField*);

	void IntOverShape_Or_NestedFor(radTField* FieldPtr, const radTlphg::iterator& Iter)
	{
		if(Iter == g3dListOfTransform.end()) IntOverShape(FieldPtr);
		else NestedFor_IntOverShape(FieldPtr, Iter);
	}
	void NestedFor_IntOverShape(radTField*, const radTlphg::iterator&);
	virtual void IntOverShape(radTField*) {}

	virtual radTg3dGraphPresent* CreateGraphPresent() { return NULL;}

	virtual int SubdivideItself(double*, radThg&, radTApplication*, radTSubdivOptions*) { return 1;}
	virtual int CutItself(TVector3d*, radThg&, radTPair_int_hg&, radTPair_int_hg&, radTApplication*, radTSubdivOptions*) { return 1;}
	virtual int SubdivideItselfByParPlanes(double*, int, radThg&, radTApplication*, radTSubdivOptions*) { return 1;}
	virtual int SubdivideItselfByOneSetOfParPlanes(TVector3d&, TVector3d*, int, radThg&, radTApplication*, radTSubdivOptions*, radTvhg*) { return 1;}
	void CheckAxesExchangeForSubdInLabFrame(double*, char&);
	int TransferSubdivisionStructToLocalFrame(TVector3d&, TVector3d*, int, radTSubdivOptions*, TVector3d&, TVector3d*&);
	virtual int ConvertToPolyhedron(radThg&, radTApplication*, char) { return 0;} // 0 is essential
	virtual int FindLowestAndUppestVertices(TVector3d&, radTSubdivOptions*, TVector3d&, TVector3d&, radTrans&, char&, char&) { return 1;}

	virtual int SubdivideItselfByEllipticCylinder(double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*) { return 1;}

	int FinishDuplication(radTg3d* g3dPtr, radThg& hg)
	{
		radTSend Send;
		if(g3dPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		g3dPtr->IsGroupMember = 0;
		radThg hgLoc(g3dPtr); hg = hgLoc; return 1;
	}

	virtual int CreateFromSym(radThg&, radTApplication*, char);
	int NestedFor_CreateFromSym(radTGroup*, radTApplication*, char, radTrans*, const radTlphg::iterator&);
	int CreateAndAddToGroupOrNestedFor(radTGroup*, radTApplication*, char, radTrans*, const radTlphg::iterator&);

	virtual int SetMaterial(radThg&, radTApplication*) { return 0;}
	virtual void SetM(TVector3d&) {}
	virtual int ScaleCurrent(double) { return 0;} //implemented in current-carying objects

	virtual int NumberOfDegOfFreedom() { return 0;}
	virtual int ItemIsNotFullyInternalAfterCut() { return 1;}

	virtual double Volume() { return 0.;}
	double VolumeWithSym();

	void Limits(radTrans*, double*);
	void DeterminePointsLimits(TVector3d* Points, int AmOfPoints, double*);
	virtual void VerticesInLocFrame(radTVectorOfVector3d&, bool) {}
	virtual void LimitsAtTransform(radTrans*, double*);

	virtual void Push_backCenterPointAndField(radTFieldKey*, radTVectPairOfVect3d*, radTrans*, radTg3d*, radTApplication*) {}
	inline void GetTrfAndCenPointInLabFrame(radTrans* pBaseTrans, radTrans& bufTrans, radTrans*& pResTrans, TVector3d& cenPointInLabFr);

	inline double TransAtans(double, double, double&);
	inline double Argument(double x, double y); 

	double Abs(double x) { return (x<0.)? -x : x;}
	double Max(double x1, double x2) { return (x1<x2)? x2 : x1;}
	double Sign(double x) { return (x<0.)? -1. : 1.;}
	double Step(double x) { return (x>0.)? 1. : 0.;}

	void FindEllipticCoordOfPoint(radTCylindricSubdivSpec*, TVector3d&, double&, double&);
	double EstimateLengthAlongEllipse(double, double, double, double);
	inline char AngleIsBetween(double, double, double);
	inline double AngularDifference(double, double);

	short CheckIfPosEven(int ii)
	{
		double x = 0.5*(double(ii) + 1.E-08); return ((ii > 0) && ((x - int(x)) < 0.1));
	}

	void Dump(std::ostream& o, int ShortSign =0) // Porting
	{
		radTg::Dump(o);
		o << "Magnetic field source object: ";
	}
	
	//void DumpBin_g3d_TreatTrfs(CAuxBinStrVect& oStr, radTmhg& mEl, vector<pair<int, int> >& vTrfKeys) 
	//void DumpBin_g3d_TreatTrfs(CAuxBinStrVect& oStr, radTmhg& mEl, radTmhg& gMapOfHandlers, int& gUniqueMapKey, vector<pair<int, int> >& vTrfKeys) 
	void DumpBin_g3d_TreatTrfs(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, vector<pair<int, int> >& vTrfKeys) 
	{//this will be called from derived classes' DumpBin
		for(radTlphg::reverse_iterator it = g3dListOfTransform.rbegin(); it != g3dListOfTransform.rend(); ++it) 
		{
			int existKey = 0;
			radThg &curTr_hg = it->Handler_g;
			//for(radTmhg::iterator mit = mEl.begin(); mit != mEl.end(); ++mit)
			//{
			//	if(mit->second == curTr_hg)
			//	{
			//		existKey = mit->first; break;
			//	}
			//}

			//if(existKey == 0)
			//{
			for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
			{
				if(mit->second == curTr_hg)
				{
					existKey = mit->first; break;
				}
			}
			if(existKey == 0)
			{
				existKey = gUniqueMapKey;
				gMapOfHandlers[gUniqueMapKey++] = curTr_hg;
			}

			//curTr_hg.rep->DumpBin(oStr, mEl, curTr_hg); //adding element to the map (mEl) should happen here
			//existKey = (int)mEl.size();
			int indExist = CAuxParse::FindElemInd(existKey, vElemKeysOut);
			if(indExist < 0)
			{
				curTr_hg.rep->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, existKey); //adding element to the map (mEl) should happen here
			}
			//}
			int mult = it->m;
			vTrfKeys.push_back(pair<int, int>(existKey, mult));
		}
	}

	void DumpBin_g3d(CAuxBinStrVect& oStr, vector<pair<int, int> >& vTrfKeys)
	{
		//radTlphg g3dListOfTransform;
		int nTrfs = (int)vTrfKeys.size();
		oStr << nTrfs;
		for(int i=0; i<nTrfs; i++)
		{
			pair<int, int> &p = vTrfKeys[i];
			oStr << p.first; oStr << p.second;
		}

		//int IsGroupMember;
		oStr << IsGroupMember;

		//char ConsiderOnlyWithTrans;
		oStr << ConsiderOnlyWithTrans;

		//TVector3d CentrPoint;
		oStr << CentrPoint;

		//char MessageChar; 
		oStr << MessageChar;

		//radTHandleAuxCompDataG3D HandleAuxCompData;
		radTAuxCompDataG3D *pAux = HandleAuxCompData.rep;
		if(pAux != 0)
		{
			oStr << (char)1;
			//radTAuxCompDataG3D
			//TVector3d Force;
			oStr << pAux->Force;

			//TVector3d Torque;
			oStr << pAux->Torque;

			//double Energy;
			oStr << pAux->Energy;

			//char SubdNeedInd;
			oStr << pAux->SubdNeedInd;
		}
		else oStr << (char)0;
	}

	void DumpBinParse_g3d(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers);
	void DumpTransApplied(std::ostream& o); // Porting

	static int Round(double dVal) 
	{
		int iBuf = int(dVal);
		return ((dVal - iBuf) < 0.499)? iBuf : (iBuf + 1);
	}

	static bool CheckIfThreePointsAreOnOneLine(const TVector3d& P1, const TVector3d& P2, const TVector3d& P3, double absTol)
	{//to test!
		TVector3d v1 = P2 - P1, v2 = P3 - P1;
		double v1e2 = v1.x*v1.x + v1.y*v1.y + v1.z*v1.z;
		double v2e2 = v2.x*v2.x + v2.y*v2.y + v2.z*v2.z;
		double absTolE2 = absTol*absTol;
		if((v1e2 < absTolE2) || (v2e2 < absTolE2)) return true;

		double maxAbsVe2 = (v1e2 > v2e2)? v1e2 : v2e2;
		TVector3d vTestVectProd = v1^v2;
		if(vTestVectProd.Abs() < absTol*sqrt(maxAbsVe2)) return true;
		else return false;
	}

	virtual void SetMessageChar(char InMessageChar) { MessageChar = InMessageChar;}
};

//-------------------------------------------------------------------------

struct radTFieldKey { 
	short B_, H_, A_, M_, J_, Phi_, PreRelax_, Ib_, Ih_, FinInt_, Force_, ForceEnr_, Torque_, Energy_, Q_;
	radTFieldKey(short InB_ =0, short InH_ =0, short InA_ =0, short InM_ =0, short InJ_ =0, short InPhi_ =0, short InPreRelax_ =0, short InIb_ =0, short InIh_ =0, short InFinInt_ =0, short InForce_ =0, short InForceEnr_ =0, short InTorque_ =0, short InEnergy_ =0, short InQ_ =0)
	{ 
		B_=InB_; H_=InH_; A_=InA_; M_=InM_; J_=InJ_; Phi_=InPhi_; PreRelax_=InPreRelax_; Ib_=InIb_; Ih_=InIh_; FinInt_=InFinInt_; Force_=InForce_; ForceEnr_=InForceEnr_; Torque_=InTorque_; Energy_=InEnergy_; Q_=InQ_;
		if(Q_)
		{
			B_= H_= PreRelax_= 1;
		}
	}
};

//-------------------------------------------------------------------------

struct radTCompCriterium {

	short BasedOnPrecLevel; // Actually this is used nowhere at the moment
	double AbsPrecB;
	double AbsPrecA;
	double AbsPrecB_int;
	double AbsPrecForce;
	double AbsPrecTorque;
	double AbsPrecEnergy;
	double AbsPrecTrjCoord;
	double AbsPrecTrjAngle;
	double MltplThresh[4]; // Threshold ratios for 4 diff. orders of multipole approx. at field computation

	double WorstRelPrec;

	char BasedOnWorstRelPrec; // Used at energy - force computation

	radTCompCriterium() 
	{
		//Default values for all the Project:
		AbsPrecB = 0.0001; // Tesla
		AbsPrecA = 0.001;  // Tesla * mm
		AbsPrecB_int = 0.001;  // Tesla * mm
		AbsPrecForce = 1.;  // Newton
		AbsPrecTorque = 10.;  // Newton * mm
		AbsPrecEnergy = 10.;
		AbsPrecTrjCoord = AbsPrecTrjAngle = -1.;

		WorstRelPrec = 0.1;  // Used for Force computation through energy
		BasedOnWorstRelPrec = 0;

		MltplThresh[0] = 0.; // No Computatation
		MltplThresh[1] = 0.; // Dipole only
		MltplThresh[2] = 0.; // Dipole + Quadr.
		MltplThresh[3] = 0.; // Dipole + Quadr. + Next order
	}
};

//-------------------------------------------------------------------------

struct radTStructForShapeInt {
	radThg HandleOfSource, HandleOfShape;
	int IntegrandLength; // Number of elements in TVector3d* to be integrated over a Shape
	TVector3d Normal;
	TVector3d* VectArray;
	char* VectTypeArray; // 'a' - axial, 'r' - regular
	void (radTg3d::*IntegrandFunPtr)(radTField*);
	short IntOverLine_, IntOverSurf_, IntOverVol_;
	double* AbsPrecArray;

	radTStructForShapeInt(const radTStructForShapeInt&);
	radTStructForShapeInt() {}
	~radTStructForShapeInt() {}
};

//-------------------------------------------------------------------------

struct radTStructForEnergyForceTorqueComp {
	radThg hSource, hDest;
	double* DestSubdivArray;
	radTApplication* radPtr;
	char AutoDestSubdivision;
	char SomethingIsWrong;

	radTStructForEnergyForceTorqueComp() { SomethingIsWrong = 0;}
};

//-------------------------------------------------------------------------

typedef radTHandle<radTStructForEnergyForceTorqueComp> radTHandleStructForEnergyForceTorqueComp;

//-------------------------------------------------------------------------

class radTField {
public:
	TVector3d P, B, H, A, M, J, Ib, Ih, NextP, Force, Torque;
	double Phi, Energy;
	int AmOfIntrctElemWithSym; // Place this into separate Relaxation dedicated structure, if more Relax specific data appear
	short PointIsInsideFrame; // This is only used with Polyhedrons (or not used at all ?)

	radTFieldKey FieldKey;
	radTCompCriterium CompCriterium;

	radTHandleStructForEnergyForceTorqueComp HandleEnergyForceTorqueCompData;
	radTStructForShapeInt* ShapeIntDataPtr;

	radTField(const radTFieldKey& InFieldKey, const radTCompCriterium& InCompCriterium, 
			  const TVector3d& InP, const TVector3d& InB, const TVector3d& InH, 
			  const TVector3d& InA, const TVector3d& InM, const TVector3d& InJ, double InPhi =0.)
	{ 
		FieldKey = InFieldKey;
		P = InP; B = InB; H = InH; A = InA; M = InM; J = InJ; Phi = InPhi;
		CompCriterium = InCompCriterium;
		PointIsInsideFrame = 0;
	}
	radTField(const radTFieldKey& InFieldKey, const radTCompCriterium& InCompCriterium, 
			  const TVector3d& InP, const TVector3d& InNextP, const TVector3d& InIb, const TVector3d& InIh)
	{// This is for Field Integral
		FieldKey = InFieldKey;
		P = InP; NextP = InNextP; Ib = InIb; Ih = InIh;
		CompCriterium = InCompCriterium;
		PointIsInsideFrame = 0;
	}
	radTField(const radTFieldKey& InFieldKey, const radTCompCriterium& InCompCriterium, 
			  radTStructForShapeInt* InShapeIntDataPtr =NULL)
	{// This is used in NestedFor_IntOverShape
		FieldKey = InFieldKey; CompCriterium = InCompCriterium; ShapeIntDataPtr = InShapeIntDataPtr;
		TVector3d ZeroVect(0.,0.,0.);
		B = H = A = M = J = Ib = Ih = Force = Torque = ZeroVect; // Add here more members should they appear
		Phi = Energy = 0.;
		PointIsInsideFrame = 0;
	}
	radTField(const radTFieldKey& InFieldKey, const radTCompCriterium& InCompCriterium, 
			  radTHandleStructForEnergyForceTorqueComp& InHandleEnergyForceTorqueCompData)
	{// This is used in NestedFor_IntOverShape
		FieldKey = InFieldKey; CompCriterium = InCompCriterium; 
		HandleEnergyForceTorqueCompData = InHandleEnergyForceTorqueCompData;
		TVector3d ZeroVect(0.,0.,0.);
		B = H = A = M = J = Ib = Ih = Force = Torque = ZeroVect; // Add here more members should they appear
		Phi = Energy = 0.;
		PointIsInsideFrame = 0;
	}
	radTField(const radTFieldKey& InFieldKey, const TVector3d& InP, const TVector3d& InB, 
			  const TVector3d& InH, const TVector3d& InA, const TVector3d& InM, const TVector3d& InJ, double InPhi =0.)
	{ 
		FieldKey = InFieldKey;
		P = InP; B = InB; H = InH; A = InA; M = InM; J = InJ; Phi = InPhi;
		CompCriterium.BasedOnPrecLevel = 0;
		PointIsInsideFrame = 0;
	}
	radTField(const radTFieldKey& InFieldKey, const TVector3d& InP, const TVector3d& InVect, double InPhi =0.)
	{ 
		P = InP;
		FieldKey = InFieldKey;
		CompCriterium.BasedOnPrecLevel = 0;
		if(FieldKey.B_) B = InVect;
		else if(FieldKey.H_) H = InVect;
		else if(FieldKey.A_) A = InVect;
		else if(FieldKey.M_) M = InVect;
		else if(FieldKey.J_) J = InVect;
		else if(FieldKey.Phi_) Phi = InPhi;
		PointIsInsideFrame = 0;
	}
	radTField(const TVector3d& InP, const TVector3d& InB) 
	{ 
		P = InP; B = InB;
		CompCriterium.BasedOnPrecLevel = 0;
		PointIsInsideFrame = 0;
	}
	radTField() { PointIsInsideFrame = 0;}

	radTField& operator +=(const radTField& AnotherField)
	{
		if(FieldKey.B_) B+=AnotherField.B;
		if(FieldKey.H_) H+=AnotherField.H;
		if(FieldKey.A_) A+=AnotherField.A;
		if(FieldKey.M_) M+=AnotherField.M;
		if(FieldKey.J_) J+=AnotherField.J;
		if(FieldKey.Phi_) Phi+=AnotherField.Phi;
		if(FieldKey.Ib_) Ib+=AnotherField.Ib;
		if(FieldKey.Ih_) Ih+=AnotherField.Ih;
		
		if(FieldKey.Q_) //matrix for H calculation //OC191005
		{
			//B+=AnotherField.B; //assumed to be done already, because Q_ enforces B_ and H_
			//H+=AnotherField.H;
			if(!FieldKey.A_) A+=AnotherField.A;
		}
		
		return *this;
	}

	inline friend radTField operator +(const radTField&, const radTField&);
	inline friend radTField operator -(const radTField&, const radTField&);
};

//-------------------------------------------------------------------------

inline radTField operator +(const radTField& F1, const radTField& F2)
{
	radTField resF(F1);
	if(F1.FieldKey.B_ && F2.FieldKey.B_) resF.B += F2.B;
	if(F1.FieldKey.H_ && F2.FieldKey.H_) resF.H += F2.H;
	if(F1.FieldKey.A_ && F2.FieldKey.A_) resF.A += F2.A;
	if(F1.FieldKey.M_ && F2.FieldKey.M_) resF.M += F2.M;
	if(F1.FieldKey.J_ && F2.FieldKey.J_) resF.J += F2.J;
	if(F1.FieldKey.Phi_ && F2.FieldKey.Phi_) resF.Phi += F2.Phi;
	if(F1.FieldKey.Ib_ && F2.FieldKey.Ib_) resF.Ib += F2.Ib;
	if(F1.FieldKey.Ih_ && F2.FieldKey.Ih_) resF.Ih += F2.Ih;
	
	if(F1.FieldKey.Q_ && F2.FieldKey.Q_) //matrix for H calculation //OC191005
	{
		//resF.B += F2.B; //assumed to be done already, because Q_ enforces B_ and H_
		//resF.H += F2.H;
		if(!resF.FieldKey.A_) resF.A += F2.A;
	}
	return resF;
}

//-------------------------------------------------------------------------

inline radTField operator -(const radTField& F1, const radTField& F2)
{
	radTField resF(F1);
	if(F1.FieldKey.B_ && F2.FieldKey.B_) resF.B -= F2.B;
	if(F1.FieldKey.H_ && F2.FieldKey.H_) resF.H -= F2.H;
	if(F1.FieldKey.A_ && F2.FieldKey.A_) resF.A -= F2.A;
	if(F1.FieldKey.M_ && F2.FieldKey.M_) resF.M -= F2.M;
	if(F1.FieldKey.J_ && F2.FieldKey.J_) resF.J -= F2.J;
	if(F1.FieldKey.Phi_ && F2.FieldKey.Phi_) resF.Phi -= F2.Phi;
	if(F1.FieldKey.Ib_ && F2.FieldKey.Ib_) resF.Ib -= F2.Ib;
	if(F1.FieldKey.Ih_ && F2.FieldKey.Ih_) resF.Ih -= F2.Ih;
	
	if(F1.FieldKey.Q_ && F2.FieldKey.Q_) //matrix for H calculation //OC191005
	{
		//resF.B -= F2.B; //assumed to be done already, because Q_ enforces B_ and H_
		//resF.H -= F2.H;
		if(!resF.FieldKey.A_) resF.A -= F2.A;
	}
	return resF;
}

//-------------------------------------------------------------------------

class radTg3dRelax : public radTg3d {
protected:
	TMatrix3d* pM_LinCoef;

public:
	radThg MaterHandle;

	//TVector3d CentrPoint; //moved to radTg3d
	TVector3d Magn;

	float AuxFloat1, AuxFloat2, AuxFloat3;

	radTg3dRelax(const TVector3d& InCPoiVect, const TVector3d& InMagnVect, const radThg& InMaterHandle)
	{
		CentrPoint=InCPoiVect; Magn=InMagnVect; MaterHandle=InMaterHandle;
		pM_LinCoef=0;
	}
	radTg3dRelax(const TVector3d& InMagnVect, const radThg& InMaterHandle)
	{
		Magn=InMagnVect; MaterHandle=InMaterHandle;
		pM_LinCoef=0;
	}
	radTg3dRelax(const TVector3d& InMagnVect)
	{
		Magn=InMagnVect;
		pM_LinCoef=0;
	}
	radTg3dRelax(const TVector3d& InMagnVect, TMatrix3d& InM_LinCoef)
	{
		Magn=InMagnVect; 
		if(!InM_LinCoef.isZero()) pM_LinCoef = new TMatrix3d(InM_LinCoef);
		else pM_LinCoef = 0;
	}
	//radTg3dRelax(const TVector3d& InMagnVect, TMatrix3d& InM_LinCoef, const radThg& InMaterHandle)
	radTg3dRelax(const TVector3d* pInMagnVect, TMatrix3d* pInM_LinCoef, const radThg& InMaterHandle)
	{
		if(pInMagnVect != 0) Magn = *pInMagnVect;

		pM_LinCoef = 0;
		if(pInM_LinCoef != 0)
		{
			if(!pInM_LinCoef->isZero()) pM_LinCoef = new TMatrix3d(*pInM_LinCoef);
		}
		MaterHandle = InMaterHandle;
	}
	radTg3dRelax() 
	{
		MaterHandle.rep = 0;
		pM_LinCoef=0;
	}
	radTg3dRelax(radTg3dRelax& aRelax) 
	{
		*this = aRelax;
		if(aRelax.pM_LinCoef != 0) pM_LinCoef = new TMatrix3d(*(aRelax.pM_LinCoef));
	}

	~radTg3dRelax()
	{//check if this is called
		if(pM_LinCoef != 0) delete pM_LinCoef;
	}

	int Type_g3d() { return 1;}
	virtual int Type_g3dRelax() { return 0;}

	void SimpleEnergyComp(radTField* FieldPtr)
	{
		const double PI = 3.14159265358979;
		const double ConstForM = -1./(4.*PI*100.);
		radTFieldKey LocFieldKey; LocFieldKey.B_ = 1;
		radTField LocField(LocFieldKey, FieldPtr->CompCriterium);
		LocField.P = CentrPoint;
		((radTg3d*)(FieldPtr->HandleEnergyForceTorqueCompData.rep->hSource.rep))->B_genComp(&LocField);
		FieldPtr->Energy += (ConstForM*Volume())*(Magn*LocField.B);
	}

	int SetMaterial(radThg& InMatHandle, radTApplication*) 
	{
		MaterHandle = InMatHandle;
		radTMaterial* MaterPtr = (radTMaterial*)(MaterHandle.rep);
		return MaterPtr->FinishSetup(Magn);
	}
	void SetM(TVector3d& M) //virtual
	{
		Magn = M;
	}	

	void Push_backCenterPointAndField(radTFieldKey*, radTVectPairOfVect3d*, radTrans*, radTg3d*, radTApplication*);

	virtual TVector3d& ReturnCentrPoint() { return CentrPoint;}
	virtual radTg3dRelax* FormalIntrctMemberPtr() { return this;}

	void CheckCenPtPositionWithRespectToPlane(TVector3d* CuttingPlane, char& CenPtPositionChar)
	{
		TVector3d& PointOnPlane = *CuttingPlane;
		TVector3d& PlaneNormal = CuttingPlane[1];
		TVector3d V = CentrPoint - PointOnPlane;
		double ScalProd = PlaneNormal*V;
		if(ScalProd >= 0.) CenPtPositionChar = 'H';
		else CenPtPositionChar = 'L';
	}

	void Dump(std::ostream& o, int ShortSign =0)
	{
		radTg3d::Dump(o);
		//o << "Relaxable: ";
	}
	
	void DumpBin_g3dRelax_TreatMat(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int& matKey)
	{
		matKey = 0;
		if(MaterHandle.rep != 0)
		{
			//MaterHandle.rep->DumpBin(oStr, mEl, MaterHandle); //adding element to the map (mEl) should happen here
			//matKey = (int)mEl.size();

			for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
			{
				if(mit->second == MaterHandle)
				{
					matKey = mit->first; break;
				}
			}
			if(matKey == 0)
			{
				matKey = gUniqueMapKey;
				gMapOfHandlers[gUniqueMapKey++] = MaterHandle;
			}

			//curTr_hg.rep->DumpBin(oStr, mEl, curTr_hg); //adding element to the map (mEl) should happen here
			//existKey = (int)mEl.size();
			int indExist = CAuxParse::FindElemInd(matKey, vElemKeysOut);
			if(indExist < 0)
			{
				MaterHandle.rep->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, matKey); //adding element to the map (mEl) should happen here
			}
		}
	}

	void DumpBin_g3dRelax(CAuxBinStrVect& oStr, int matKey)
	{
		//radThg MaterHandle;
		oStr << matKey;

		//TMatrix3d* pM_LinCoef;
		if(pM_LinCoef == 0) oStr << (char)0;
		else
		{
			oStr << (char)0;
			oStr << (*pM_LinCoef);
		}

		//TVector3d Magn;
		oStr << Magn;

		//float AuxFloat1, AuxFloat2, AuxFloat3;
		oStr << AuxFloat1 << AuxFloat2 << AuxFloat3;
	}

	void DumpBinParse_g3dRelax(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers);
	void DumpMaterApplied(std::ostream&);
};

//-------------------------------------------------------------------------

inline void radTg3d::AddTransform(int Multiplicity, const radThg& hg)
{
	radTPair_int_hg aPair(Multiplicity, hg);
	g3dListOfTransform.push_front(aPair);
}

//-------------------------------------------------------------------------

inline void radTg3d::AddTransform_OtherSide(int Multiplicity, const radThg& hg)
{
	radTPair_int_hg aPair(Multiplicity, hg);
	g3dListOfTransform.push_back(aPair);
}

//-------------------------------------------------------------------------

inline void radTg3d::EraseOuterTransform()
{
	g3dListOfTransform.pop_front();
}

//-------------------------------------------------------------------------

inline void radTg3d::EraseInnerTransform()
{
	g3dListOfTransform.pop_back();
}

//-------------------------------------------------------------------------

inline void radTg3d::EraseAllTransformations()
{
	g3dListOfTransform.erase(g3dListOfTransform.begin(), g3dListOfTransform.end());
}

//-------------------------------------------------------------------------

inline void radTg3d::B_genComp(radTField* FieldPtr)
{
	radTFieldKey& FieldKey = FieldPtr->FieldKey;
	if(g3dListOfTransform.empty()) 
	{
		if(FieldKey.Ib_ || FieldKey.Ih_) B_intComp(FieldPtr);
		else if(FieldKey.Force_) IntOverShape(FieldPtr);
		else B_comp(FieldPtr);
	}
	else
	{
		if(FieldKey.Force_) NestedFor_IntOverShape(FieldPtr, g3dListOfTransform.begin());
		else NestedFor_B(FieldPtr, g3dListOfTransform.begin());
	}
}

//-------------------------------------------------------------------------

inline void radTg3d::B_comp_Or_NestedFor(radTField* FieldPtr, const radTlphg::iterator& Iter)
{
	if(Iter == g3dListOfTransform.end())
	{
		if(FieldPtr->FieldKey.Ib_ || FieldPtr->FieldKey.Ih_) B_intComp(FieldPtr);
		else B_comp(FieldPtr);
	}
	else NestedFor_B(FieldPtr, Iter);
}

//-------------------------------------------------------------------------

inline double radTg3d::TransAtans(double x, double y, double& PiMult) 
{// To optimally compute sums of atans in derived classes
	double Buf = 1.-x*y;

	if(Buf == 0.) Buf = 1.e-50; //OC040504

	PiMult = (((Buf > 0)? 0.:1.) * ((x < 0)? -1.:1.));
	return (x+y)/Buf;
}

//-------------------------------------------------------------------------
/** Calculates principal value of argument of a complex number (-Pi < Phi <= Pi)  
	@param [out] x real part
	@param [out] y imaginary part
 	@return	calculated argument value
 	@see		... */
inline double radTg3d::Argument(double x, double y)
{

	const double Pi = 3.1415926535897932;
	if(x == 0)
	{
		if(y < 0) return -0.5*Pi;
		else if(y == 0) return 0;
		else return 0.5*Pi;
	}
	if(y == 0)
	{
		if(x >= 0) return 0.;
		else return Pi;
	}
	if(y < 0)
	{
		if(x < 0) return -Pi + atan(y/x);
		else return atan(y/x); // x > 0
	}
	else // y > 0
	{
		if(x < 0) return Pi + atan(y/x);
		else return atan(y/x); // x > 0
	}
}

//-------------------------------------------------------------------------

inline double radTg3d::AngularDifference(double Phi1, double Phi2)
{// This assumes Phi1 and Phi2 between 0 and TwoPI
	if(Phi1 > Phi2) return Phi1 - Phi2;
	else return Phi1 - (Phi2 - 6.28318530717959);
}

//-------------------------------------------------------------------------

inline char radTg3d::AngleIsBetween(double Phi, double PhiSt, double PhiFi)
{// This assumes Phi, PhiSt, PhiFi between 0 and TwoPI
	if(PhiSt > PhiFi)
	{
		if((Phi > PhiSt) && (Phi < 6.28318530717959)) return 1;
		if(Phi < PhiFi) return 1;
		return 0;
	}
	else
	{
		if((Phi > PhiSt) && (Phi < PhiFi)) return 1;
		return 0;
	}
}

//-------------------------------------------------------------------------

#include "radtrans.h" //to allow making subsequent inline

//-------------------------------------------------------------------------

inline void radTg3d::GetTrfAndCenPointInLabFrame(radTrans* pInBaseTrans, radTrans& bufTrans, radTrans*& pOutResTrans, TVector3d& outCenPointInLabFr)
{//assumes that pInBaseTrans and pResTrans were allocated by calling function(s)
 //in any case, it doesn't (re-)allocate pOutResTrans
	pOutResTrans = 0;
	radTrans* pTrans = (g3dListOfTransform.empty())? 0 : (radTrans*)((*(g3dListOfTransform.begin())).Handler_g.rep);
	if(pTrans != 0)
	{
		if(pInBaseTrans != 0) 
		{
			TrProduct(pInBaseTrans, pTrans, bufTrans);
			pOutResTrans = &bufTrans;
		}
		else //OC04082010
		{
			pOutResTrans = pTrans; //?
		}
	}
	else
	{
		if(pInBaseTrans != 0) pOutResTrans = pInBaseTrans;
	}
	outCenPointInLabFr = CentrPoint;
	if(pOutResTrans != 0) outCenPointInLabFr = pOutResTrans->TrPoint(outCenPointInLabFr);
}

//-------------------------------------------------------------------------

inline void radTAuxCompDataG3D::StoreDataFromField(radTField* FieldPtr)
{
	radTFieldKey &FieldKey = FieldPtr->FieldKey;
	if(FieldKey.Energy_) Energy = FieldPtr->Energy;
	if(FieldKey.ForceEnr_) Force = FieldPtr->Force;
	if(FieldKey.Torque_) Torque = FieldPtr->Torque;
}

//-------------------------------------------------------------------------

inline void radTAuxCompDataG3D::PutStoredDataToField(radTField* FieldPtr)
{
	radTFieldKey &FieldKey = FieldPtr->FieldKey;
	if(FieldKey.Energy_) FieldPtr->Energy = Energy;
	if(FieldKey.ForceEnr_) FieldPtr->Force = Force;
	if(FieldKey.Torque_) FieldPtr->Torque = Torque;
}

//-------------------------------------------------------------------------

#endif
