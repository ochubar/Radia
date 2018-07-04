/*-------------------------------------------------------------------------
*
* File name:      radarccu.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 rectangular cross-section arc with azimuthal current
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADARCCU_H
#define __RADARCCU_H

#include "radg3d.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTArcCur : public radTg3d {
public:
	//TVector3d CentrPoint; //moved to radTg3d
	TVector3d CircleCentrPoint;
	double R_min, R_max;
	double Phi_min, Phi_max;
	double Height;
	double J_azim;
	int NumberOfSectors;

	short BasedOnPrecLevel;

	short InternalFacesAfterCut;
	char J_IsNotZero;

	radTArcCur(const TVector3d& InCPoiVect, const double* InRadiiPtr, 
			   const double* InAnglesPtr, double InHeight, 
			   double InJ_azim, int InNumberOfSectors =0, short InBasedOnPrecLevel =0)
	{
		CircleCentrPoint = InCPoiVect;
		R_min = *InRadiiPtr; R_max = *(InRadiiPtr+1);
		Phi_min = *InAnglesPtr; Phi_max = *(InAnglesPtr+1);
		Height = InHeight;
		J_azim = InJ_azim;
		NumberOfSectors = InNumberOfSectors;
		BasedOnPrecLevel = InBasedOnPrecLevel;

		ComputeCentrPoint();

		if(J_azim != 0.) J_IsNotZero = 1;
		InternalFacesAfterCut = 0;
		ConsiderOnlyWithTrans = 0;
	}
	radTArcCur(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{//Instantiates from string according to DumpBin
		DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);

		//TVector3d CircleCentrPoint;
		inStr >> CircleCentrPoint;

		//double R_min, R_max;
		inStr >> R_min;
		inStr >> R_max;

		//double Phi_min, Phi_max;
		inStr >> Phi_min;
		inStr >> Phi_max;

		//double Height;
		inStr >> Height;

		//double J_azim;
		inStr >> J_azim;

		//int NumberOfSectors;
		inStr >> NumberOfSectors;

		//short BasedOnPrecLevel;
		inStr >> BasedOnPrecLevel;

		//short InternalFacesAfterCut;
		inStr >> InternalFacesAfterCut;

		//char J_IsNotZero;
		inStr >> J_IsNotZero;
	}

	radTArcCur() 
	{ 
		J_IsNotZero = 0; InternalFacesAfterCut = 0; ConsiderOnlyWithTrans = 0;
	}

	int Type_g3d() { return 3;}

	void ComputeCentrPoint()
	{
		double Phic = 0.5*(Phi_min + Phi_max);
		double rc = (2./3.)*(R_max*R_max + R_max*R_min + R_min*R_min)/(R_max + R_min);
		CentrPoint = TVector3d(rc*cos(Phic) + CircleCentrPoint.x, rc*sin(Phic) + CircleCentrPoint.y, CircleCentrPoint.z);
	}

	void B_comp(radTField* FieldPtr)
	{
		if(FieldPtr->FieldKey.J_) J_comp(FieldPtr);
		if(!(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_ || FieldPtr->FieldKey.A_ || FieldPtr->FieldKey.PreRelax_)) return;

		if(BasedOnPrecLevel) B_compWithNewtonCotes4(FieldPtr);
		else B_compWithTrapeth(FieldPtr);
	}
	void B_compWithNewtonCotes4(radTField*);
	void B_compWithTrapeth(radTField*);

	void J_comp(radTField* FieldPtr)
	{
		const double twoPi = 2.*3.141592653589793238;
		const int maxNumTwoPi = 10;

		TVector3d vR = FieldPtr->P - CircleCentrPoint;
		double r = vR.Abs();
		if((r < R_min) || (r > R_max)) return; //point is outside
		if(::fabs(vR.z) > 0.5*Height) return; //point is outside

		double phi = Argument(vR.x, vR.y); // -Pi < phi <= Pi
		if(phi < Phi_min)
		{
			for(int i=0; i<maxNumTwoPi; i++)
			{
				phi += twoPi; if(phi >= Phi_min) break;
			}
		}
		else if(phi > Phi_max)
		{
			for(int i=0; i<maxNumTwoPi; i++)
			{
				phi -= twoPi; if(phi <= Phi_max) break;
			}
		}
		if((phi < Phi_min) || (phi > Phi_max)) return; //point is outside

		TVector3d locJ(-J_azim*sin(phi), J_azim*cos(phi), 0);
		FieldPtr->J += locJ;
	}

	void B_intComp(radTField* FieldPtr)
	{
		if(FieldPtr->FieldKey.FinInt_) B_intCompFinNum(FieldPtr);
		else if(BasedOnPrecLevel) B_intCompWithNewton3(FieldPtr);
		else B_intCompWithTrapeth(FieldPtr);
	}
	void B_intCompWithNewton3(radTField*);
	void B_intCompWithTrapeth(radTField*);
	void B_intUtilSpecCaseZeroVxVy(double, double, double, double, double, double, double, double&);

	double Volume() { return 0.5*(Phi_max - Phi_min)*(R_max*R_max - R_min*R_min)*Height;}
	void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique);

	void SimpleEnergyComp(radTField* FieldPtr)
	{
		const double ConstForJ = -1.E-06;
		radTFieldKey LocFieldKey; LocFieldKey.A_ = 1;
		radTField LocField(LocFieldKey, FieldPtr->CompCriterium); LocField.P = CentrPoint;
		((radTg3d*)(FieldPtr->HandleEnergyForceTorqueCompData.rep->hSource.rep))->B_genComp(&LocField);
		double Phic = 0.5*(Phi_min + Phi_max);
		TVector3d J(-J_azim*sin(Phic), J_azim*cos(Phic), 0.);
		FieldPtr->Energy += (ConstForJ*Volume())*(J*LocField.A);
	}

	void Push_backCenterPointAndField(radTFieldKey*, radTVectPairOfVect3d*, radTrans* pBaseTrans, radTg3d*, radTApplication*);

	void Dump(std::ostream&, int ShortSign =0);
	void DumpPureObjInfo(std::ostream&, int);
	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey);

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
		return FinishDuplication(new radTArcCur(*this), hg);
	}
	int SizeOfThis() { return sizeof(radTArcCur);}

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
	void ListFacesInternalAfterCut(short* FacesState)
	{
		short BufNum = InternalFacesAfterCut;
		for(int k=0; k<6; k++) { *(FacesState++) = BufNum & 1; BufNum >>= 1;}
	}
	void SetFacesInternalAfterCut(short* FacesState)
	{
		for(int k=0; k<6; k++) if(*(FacesState++)) MapFaceAsInternalAfterCut(k+1);
	}

	int ScaleCurrent(double scaleCoef) //virtual in g3d
	{
		J_azim *= scaleCoef; 
		return 1;
	}
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTBackgroundFieldSource : public radTg3d {
public:
	TVector3d BackgrB;

	radTBackgroundFieldSource(const TVector3d& InBackgrB)
	{
		BackgrB = InBackgrB;
	}

	radTBackgroundFieldSource(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{//Instantiates from string according to DumpBin
		DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);

		//TVector3d BackgrB;
		inStr >> BackgrB;
	}

	radTBackgroundFieldSource()
	{
		BackgrB.Zero();
	}

	void B_comp(radTField* FieldPtr)
	{
		if(FieldPtr->FieldKey.B_) FieldPtr->B += BackgrB;
		if(FieldPtr->FieldKey.H_) FieldPtr->H += BackgrB;
		if(FieldPtr->FieldKey.A_) 
		{
			TVector3d BufA(FieldPtr->P.z*BackgrB.y, FieldPtr->P.x*BackgrB.z, FieldPtr->P.y*BackgrB.x);
			FieldPtr->A += BufA;
		}
	}
	void B_intComp(radTField* FieldPtr)
	{
		if(FieldPtr->FieldKey.FinInt_)
		{
			TVector3d D = FieldPtr->NextP - FieldPtr->P;
			TVector3d BufIb = sqrt(D.x*D.x + D.y*D.y + D.z*D.z) * BackgrB;
			if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIb;
			if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIb;
		}
		// Infinite integral is set to zero (though, formally, its contribution is infinite)
	}
	radTg3dGraphPresent* CreateGraphPresent();

	void Dump(std::ostream& o, int ShortSign =0) // Porting
	{
		radTg3d::Dump(o);
		o << "Uniform background field source";
		if(ShortSign==1) return;
		o << endl;
		o << "   {bx,by,bz}= {" << BackgrB.x << ',' << BackgrB.y << ',' << BackgrB.z << "}";
		o << endl;
		o << "   Memory occupied: " << SizeOfThis() << " bytes";
	}

	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
	{
		//Dumping objects that may be used by this object
		vector<pair<int, int> > vTrfKeys;
		DumpBin_g3d_TreatTrfs(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vTrfKeys);

		vElemKeysOut.push_back(elemKey);
		oStr << elemKey;

		//Next 5 bytes define/encode element type:
		oStr << (char)Type_g();
		oStr << (char)Type_g3d();
		oStr << (char)0;
		oStr << (char)0;
		oStr << (char)0;

		//Members of radTg3d
		DumpBin_g3d(oStr, vTrfKeys);

		//TVector3d BackgrB;
		oStr << BackgrB;
	}

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
		return FinishDuplication(new radTBackgroundFieldSource(*this), hg);
	}
	int SizeOfThis() { return sizeof(radTBackgroundFieldSource);}
};

//-------------------------------------------------------------------------

#endif
