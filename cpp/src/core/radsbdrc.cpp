/*-------------------------------------------------------------------------
*
* File name:      radsbdrc.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 subdivided rectangular parallelepiped with constant magnetization
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radsbdrc.h"
#include "radg3dgr.h"
#include "radmamet.h"
#include "radcast.h"
#include "radappl.h"

#include <math.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTSubdividedRecMag::radTSubdividedRecMag(const TVector3d& InCPoiVect, const TVector3d& InDimsVect, const TVector3d& InMagnVect, 
										   const radThg& InMaterHandle, const double* InSubdivArray) 
{// This seems not to be used
	radTRecMag::CentrPoint = radTGroup::CentrPoint = InCPoiVect; //OC101008
	Dimensions=InDimsVect;
	Magn=InMagnVect;
	MaterHandle = InMaterHandle;

	double Small = 1.E-09;
	kx = int(InSubdivArray[0]+Small); ky = int(InSubdivArray[2]+Small); kz = int(InSubdivArray[4]+Small);
	qx = InSubdivArray[1]; qy = InSubdivArray[3]; qz = InSubdivArray[5];

	InternalFacesAfterCut = 0;

	FldCmpMeth = 0;
	Q_forM = NULL;
	FormCenPoPtrArray = NULL;
	AmOfSubElem = int(kx)*int(ky)*int(kz);
	CenPoLoopCounter = IntrcMatrConstrCounter = FormIntrctMembCounter = -1;

	AlgsBasedOnKsQsMayNotWork = false;
}

//-------------------------------------------------------------------------

radTSubdividedRecMag::radTSubdividedRecMag(const TVector3d& InCPoiVect, const TVector3d& InDimsVect, const double* InSubdivArray)
{// This seems not to be used
	radTRecMag::CentrPoint = radTGroup::CentrPoint = InCPoiVect; //OC101008
	Dimensions=InDimsVect;

	double Small = 1.E-09;
	kx = int(InSubdivArray[0]+Small); ky = int(InSubdivArray[2]+Small); kz = int(InSubdivArray[4]+Small);
	qx = InSubdivArray[1]; qy = InSubdivArray[3]; qz = InSubdivArray[5];

	InternalFacesAfterCut = 0;

	FldCmpMeth = 0;
	Q_forM = NULL;
	FormCenPoPtrArray = NULL;
	AmOfSubElem = int(kx)*int(ky)*int(kz);
	CenPoLoopCounter = IntrcMatrConstrCounter = FormIntrctMembCounter = -1;

	AlgsBasedOnKsQsMayNotWork = false;
}

//-------------------------------------------------------------------------

radTSubdividedRecMag::radTSubdividedRecMag(const radTRecMag* RecMagPtr, const double* InSubdivArray)
{
	radTRecMag::CentrPoint = radTGroup::CentrPoint = RecMagPtr->CentrPoint; //OC101008
	Dimensions = RecMagPtr->Dimensions;
	Magn = RecMagPtr->Magn;
	MaterHandle = RecMagPtr->MaterHandle;

	InternalFacesAfterCut = RecMagPtr->InternalFacesAfterCut;
	J_IsNotZero = RecMagPtr->J_IsNotZero;

	radTGroup::IsGroupMember = RecMagPtr->IsGroupMember;
	radTGroup::g3dListOfTransform = RecMagPtr->g3dListOfTransform;
	radTGroup::HandleAuxCompData = RecMagPtr->HandleAuxCompData;
	radTGroup::ConsiderOnlyWithTrans = RecMagPtr->ConsiderOnlyWithTrans;

	radTGroup::MessageChar = RecMagPtr->MessageChar;

	double Small = 1.E-09;
	kx = int(InSubdivArray[0]+Small); ky = int(InSubdivArray[2]+Small); kz = int(InSubdivArray[4]+Small);
	qx = InSubdivArray[1]; qy = InSubdivArray[3]; qz = InSubdivArray[5];

	FldCmpMeth = 0;
	Q_forM = NULL;
	FormCenPoPtrArray = NULL;
	AmOfSubElem = int(kx)*int(ky)*int(kz);
	CenPoLoopCounter = IntrcMatrConstrCounter = FormIntrctMembCounter = -1;

	AlgsBasedOnKsQsMayNotWork = false;
}

//-------------------------------------------------------------------------

radTSubdividedRecMag::radTSubdividedRecMag(const radTRecMag* RecMagPtr)
{
	radTRecMag::CentrPoint = radTGroup::CentrPoint = RecMagPtr->CentrPoint; //OC101008
	Dimensions = RecMagPtr->Dimensions;
	Magn = RecMagPtr->Magn;
	MaterHandle = RecMagPtr->MaterHandle;

	InternalFacesAfterCut = RecMagPtr->InternalFacesAfterCut;
	J_IsNotZero = RecMagPtr->J_IsNotZero;

	radTGroup::IsGroupMember = RecMagPtr->IsGroupMember;
	radTGroup::g3dListOfTransform = RecMagPtr->g3dListOfTransform;
	radTGroup::HandleAuxCompData = RecMagPtr->HandleAuxCompData;
	radTGroup::ConsiderOnlyWithTrans = RecMagPtr->ConsiderOnlyWithTrans;

	radTGroup::MessageChar = RecMagPtr->MessageChar;
	AlgsBasedOnKsQsMayNotWork = true;

	FldCmpMeth = 0;
	Q_forM = NULL;
	FormCenPoPtrArray = NULL;
	CenPoLoopCounter = IntrcMatrConstrCounter = FormIntrctMembCounter = -1;
}

//-------------------------------------------------------------------------

radTSubdividedRecMag::radTSubdividedRecMag() 
{
	InternalFacesAfterCut = 0;

	FldCmpMeth = 0;
	Q_forM = NULL;
	FormCenPoPtrArray = NULL;

	CenPoLoopCounter = IntrcMatrConstrCounter = FormIntrctMembCounter = -1;
	AlgsBasedOnKsQsMayNotWork = true;
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTSubdividedRecMag::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = 0;
	if(!AlgsBasedOnKsQsMayNotWork) g3dGraphPresentPtr = new radTSubdivRecMagGraphPresent((radTGroup*)this);
	else g3dGraphPresentPtr = new radTGroupGraphPresent((radTGroup*)this);
	g3dGraphPresentPtr->ShowInternalFacesAfterCut = false;
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------

void radTSubdividedRecMag::B_comp(radTField* FieldPtr)
{
	if(FldCmpMeth==0 || AlgsBasedOnKsQsMayNotWork) { radTGroup::B_comp(FieldPtr); return;}
	else if(FldCmpMeth==1)
	{
		if(FieldPtr->FieldKey.PreRelax_)
		{
			IntrcMatrConstrCounter++;
			int SubElemNo = int(double(IntrcMatrConstrCounter)/double(FieldPtr->AmOfIntrctElemWithSym));
			
			int IndxInCont = -1;
			if(IntrcMatrConstrCounter >= FieldPtr->AmOfIntrctElemWithSym)
				IndxInCont = IntrcMatrConstrCounter - SubElemNo*FieldPtr->AmOfIntrctElemWithSym;

			if(IndxInCont==-1)
			{
				TVector3d* NewArrayOfT = new TVector3d[AmOfSubElem]; // Wrap with try-catch?
				TVector3d* NewArrayOfS = new TVector3d[AmOfSubElem];
				VectOfTsComputed.push_back(NewArrayOfT);
				VectOfSsComputed.push_back(NewArrayOfS);

				B_compPolynomial(FieldPtr); // This fills-in last arrays in VectOfTs(Ss)Computed
				IndxInCont = IntrcMatrConstrCounter;
			}
			
			TVector3d& RefB = FieldPtr->B;
			TVector3d& RefH = FieldPtr->H;
			TVector3d& RefA = FieldPtr->A;
			TVector3d& T = (VectOfTsComputed[IndxInCont])[SubElemNo];
			TVector3d& S = (VectOfSsComputed[IndxInCont])[SubElemNo];
			RefB.x = T.x; RefB.y = -S.z; RefB.z = -S.y;
			RefH.x = -S.z; RefH.y = T.y; RefH.z = -S.x;
			RefA.x = -S.y; RefA.y = -S.x; RefA.z = T.z;

			if(IntrcMatrConstrCounter == (FieldPtr->AmOfIntrctElemWithSym)*AmOfSubElem-1)
			{
				IntrcMatrConstrCounter = -1;
				for(int i=0; i<(FieldPtr->AmOfIntrctElemWithSym); i++)
				{
					delete[] (VectOfTsComputed[i]);
					delete[] (VectOfSsComputed[i]);
				}
				VectOfTsComputed.erase(VectOfTsComputed.begin(), VectOfTsComputed.end());
				VectOfSsComputed.erase(VectOfSsComputed.begin(), VectOfSsComputed.end());
			}
			return;
		}
		B_compPolynomial(FieldPtr);
	}
}

//-------------------------------------------------------------------------

void radTSubdividedRecMag::B_compPolynomial(radTField* FieldPtr)
{
	TVector3d CenPo_mi_P = radTRecMag::CentrPoint - FieldPtr->P;
	TVector3d HalfDim = 0.5*Dimensions;
	const double Pi = 3.141592653589793238;
	const double ConForH = 1./4./Pi;

	double x1 = CenPo_mi_P.x - HalfDim.x;
	double x2 = CenPo_mi_P.x + HalfDim.x;
	double y1 = CenPo_mi_P.y - HalfDim.y;
	double y2 = CenPo_mi_P.y + HalfDim.y;
	double z1 = CenPo_mi_P.z - HalfDim.z;
	double z2 = CenPo_mi_P.z + HalfDim.z;
	// Check & repair trapping to borders!!!

	TVector3d T[3][3][3], S[3][3][3];  // Make larger or create in heap later

	double InvDelX = 1.;
	double InvDelY = 1.;
	double InvDelZ = 1.;

	short FieldCompNeeded = FieldPtr->FieldKey.PreRelax_ || FieldPtr->FieldKey.H_ || FieldPtr->FieldKey.B_;
	short B_orH_CompNeeded = FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_;
	short InsideBlock = (x1*x2<0) && (y1*y2<0) && (z1*z2<0);
	short MagnCompNeeded = (FieldPtr->FieldKey.M_ || FieldPtr->FieldKey.B_) && InsideBlock;

	if(FieldCompNeeded)
	{
		double x1x1 = x1*x1;
		double x2x2 = x2*x2;
		double y1y1 = y1*y1;
		double y2y2 = y2*y2;
		double z1z1 = z1*z1;
		double z2z2 = z2*z2;
		double R111 = sqrt(x1x1+y1y1+z1z1);
		double R211 = sqrt(x2x2+y1y1+z1z1);
		double R121 = sqrt(x1x1+y2y2+z1z1);
		double R221 = sqrt(x2x2+y2y2+z1z1);
		double R112 = sqrt(x1x1+y1y1+z2z2);
		double R212 = sqrt(x2x2+y1y1+z2z2);
		double R122 = sqrt(x1x1+y2y2+z2z2);
		double R222 = sqrt(x2x2+y2y2+z2z2);

		double PiMult1, PiMult2, PiMult3;
		PiMult1 = PiMult2 = PiMult3 = 0.;

		double SumAtansX1 = atan(radTGroup::TransAtans(radTGroup::TransAtans(y1*z1/x1/R111, -y1*z2/x1/R112, PiMult1), radTGroup::TransAtans(-y2*z1/x1/R121, y2*z2/x1/R122, PiMult2), PiMult3)) + Pi*(PiMult1+PiMult2+PiMult3);
		double SumAtansX2 = atan(radTGroup::TransAtans(radTGroup::TransAtans(y1*z1/x2/R211, -y1*z2/x2/R212, PiMult1), radTGroup::TransAtans(-y2*z1/x2/R221, y2*z2/x2/R222, PiMult2), PiMult3)) + Pi*(PiMult1+PiMult2+PiMult3);
		double SumAtansY1 = atan(radTGroup::TransAtans(radTGroup::TransAtans(x1*z1/y1/R111, -x1*z2/y1/R112, PiMult1), radTGroup::TransAtans(-x2*z1/y1/R211, x2*z2/y1/R212, PiMult2), PiMult3)) + Pi*(PiMult1+PiMult2+PiMult3);
		double SumAtansY2 = atan(radTGroup::TransAtans(radTGroup::TransAtans(x1*z1/y2/R121, -x1*z2/y2/R122, PiMult1), radTGroup::TransAtans(-x2*z1/y2/R221, x2*z2/y2/R222, PiMult2), PiMult3)) + Pi*(PiMult1+PiMult2+PiMult3);
		double SumAtansZ1 = atan(radTGroup::TransAtans(radTGroup::TransAtans(x1*y1/z1/R111, -x1*y2/z1/R121, PiMult1), radTGroup::TransAtans(-x2*y1/z1/R211, x2*y2/z1/R221, PiMult2), PiMult3)) + Pi*(PiMult1+PiMult2+PiMult3);
		double SumAtansZ2 = atan(radTGroup::TransAtans(radTGroup::TransAtans(x1*y1/z2/R112, -x1*y2/z2/R122, PiMult1), radTGroup::TransAtans(-x2*y1/z2/R212, x2*y2/z2/R222, PiMult2), PiMult3)) + Pi*(PiMult1+PiMult2+PiMult3);
		double SumLogsX1Z1 = log((y1+R111)/(y2+R121));
		double SumLogsX1Z2 = log((y1+R112)/(y2+R122));
		double SumLogsX2Z1 = log((y1+R211)/(y2+R221));
		double SumLogsX2Z2 = log((y1+R212)/(y2+R222));
		double SumLogsX1Y1 = log((z1+R111)/(z2+R112));
		double SumLogsX1Y2 = log((z1+R121)/(z2+R122));
		double SumLogsX2Y1 = log((z1+R211)/(z2+R212));
		double SumLogsX2Y2 = log((z1+R221)/(z2+R222));
		double SumLogsY1Z1 = log((x1+R111)/(x2+R211));
		double SumLogsY1Z2 = log((x1+R112)/(x2+R212));
		double SumLogsY2Z1 = log((x1+R121)/(x2+R221));
		double SumLogsY2Z2 = log((x1+R122)/(x2+R222));

		double InvDelZe2, InvDelYInvDelZe2, InvDelXInvDelZe2, InvDelXInvDelYInvDelZ, InvDelXInvDelYInvDelZe2,
			   InvDelXInvDelY, InvDelXInvDelZ, InvDelYInvDelZ, InvDelYe2, InvDelYe2InvDelZ, InvDelYe2InvDelZe2, 
			   InvDelXInvDelYe2, InvDelXInvDelYe2InvDelZ, InvDelXInvDelYe2InvDelZe2, InvDelXe2, 
			   InvDelXe2InvDelZ, InvDelXe2InvDelZe2, InvDelXe2InvDelY, InvDelXe2InvDelYInvDelZ, 
			   InvDelXe2InvDelYInvDelZe2, InvDelXe2InvDelYe2, InvDelXe2InvDelYe2InvDelZ, InvDelXe2InvDelYe2InvDelZe2, 
			   x1e2=x1x1, x2e2=x2x2, y1e2=y1y1, y2e2=y2y2, z1e2=z1z1, z2e2=z2z2,
			   x1e3, x2e3, y1e3, y2e3, z1e3, z2e3,
			   x1e4, x2e4, y1e4, y2e4, z1e4, z2e4,
			   x1e5, x2e5, y1e5, y2e5, z1e5, z2e5,
			   x1e6, x2e6, y1e6, y2e6, z1e6, z2e6,
			   x1e7, x2e7, y1e7, y2e7, z1e7, z2e7,
			   x1e8, x2e8, y1e8, y2e8, z1e8, z2e8,
			   x1e9, x2e9, y1e9, y2e9, z1e9, z2e9;

		int MaxOrd = int(kx) + int(ky) + int(kz) - 3;
		if(MaxOrd>2) { x1e3=x1e2*x1; x2e3=x2e2*x2; y1e3=y1e2*y1; y2e3=y2e2*y2; z1e3=z1e2*z1; z2e3=z2e2*z2;}
		if(MaxOrd>3) { x1e4=x1e3*x1; x2e4=x2e3*x2; y1e4=y1e3*y1; y2e4=y2e3*y2; z1e4=z1e3*z1; z2e4=z2e3*z2;}
		if(MaxOrd>4) { x1e5=x1e4*x1; x2e5=x2e4*x2; y1e5=y1e4*y1; y2e5=y2e4*y2; z1e5=z1e4*z1; z2e5=z2e4*z2;}
		if(MaxOrd>5) { x1e6=x1e5*x1; x2e6=x2e5*x2; y1e6=y1e5*y1; y2e6=y2e5*y2; z1e6=z1e5*z1; z2e6=z2e5*z2;}
		if(MaxOrd>6) { x1e7=x1e6*x1; x2e7=x2e6*x2; y1e7=y1e6*y1; y2e7=y2e6*y2; z1e7=z1e6*z1; z2e7=z2e6*z2;}
		if(MaxOrd>7) { x1e8=x1e7*x1; x2e8=x2e7*x2; y1e8=y1e7*y1; y2e8=y2e7*y2; z1e8=z1e7*z1; z2e8=z2e7*z2;}
		if(MaxOrd>8) { x1e9=x1e8*x1; x2e9=x2e8*x2; y1e9=y1e8*y1; y2e9=y2e8*y2; z1e9=z1e8*z1; z2e9=z2e8*z2;}

		if(int(kx)>0)
		{
			if(int(ky)>0)
			{
				if(int(kz)>0)
				{
					TVector3d& T000 = T[0][0][0];
					TVector3d& S000 = S[0][0][0];

					T000.x = SumAtansX1-SumAtansX2;
					S000.x = -SumLogsX1Y1+SumLogsX1Y2+SumLogsX2Y1-SumLogsX2Y2;
					S000.y = -SumLogsX1Z1+SumLogsX1Z2+SumLogsX2Z1-SumLogsX2Z2;
					T000.y = SumAtansY1-SumAtansY2;
					S000.z = -SumLogsY1Z1+SumLogsY1Z2+SumLogsY2Z1-SumLogsY2Z2;
					T000.z = SumAtansZ1-SumAtansZ2;
				}
				if(int(kz)>1)
				{
					TVector3d& T001 = T[0][0][1];
					TVector3d& S001 = S[0][0][1];

					double CmnPrt1 = x1*(SumLogsX1Z1-SumLogsX1Z2) - x2*(SumLogsX2Z1-SumLogsX2Z2);
					double CmnPrt2 = y1*(SumLogsY1Z1-SumLogsY1Z2) - y2*(SumLogsY2Z1-SumLogsY2Z2);
					T001.x = -CmnPrt1;
					S001.x = -R111+R121+R211-R221+R112-R122-R212+R222;
					S001.y = -x1*SumAtansX1 + x2*SumAtansX2 + y1*(SumLogsX1Y1-SumLogsX2Y1) - y2*(SumLogsX1Y2-SumLogsX2Y2);
					T001.y = -CmnPrt2;
					S001.z = -y1*SumAtansY1 + y2*SumAtansY2 + x1*(SumLogsX1Y1-SumLogsX1Y2) - x2*(SumLogsX2Y1-SumLogsX2Y2);
					T001.z = CmnPrt2 + CmnPrt1;

					T001 = InvDelZ*T001; S001 = InvDelZ*S001;
				}
				if(int(kz)>2)
				{
					TVector3d& T002 = T[0][0][2];
					TVector3d& S002 = S[0][0][2];

					double CmnPrt1 = x1*(SumLogsX1Y1*y1 - SumLogsX1Y2*y2) + x2*(-SumLogsX2Y1*y1 + SumLogsX2Y2*y2);
					T002.x = -SumAtansX1*x1x1 + SumAtansX2*x2x2 + CmnPrt1;
					S002.x = 0.5*((SumLogsX1Y1 - SumLogsX1Y2)*x1x1 + (-SumLogsX2Y1 + SumLogsX2Y2)*x2x2 + (SumLogsX1Y1 - SumLogsX2Y1)*y1y1 + (-SumLogsX1Y2 + SumLogsX2Y2)*y2y2 + (-R111 + R121 + R211 - R221)*z1 + (R112 - R122 - R212 + R222)*z2);
					S002.y = (SumLogsX1Z1 - SumLogsX1Z2)*x1x1 + (-SumLogsX2Z1 + SumLogsX2Z2)*x2x2 + (R111 - R112 - R211 + R212)*y1 + (-R121 + R122 + R221 - R222)*y2;
					T002.y = -SumAtansY1*y1y1 + SumAtansY2*y2y2 + CmnPrt1;
					S002.z = (SumLogsY1Z1 - SumLogsY1Z2)*y1y1 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2y2 + (R111 - R112 - R121 + R122)*x1 + (-R211 + R212 + R221 - R222)*x2;
					T002.z = SumAtansX1*x1x1 - SumAtansX2*x2x2 + SumAtansY1*y1y1 - SumAtansY2*y2y2 - 2.*CmnPrt1;

					InvDelZe2 = InvDelZ*InvDelZ;
					T002 = InvDelZe2*T002; S002 = InvDelZe2*S002;
				}
			}
			if(int(ky)>1)
			{
				if(int(kz)>0)
				{
					TVector3d& T010 = T[0][1][0];
					TVector3d& S010 = S[0][1][0];

					double CmnPrt1 = x1*(SumLogsX1Y1-SumLogsX1Y2) - x2*(SumLogsX2Y1-SumLogsX2Y2);
					double CmnPrt2 = z1*(SumLogsY1Z1-SumLogsY2Z1) - z2*(SumLogsY1Z2-SumLogsY2Z2);
					T010.x = -CmnPrt1;
					S010.x = -x1*SumAtansX1 + x2*SumAtansX2 + z1*(SumLogsX1Z1-SumLogsX2Z1) - z2*(SumLogsX1Z2-SumLogsX2Z2);
					S010.y = -R111+R121+R211-R221+R112-R122-R212+R222;
					T010.y =  CmnPrt2 + CmnPrt1;
					S010.z = -z1*SumAtansZ1 + z2*SumAtansZ2 + x1*(SumLogsX1Z1-SumLogsX1Z2) - x2*(SumLogsX2Z1-SumLogsX2Z2);
					T010.z = -CmnPrt2;

					T010 = InvDelY*T010; S010 = InvDelY*S010;
				}
				if(int(kz)>1)
				{
					TVector3d& T011 = T[0][1][1];
					TVector3d& S011 = S[0][1][1];

					double CmnPrt1 = 0.5*(x1*(R111-R121-R112+R122) - x2*(R211-R221-R212+R222));
					double CmnPrt2 = 0.5*(-(y1y1-z1z1)*SumLogsY1Z1 + (y1y1-z2z2)*SumLogsY1Z2 + (y2y2-z1z1)*SumLogsY2Z1 - (y2y2-z2z2)*SumLogsY2Z2);
					T011.x = -2.*CmnPrt1;
					S011.x = 0.5*(y1*(-R111+R211+R112-R212) - y2*(-R121+R221+R122-R222) + (x1x1+z1z1)*SumLogsX1Z1 - (x1x1+z2z2)*SumLogsX1Z2 - (x2x2+z1z1)*SumLogsX2Z1 + (x2x2+z2z2)*SumLogsX2Z2);
					S011.y = 0.5*(z1*(-R111+R211+R121-R221) - z2*(-R112+R212+R122-R222) + (x1x1+y1y1)*SumLogsX1Y1 - (x1x1+y2y2)*SumLogsX1Y2 - (x2x2+y1y1)*SumLogsX2Y1 + (x2x2+y2y2)*SumLogsX2Y2);
					T011.y = CmnPrt1 + CmnPrt2;
					S011.z = 0.5*(x1x1*SumAtansX1 - x2x2*SumAtansX2 - y1y1*SumAtansY1 + y2y2*SumAtansY2 - z1z1*SumAtansZ1 + z2z2*SumAtansZ2);
					T011.z = CmnPrt1 - CmnPrt2;

					InvDelYInvDelZ = InvDelY*InvDelZ;
					T011 = InvDelYInvDelZ*T011; S011 = InvDelYInvDelZ*S011;
				}
				if(int(kz)>2)
				{
					TVector3d& T012 = T[0][1][2];
					TVector3d& S012 = S[0][1][2];

					double CmnPrt1 = (SumLogsX1Y1 - SumLogsX1Y2)*x1e3 + (-SumLogsX2Y1 + SumLogsX2Y2)*x2e3;
					double CmnPrt2 = (R121-R111)*z1 + (R112-R122)*z2;
					double CmnPrt3 = (R211-R221)*z1 + (R222-R212)*z2;
					double CmnPrt4 = (R211-R111)*z1 + (R112-R212)*z2;
					double CmnPrt5 = (R121-R221)*z1 + (R222-R122)*z2;
					T012.x = 0.5*(CmnPrt1 + x1*(SumLogsX1Y1*y1y1 - SumLogsX1Y2*y2y2 + CmnPrt2) + x2*(-SumLogsX2Y1*y1y1 + SumLogsX2Y2*y2y2 + CmnPrt3));
					S012.x = (SumAtansX1*x1e3 - SumAtansX2*x2e3 + (SumLogsX1Y1 - SumLogsX2Y1)*y1e3 + (-SumLogsX1Y2 + SumLogsX2Y2)*y2e3 + (SumLogsX1Z1 - SumLogsX2Z1)*z1e3 + (-SumLogsX1Z2 + SumLogsX2Z2)*z2e3 + y1*CmnPrt4 + y2*CmnPrt5)/3.;
					S012.y = (2.*((R111 - R112 - R121 + R122)*x1x1 + (-R211 + R212 + R221 - R222)*x2x2 + (R111 - R112 - R211 + R212)*y1y1 + (-R121 + R122 + R221 - R222)*y2y2) + (-R111 + R121 + R211 - R221)*z1z1 + (R112 - R122 - R212 + R222)*z2z2)/3.;
					T012.y = (-CmnPrt1 + 4.*(-SumAtansY1*y1e3 + SumAtansY2*y2e3) + 2.*((SumLogsY1Z1 - SumLogsY2Z1)*z1e3 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e3) + x1*(3.*(SumLogsX1Y1*y1y1 - SumLogsX1Y2*y2y2) - CmnPrt2) + x2*(3.*(-SumLogsX2Y1*y1y1 + SumLogsX2Y2*y2y2) - CmnPrt3))/6.;
					S012.z = ((-SumLogsX1Z1 + SumLogsX1Z2)*x1e3 + (SumLogsX2Z1 - SumLogsX2Z2)*x2e3 + 2.*((SumLogsY1Z1 - SumLogsY1Z2)*y1e3 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2e3) + x1*((R111 - R112)*y1 + (-R121 + R122)*y2) + x2*((-R211 + R212)*y1 + (R221 - R222)*y2) - SumAtansZ1*z1e3 + SumAtansZ2*z2e3)/3.;
					T012.z = (-CmnPrt1 + 2.*(SumAtansY1*y1e3 - SumAtansY2*y2e3) + (-SumLogsY1Z1 + SumLogsY2Z1)*z1e3 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e3 + x1*(3.*(-SumLogsX1Y1*y1y1 + SumLogsX1Y2*y2y2) - CmnPrt2) + x2*(3.*(SumLogsX2Y1*y1y1 - SumLogsX2Y2*y2y2) - CmnPrt3))/3.;

					InvDelYInvDelZe2 = InvDelY*InvDelZe2;
					T012 = InvDelYInvDelZe2*T012; S012 = InvDelYInvDelZe2*S012;
				}
			}
			if(int(ky)>2)
			{
				if(int(kz)>0)
				{
					TVector3d& T020 = T[0][2][0];
					TVector3d& S020 = S[0][2][0];

					double CmnPrt1 = -SumAtansX1*x1e2 + SumAtansX2*x2e2;
					double CmnPrt2 = x1*(SumLogsX1Z1*z1 - SumLogsX1Z2*z2) + x2*(-SumLogsX2Z1*z1 + SumLogsX2Z2*z2);
					T020.x = CmnPrt1 + CmnPrt2;
					S020.x = (SumLogsX1Y1 - SumLogsX1Y2)*x1e2 + (-SumLogsX2Y1 + SumLogsX2Y2)*x2e2 + (R111 - R121 - R211 + R221)*z1 + (-R112 + R122 + R212 - R222)*z2;
					S020.y = ((SumLogsX1Z1 - SumLogsX1Z2)*x1e2 + (-SumLogsX2Z1 + SumLogsX2Z2)*x2e2 + (-R111 + R112 + R211 - R212)*y1 + (R121 - R122 - R221 + R222)*y2 + (SumLogsX1Z1 - SumLogsX2Z1)*z1e2 + (-SumLogsX1Z2 + SumLogsX2Z2)*z2e2)/2.;
					T020.y = -CmnPrt1 + SumAtansZ1*z1e2 - SumAtansZ2*z2e2 - 2.*CmnPrt2;
					S020.z = (R111 - R112 - R121 + R122)*x1 + (-R211 + R212 + R221 - R222)*x2 + (SumLogsY1Z1 - SumLogsY2Z1)*z1e2 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e2;
					T020.z = -SumAtansZ1*z1e2 + SumAtansZ2*z2e2 + CmnPrt2;

					InvDelYe2 = InvDelY*InvDelY;
					T020 = InvDelYe2*T020; S020 = InvDelYe2*S020;
				}
				if(int(kz)>1)
				{
					TVector3d& T021 = T[0][2][1];
					TVector3d& S021 = S[0][2][1];

					double CmnPrt1 = (-SumLogsX1Z1 + SumLogsX1Z2)*x1e3 + (SumLogsX2Z1 - SumLogsX2Z2)*x2e3;
					double CmnPrt2 = 3.*(SumLogsX2Z1*z1e2 - SumLogsX2Z2*z2e2);
					T021.x = (-CmnPrt1 + x1*((-R111 + R112)*y1 + (R121 - R122)*y2 + SumLogsX1Z1*z1e2 - SumLogsX1Z2*z2e2) + x2*((R211 - R212)*y1 + (-R221 + R222)*y2 - SumLogsX2Z1*z1e2 + SumLogsX2Z2*z2e2))/2.;
					S021.x = (2.*((R111 - R112 - R121 + R122)*x1e2 + (-R211 + R212 + R221 - R222)*x2e2 + (R111 - R121 - R211 + R221)*z1e2 + (-R112 + R122 + R212 - R222)*z2e2) + (-R111 + R112 + R211 - R212)*y1e2 + (R121 - R122 - R221 + R222)*y2e2)/3.;
					S021.y = (SumAtansX1*x1e3 - SumAtansX2*x2e3 + (SumLogsX1Y1 - SumLogsX2Y1)*y1e3 + (-SumLogsX1Y2 + SumLogsX2Y2)*y2e3 + (SumLogsX1Z1 - SumLogsX2Z1)*z1e3 + (-SumLogsX1Z2 + SumLogsX2Z2)*z2e3 + y1*((-R111 + R211)*z1 + (R112 - R212)*z2) + y2*((R121 - R221)*z1 + (-R122 + R222)*z2))/3.;
					T021.y = (CmnPrt1 + (-SumLogsY1Z1 + SumLogsY1Z2)*y1e3 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e3 + 2.*(SumAtansZ1*z1e3 - SumAtansZ2*z2e3) + x1*((R111 - R112)*y1 + (-R121 + R122)*y2 + 3.*(-SumLogsX1Z1*z1e2 + SumLogsX1Z2*z2e2)) + x2*((-R211 + R212)*y1 + (R221 - R222)*y2 + CmnPrt2))/3.;
					S021.z = ((-SumLogsX1Y1 + SumLogsX1Y2)*x1e3 + (SumLogsX2Y1 - SumLogsX2Y2)*x2e3 - SumAtansY1*y1e3 + SumAtansY2*y2e3 + 2.*((SumLogsY1Z1 - SumLogsY2Z1)*z1e3 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e3) + x1*((R111 - R121)*z1 + (-R112 + R122)*z2) + x2*((-R211 + R221)*z1 + (R212 - R222)*z2))/3.;
					T021.z = (CmnPrt1 + 2.*((SumLogsY1Z1 - SumLogsY1Z2)*y1e3 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2e3) + 4.*(-SumAtansZ1*z1e3 + SumAtansZ2*z2e3) + x1*((R111 - R112)*y1 + (-R121 + R122)*y2 + 3.*(SumLogsX1Z1*z1e2 - SumLogsX1Z2*z2e2)) + x2*((-R211 + R212)*y1 + (R221 - R222)*y2 - CmnPrt2))/6.;

					InvDelYe2InvDelZ = InvDelYe2*InvDelZ;
					T021 = InvDelYe2InvDelZ*T021; S021 = InvDelYe2InvDelZ*S021;
				}
				if(int(kz)>2)
				{
					TVector3d& T022 = T[0][2][2];
					TVector3d& S022 = S[0][2][2];

					double CmnPrt1 = -SumAtansX1*x1e4 + SumAtansX2*x2e4;
					double CmnPrt2 = y1*(R111*z1 - R112*z2) + y2*(-R121*z1 + R122*z2);
					double CmnPrt3 = SumLogsX2Y1*y1e3 - SumLogsX2Y2*y2e3;
					double CmnPrt4 = SumLogsX2Z1*z1e3 - SumLogsX2Z2*z2e3;
					double CmnPrt5 = y1*(-R211*z1 + R212*z2) + y2*(R221*z1 - R222*z2);
					T022.x = (-CmnPrt1 + x1*(SumLogsX1Y1*y1e3 - SumLogsX1Y2*y2e3 + SumLogsX1Z1*z1e3 - SumLogsX1Z2*z2e3 + y1*(-R111*z1 + R112*z2) + y2*(R121*z1 - R122*z2)) + x2*(-SumLogsX2Y1*y1e3 + SumLogsX2Y2*y2e3 - SumLogsX2Z1*z1e3 + SumLogsX2Z2*z2e3 + y1*(R211*z1 - R212*z2) + y2*(-R221*z1 + R222*z2)))/3.;
					S022.x = ((-SumLogsX1Y1 + SumLogsX1Y2)*x1e4 + (SumLogsX2Y1 - SumLogsX2Y2)*x2e4 + (SumLogsX1Y1 - SumLogsX2Y1)*y1e4 + (-SumLogsX1Y2 + SumLogsX2Y2)*y2e4 + 2.*((R111 - R121 - R211 + R221)*z1e3 + (-R112 + R122 + R212 - R222)*z2e3) + x1e2*((R111 - R121)*z1 + (-R112 + R122)*z2) + y1e2*((-R111 + R211)*z1 + (R112 - R212)*z2) + x2e2*((-R211 + R221)*z1 + (R212 - R222)*z2) + y2e2*((R121 - R221)*z1 + (-R122 + R222)*z2))/4.;
					S022.y = ((-SumLogsX1Z1 + SumLogsX1Z2)*x1e4 + (SumLogsX2Z1 - SumLogsX2Z2)*x2e4 + 2.*((R111 - R112 - R211 + R212)*y1e3 + (-R121 + R122 + R221 - R222)*y2e3) + x1e2*((R111 - R112)*y1 + (-R121 + R122)*y2) + x2e2*((-R211 + R212)*y1 + (R221 - R222)*y2) + (SumLogsX1Z1 - SumLogsX2Z1)*z1e4 + (-SumLogsX1Z2 + SumLogsX2Z2)*z2e4 + y1*((-R111 + R211)*z1e2 + (R112 - R212)*z2e2) + y2*((R121 - R221)*z1e2 + (-R122 + R222)*z2e2))/4.;
					T022.y = (CmnPrt1 + 3.*(-SumAtansY1*y1e4 + SumAtansY2*y2e4 + SumAtansZ1*z1e4 - SumAtansZ2*z2e4) + x1*(2.*(SumLogsX1Y1*y1e3 - SumLogsX1Y2*y2e3) + 4.*(-SumLogsX1Z1*z1e3 + SumLogsX1Z2*z2e3) + CmnPrt2) + x2*(-2.*CmnPrt3 + 4.*CmnPrt4 + CmnPrt5))/6.;
					S022.z = (2.*((-R111 + R112 + R121 - R122)*x1e3 + (R211 - R212 - R221 + R222)*x2e3) + 3.*((SumLogsY1Z1 - SumLogsY1Z2)*y1e4 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2e4 + (SumLogsY1Z1 - SumLogsY2Z1)*z1e4 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e4) + x1*((R111 - R112)*y1e2 + (-R121 + R122)*y2e2 + (R111 - R121)*z1e2 + (-R112 + R122)*z2e2) + x2*((-R211 + R212)*y1e2 + (R221 - R222)*y2e2 + (-R211 + R221)*z1e2 + (R212 - R222)*z2e2))/6.;
					T022.z = (CmnPrt1 + 3.*(SumAtansY1*y1e4 - SumAtansY2*y2e4 - SumAtansZ1*z1e4 + SumAtansZ2*z2e4) + x1*(4.*(-SumLogsX1Y1*y1e3 + SumLogsX1Y2*y2e3) + 2.*(SumLogsX1Z1*z1e3 - SumLogsX1Z2*z2e3) + CmnPrt2) + x2*(4.*CmnPrt3 - 2.*CmnPrt4 + CmnPrt5))/6.;

					InvDelYe2InvDelZe2 = InvDelYe2*InvDelZe2;
					T022 = InvDelYe2InvDelZe2*T022; S022 = InvDelYe2InvDelZe2*S022;
				}
			}
		}
		if(int(kx)>1)
		{
			if(int(ky)>0)
			{
				if(int(kz)>0)
				{
					TVector3d& T100 = T[1][0][0];
					TVector3d& S100 = S[1][0][0];

					double CmnPrt1 = y1*(SumLogsX1Y1-SumLogsX2Y1) - y2*(SumLogsX1Y2-SumLogsX2Y2);
					double CmnPrt2 = z1*(SumLogsX1Z1-SumLogsX2Z1) - z2*(SumLogsX1Z2-SumLogsX2Z2);
					T100.x = CmnPrt2 + CmnPrt1;
					S100.x = -y1*SumAtansY1 + y2*SumAtansY2 + z1*(SumLogsY1Z1-SumLogsY2Z1) - z2*(SumLogsY1Z2-SumLogsY2Z2);
					S100.y = -z1*SumAtansZ1 + z2*SumAtansZ2 + y1*(SumLogsY1Z1-SumLogsY1Z2) - y2*(SumLogsY2Z1-SumLogsY2Z2);
					T100.y = -CmnPrt1;
					S100.z = -R111+R121+R211-R221+R112-R122-R212+R222;
					T100.z = -CmnPrt2;

					T100 = InvDelX*T100; S100 = InvDelX*S100;
				}
				if(int(kz)>1)
				{
					TVector3d& T101 = T[1][0][1];
					TVector3d& S101 = S[1][0][1];

					double CmnPrt1 = 0.5*(y1*(R111-R211-R112+R212) - y2*(R121-R221-R122+R222));
					double CmnPrt2 = 0.5*(-(x1x1-z1z1)*SumLogsX1Z1 + (x1x1-z2z2)*SumLogsX1Z2 + (x2x2-z1z1)*SumLogsX2Z1 - (x2x2-z2z2)*SumLogsX2Z2);
					T101.x = CmnPrt1 + CmnPrt2;
					S101.x = 0.5*(x1*(-R111+R121+R112-R122) - x2*(-R211+R221+R212-R222) + (y1y1+z1z1)*SumLogsY1Z1 - (y1y1+z2z2)*SumLogsY1Z2 - (y2y2+z1z1)*SumLogsY2Z1 + (y2y2+z2z2)*SumLogsY2Z2);
					S101.y = 0.5*(-x1x1*SumAtansX1 + x2x2*SumAtansX2 + y1y1*SumAtansY1 - y2y2*SumAtansY2 - z1z1*SumAtansZ1 + z2z2*SumAtansZ2);
					T101.y = -2.*CmnPrt1;
					S101.z = 0.5*(z1*(-R111+R211+R121-R221) - z2*(-R112+R212+R122-R222) + (x1x1+y1y1)*SumLogsX1Y1 - (x1x1+y2y2)*SumLogsX1Y2 - (x2x2+y1y1)*SumLogsX2Y1 + (x2x2+y2y2)*SumLogsX2Y2);
					T101.z = CmnPrt1 - CmnPrt2;

					InvDelXInvDelZ = InvDelX*InvDelZ;
					T101 = InvDelXInvDelZ*T101; S101 = InvDelXInvDelZ*S101;
				}
				if(int(kz)>2)
				{
					TVector3d& T102 = T[1][0][2];
					TVector3d& S102 = S[1][0][2];

					double CmnPrt1 = 2.*(-SumAtansX1*x1e3 + SumAtansX2*x2e3);
					double CmnPrt2 = (SumLogsX2Y1 - SumLogsX1Y1)*y1e3 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e3;
					double CmnPrt3 = (SumLogsX1Z1 - SumLogsX2Z1)*z1e3 + (SumLogsX2Z2 - SumLogsX1Z2)*z2e3;
					double CmnPrt4 = y1*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2*((-R121 + R221)*z1 + (R122 - R222)*z2);
					T102.x = CmnPrt1/3. + CmnPrt2/6. + 0.5*(x1x1*(SumLogsX1Y1*y1 - SumLogsX1Y2*y2) + x2x2*(-SumLogsX2Y1*y1 + SumLogsX2Y2*y2)) + CmnPrt3/3. + CmnPrt4/6.;
					S102.x = ((SumLogsX1Y1 - SumLogsX1Y2)*x1e3 + (-SumLogsX2Y1 + SumLogsX2Y2)*x2e3 + SumAtansY1*y1e3 - SumAtansY2*y2e3 + (SumLogsY1Z1 - SumLogsY2Z1)*z1e3 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e3 + x1*((-R111 + R121)*z1 + (R112 - R122)*z2) + x2*((R211 - R221)*z1 + (-R212 + R222)*z2))/3.;
					S102.y = (2.*((SumLogsX1Z1 - SumLogsX1Z2)*x1e3 + (-SumLogsX2Z1 + SumLogsX2Z2)*x2e3) + (-SumLogsY1Z1 + SumLogsY1Z2)*y1e3 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e3 + x1*((R111 - R112)*y1 + (-R121 + R122)*y2) + x2*((-R211 + R212)*y1 + (R221 - R222)*y2) - SumAtansZ1*z1e3 + SumAtansZ2*z2e3)/3.;
					T102.y = 0.5*(-CmnPrt2 + x1x1*(SumLogsX1Y1*y1 - SumLogsX1Y2*y2) + x2x2*(-SumLogsX2Y1*y1 + SumLogsX2Y2*y2) - CmnPrt4);
					S102.z = (2.*((R111 - R112 - R121 + R122)*x1x1 + (-R211 + R212 + R221 - R222)*x2x2 + (R111 - R112 - R211 + R212)*y1y1 + (-R121 + R122 + R221 - R222)*y2y2) + (-R111 + R121 + R211 - R221)*z1z1 + (R112 - R122 - R212 + R222)*z2z2)/3.;
					T102.z = (-CmnPrt1 + CmnPrt2 + 3.*(x1x1*(-SumLogsX1Y1*y1 + SumLogsX1Y2*y2) + x2x2*(SumLogsX2Y1*y1 - SumLogsX2Y2*y2)) + (-SumLogsX1Z1 + SumLogsX2Z1)*z1e3 + (SumLogsX1Z2 - SumLogsX2Z2)*z2e3 + CmnPrt4)/3.;

					InvDelXInvDelZe2 = InvDelX*InvDelZe2;
					T102 = InvDelXInvDelZe2*T102; S102 = InvDelXInvDelZe2*S102;
				}
			}
			if(int(ky)>1)
			{
				if(int(kz)>0)
				{
					TVector3d& T110 = T[1][1][0];
					TVector3d& S110 = S[1][1][0];

					double CmnPrt1 = 0.5*(z1*(R111-R121-R211+R221) - z2*(R112-R122-R212+R222));
					double CmnPrt2 = 0.5*(-(x1x1-y1y1)*SumLogsX1Y1 + (x1x1-y2y2)*SumLogsX1Y2 + (x2x2-y1y1)*SumLogsX2Y1 - (x2x2-y2y2)*SumLogsX2Y2);
					T110.x = CmnPrt1 + CmnPrt2;
					S110.x = 0.5*(-x1x1*SumAtansX1 + x2x2*SumAtansX2 - y1y1*SumAtansY1 + y2y2*SumAtansY2 + z1z1*SumAtansZ1 - z2z2*SumAtansZ2);
					S110.y = 0.5*(x1*(-R111+R121+R112-R122) - x2*(-R211+R221+R212-R222) + (y1y1+z1z1)*SumLogsY1Z1 - (y1y1+z2z2)*SumLogsY1Z2 - (y2y2+z1z1)*SumLogsY2Z1 + (y2y2+z2z2)*SumLogsY2Z2);
					T110.y = CmnPrt1 - CmnPrt2;
					S110.z = 0.5*(y1*(-R111+R211+R112-R212) - y2*(-R121+R221+R122-R222) + (x1x1+z1z1)*SumLogsX1Z1 - (x1x1+z2z2)*SumLogsX1Z2 - (x2x2+z1z1)*SumLogsX2Z1 + (x2x2+z2z2)*SumLogsX2Z2);
					T110.z = -2.*CmnPrt1;
	
					InvDelXInvDelY = InvDelX*InvDelY;
					T110 = InvDelXInvDelY*T110; S110 = InvDelXInvDelY*S110;
				}
				if(int(kz)>1)
				{
					TVector3d& T111 = T[1][1][1];
					TVector3d& S111 = S[1][1][1];

					T111.x = (-(2.*x1x1-y1y1-z1z1)*R111 + (2.*x1x1-y2y2-z1z1)*R121 + (2.*x2x2-y1y1-z1z1)*R211 - (2.*x2x2-y2y2-z1z1)*R221 + (2.*x1x1-y1y1-z2z2)*R112 - (2.*x1x1-y2y2-z2z2)*R122 - (2.*x2x2-y1y1-z2z2)*R212 + (2.*x2x2-y2y2-z2z2)*R222)/3.;
					S111.x = (-x1*y1*(R111-R112) + x1*y2*(R121-R122) + x2*y1*(R211-R212) - x2*y2*(R221-R222) + z1e3*SumAtansZ1 - z2e3*SumAtansZ2 + x1e3*(SumLogsX1Z1-SumLogsX1Z2) - x2e3*(SumLogsX2Z1-SumLogsX2Z2) + y1e3*(SumLogsY1Z1-SumLogsY1Z2) - y2e3*(SumLogsY2Z1-SumLogsY2Z2))/3.;
					S111.y = (-x1*z1*(R111-R121) + x1*z2*(R112-R122) + x2*z1*(R211-R221) - x2*z2*(R212-R222) + y1e3*SumAtansY1 - y2e3*SumAtansY2 + x1e3*(SumLogsX1Y1-SumLogsX1Y2) - x2e3*(SumLogsX2Y1-SumLogsX2Y2) + z1e3*(SumLogsY1Z1-SumLogsY2Z1) - z2e3*(SumLogsY1Z2-SumLogsY2Z2))/3.;
					T111.y = (-(2.*y1y1-x1x1-z1z1)*R111 + (2.*y1y1-x2x2-z1z1)*R211 + (2.*y2y2-x1x1-z1z1)*R121 - (2.*y2y2-x2x2-z1z1)*R221 + (2.*y1y1-x1x1-z2z2)*R112 - (2.*y1y1-x2x2-z2z2)*R212 - (2.*y2y2-x1x1-z2z2)*R122 + (2.*y2y2-x2x2-z2z2)*R222)/3.;
					S111.z = (-y1*z1*(R111-R211) + y1*z2*(R112-R212) + y2*z1*(R121-R221) - y2*z2*(R122-R222) + x1e3*SumAtansX1 - x2e3*SumAtansX2 + y1e3*(SumLogsX1Y1-SumLogsX2Y1) - y2e3*(SumLogsX1Y2-SumLogsX2Y2) + z1e3*(SumLogsX1Z1-SumLogsX2Z1) - z2e3*(SumLogsX1Z2-SumLogsX2Z2))/3.;
					T111.z = (-(2.*z1z1-x1x1-y1y1)*R111 + (2.*z1z1-x2x2-y1y1)*R211 + (2.*z2z2-x1x1-y1y1)*R112 - (2.*z2z2-x2x2-y1y1)*R212 + (2.*z1z1-x1x1-y2y2)*R121 - (2.*z1z1-x2x2-y2y2)*R221 - (2.*z2z2-x1x1-y2y2)*R122 + (2.*z2z2-x2x2-y2y2)*R222)/3.;

					InvDelXInvDelYInvDelZ = InvDelX*InvDelY*InvDelZ;
					T111 = InvDelXInvDelYInvDelZ*T111; S111 = InvDelXInvDelYInvDelZ*S111;
				}
				if(int(kz)>2)
				{
					TVector3d& T112 = T[1][1][2];
					TVector3d& S112 = S[1][1][2];

					double CmnPrt1 = (-SumLogsX1Y1 + SumLogsX1Y2)*x1e4 + (SumLogsX2Y1 - SumLogsX2Y2)*x2e4;
					double CmnPrt2 = (-SumLogsX1Z1 + SumLogsX1Z2)*x1e4 + (SumLogsX2Z1 - SumLogsX2Z2)*x2e4;
					double CmnPrt3 = 2.*((R111 - R121 - R211 + R221)*z1e3 + (-R112 + R122 + R212 - R222)*z2e3);
					double CmnPrt4 = (-R111 + R121)*z1 + (R112 - R122)*z2;
					double CmnPrt5 = y1y1*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2y2*((-R121 + R221)*z1 + (R122 - R222)*z2);
					double CmnPrt6 = 2.*(-SumLogsX2Y1*y1y1 + SumLogsX2Y2*y2y2);
					double CmnPrt7 = (R211 - R221)*z1 + (-R212 + R222)*z2;
					T112.x = (3.*(-CmnPrt1) + (-SumLogsX1Y1 + SumLogsX2Y1)*y1e4 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e4 + CmnPrt3 + x1x1*(2.*(SumLogsX1Y1*y1y1 - SumLogsX1Y2*y2y2) + 3.*CmnPrt4) + CmnPrt5 + x2x2*(CmnPrt6 + 3.*CmnPrt7))/8.;
					S112.x = (SumAtansX1*x1e4 - SumAtansX2*x2e4 + SumAtansY1*y1e4 - SumAtansY2*y2e4 + SumAtansZ1*z1e4 - SumAtansZ2*z2e4 + x1*(y1*(-R111*z1 + R112*z2) + y2*(R121*z1 - R122*z2)) + x2*(y1*(R211*z1 - R212*z2) + y2*(-R221*z1 + R222*z2)))/4.;
					S112.y = (2.*((R111 - R112 - R121 + R122)*x1e3 + (-R211 + R212 + R221 - R222)*x2e3) + (-SumLogsY1Z1 + SumLogsY1Z2)*y1e4 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e4 + (SumLogsY1Z1 - SumLogsY2Z1)*z1e4 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e4 + x1*((R111 - R112)*y1y1 + (-R121 + R122)*y2y2 + (-R111 + R121)*z1z1 + (R112 - R122)*z2z2) + x2*((-R211 + R212)*y1y1 + (R221 - R222)*y2y2 + (R211 - R221)*z1z1 + (-R212 + R222)*z2z2))/4.;
					T112.y = (CmnPrt1 + 3.*(SumLogsX1Y1 - SumLogsX2Y1)*y1e4 + 3.*(-SumLogsX1Y2 + SumLogsX2Y2)*y2e4 + CmnPrt3 + x1x1*(2.*(SumLogsX1Y1*y1y1 - SumLogsX1Y2*y2y2) - CmnPrt4) - 3.*CmnPrt5 + x2x2*(CmnPrt6 - CmnPrt7))/8.;
					S112.z = (CmnPrt2 + 2.*((R111 - R112 - R211 + R212)*y1e3 + (-R121 + R122 + R221 - R222)*y2e3) + x1e2*((R111 - R112)*y1 + (-R121 + R122)*y2) + x2e2*((-R211 + R212)*y1 + (R221 - R222)*y2) + (SumLogsX1Z1 - SumLogsX2Z1)*z1e4 + (-SumLogsX1Z2 + SumLogsX2Z2)*z2e4 + y1*((-R111 + R211)*z1e2 + (R112 - R212)*z2e2) + y2*((R121 - R221)*z1e2 + (-R122 + R222)*z2e2))/4.;
					T112.z = (CmnPrt1 + (-SumLogsX1Y1 + SumLogsX2Y1)*y1e4 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e4 - CmnPrt3 + x1e2*(2.*(-SumLogsX1Y1*y1e2 + SumLogsX1Y2*y2e2) - CmnPrt4) + y1e2*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2e2*((-R121 + R221)*z1 + (R122 - R222)*z2) + x2e2*(2.*(SumLogsX2Y1*y1e2 - SumLogsX2Y2*y2e2) + (-R211 + R221)*z1 + (R212 - R222)*z2))/4.;

					InvDelXInvDelYInvDelZe2 = InvDelXInvDelYInvDelZ*InvDelZ;
					T112 = InvDelXInvDelYInvDelZe2*T112; S112 = InvDelXInvDelYInvDelZe2*S112;
				}
			}
			if(int(ky)>2)
			{
				if(int(kz)>0)
				{
					TVector3d& T120 = T[1][2][0];
					TVector3d& S120 = S[1][2][0];

					double CmnPrt1 = 2.*(SumAtansX1*x1e3 - SumAtansX2*x2e3) + (-SumLogsX1Y1 + SumLogsX2Y1)*y1e3 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e3;
					double CmnPrt2 = x1e2*(SumLogsX1Z1*z1 - SumLogsX1Z2*z2) + x2e2*(-SumLogsX2Z1*z1 + SumLogsX2Z2*z2);
					double CmnPrt3 = (-SumLogsX1Z1 + SumLogsX2Z1)*z1e3 + (SumLogsX1Z2 - SumLogsX2Z2)*z2e3 + y1*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2*((-R121 + R221)*z1 + (R122 - R222)*z2);
					T120.x = (-2.*CmnPrt1 + CmnPrt3 + 3.*CmnPrt2)/6.;
					S120.x = (2.*((SumLogsX1Y1 - SumLogsX1Y2)*x1e3 + (-SumLogsX2Y1 + SumLogsX2Y2)*x2e3) - SumAtansY1*y1e3 + SumAtansY2*y2e3 + (-SumLogsY1Z1 + SumLogsY2Z1)*z1e3 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e3 + x1*((R111 - R121)*z1 + (-R112 + R122)*z2) + x2*((-R211 + R221)*z1 + (R212 - R222)*z2))/3.;
					S120.y = ((SumLogsX1Z1 - SumLogsX1Z2)*x1e3 + (-SumLogsX2Z1 + SumLogsX2Z2)*x2e3 + (SumLogsY1Z1 - SumLogsY1Z2)*y1e3 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2e3 + x1*((-R111 + R112)*y1 + (R121 - R122)*y2) + x2*((R211 - R212)*y1 + (-R221 + R222)*y2) + SumAtansZ1*z1e3 - SumAtansZ2*z2e3)/3.;
					T120.y = (CmnPrt1 + CmnPrt3 - 3.*CmnPrt2)/3.;
					S120.z = (2.*((R111 - R112 - R121 + R122)*x1e2 + (-R211 + R212 + R221 - R222)*x2e2 + (R111 - R121 - R211 + R221)*z1e2 + (-R112 + R122 + R212 - R222)*z2e2) + (-R111 + R112 + R211 - R212)*y1e2 + (R121 - R122 - R221 + R222)*y2e2)/3.;
					T120.z = ((SumLogsX1Z1 - SumLogsX2Z1)*z1e3 + (-SumLogsX1Z2 + SumLogsX2Z2)*z2e3 + y1*((-R111 + R211)*z1 + (R112 - R212)*z2) + y2*((R121 - R221)*z1 + (-R122 + R222)*z2) + CmnPrt2)/2.;
	
					InvDelXInvDelYe2 = InvDelXInvDelY*InvDelY;
					T120 = InvDelXInvDelYe2*T120; S120 = InvDelXInvDelYe2*S120;
				}
				if(int(kz)>1)
				{
					TVector3d& T121 = T[1][2][1];
					TVector3d& S121 = S[1][2][1];

					double CmnPrt1 = (-SumLogsX1Z1 + SumLogsX1Z2)*x1e4 + (SumLogsX2Z1 - SumLogsX2Z2)*x2e4;
					double CmnPrt2 = 2.*((R111 - R112 - R211 + R212)*y1e3 + (-R121 + R122 + R221 - R222)*y2e3);
					double CmnPrt3 = (SumLogsX1Z1 - SumLogsX2Z1)*z1e4 + (-SumLogsX1Z2 + SumLogsX2Z2)*z2e4 + y1*((-R111 + R211)*z1e2 + (R112 - R212)*z2e2) + y2*((R121 - R221)*z1e2 + (-R122 + R222)*z2e2);
					double CmnPrt4 = (R111 - R112)*y1 + (-R121 + R122)*y2;
					double CmnPrt5 = 2.*(SumLogsX1Z1*z1e2 - SumLogsX1Z2*z2e2);
					double CmnPrt6 = (-R211 + R212)*y1 + (R221 - R222)*y2;
					double CmnPrt7 = 2.*(-SumLogsX2Z1*z1e2 + SumLogsX2Z2*z2e2);
					T121.x = (3.*(-CmnPrt1) + CmnPrt2 - CmnPrt3 + x1e2*(3.*(-CmnPrt4) + CmnPrt5) + x2e2*(3.*(-CmnPrt6) + CmnPrt7))/8.;
					S121.x = (2.*((R111 - R112 - R121 + R122)*x1e3 + (-R211 + R212 + R221 - R222)*x2e3) + (SumLogsY1Z1 - SumLogsY1Z2)*y1e4 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2e4 + (-SumLogsY1Z1 + SumLogsY2Z1)*z1e4 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e4 + x1*((-R111 + R112)*y1e2 + (R121 - R122)*y2e2 + (R111 - R121)*z1e2 + (-R112 + R122)*z2e2) + x2*((R211 - R212)*y1e2 + (-R221 + R222)*y2e2 + (-R211 + R221)*z1e2 + (R212 - R222)*z2e2))/4.;
					S121.y = (SumAtansX1*x1e4 - SumAtansX2*x2e4 + SumAtansY1*y1e4 - SumAtansY2*y2e4 + SumAtansZ1*z1e4 - SumAtansZ2*z2e4 + x1*(y1*(-(R111*z1) + R112*z2) + y2*(R121*z1 - R122*z2)) + x2*(y1*(R211*z1 - R212*z2) + y2*(-(R221*z1) + R222*z2)))/4.;
					T121.y = (CmnPrt1 - CmnPrt2 - CmnPrt3 + x1e2*(CmnPrt4 - CmnPrt5) + x2e2*(CmnPrt6 - CmnPrt7))/4.;
					S121.z = ((-SumLogsX1Y1 + SumLogsX1Y2)*x1e4 + (SumLogsX2Y1 - SumLogsX2Y2)*x2e4 + (SumLogsX1Y1 - SumLogsX2Y1)*y1e4 + (-SumLogsX1Y2 + SumLogsX2Y2)*y2e4 + 2.*((R111 - R121 - R211 + R221)*z1e3 + (-R112 + R122 + R212 - R222)*z2e3) + x1e2*((R111 - R121)*z1 + (-R112 + R122)*z2) + y1e2*((-R111 + R211)*z1 + (R112 - R212)*z2) + x2e2*((-R211 + R221)*z1 + (R212 - R222)*z2) + y2e2*((R121 - R221)*z1 + (-R122 + R222)*z2))/4.;
					T121.z = (CmnPrt1 + CmnPrt2 + 3.*CmnPrt3 + x1e2*(CmnPrt4 + CmnPrt5) + x2e2*(CmnPrt6 + CmnPrt7))/8.;

					InvDelXInvDelYe2InvDelZ = InvDelXInvDelYe2*InvDelZ;
					T121 = InvDelXInvDelYe2InvDelZ*T121; S121 = InvDelXInvDelYe2InvDelZ*S121;
				}
				if(int(kz)>2)
				{
					TVector3d& T122 = T[1][2][2];
					TVector3d& S122 = S[1][2][2];

					T122.x = (8.*(SumAtansX1*x1e5 - SumAtansX2*x2e5) + 3.*((-SumLogsX1Y1 + SumLogsX2Y1)*y1e5 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e5 + (-SumLogsX1Z1 + SumLogsX2Z1)*z1e5 + (SumLogsX1Z2 - SumLogsX2Z2)*z2e5 + y1e3*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2e3*((-R121 + R221)*z1 + (R122 - R222)*z2) + y1*((R111 - R211)*z1e3 + (-R112 + R212)*z2e3) + y2*((-R121 + R221)*z1e3 + (R122 - R222)*z2e3)) + x1e2*(5.*(SumLogsX1Y1*y1e3 - SumLogsX1Y2*y2e3 + SumLogsX1Z1*z1e3 - SumLogsX1Z2*z2e3) + 8.*(y1*(-R111*z1 + R112*z2) + y2*(R121*z1 - R122*z2))) + x2e2*(5.*(-SumLogsX2Y1*y1e3 + SumLogsX2Y2*y2e3 - SumLogsX2Z1*z1e3 + SumLogsX2Z2*z2e3) + 8.*(y1*(R211*z1 - R212*z2) + y2*(-R221*z1 + R222*z2))))/30.;
					S122.x = ((-SumLogsX1Y1 + SumLogsX1Y2)*x1e5 + (SumLogsX2Y1 - SumLogsX2Y2)*x2e5 + SumAtansY1*y1e5 - SumAtansY2*y2e5 + (-SumLogsY1Z1 + SumLogsY2Z1)*z1e5 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e5 + x1e3*((R111 - R121)*z1 + (-R112 + R122)*z2) + x2e3*((-R211 + R221)*z1 + (R212 - R222)*z2) + x1*((R111 - R121)*z1e3 + (-R112 + R122)*z2e3 + y1e2*(-(R111*z1) + R112*z2) + y2e2*(R121*z1 - R122*z2)) + x2*((-R211 + R221)*z1e3 + (R212 - R222)*z2e3 + y1e2*(R211*z1 - R212*z2) + y2e2*(-(R221*z1) + R222*z2)))/5.;
					S122.y = ((-SumLogsX1Z1 + SumLogsX1Z2)*x1e5 + (SumLogsX2Z1 - SumLogsX2Z2)*x2e5 + (-SumLogsY1Z1 + SumLogsY1Z2)*y1e5 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e5 + x1e3*((R111 - R112)*y1 + (-R121 + R122)*y2) + x2e3*((-R211 + R212)*y1 + (R221 - R222)*y2) + SumAtansZ1*z1e5 - SumAtansZ2*z2e5 + x1*((R111 - R112)*y1e3 + (-R121 + R122)*y2e3 + y1*(-(R111*z1e2) + R112*z2e2) + y2*(R121*z1e2 - R122*z2e2)) + x2*((-R211 + R212)*y1e3 + (R221 - R222)*y2e3 + y1*(R211*z1e2 - R212*z2e2) + y2*(-(R221*z1e2) + R222*z2e2)))/5.;
					T122.y = (4.*(-SumAtansX1*x1e5 + SumAtansX2*x2e5) + 9.*((SumLogsX1Y1 - SumLogsX2Y1)*y1e5 + (-SumLogsX1Y2 + SumLogsX2Y2)*y2e5) + 6.*((-SumLogsX1Z1 + SumLogsX2Z1)*z1e5 + (SumLogsX1Z2 - SumLogsX2Z2)*z2e5) + 9.*(y1e3*((-R111 + R211)*z1 + (R112 - R212)*z2) + y2e3*((R121 - R221)*z1 + (-R122 + R222)*z2)) + 6.*(y1*((R111 - R211)*z1e3 + (-R112 + R212)*z2e3) + y2*((-R121 + R221)*z1e3 + (R122 - R222)*z2e3)) + x1e2*(5.*(SumLogsX1Y1*y1e3 - SumLogsX1Y2*y2e3) + 10.*(-SumLogsX1Z1*z1e3 + SumLogsX1Z2*z2e3) + 4.*(y1*(R111*z1 - R112*z2) + y2*(-R121*z1 + R122*z2))) + x2e2*(5.*(-SumLogsX2Y1*y1e3 + SumLogsX2Y2*y2e3) + 10.*(SumLogsX2Z1*z1e3 - SumLogsX2Z2*z2e3) + 4.*(y1*(-R211*z1 + R212*z2) + y2*(R221*z1 - R222*z2))))/30.;
					S122.z = (4.*((-R111 + R112 + R121 - R122)*x1e4 + (R211 - R212 - R221 + R222)*x2e4) + 6.*((R111 - R112 - R211 + R212)*y1e4 + (-R121 + R122 + R221 - R222)*y2e4 + (R111 - R121 - R211 + R221)*z1e4 + (-R112 + R122 + R212 - R222)*z2e4) + 2.*(x1e2*((R111 - R112)*y1e2 + (-R121 + R122)*y2e2 + (R111 - R121)*z1e2 + (-R112 + R122)*z2e2)) + 3.*y1e2*((-R111 + R211)*z1e2 + (R112 - R212)*z2e2) + 2.*x2e2*((-R211 + R212)*y1e2 + (R221 - R222)*y2e2 + (-R211 + R221)*z1e2 + (R212 - R222)*z2e2) + 3.*y2e2*((R121 - R221)*z1e2 + (-R122 + R222)*z2e2))/15.;
					T122.z = (4.*(-SumAtansX1*x1e5 + SumAtansX2*x2e5) + 6.*((-SumLogsX1Y1 + SumLogsX2Y1)*y1e5 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e5) + 9.*((SumLogsX1Z1 - SumLogsX2Z1)*z1e5 + (-SumLogsX1Z2 + SumLogsX2Z2)*z2e5) + 6.*(y1e3*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2e3*((-R121 + R221)*z1 + (R122 - R222)*z2)) + 9.*(y1*((-R111 + R211)*z1e3 + (R112 - R212)*z2e3) + y2*((R121 - R221)*z1e3 + (-R122 + R222)*z2e3)) + x1e2*(10.*(-SumLogsX1Y1*y1e3 + SumLogsX1Y2*y2e3) + 5.*(SumLogsX1Z1*z1e3 - SumLogsX1Z2*z2e3) + 4.*(y1*(R111*z1 - R112*z2) + y2*(-R121*z1 + R122*z2))) + x2e2*(10.*(SumLogsX2Y1*y1e3 - SumLogsX2Y2*y2e3) + 5.*(-SumLogsX2Z1*z1e3 + SumLogsX2Z2*z2e3) + 4.*(y1*(-R211*z1 + R212*z2) + y2*(R221*z1 - R222*z2))))/30.;

					InvDelXInvDelYe2InvDelZe2 = InvDelXInvDelYe2InvDelZ*InvDelZ;
					T122 = InvDelXInvDelYe2InvDelZe2*T122; S122 = InvDelXInvDelYe2InvDelZe2*S122;
				}
			}
		}
		if(int(kx)>2)
		{
			if(int(ky)>0)
			{
				if(int(kz)>0)
				{
					TVector3d& T200 = T[2][0][0];
					TVector3d& S200 = S[2][0][0];

					T200.x = SumAtansY1*y1e2 - SumAtansY2*y2e2 + SumAtansZ1*z1e2 - SumAtansZ2*z2e2 + 2.*(y1*(-SumLogsY1Z1*z1 + SumLogsY1Z2*z2) + y2*(SumLogsY2Z1*z1 - SumLogsY2Z2*z2));
					S200.x = (SumLogsX1Y1 - SumLogsX2Y1)*y1e2 + (-SumLogsX1Y2 + SumLogsX2Y2)*y2e2 + (R111 - R121 - R211 + R221)*z1 + (-R112 + R122 + R212 - R222)*z2;
					S200.y = (R111 - R112 - R211 + R212)*y1 + (-R121 + R122 + R221 - R222)*y2 + (SumLogsX1Z1 - SumLogsX2Z1)*z1e2 + (-SumLogsX1Z2 + SumLogsX2Z2)*z2e2;
					T200.y = -(SumAtansY1*y1e2) + SumAtansY2*y2e2 + y1*(SumLogsY1Z1*z1 - SumLogsY1Z2*z2) + y2*(-(SumLogsY2Z1*z1) + SumLogsY2Z2*z2);
					S200.z = ((-R111 + R112 + R121 - R122)*x1 + (R211 - R212 - R221 + R222)*x2 + (SumLogsY1Z1 - SumLogsY1Z2)*y1e2 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2e2 + (SumLogsY1Z1 - SumLogsY2Z1)*z1e2 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e2)/2.;
					T200.z = -(SumAtansZ1*z1e2) + SumAtansZ2*z2e2 + y1*(SumLogsY1Z1*z1 - SumLogsY1Z2*z2) + y2*(-(SumLogsY2Z1*z1) + SumLogsY2Z2*z2);

					InvDelXe2 = InvDelX*InvDelX;
					T200 = InvDelXe2*T200; S200 = InvDelXe2*S200;
				}
				if(int(kz)>1)
				{
					TVector3d& T201 = T[2][0][1];
					TVector3d& S201 = S[2][0][1];

					T201.x = ((-SumLogsX1Z1 + SumLogsX1Z2)*x1e3 + (SumLogsX2Z1 - SumLogsX2Z2)*x2e3 + (-SumLogsY1Z1 + SumLogsY1Z2)*y1e3 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e3 + x1*((R111 - R112)*y1 + (-R121 + R122)*y2) + x2*((-R211 + R212)*y1 + (R221 - R222)*y2) + 2.*(SumAtansZ1*z1e3 - SumAtansZ2*z2e3) + 3.*(y1*(-SumLogsY1Z1*z1e2 + SumLogsY1Z2*z2e2) + y2*(SumLogsY2Z1*z1e2 - SumLogsY2Z2*z2e2)))/3.;
					S201.x = ((-R111 + R112 + R121 - R122)*x1e2 + (R211 - R212 - R221 + R222)*x2e2 + 2.*((R111 - R112 - R211 + R212)*y1e2 + (-R121 + R122 + R221 - R222)*y2e2 + (R111 - R121 - R211 + R221)*z1e2 + (-R112 + R122 + R212 - R222)*z2e2))/3.;
					S201.y = (-(SumAtansX1*x1e3) + SumAtansX2*x2e3 + (-SumLogsX1Y1 + SumLogsX2Y1)*y1e3 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e3 + 2.*((SumLogsX1Z1 - SumLogsX2Z1)*z1e3 + (-SumLogsX1Z2 + SumLogsX2Z2)*z2e3) + y1*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2*((-R121 + R221)*z1 + (R122 - R222)*z2))/3.;
					T201.y = ((SumLogsY1Z1 - SumLogsY1Z2)*y1e3 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2e3 + x1*((-R111 + R112)*y1 + (R121 - R122)*y2) + x2*((R211 - R212)*y1 + (-R221 + R222)*y2) + y1*(SumLogsY1Z1*z1e2 - SumLogsY1Z2*z2e2) + y2*(-(SumLogsY2Z1*z1e2) + SumLogsY2Z2*z2e2))/2.;
					S201.z = ((SumLogsX1Y1 - SumLogsX1Y2)*x1e3 + (-SumLogsX2Y1 + SumLogsX2Y2)*x2e3 + SumAtansY1*y1e3 - SumAtansY2*y2e3 + (SumLogsY1Z1 - SumLogsY2Z1)*z1e3 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e3 + x1*((-R111 + R121)*z1 + (R112 - R122)*z2) + x2*((R211 - R221)*z1 + (-R212 + R222)*z2))/3.;
					T201.z = (2.*((SumLogsX1Z1 - SumLogsX1Z2)*x1e3 + (-SumLogsX2Z1 + SumLogsX2Z2)*x2e3) + (-SumLogsY1Z1 + SumLogsY1Z2)*y1e3 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e3 + x1*((R111 - R112)*y1 + (-R121 + R122)*y2) + x2*((-R211 + R212)*y1 + (R221 - R222)*y2) + 4.*(-SumAtansZ1*z1e3 + SumAtansZ2*z2e3) + 3.*(y1*(SumLogsY1Z1*z1e2 - SumLogsY1Z2*z2e2) + y2*(-SumLogsY2Z1*z1e2 + SumLogsY2Z2*z2e2)))/6.;

					InvDelXe2InvDelZ = InvDelXe2*InvDelZ;
					T201 = InvDelXe2InvDelZ*T201; S201 = InvDelXe2InvDelZ*S201;
				}
				if(int(kz)>2)
				{
					TVector3d& T202 = T[2][0][2];
					TVector3d& S202 = S[2][0][2];

					T202.x = (3.*(-SumAtansX1*x1e4 + SumAtansX2*x2e4) - SumAtansY1*y1e4 + SumAtansY2*y2e4 + 2.*(x1e3*(SumLogsX1Y1*y1 - SumLogsX1Y2*y2) + x2e3*(-SumLogsX2Y1*y1 + SumLogsX2Y2*y2)) + 3.*(SumAtansZ1*z1e4 - SumAtansZ2*z2e4) + 4.*(y1*(-SumLogsY1Z1*z1e3 + SumLogsY1Z2*z2e3) + y2*(SumLogsY2Z1*z1e3 - SumLogsY2Z2*z2e3)) + x1*(y1*(R111*z1 - R112*z2) + y2*(-(R121*z1) + R122*z2)) + x2*(y1*(-(R211*z1) + R212*z2) + y2*(R221*z1 - R222*z2)))/6.;
					S202.x = ((SumLogsX1Y1 - SumLogsX1Y2)*x1e4 + (-SumLogsX2Y1 + SumLogsX2Y2)*x2e4 + (-SumLogsX1Y1 + SumLogsX2Y1)*y1e4 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e4 + 2.*((R111 - R121 - R211 + R221)*z1e3 + (-R112 + R122 + R212 - R222)*z2e3) + x1e2*((-R111 + R121)*z1 + (R112 - R122)*z2) + y1e2*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2e2*((-R121 + R221)*z1 + (R122 - R222)*z2) + x2e2*((R211 - R221)*z1 + (-R212 + R222)*z2))/4.;
					S202.y = (3.*((SumLogsX1Z1 - SumLogsX1Z2)*x1e4 + (-SumLogsX2Z1 + SumLogsX2Z2)*x2e4) + 2.*((-R111 + R112 + R211 - R212)*y1e3 + (R121 - R122 - R221 + R222)*y2e3) + x1e2*((R111 - R112)*y1 + (-R121 + R122)*y2) + x2e2*((-R211 + R212)*y1 + (R221 - R222)*y2) + 3.*((SumLogsX1Z1 - SumLogsX2Z1)*z1e4 + (-SumLogsX1Z2 + SumLogsX2Z2)*z2e4) + y1*((R111 - R211)*z1e2 + (-R112 + R212)*z2e2) + y2*((-R121 + R221)*z1e2 + (R122 - R222)*z2e2))/6.;
					T202.y = (SumAtansY1*y1e4 - SumAtansY2*y2e4 + x1e3*(SumLogsX1Y1*y1 - SumLogsX1Y2*y2) + x2e3*(-(SumLogsX2Y1*y1) + SumLogsX2Y2*y2) + y1*(SumLogsY1Z1*z1e3 - SumLogsY1Z2*z2e3) + y2*(-(SumLogsY2Z1*z1e3) + SumLogsY2Z2*z2e3) + x1*(y1*(-(R111*z1) + R112*z2) + y2*(R121*z1 - R122*z2)) + x2*(y1*(R211*z1 - R212*z2) + y2*(-(R221*z1) + R222*z2)))/3.;
					S202.z = (2.*((R111 - R112 - R121 + R122)*x1e3 + (-R211 + R212 + R221 - R222)*x2e3) + (-SumLogsY1Z1 + SumLogsY1Z2)*y1e4 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e4 + (SumLogsY1Z1 - SumLogsY2Z1)*z1e4 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e4 + x1*((R111 - R112)*y1e2 + (-R121 + R122)*y2e2 + (-R111 + R121)*z1e2 + (R112 - R122)*z2e2) + x2*((-R211 + R212)*y1e2 + (R221 - R222)*y2e2 + (R211 - R221)*z1e2 + (-R212 + R222)*z2e2))/4.;
					T202.z = (3.*(SumAtansX1*x1e4 - SumAtansX2*x2e4) - SumAtansY1*y1e4 + SumAtansY2*y2e4 + 4.*(x1e3*(-SumLogsX1Y1*y1 + SumLogsX1Y2*y2) + x2e3*(SumLogsX2Y1*y1 - SumLogsX2Y2*y2)) + 3.*(-SumAtansZ1*z1e4 + SumAtansZ2*z2e4) + 2.*(y1*(SumLogsY1Z1*z1e3 - SumLogsY1Z2*z2e3) + y2*(-SumLogsY2Z1*z1e3 + SumLogsY2Z2*z2e3)) + x1*(y1*(R111*z1 - R112*z2) + y2*(-(R121*z1) + R122*z2)) + x2*(y1*(-(R211*z1) + R212*z2) + y2*(R221*z1 - R222*z2)))/6.;

					InvDelXe2InvDelZe2 = InvDelXe2InvDelZ*InvDelZ;
					T202 = InvDelXe2InvDelZe2*T202; S202 = InvDelXe2InvDelZe2*S202;
				}
			}
			if(int(ky)>1)
			{
				if(int(kz)>0)
				{
					TVector3d& T210 = T[2][1][0];
					TVector3d& S210 = S[2][1][0];

					T210.x = ((-SumLogsX1Y1 + SumLogsX1Y2)*x1e3 + (SumLogsX2Y1 - SumLogsX2Y2)*x2e3 + 2.*(SumAtansY1*y1e3 - SumAtansY2*y2e3) + (-SumLogsY1Z1 + SumLogsY2Z1)*z1e3 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e3 + x1*((R111 - R121)*z1 + (-R112 + R122)*z2) + x2*((-R211 + R221)*z1 + (R212 - R222)*z2) + 3.*(y1e2*(-SumLogsY1Z1*z1 + SumLogsY1Z2*z2) + y2e2*(SumLogsY2Z1*z1 - SumLogsY2Z2*z2)))/3.;
					S210.x = (-(SumAtansX1*x1e3) + SumAtansX2*x2e3 + 2.*((SumLogsX1Y1 - SumLogsX2Y1)*y1e3 + (-SumLogsX1Y2 + SumLogsX2Y2)*y2e3) + (-SumLogsX1Z1 + SumLogsX2Z1)*z1e3 + (SumLogsX1Z2 - SumLogsX2Z2)*z2e3 + y1*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2*((-R121 + R221)*z1 + (R122 - R222)*z2))/3.;
					S210.y = ((-R111 + R112 + R121 - R122)*x1e2 + (R211 - R212 - R221 + R222)*x2e2 + 2.*((R111 - R112 - R211 + R212)*y1e2 + (-R121 + R122 + R221 - R222)*y2e2 + (R111 - R121 - R211 + R221)*z1e2 + (-R112 + R122 + R212 - R222)*z2e2))/3.;
					T210.y = (2.*((SumLogsX1Y1 - SumLogsX1Y2)*x1e3 + (-SumLogsX2Y1 + SumLogsX2Y2)*x2e3) + 4.*(-SumAtansY1*y1e3 + SumAtansY2*y2e3) + (-SumLogsY1Z1 + SumLogsY2Z1)*z1e3 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e3 + x1*((R111 - R121)*z1 + (-R112 + R122)*z2) + x2*((-R211 + R221)*z1 + (R212 - R222)*z2) + 3.*(y1e2*(SumLogsY1Z1*z1 - SumLogsY1Z2*z2) + y2e2*(-SumLogsY2Z1*z1 + SumLogsY2Z2*z2)))/6.;
					S210.z = ((SumLogsX1Z1 - SumLogsX1Z2)*x1e3 + (-SumLogsX2Z1 + SumLogsX2Z2)*x2e3 + (SumLogsY1Z1 - SumLogsY1Z2)*y1e3 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2e3 + x1*((-R111 + R112)*y1 + (R121 - R122)*y2) + x2*((R211 - R212)*y1 + (-R221 + R222)*y2) + SumAtansZ1*z1e3 - SumAtansZ2*z2e3)/3.;
					T210.z = ((SumLogsY1Z1 - SumLogsY2Z1)*z1e3 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e3 + x1*((-R111 + R121)*z1 + (R112 - R122)*z2) + x2*((R211 - R221)*z1 + (-R212 + R222)*z2) + y1e2*(SumLogsY1Z1*z1 - SumLogsY1Z2*z2) + y2e2*(-(SumLogsY2Z1*z1) + SumLogsY2Z2*z2))/2.;

					InvDelXe2InvDelY = InvDelXe2*InvDelY;
					T210 = InvDelXe2InvDelY*T210; S210 = InvDelXe2InvDelY*S210;
				}
				if(int(kz)>1)
				{
					TVector3d& T211 = T[2][1][1];
					TVector3d& S211 = S[2][1][1];

					T211.x = (2.*((-R111 + R112 + R121 - R122)*x1e3 + (R211 - R212 - R221 + R222)*x2e3) + (-SumLogsY1Z1 + SumLogsY1Z2)*y1e4 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e4 + (-SumLogsY1Z1 + SumLogsY2Z1)*z1e4 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e4 + x1*((R111 - R112)*y1e2 + (-R121 + R122)*y2e2 + (R111 - R121)*z1e2 + (-R112 + R122)*z2e2) + x2*((-R211 + R212)*y1e2 + (R221 - R222)*y2e2 + (-R211 + R221)*z1e2 + (R212 - R222)*z2e2) + 2.*(y1e2*(-SumLogsY1Z1*z1e2 + SumLogsY1Z2*z2e2) + y2e2*(SumLogsY2Z1*z1e2 - SumLogsY2Z2*z2e2)))/4.;
					S211.x = ((SumLogsX1Z1 - SumLogsX1Z2)*x1e4 + (-SumLogsX2Z1 + SumLogsX2Z2)*x2e4 + 2.*((R111 - R112 - R211 + R212)*y1e3 + (-R121 + R122 + R221 - R222)*y2e3) + x1e2*((-R111 + R112)*y1 + (R121 - R122)*y2) + x2e2*((R211 - R212)*y1 + (-R221 + R222)*y2) + (-SumLogsX1Z1 + SumLogsX2Z1)*z1e4 + (SumLogsX1Z2 - SumLogsX2Z2)*z2e4 + y1*((R111 - R211)*z1e2 + (-R112 + R212)*z2e2) + y2*((-R121 + R221)*z1e2 + (R122 - R222)*z2e2))/4.;
					S211.y = ((SumLogsX1Y1 - SumLogsX1Y2)*x1e4 + (-SumLogsX2Y1 + SumLogsX2Y2)*x2e4 + (-SumLogsX1Y1 + SumLogsX2Y1)*y1e4 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e4 + 2.*((R111 - R121 - R211 + R221)*z1e3 + (-R112 + R122 + R212 - R222)*z2e3) + x1e2*((-R111 + R121)*z1 + (R112 - R122)*z2) + y1e2*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2e2*((-R121 + R221)*z1 + (R122 - R222)*z2) + x2e2*((R211 - R221)*z1 + (-R212 + R222)*z2))/4.;
					T211.y = (2.*((R111 - R112 - R121 + R122)*x1e3 + (-R211 + R212 + R221 - R222)*x2e3) + 3.*((SumLogsY1Z1 - SumLogsY1Z2)*y1e4 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2e4) + (-SumLogsY1Z1 + SumLogsY2Z1)*z1e4 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e4 + x1*(3.*((-R111 + R112)*y1e2 + (R121 - R122)*y2e2) + (R111 - R121)*z1e2 + (-R112 + R122)*z2e2) + x2*(3.*((R211 - R212)*y1e2 + (-R221 + R222)*y2e2) + (-R211 + R221)*z1e2 + (R212 - R222)*z2e2) + 2.*(y1e2*(SumLogsY1Z1*z1e2 - SumLogsY1Z2*z2e2) + y2e2*(-SumLogsY2Z1*z1e2 + SumLogsY2Z2*z2e2)))/8.;
					S211.z = (SumAtansX1*x1e4 - SumAtansX2*x2e4 + SumAtansY1*y1e4 - SumAtansY2*y2e4 + SumAtansZ1*z1e4 - SumAtansZ2*z2e4 + x1*(y1*(-(R111*z1) + R112*z2) + y2*(R121*z1 - R122*z2)) + x2*(y1*(R211*z1 - R212*z2) + y2*(-(R221*z1) + R222*z2)))/4.;
					T211.z = (2.*((R111 - R112 - R121 + R122)*x1e3 + (-R211 + R212 + R221 - R222)*x2e3) + (-SumLogsY1Z1 + SumLogsY1Z2)*y1e4 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e4 + 3.*((SumLogsY1Z1 - SumLogsY2Z1)*z1e4 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e4) + x1*((R111 - R112)*y1e2 + (-R121 + R122)*y2e2 + 3.*((-R111 + R121)*z1e2 + (R112 - R122)*z2e2)) + x2*((-R211 + R212)*y1e2 + (R221 - R222)*y2e2 + 3.*((R211 - R221)*z1e2 + (-R212 + R222)*z2e2)) + 2.*(y1e2*(SumLogsY1Z1*z1e2 - SumLogsY1Z2*z2e2) + y2e2*(-SumLogsY2Z1*z1e2 + SumLogsY2Z2*z2e2)))/8.;

					InvDelXe2InvDelYInvDelZ = InvDelXe2InvDelY*InvDelZ;
					T211 = InvDelXe2InvDelYInvDelZ*T211; S211 = InvDelXe2InvDelYInvDelZ*S211;
				}
				if(int(kz)>2)
				{
					TVector3d& T212 = T[2][1][2];
					TVector3d& S212 = S[2][1][2];

					T212.x = (9.*((SumLogsX1Y1 - SumLogsX1Y2)*x1e5 + (-SumLogsX2Y1 + SumLogsX2Y2)*x2e5) + 4.*(-SumAtansY1*y1e5 + SumAtansY2*y2e5) + 6.*((-SumLogsY1Z1 + SumLogsY2Z1)*z1e5 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e5) + x1e3*(5.*(SumLogsX1Y1*y1e2 - SumLogsX1Y2*y2e2) + 9.*((-R111 + R121)*z1 + (R112 - R122)*z2)) + x2e3*(5.*(-SumLogsX2Y1*y1e2 + SumLogsX2Y2*y2e2) + 9.*((R211 - R221)*z1 + (-R212 + R222)*z2)) + 10.*(y1e2*(-SumLogsY1Z1*z1e3 + SumLogsY1Z2*z2e3) + y2e2*(SumLogsY2Z1*z1e3 - SumLogsY2Z2*z2e3)) + x1*(6.*((R111 - R121)*z1e3 + (-R112 + R122)*z2e3) + 4.*(y1e2*(R111*z1 - R112*z2) + y2e2*(-R121*z1 + R122*z2))) + x2*(6.*((-R211 + R221)*z1e3 + (R212 - R222)*z2e3) + 4.*(y1e2*(-R211*z1 + R212*z2) + y2e2*(R221*z1 - R222*z2))))/30.;
					S212.x = (SumAtansX1*x1e5 - SumAtansX2*x2e5 + (-SumLogsX1Y1 + SumLogsX2Y1)*y1e5 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e5 + (-SumLogsX1Z1 + SumLogsX2Z1)*z1e5 + (SumLogsX1Z2 - SumLogsX2Z2)*z2e5 + y1e3*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2e3*((-R121 + R221)*z1 + (R122 - R222)*z2) + y1*((R111 - R211)*z1e3 + (-R112 + R212)*z2e3) + y2*((-R121 + R221)*z1e3 + (R122 - R222)*z2e3) + x1e2*(y1*(-(R111*z1) + R112*z2) + y2*(R121*z1 - R122*z2)) + x2e2*(y1*(R211*z1 - R212*z2) + y2*(-(R221*z1) + R222*z2)))/5.;
					S212.y = (6.*((R111 - R112 - R121 + R122)*x1e4 + (-R211 + R212 + R221 - R222)*x2e4) + 4.*((-R111 + R112 + R211 - R212)*y1e4 + (R121 - R122 - R221 + R222)*y2e4) + 6.*((R111 - R121 - R211 + R221)*z1e4 + (-R112 + R122 + R212 - R222)*z2e4) + x1e2*(2.*((R111 - R112)*y1e2 + (-R121 + R122)*y2e2) + 3.*((-R111 + R121)*z1e2 + (R112 - R122)*z2e2)) + 2.*(y1e2*((R111 - R211)*z1e2 + (-R112 + R212)*z2e2) + y2e2*((-R121 + R221)*z1e2 + (R122 - R222)*z2e2)) + x2e2*(2.*((-R211 + R212)*y1e2 + (R221 - R222)*y2e2) + 3.*((R211 - R221)*z1e2 + (-R212 + R222)*z2e2)))/15.;
					T212.y = (3.*((-SumLogsX1Y1 + SumLogsX1Y2)*x1e5 + (SumLogsX2Y1 - SumLogsX2Y2)*x2e5) + 8.*(SumAtansY1*y1e5 - SumAtansY2*y2e5) + 3.*((-SumLogsY1Z1 + SumLogsY2Z1)*z1e5 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e5) + x1e3*(5.*(SumLogsX1Y1*y1e2 - SumLogsX1Y2*y2e2) + 3.*((R111 - R121)*z1 + (-R112 + R122)*z2)) + x2e3*(5.*(-SumLogsX2Y1*y1e2 + SumLogsX2Y2*y2e2) + 3.*((-R211 + R221)*z1 + (R212 - R222)*z2)) + 5.*(y1e2*(SumLogsY1Z1*z1e3 - SumLogsY1Z2*z2e3) + y2e2*(-SumLogsY2Z1*z1e3 + SumLogsY2Z2*z2e3)) + x1*(3.*((R111 - R121)*z1e3 + (-R112 + R122)*z2e3) + 8.*(y1e2*(-R111*z1 + R112*z2) + y2e2*(R121*z1 - R122*z2))) + x2*(3.*((-R211 + R221)*z1e3 + (R212 - R222)*z2e3) + 8.*(y1e2*(R211*z1 - R212*z2) + y2e2*(-R221*z1 + R222*z2))))/30.;
					S212.z = ((-SumLogsX1Z1 + SumLogsX1Z2)*x1e5 + (SumLogsX2Z1 - SumLogsX2Z2)*x2e5 + (-SumLogsY1Z1 + SumLogsY1Z2)*y1e5 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e5 + x1e3*((R111 - R112)*y1 + (-R121 + R122)*y2) + x2e3*((-R211 + R212)*y1 + (R221 - R222)*y2) + SumAtansZ1*z1e5 - SumAtansZ2*z2e5 + x1*((R111 - R112)*y1e3 + (-R121 + R122)*y2e3 + y1*(-(R111*z1e2) + R112*z2e2) + y2*(R121*z1e2 - R122*z2e2)) + x2*((-R211 + R212)*y1e3 + (R221 - R222)*y2e3 + y1*(R211*z1e2 - R212*z2e2) + y2*(-(R221*z1e2) + R222*z2e2)))/5.;
					T212.z = (6.*((-SumLogsX1Y1 + SumLogsX1Y2)*x1e5 + (SumLogsX2Y1 - SumLogsX2Y2)*x2e5) + 4.*(-SumAtansY1*y1e5 + SumAtansY2*y2e5) + 9.*((SumLogsY1Z1 - SumLogsY2Z1)*z1e5 + (-SumLogsY1Z2 + SumLogsY2Z2)*z2e5) + x1e3*(10.*(-SumLogsX1Y1*y1e2 + SumLogsX1Y2*y2e2) + 6.*((R111 - R121)*z1 + (-R112 + R122)*z2)) + x2e3*(10.*(SumLogsX2Y1*y1e2 - SumLogsX2Y2*y2e2) + 6.*((-R211 + R221)*z1 + (R212 - R222)*z2)) + 5.*(y1e2*(SumLogsY1Z1*z1e3 - SumLogsY1Z2*z2e3) + y2e2*(-SumLogsY2Z1*z1e3 + SumLogsY2Z2*z2e3)) + x1*(9.*((-R111 + R121)*z1e3 + (R112 - R122)*z2e3) + 4.*(y1e2*(R111*z1 - R112*z2) + y2e2*(-R121*z1 + R122*z2))) + x2*(9.*((R211 - R221)*z1e3 + (-R212 + R222)*z2e3) + 4.*(y1e2*(-R211*z1 + R212*z2) + y2e2*(R221*z1 - R222*z2))))/30.;

					InvDelXe2InvDelYInvDelZe2 = InvDelXe2InvDelYInvDelZ*InvDelZ;
					T212 = InvDelXe2InvDelYInvDelZe2*T212; S212 = InvDelXe2InvDelYInvDelZe2*S212;
				}
			}
			if(int(ky)>2)
			{
				if(int(kz)>0)
				{
					TVector3d& T220 = T[2][2][0];
					TVector3d& S220 = S[2][2][0];

					T220.x = (3.*(-SumAtansX1*x1e4 + SumAtansX2*x2e4 + SumAtansY1*y1e4 - SumAtansY2*y2e4) - SumAtansZ1*z1e4 + SumAtansZ2*z2e4 + 2.*(x1e3*(SumLogsX1Z1*z1 - SumLogsX1Z2*z2) + x2e3*(-SumLogsX2Z1*z1 + SumLogsX2Z2*z2)) + 4.*(y1e3*(-SumLogsY1Z1*z1 + SumLogsY1Z2*z2) + y2e3*(SumLogsY2Z1*z1 - SumLogsY2Z2*z2)) + x1*(y1*(R111*z1 - R112*z2) + y2*(-(R121*z1) + R122*z2)) + x2*(y1*(-(R211*z1) + R212*z2) + y2*(R221*z1 - R222*z2)))/6.;
					S220.x = (3.*((SumLogsX1Y1 - SumLogsX1Y2)*x1e4 + (-SumLogsX2Y1 + SumLogsX2Y2)*x2e4 + (SumLogsX1Y1 - SumLogsX2Y1)*y1e4 + (-SumLogsX1Y2 + SumLogsX2Y2)*y2e4) + 2.*((-R111 + R121 + R211 - R221)*z1e3 + (R112 - R122 - R212 + R222)*z2e3) + x1e2*((R111 - R121)*z1 + (-R112 + R122)*z2) + y1e2*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2e2*((-R121 + R221)*z1 + (R122 - R222)*z2) + x2e2*((-R211 + R221)*z1 + (R212 - R222)*z2))/6.;
					S220.y = ((SumLogsX1Z1 - SumLogsX1Z2)*x1e4 + (-SumLogsX2Z1 + SumLogsX2Z2)*x2e4 + 2.*((R111 - R112 - R211 + R212)*y1e3 + (-R121 + R122 + R221 - R222)*y2e3) + x1e2*((-R111 + R112)*y1 + (R121 - R122)*y2) + x2e2*((R211 - R212)*y1 + (-R221 + R222)*y2) + (-SumLogsX1Z1 + SumLogsX2Z1)*z1e4 + (SumLogsX1Z2 - SumLogsX2Z2)*z2e4 + y1*((R111 - R211)*z1e2 + (-R112 + R212)*z2e2) + y2*((-R121 + R221)*z1e2 + (R122 - R222)*z2e2))/4.;
					T220.y = (3.*(SumAtansX1*x1e4 - SumAtansX2*x2e4 - SumAtansY1*y1e4 + SumAtansY2*y2e4) - SumAtansZ1*z1e4 + SumAtansZ2*z2e4 + 4.*(x1e3*(-SumLogsX1Z1*z1 + SumLogsX1Z2*z2) + x2e3*(SumLogsX2Z1*z1 - SumLogsX2Z2*z2)) + 2.*(y1e3*(SumLogsY1Z1*z1 - SumLogsY1Z2*z2) + y2e3*(-SumLogsY2Z1*z1 + SumLogsY2Z2*z2)) + x1*(y1*(R111*z1 - R112*z2) + y2*(-(R121*z1) + R122*z2)) + x2*(y1*(-(R211*z1) + R212*z2) + y2*(R221*z1 - R222*z2)))/6.;
					S220.z = (2.*((R111 - R112 - R121 + R122)*x1e3 + (-R211 + R212 + R221 - R222)*x2e3) + (SumLogsY1Z1 - SumLogsY1Z2)*y1e4 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2e4 + (-SumLogsY1Z1 + SumLogsY2Z1)*z1e4 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e4 + x1*((-R111 + R112)*y1e2 + (R121 - R122)*y2e2 + (R111 - R121)*z1e2 + (-R112 + R122)*z2e2) + x2*((R211 - R212)*y1e2 + (-R221 + R222)*y2e2 + (-R211 + R221)*z1e2 + (R212 - R222)*z2e2))/4.;
					T220.z = (SumAtansZ1*z1e4 - SumAtansZ2*z2e4 + x1e3*(SumLogsX1Z1*z1 - SumLogsX1Z2*z2) + x2e3*(-(SumLogsX2Z1*z1) + SumLogsX2Z2*z2) + y1e3*(SumLogsY1Z1*z1 - SumLogsY1Z2*z2) + y2e3*(-(SumLogsY2Z1*z1) + SumLogsY2Z2*z2) + x1*(y1*(-(R111*z1) + R112*z2) + y2*(R121*z1 - R122*z2)) + x2*(y1*(R211*z1 - R212*z2) + y2*(-(R221*z1) + R222*z2)))/3.;

					InvDelXe2InvDelYe2 = InvDelXe2InvDelY*InvDelY;
					T220 = InvDelXe2InvDelYe2*T220; S220 = InvDelXe2InvDelYe2*S220;
				}
				if(int(kz)>1)
				{
					TVector3d& T221 = T[2][2][1];
					TVector3d& S221 = S[2][2][1];

					T221.x = (9.*((SumLogsX1Z1 - SumLogsX1Z2)*x1e5 + (-SumLogsX2Z1 + SumLogsX2Z2)*x2e5) + 6.*((-SumLogsY1Z1 + SumLogsY1Z2)*y1e5 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e5) + 4.*(-SumAtansZ1*z1e5 + SumAtansZ2*z2e5) + x1e3*(9.*((-R111 + R112)*y1 + (R121 - R122)*y2) + 5.*(SumLogsX1Z1*z1e2 - SumLogsX1Z2*z2e2)) + x2e3*(9.*((R211 - R212)*y1 + (-R221 + R222)*y2) + 5.*(-SumLogsX2Z1*z1e2 + SumLogsX2Z2*z2e2)) + 10.*(y1e3*(-SumLogsY1Z1*z1e2 + SumLogsY1Z2*z2e2) + y2e3*(SumLogsY2Z1*z1e2 - SumLogsY2Z2*z2e2)) + x1*(6.*((R111 - R112)*y1e3 + (-R121 + R122)*y2e3) + 4.*(y1*(R111*z1e2 - R112*z2e2) + y2*(-R121*z1e2 + R122*z2e2))) + x2*(6.*((-R211 + R212)*y1e3 + (R221 - R222)*y2e3) + 4.*(y1*(-R211*z1e2 + R212*z2e2) + y2*(R221*z1e2 - R222*z2e2))))/30.;
					S221.x = (6.*((R111 - R112 - R121 + R122)*x1e4 + (-R211 + R212 + R221 - R222)*x2e4 + (R111 - R112 - R211 + R212)*y1e4 + (-R121 + R122 + R221 - R222)*y2e4) + 4.*((-R111 + R121 + R211 - R221)*z1e4 + (R112 - R122 - R212 + R222)*z2e4) + x1e2*(3.*((-R111 + R112)*y1e2 + (R121 - R122)*y2e2) + 2.*((R111 - R121)*z1e2 + (-R112 + R122)*z2e2)) + 2.*(y1e2*((R111 - R211)*z1e2 + (-R112 + R212)*z2e2) + y2e2*((-R121 + R221)*z1e2 + (R122 - R222)*z2e2)) + x2e2*(3.*((R211 - R212)*y1e2 + (-R221 + R222)*y2e2) + 2.*((-R211 + R221)*z1e2 + (R212 - R222)*z2e2)))/15.;
					S221.y = (SumAtansX1*x1e5 - SumAtansX2*x2e5 + (-SumLogsX1Y1 + SumLogsX2Y1)*y1e5 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e5 + (-SumLogsX1Z1 + SumLogsX2Z1)*z1e5 + (SumLogsX1Z2 - SumLogsX2Z2)*z2e5 + y1e3*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2e3*((-R121 + R221)*z1 + (R122 - R222)*z2) + y1*((R111 - R211)*z1e3 + (-R112 + R212)*z2e3) + y2*((-R121 + R221)*z1e3 + (R122 - R222)*z2e3) + x1e2*(y1*(-(R111*z1) + R112*z2) + y2*(R121*z1 - R122*z2)) + x2e2*(y1*(R211*z1 - R212*z2) + y2*(-(R221*z1) + R222*z2)))/5.;
					T221.y = (6.*((-SumLogsX1Z1 + SumLogsX1Z2)*x1e5 + (SumLogsX2Z1 - SumLogsX2Z2)*x2e5) + 9.*((SumLogsY1Z1 - SumLogsY1Z2)*y1e5 + (-SumLogsY2Z1 + SumLogsY2Z2)*y2e5) + 4.*(-SumAtansZ1*z1e5 + SumAtansZ2*z2e5) + x1e3*(6.*((R111 - R112)*y1 + (-R121 + R122)*y2) + 10.*(-SumLogsX1Z1*z1e2 + SumLogsX1Z2*z2e2)) + x2e3*(6.*((-R211 + R212)*y1 + (R221 - R222)*y2) + 10.*(SumLogsX2Z1*z1e2 - SumLogsX2Z2*z2e2)) + 5.*(y1e3*(SumLogsY1Z1*z1e2 - SumLogsY1Z2*z2e2) + y2e3*(-SumLogsY2Z1*z1e2 + SumLogsY2Z2*z2e2)) + x1*(9.*((-R111 + R112)*y1e3 + (R121 - R122)*y2e3) + 4.*(y1*(R111*z1e2 - R112*z2e2) + y2*(-R121*z1e2 + R122*z2e2))) + x2*(9.*((R211 - R212)*y1e3 + (-R221 + R222)*y2e3) + 4.*(y1*(-R211*z1e2 + R212*z2e2) + y2*(R221*z1e2 - R222*z2e2))))/30.;
					S221.z = ((-SumLogsX1Y1 + SumLogsX1Y2)*x1e5 + (SumLogsX2Y1 - SumLogsX2Y2)*x2e5 + SumAtansY1*y1e5 - SumAtansY2*y2e5 + (-SumLogsY1Z1 + SumLogsY2Z1)*z1e5 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e5 + x1e3*((R111 - R121)*z1 + (-R112 + R122)*z2) + x2e3*((-R211 + R221)*z1 + (R212 - R222)*z2) + x1*((R111 - R121)*z1e3 + (-R112 + R122)*z2e3 + y1e2*(-(R111*z1) + R112*z2) + y2e2*(R121*z1 - R122*z2)) + x2*((-R211 + R221)*z1e3 + (R212 - R222)*z2e3 + y1e2*(R211*z1 - R212*z2) + y2e2*(-(R221*z1) + R222*z2)))/5.;
					T221.z = (3.*((-SumLogsX1Z1 + SumLogsX1Z2)*x1e5 + (SumLogsX2Z1 - SumLogsX2Z2)*x2e5 + (-SumLogsY1Z1 + SumLogsY1Z2)*y1e5 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e5) + 8.*(SumAtansZ1*z1e5 - SumAtansZ2*z2e5) + x1e3*(3.*((R111 - R112)*y1 + (-R121 + R122)*y2) + 5.*(SumLogsX1Z1*z1e2 - SumLogsX1Z2*z2e2)) + x2e3*(3.*((-R211 + R212)*y1 + (R221 - R222)*y2) + 5.*(-SumLogsX2Z1*z1e2 + SumLogsX2Z2*z2e2)) + 5.*(y1e3*(SumLogsY1Z1*z1e2 - SumLogsY1Z2*z2e2) + y2e3*(-SumLogsY2Z1*z1e2 + SumLogsY2Z2*z2e2)) + x1*(3.*((R111 - R112)*y1e3 + (-R121 + R122)*y2e3) + 8.*(y1*(-R111*z1e2 + R112*z2e2) + y2*(R121*z1e2 - R122*z2e2))) + x2*(3.*((-R211 + R212)*y1e3 + (R221 - R222)*y2e3) + 8.*(y1*(R211*z1e2 - R212*z2e2) + y2*(-R221*z1e2 + R222*z2e2))))/30.;

					InvDelXe2InvDelYe2InvDelZ = InvDelXe2InvDelYe2*InvDelZ;
					T221 = InvDelXe2InvDelYe2InvDelZ*T221; S221 = InvDelXe2InvDelYe2InvDelZ*S221;
				}
				if(int(kz)>2)
				{
					TVector3d& T222 = T[2][2][2];
					TVector3d& S222 = S[2][2][2];

					T222.x = (2.*(SumAtansX1*x1e6 - SumAtansX2*x2e6) - SumAtansY1*y1e6 + SumAtansY2*y2e6 - SumAtansZ1*z1e6 + SumAtansZ2*z2e6 + 2.*(y1e3*(-SumLogsY1Z1*z1e3 + SumLogsY1Z2*z2e3) + y2e3*(SumLogsY2Z1*z1e3 - SumLogsY2Z2*z2e3)) + x1e3*(SumLogsX1Y1*y1e3 - SumLogsX1Y2*y2e3 + SumLogsX1Z1*z1e3 - SumLogsX1Z2*z2e3 + 2.*(y1*(-R111*z1 + R112*z2) + y2*(R121*z1 - R122*z2))) + x2e3*(-(SumLogsX2Y1*y1e3) + SumLogsX2Y2*y2e3 - SumLogsX2Z1*z1e3 + SumLogsX2Z2*z2e3 + 2.*(y1*(R211*z1 - R212*z2) + y2*(-R221*z1 + R222*z2))) + x1*(y1e3*(R111*z1 - R112*z2) + y2e3*(-(R121*z1) + R122*z2) + y1*(R111*z1e3 - R112*z2e3) + y2*(-(R121*z1e3) + R122*z2e3)) + x2*(y1e3*(-(R211*z1) + R212*z2) + y2e3*(R221*z1 - R222*z2) + y1*(-(R211*z1e3) + R212*z2e3) + y2*(R221*z1e3 - R222*z2e3)))/9.;
					S222.x = (3.*((-SumLogsX1Y1 + SumLogsX1Y2)*x1e6 + (SumLogsX2Y1 - SumLogsX2Y2)*x2e6 + (-SumLogsX1Y1 + SumLogsX2Y1)*y1e6 + (SumLogsX1Y2 - SumLogsX2Y2)*y2e6) + 4.*((-R111 + R121 + R211 - R221)*z1e5 + (R112 - R122 - R212 + R222)*z2e5) + 3.*(x1e4*((R111 - R121)*z1 + (-R112 + R122)*z2) + y1e4*((R111 - R211)*z1 + (-R112 + R212)*z2) + y2e4*((-R121 + R221)*z1 + (R122 - R222)*z2) + x2e4*((-R211 + R221)*z1 + (R212 - R222)*z2)) + 2.*(y1e2*((R111 - R211)*z1e3 + (-R112 + R212)*z2e3) + y2e2*((-R121 + R221)*z1e3 + (R122 - R222)*z2e3)) + x1e2*(2.*((R111 - R121)*z1e3 + (-R112 + R122)*z2e3) + 3.*(y1e2*(-R111*z1 + R112*z2) + y2e2*(R121*z1 - R122*z2))) + x2e2*(2.*((-R211 + R221)*z1e3 + (R212 - R222)*z2e3) + 3.*(y1e2*(R211*z1 - R212*z2) + y2e2*(-R221*z1 + R222*z2))))/18.;
					S222.y = (3.*((-SumLogsX1Z1 + SumLogsX1Z2)*x1e6 + (SumLogsX2Z1 - SumLogsX2Z2)*x2e6) + 4.*((-R111 + R112 + R211 - R212)*y1e5 + (R121 - R122 - R221 + R222)*y2e5) + 3.*(x1e4*((R111 - R112)*y1 + (-R121 + R122)*y2) + x2e4*((-R211 + R212)*y1 + (R221 - R222)*y2) + (-SumLogsX1Z1 + SumLogsX2Z1)*z1e6 + (SumLogsX1Z2 - SumLogsX2Z2)*z2e6) + 2.*(y1e3*((R111 - R211)*z1e2 + (-R112 + R212)*z2e2) + y2e3*((-R121 + R221)*z1e2 + (R122 - R222)*z2e2)) + 3.*(y1*((R111 - R211)*z1e4 + (-R112 + R212)*z2e4) + y2*((-R121 + R221)*z1e4 + (R122 - R222)*z2e4)) + x1e2*(2.*((R111 - R112)*y1e3 + (-R121 + R122)*y2e3) + 3.*(y1*(-R111*z1e2 + R112*z2e2) + y2*(R121*z1e2 - R122*z2e2))) + x2e2*(2.*((-R211 + R212)*y1e3 + (R221 - R222)*y2e3) + 3.*(y1*(R211*z1e2 - R212*z2e2) + y2*(-R221*z1e2 + R222*z2e2))))/18.;
					T222.y = (-(SumAtansX1*x1e6) + SumAtansX2*x2e6 + 2.*(SumAtansY1*y1e6 - SumAtansY2*y2e6) - SumAtansZ1*z1e6 + SumAtansZ2*z2e6 + y1e3*(SumLogsY1Z1*z1e3 - SumLogsY1Z2*z2e3) + y2e3*(-(SumLogsY2Z1*z1e3) + SumLogsY2Z2*z2e3) + x1e3*(SumLogsX1Y1*y1e3 - SumLogsX1Y2*y2e3 + 2.*(-SumLogsX1Z1*z1e3 + SumLogsX1Z2*z2e3) + y1*(R111*z1 - R112*z2) + y2*(-(R121*z1) + R122*z2)) + x2e3*(-(SumLogsX2Y1*y1e3) + SumLogsX2Y2*y2e3 + 2.*(SumLogsX2Z1*z1e3 - SumLogsX2Z2*z2e3) + y1*(-(R211*z1) + R212*z2) + y2*(R221*z1 - R222*z2)) + x1*(2.*(y1e3*(-R111*z1 + R112*z2) + y2e3*(R121*z1 - R122*z2)) + y1*(R111*z1e3 - R112*z2e3) + y2*(-(R121*z1e3) + R122*z2e3)) + x2*(2.*(y1e3*(R211*z1 - R212*z2) + y2e3*(-R221*z1 + R222*z2)) + y1*(-(R211*z1e3) + R212*z2e3) + y2*(R221*z1e3 - R222*z2e3)))/9.;
					S222.z = (4.*((-R111 + R112 + R121 - R122)*x1e5 + (R211 - R212 - R221 + R222)*x2e5) + 3.*((-SumLogsY1Z1 + SumLogsY1Z2)*y1e6 + (SumLogsY2Z1 - SumLogsY2Z2)*y2e6 + (-SumLogsY1Z1 + SumLogsY2Z1)*z1e6 + (SumLogsY1Z2 - SumLogsY2Z2)*z2e6) + 2.*(x1e3*((R111 - R112)*y1e2 + (-R121 + R122)*y2e2 + (R111 - R121)*z1e2 + (-R112 + R122)*z2e2) + x2e3*((-R211 + R212)*y1e2 + (R221 - R222)*y2e2 + (-R211 + R221)*z1e2 + (R212 - R222)*z2e2)) + 3.*(x1*((R111 - R112)*y1e4 + (-R121 + R122)*y2e4 + (R111 - R121)*z1e4 + (-R112 + R122)*z2e4 + y1e2*(-R111*z1e2 + R112*z2e2) + y2e2*(R121*z1e2 - R122*z2e2)) + x2*((-R211 + R212)*y1e4 + (R221 - R222)*y2e4 + (-R211 + R221)*z1e4 + (R212 - R222)*z2e4 + y1e2*(R211*z1e2 - R212*z2e2) + y2e2*(-R221*z1e2 + R222*z2e2))))/18.;
					T222.z = (-(SumAtansX1*x1e6) + SumAtansX2*x2e6 - SumAtansY1*y1e6 + SumAtansY2*y2e6 + 2.*(SumAtansZ1*z1e6 - SumAtansZ2*z2e6) + y1e3*(SumLogsY1Z1*z1e3 - SumLogsY1Z2*z2e3) + y2e3*(-(SumLogsY2Z1*z1e3) + SumLogsY2Z2*z2e3) + x1e3*(2.*(-SumLogsX1Y1*y1e3 + SumLogsX1Y2*y2e3) + SumLogsX1Z1*z1e3 - SumLogsX1Z2*z2e3 + y1*(R111*z1 - R112*z2) + y2*(-(R121*z1) + R122*z2)) + x2e3*(2.*(SumLogsX2Y1*y1e3 - SumLogsX2Y2*y2e3) - SumLogsX2Z1*z1e3 + SumLogsX2Z2*z2e3 + y1*(-(R211*z1) + R212*z2) + y2*(R221*z1 - R222*z2)) + x1*(y1e3*(R111*z1 - R112*z2) + y2e3*(-(R121*z1) + R122*z2) + 2.*(y1*(-R111*z1e3 + R112*z2e3) + y2*(R121*z1e3 - R122*z2e3))) + x2*(y1e3*(-(R211*z1) + R212*z2) + y2e3*(R221*z1 - R222*z2) + 2.*(y1*(R211*z1e3 - R212*z2e3) + y2*(-R221*z1e3 + R222*z2e3))))/9.;

					InvDelXe2InvDelYe2InvDelZe2 = InvDelXe2InvDelYe2InvDelZ*InvDelZ;
					T222 = InvDelXe2InvDelYe2InvDelZe2*T222; S222 = InvDelXe2InvDelYe2InvDelZe2*S222;
				}
			}
		}
	}

	double Xpn = -CenPo_mi_P.x*InvDelX;
	double Ypn = -CenPo_mi_P.y*InvDelY;
	double Zpn = -CenPo_mi_P.z*InvDelZ;
	
	double XpnMlt[3], YpnMlt[3], ZpnMlt[3]; // Make larger or create in heap later
	XpnMlt[0] = YpnMlt[0] = ZpnMlt[0] = 1;
	for(int iix=1; iix<kx; iix++) XpnMlt[iix] = XpnMlt[iix-1]*Xpn;
	for(int iiy=1; iiy<ky; iiy++) YpnMlt[iiy] = YpnMlt[iiy-1]*Ypn;
	for(int iiz=1; iiz<kz; iiz++) ZpnMlt[iiz] = ZpnMlt[iiz-1]*Zpn;

	const double C0[] = {1., 0., 0., 0.}; // First String of matrix of Binomial coefficients C[n][m]; C[0][1] not determined
	const double C1[] = {1., 1., 0., 0.}; // C(1,m)
	const double C2[] = {1., 2., 1., 0.}; // C(2,m)
	const double C3[] = {1., 3., 3., 1.}; // C(3,m)
	const double* C[] = {C0, C1, C2, C3}; // C[n][m] - Binomial coefficient (n,m)

	TVector3d TmQ, SmQ, BufH(0.,0.,0.), BufM(0.,0.,0.), ZeroVect(0.,0.,0.);
	double Cm;

	//radTCast Cast;

	radTmhg::const_iterator MapIter = GroupMapOfHandlers.begin();

	int kykz = int(ky)*int(kz);
	for(int m=0; m<AmOfSubElem; m++)
	{
		TmQ = SmQ = ZeroVect;
		Cm = 0.; int XqStNoForM, YqStNoForM;
		for(int ix=0; ix<int(kx); ix++)
		{
			if(MagnCompNeeded) XqStNoForM = ix*kykz;
			for(int iy=0; iy<int(ky); iy++)
			{
				if(MagnCompNeeded) YqStNoForM = iy*int(kz);
				for(int iz=0; iz<int(kz); iz++)
				{
					if(FieldCompNeeded)
					{
						double rm = 0., xBuf_rm, yBuf_rm, zBuf_rm;
						for(int lx=ix; lx<int(kx); lx++)
						{
							xBuf_rm = C[lx][ix]*XpnMlt[lx-ix];
							int qStNoX = lx*kykz;
							for(int ly=iy; ly<int(ky); ly++)
							{
								yBuf_rm = C[ly][iy]*YpnMlt[ly-iy];
								int qStNoY = ly*int(kz);
								for(int lz=iz; lz<int(kz); lz++)
								{
									zBuf_rm = C[lz][iz]*ZpnMlt[lz-iz];

									double TestBuf = Q_forM[qStNoX+qStNoY+lz][m];
									rm += xBuf_rm*yBuf_rm*zBuf_rm * Q_forM[qStNoX+qStNoY+lz][m]; 
								}
							}
						}
						TmQ += rm * T[ix][iy][iz];
						SmQ += rm * S[ix][iy][iz];
					}
					if(MagnCompNeeded)
					{
						Cm += XpnMlt[ix]*YpnMlt[iy]*ZpnMlt[iz] * Q_forM[XqStNoForM+YqStNoForM+iz][m];
					}
				}
			}
		}

		if(FieldPtr->FieldKey.PreRelax_)
		{
			radTVectPtrVect3d::iterator IterVectOfTs = VectOfTsComputed.end();
			radTVectPtrVect3d::iterator IterVectOfSs = VectOfSsComputed.end();

			(*(--IterVectOfTs))[m] = ConForH*TmQ;
			((*(--IterVectOfSs))[m]).x = -ConForH*SmQ.z;
			((*(IterVectOfSs))[m]).y = -ConForH*SmQ.y;
			((*(IterVectOfSs))[m]).z = -ConForH*SmQ.x;
		}
		else if(B_orH_CompNeeded || MagnCompNeeded)
		{
			radTg3d* g3dPtr = (radTg3d*)(((*MapIter).second).rep);
			radTGroup* GroupPtr = radTCast::GroupCast(g3dPtr);
			TVector3d& LocMagn = (GroupPtr == 0)? ((radTg3dRelax*)g3dPtr)->Magn : ((radTSubdividedRecMag*)GroupPtr)->Magn;

			++MapIter;

			if(MagnCompNeeded) BufM += Cm*LocMagn;
			if(B_orH_CompNeeded)
			{
				TVector3d Str0(TmQ.x,SmQ.x,SmQ.y), Str1(SmQ.x,TmQ.y,SmQ.z), Str2(SmQ.y,SmQ.z,TmQ.z);
				TMatrix3d Qm(Str0, Str1, Str2);
				BufH += Qm*LocMagn;
			}
		}
	}

	if(B_orH_CompNeeded)
	{
		BufH = ConForH*BufH;
		if(FieldPtr->FieldKey.B_)
		{
			if(InsideBlock) FieldPtr->B += BufH + BufM;
			else FieldPtr->B += BufH;
		}
		if(FieldPtr->FieldKey.H_) FieldPtr->H += BufH;
	}
	if(FieldPtr->FieldKey.M_) if(InsideBlock) FieldPtr->M += BufM;
}

//-------------------------------------------------------------------------

int radTSubdividedRecMag::SetupFldCmpData(short InFldCmpMeth, int SubLevel)
{
// Modify this: take into account different element sizes
	//radTCast Cast;
	int Probe = 0;
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
	{
		radTGroup* GroupPtr = radTCast::GroupCast((radTg3d*)(((*iter).second).rep));
		if(GroupPtr!=0)
		{
			radTSubdividedRecMag* SubdividedRecMagPtr = radTCast::SubdividedRecMagCast(GroupPtr);
			if(SubdividedRecMagPtr!=0)
			{
				int LocProbe = SubdividedRecMagPtr->SetupFldCmpData(InFldCmpMeth, SubLevel);
				if(LocProbe<0) return LocProbe;
				Probe = LocProbe;
			}
		}
	}
	if(Probe!=0) return Probe;

	FldCmpMeth = InFldCmpMeth;

	if((FldCmpMeth == 1) && (int(kx) > 3 || int(ky) > 3 || int(kz) > 3)) return -38;

	if(FldCmpMeth==0)
	{
		if(Q_forM!=NULL)
		{
			for(int j=0; j<AmOfSubElem; j++) delete[] (Q_forM[j]);
			delete[] Q_forM;
		}
		if(FormCenPoPtrArray!=NULL) delete[] FormCenPoPtrArray;

		return 1;
	}

	Q_forM = new double*[AmOfSubElem];
	double** DirQ = new double*[AmOfSubElem];
	for(int j=0; j<AmOfSubElem; j++)
	{
		Q_forM[j] = new double[AmOfSubElem];
		DirQ[j] = new double[AmOfSubElem];
	}
	FormCenPoPtrArray = new TVector3d*[AmOfSubElem];

	if(FldCmpMeth==1)
	{
		//Defining direct matrix
		TVector3d HalfDim = 0.5*Dimensions;

		const double AbsZeroTol = 5.E-13;
		double q0x = pow(qx, 1./(kx-1.)), q0y = pow(qy, 1./(ky-1.)), q0z = pow(qz, 1./(kz-1.));
		double BufX = qx*q0x - 1., BufY = qy*q0y - 1., BufZ = qz*q0z - 1.;
		double ax = (fabs(BufX) > AbsZeroTol)? Dimensions.x*(q0x - 1.)/BufX : Dimensions.x/kx;
		double ay = (fabs(BufY) > AbsZeroTol)? Dimensions.y*(q0y - 1.)/BufY : Dimensions.y/ky;
		double az = (fabs(BufZ) > AbsZeroTol)? Dimensions.z*(q0z - 1.)/BufZ : Dimensions.z/kz;

		double xBeg = -HalfDim.x + 0.5*ax;
		double yBeg = -HalfDim.y + 0.5*ay;
		double zBeg = -HalfDim.z + 0.5*az;

		double Xj, Yj, Zj, MultCx, MultCy, MultCz;

		double* Cx = new double[int(kx)];
		double* Cy = new double[int(ky)];
		double* Cz = new double[int(kz)];

		int ColNo, StrNo, ix, iy, iz;
		StrNo = 0;
		Xj = xBeg;
		for(int jx=0; jx<kx; jx++)
		{
			MultCx = Xj; Cx[0] = 1.;

			for(ix=1; ix<kx; ix++) Cx[ix] = Cx[ix-1]*MultCx;

			Yj = yBeg;
			for(int jy=0; jy<ky; jy++)
			{
				MultCy = Yj; Cy[0] = 1.;

				for(iy=1; iy<ky; iy++) Cy[iy] = Cy[iy-1]*MultCy;

				Zj = zBeg;
				for(int jz=0; jz<kz; jz++)
				{
					MultCz = Zj; Cz[0] = 1.;

					for(iz=1; iz<kz; iz++) Cz[iz] = Cz[iz-1]*MultCz;

					ColNo = 0;
					for(ix=0; ix<kx; ix++)
						for(iy=0; iy<ky; iy++)
							for(iz=0; iz<kz; iz++)
							{
								double TestBuf = Cx[ix]*Cy[iy]*Cz[iz];
								DirQ[StrNo][ColNo++] = TestBuf;
							}
					StrNo++;

					Zj += 0.5*az;
					az *= q0z;
					Zj += 0.5*az;
				}
				Yj += 0.5*ay;
				ay *= q0y;
				Yj += 0.5*ay;
			}
			Xj += 0.5*ax;
			ax *= q0x;
			Xj += 0.5*ax;
		}
		delete[] Cx; delete[] Cy; delete[] Cz;

		radTMathLinAlgEq MathMet(AmOfSubElem);
		MathMet.InverseMatrix(DirQ, AmOfSubElem, Q_forM);

		StrNo = 0;
		for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
			iter != GroupMapOfHandlers.end(); ++iter)
		{
			radTg3d* g3dPtr = (radTg3d*)(((*iter).second).rep);
			radTGroup* GroupPtr = radTCast::GroupCast(g3dPtr);
			//FormCenPoPtrArray[StrNo++] = (GroupPtr == 0)? &(((radTg3dRelax*)g3dPtr)->CentrPoint) : &(((radTSubdividedRecMag*)GroupPtr)->CentrPoint);
			FormCenPoPtrArray[StrNo++] = (GroupPtr == 0)? &(((radTg3dRelax*)g3dPtr)->CentrPoint) : &(GroupPtr->CentrPoint); //OC061008
		}
	}
	
	for(int jj=0; jj<AmOfSubElem; jj++) delete[] (DirQ[jj]);
	delete[] DirQ;

	return 1;
}

//-------------------------------------------------------------------------

radTg3dRelax* radTSubdividedRecMag::FormalIntrctMemberPtr()
{
	if(FldCmpMeth==1)
	{
		if((FormIntrctMembCounter == AmOfSubElem) || (FormIntrctMembCounter == -1))
		{
			FormIntrctMembIter = GroupMapOfHandlers.begin();
			FormIntrctMembCounter = -1; // ??
		}

		radTg3dRelax* FormIntrctMembPtr = NULL;
		radTg3d* g3dPtr = (radTg3d*)(((*FormIntrctMembIter).second).rep);

		//radTCast Cast;
		radTGroup* GroupPtr = radTCast::GroupCast(g3dPtr);
		if(GroupPtr==0) FormIntrctMembPtr = (radTg3dRelax*)g3dPtr;
		else FormIntrctMembPtr = (radTg3dRelax*)((radTSubdividedRecMag*)GroupPtr);

		++FormIntrctMembIter;
		++FormIntrctMembCounter;

		return FormIntrctMembPtr;
	}
	else return this;
}

//-------------------------------------------------------------------------

void radTSubdividedRecMag::Dump(std::ostream& o, int ShortSign)
{
	((radTg3d*)((radTGroup*)this))->radTg3d::Dump(o, ShortSign);

	o << "Subdivided RecMag";

	if(ShortSign==1) return;

	o << endl;
	radTRecMag::DumpPureObjInfo(o, ShortSign);
	DumpMaterApplied(o);

	o << endl;
	radTGroup::DumpPureObjInfo(o, ShortSign);

	o << endl;
	((radTg3d*)((radTGroup*)this))->radTg3d::DumpTransApplied(o);

	o << endl;
	o << "   Memory occupied (incl. the content): " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

void radTSubdividedRecMag::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	//Dumping objects that may be used by this object
	vector<pair<int, int> > vTrfKeys;
	//DumpBin_g3d_TreatTrfs(oStr, mEl, vTrfKeys);
	radTGroup::DumpBin_g3d_TreatTrfs(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vTrfKeys); //all G3D members should be accessed via radTGroup
	
	vector<int> vGroupMemKeys;
	DumpBin_Group_TreatMembers(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, vGroupMemKeys);

	int matKey=0;
	DumpBin_g3dRelax_TreatMat(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, matKey);
	
	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)radTGroup::Type_g();
	oStr << (char)radTGroup::Type_g3d();
	oStr << (char)Type_Group();
	oStr << (char)0;
	oStr << (char)0;

	//Members of radTg3d
	//((radTRecMag*)this)->DumpBin_g3d(oStr, vTrfKeys);
	radTGroup::DumpBin_g3d(oStr, vTrfKeys);

	//Members of radTGroup
	DumpBin_Group_OutMemKeys(oStr, vGroupMemKeys);

	//Members of radTg3dRelax
	DumpBin_g3dRelax(oStr, matKey);

	//Members of radTRecMag
	DumpBin_RecMag(oStr);

	//Members of radTSubdividedRecMag
	//double** Q_forM; - not included
	//TVector3d** FormCenPoPtrArray; - not included
	
	//int AmOfSubElem;
	oStr << AmOfSubElem;

	//int CenPoLoopCounter, IntrcMatrConstrCounter, FormIntrctMembCounter;
	oStr << CenPoLoopCounter << IntrcMatrConstrCounter << FormIntrctMembCounter;

	//radTmhg::const_iterator FormIntrctMembIter; - not included
	//radTVectPtrVect3d VectOfTsComputed, VectOfSsComputed; - not included

	//short GroupMembersArePureRecMags;
	oStr << GroupMembersArePureRecMags;

	//int kx, ky, kz;
	oStr << kx << ky << kz;

	//double qx, qy, qz;
	oStr << qx << qy << qz;

	//short FldCmpMeth;
	oStr << FldCmpMeth;

	//bool AlgsBasedOnKsQsMayNotWork;
	oStr << AlgsBasedOnKsQsMayNotWork;
}

//-------------------------------------------------------------------------

radTSubdividedRecMag::radTSubdividedRecMag(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
{
	//Members of radTg3d 
	//((radTRecMag*)this)->DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);
	radTGroup::DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers); //all G3D members should be accessed via radTGroup

	//Members of radTGroup 
	DumpBinParse_Group(inStr, mKeysOldNew, gMapOfHandlers);

	//Members of radTg3dRelax 
	DumpBinParse_g3dRelax(inStr, mKeysOldNew, gMapOfHandlers);

	//Members of radTRecMag 
	DumpBinParse_RecMag(inStr);

	//Members of radTSubdividedRecMag
	//double** Q_forM; - not included
	//TVector3d** FormCenPoPtrArray; - not included

	//int AmOfSubElem;
	inStr >> AmOfSubElem;

	//int CenPoLoopCounter, IntrcMatrConstrCounter, FormIntrctMembCounter;
	inStr >> CenPoLoopCounter;
	inStr >> IntrcMatrConstrCounter;
	inStr >> FormIntrctMembCounter;

	//radTmhg::const_iterator FormIntrctMembIter; - not included
	//radTVectPtrVect3d VectOfTsComputed, VectOfSsComputed; - not included

	//short GroupMembersArePureRecMags;
	inStr >> GroupMembersArePureRecMags;

	//int kx, ky, kz;
	inStr >> kx;
	inStr >> ky;
	inStr >> kz;

	//double qx, qy, qz;
	inStr >> qx;
	inStr >> qy;
	inStr >> qz;

	//short FldCmpMeth;
	inStr >> FldCmpMeth;

	//bool AlgsBasedOnKsQsMayNotWork;
	inStr >> AlgsBasedOnKsQsMayNotWork;
}

//-------------------------------------------------------------------------
