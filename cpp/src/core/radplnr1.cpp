/*-------------------------------------------------------------------------
*
* File name:      radplnr1.cpp
*
* Project:        RADIA
*
* Description:    Auxiliary 2D objects
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radplnr.h"
#include "radg3dgr.h"
#include "radappl.h"

//-------------------------------------------------------------------------

extern radTConvergRepair& radCR;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRectangle::IntOverSurf(radTField* FieldPtr)
{
	int LenVal = FieldPtr->ShapeIntDataPtr->IntegrandLength;
	int LenValp2 = LenVal+2;
	TVector3d ZeroVect(0.,0.,0.);

	TVector3d* InnerIntegVal[6]; 
	TVector3d* OuterIntegVal[6];

	short* InnerElemCompNotFinished;
	short* OuterElemCompNotFinished;
	double* InnerAbsPrecAndLimitsArray;
	double* OuterAbsPrecAndLimitsArray;
	TVector3d* LocalVectArray;
	//radTSend Send;
	//try
	//{
		int j;
		for(j=0; j<6; j++)
		{
			InnerIntegVal[j] = new TVector3d[LenVal];
			OuterIntegVal[j] = new TVector3d[LenVal];
		}
		InnerElemCompNotFinished = new short[LenVal];
		OuterElemCompNotFinished = new short[LenVal];
		InnerAbsPrecAndLimitsArray = new double[LenValp2];
		OuterAbsPrecAndLimitsArray = new double[LenValp2];
		LocalVectArray = new TVector3d[LenVal];

		SurfIntDataPtr = new radTRectangleSurfIntData();
	//}
	//catch (radTException* radExceptionPtr)
	//{
	//	Send.ErrorMessage(radExceptionPtr->what()); return;
	//}
	//catch (...)
	//{
	//	Send.ErrorMessage("Radia::Error999"); return;
	//}

	SurfIntDataPtr->IntegrandLen = LenVal;
	SurfIntDataPtr->IntegrandFunPtr = FieldPtr->ShapeIntDataPtr->IntegrandFunPtr;
	SurfIntDataPtr->InnerAbsPrecAndLimitsArray = InnerAbsPrecAndLimitsArray;
	SurfIntDataPtr->InnerElemCompNotFinished = InnerElemCompNotFinished;
	SurfIntDataPtr->InnerIntegVal = InnerIntegVal;
	
	SurfIntDataPtr->Field = *FieldPtr;
	radTStructForShapeInt LocShapeIntData = *(FieldPtr->ShapeIntDataPtr);

	TVector3d* OutVectArray = FieldPtr->ShapeIntDataPtr->VectArray;

	LocShapeIntData.VectArray = LocalVectArray;
	SurfIntDataPtr->Field.ShapeIntDataPtr = &LocShapeIntData;

	double SmallPositive = 1.E-10;

	int i;
	for(i=0; i<LenVal; i++)
		OuterAbsPrecAndLimitsArray[i] = (FieldPtr->ShapeIntDataPtr->AbsPrecArray)[i];

	OuterAbsPrecAndLimitsArray[LenVal] = CentrPoint.y - 0.5*Dimensions.y + SmallPositive;
	OuterAbsPrecAndLimitsArray[LenVal+1] = CentrPoint.y + 0.5*Dimensions.y;

	SurfIntDataPtr->PointOnSurface.z = CentrPoint.z + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(0.,0.,1.);
	FormalOneFoldInteg(this, &radTRectangle::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];

// Delete everything created above:
	delete[] InnerElemCompNotFinished;
	delete[] OuterElemCompNotFinished;
	delete[] InnerAbsPrecAndLimitsArray;
	delete[] OuterAbsPrecAndLimitsArray;
	for(j=0; j<6; j++)
	{
		delete[] InnerIntegVal[j];
		delete[] OuterIntegVal[j];
	}
	delete[] LocalVectArray;
	delete SurfIntDataPtr;
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTRectangle::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTRectangleGraphPresent(this);
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTPolygon::radTPolygon(TVector2d* InEdgePointsArray, int InAmOfEdgePoints)
{
	CoordZ =0.; 
	AmOfEdgePoints = InAmOfEdgePoints;

	CheckAndRearrangeEdgePoints(InEdgePointsArray, InAmOfEdgePoints); // This seems to be needed only here

	for(int i=0; i<AmOfEdgePoints; i++) EdgePointsVector.push_back(InEdgePointsArray[i]);

	radTSend Send;
	if(!CheckIfNotSelfIntersecting()) SomethingIsWrong = 1;
	else SomethingIsWrong = 0;
	if(SomethingIsWrong) return;

	IsConvex = CheckIfConvex();
	if(!IsConvex) Send.WarningMessage("Radia::Warning011");

	short InsideBlock=1;
	char PointsAreDifferent=1;
	//SimpleComputeCentrPoint(InsideBlock); //OC 220902
	SimpleComputeCentrPoint(InsideBlock, PointsAreDifferent); // Test!

	if(!InsideBlock) 
	{
		Send.WarningMessage("Radia::Warning010");
	}
	//if(!PointsAreDifferent) Send.WarningMessage("Radia::Warning016"); //OC 220902
}

//-------------------------------------------------------------------------

radTPolygon::radTPolygon(const radTVect2dVect& InEdgePointsVector)
{
	CoordZ =0.; 
	AmOfEdgePoints = (int)(InEdgePointsVector.size());
	for(int i=0; i<AmOfEdgePoints; i++) EdgePointsVector.push_back(InEdgePointsVector[i]);

	radTSend Send;
	IsConvex = CheckIfConvex();
	if(!IsConvex) Send.WarningMessage("Radia::Warning011");

	short InsideBlock=1;
	char PointsAreDifferent=1;
	//SimpleComputeCentrPoint(InsideBlock); //OC 220902
	SimpleComputeCentrPoint(InsideBlock, PointsAreDifferent); // Test!

	if(!InsideBlock) Send.WarningMessage("Radia::Warning010");
	//if(!PointsAreDifferent) Send.WarningMessage("Radia::Warning016"); //OC 220902

	SomethingIsWrong = 0;
}

//-------------------------------------------------------------------------

radTPolygon::radTPolygon(double InCoordZ, TVector2d* InEdgePointsArray, int InAmOfEdgePoints, const TVector3d& InMagn)
{
	CoordZ = InCoordZ; Magn = InMagn;
	AmOfEdgePoints = InAmOfEdgePoints;

	CheckAndRearrangeEdgePoints(InEdgePointsArray, InAmOfEdgePoints);

	for(int i=0; i<AmOfEdgePoints; i++) EdgePointsVector.push_back(InEdgePointsArray[i]);

	radTSend Send;
	if(!CheckIfNotSelfIntersecting()) SomethingIsWrong = 1;
	else SomethingIsWrong = 0;
	if(SomethingIsWrong) return;

	IsConvex = CheckIfConvex();
	if(!IsConvex) Send.WarningMessage("Radia::Warning011");

	short InsideBlock=1;
	char PointsAreDifferent=1;
	//SimpleComputeCentrPoint(InsideBlock); //OC 220902
	SimpleComputeCentrPoint(InsideBlock, PointsAreDifferent); // Test!

	if(!InsideBlock) Send.WarningMessage("Radia::Warning010");
	//if(!PointsAreDifferent) Send.WarningMessage("Radia::Warning016"); //OC 220902
}

//-------------------------------------------------------------------------

radTPolygon::radTPolygon(radTVect2dVect& InEdgePointsVector, double InZ, const TVector3d& InMagn)
{
	CoordZ = InZ; 
	AmOfEdgePoints = (int)(InEdgePointsVector.size());

	CheckAndRearrangeEdgePoints(InEdgePointsVector);

	for(int i=0; i<AmOfEdgePoints; i++) EdgePointsVector.push_back(InEdgePointsVector[i]);

	Magn = InMagn;

	radTSend Send;
	if(!CheckIfNotSelfIntersecting()) SomethingIsWrong = 1;
	else SomethingIsWrong = 0;
	if(SomethingIsWrong) return;

	IsConvex = CheckIfConvex();
	if(!IsConvex) Send.WarningMessage("Radia::Warning011");

	short InsideBlock=1;
	char PointsAreDifferent=1;
	//SimpleComputeCentrPoint(InsideBlock); //OC 220902
	SimpleComputeCentrPoint(InsideBlock, PointsAreDifferent); // Test!

	if(!InsideBlock) Send.WarningMessage("Radia::Warning010");
	//if(!PointsAreDifferent) Send.WarningMessage("Radia::Warning016"); //OC 220902

	SomethingIsWrong = 0;
}

//-------------------------------------------------------------------------

radTPolygon::radTPolygon(CAuxBinStrVect& inStr)
{
	//int AmOfEdgePoints;
	inStr >> AmOfEdgePoints;

	//radTVect2dVect EdgePointsVector;
	TVector2d p2d;
	for(int i=0; i<AmOfEdgePoints; i++)
	{
		inStr >> p2d;
		EdgePointsVector.push_back(p2d);
	}

	//TVector2d CentrPoint;
	inStr >> radTPolygon::CentrPoint;

	//double CoordZ;
	inStr >> CoordZ;

	//TVector3d Magn;
	inStr >> Magn;

	//short SomethingIsWrong;
	inStr >> SomethingIsWrong;

	//char IsConvex;
	inStr >> IsConvex;
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTPolygon::CreateGraphPresent()
{
	radTg3dGraphPresent* g3dGraphPresentPtr = new radTPolygonGraphPresent(this);
	return g3dGraphPresentPtr;
}

//-------------------------------------------------------------------------

void radTPolygon::IntrsctOfTwoLines(const TVector2d& V1, const TVector2d& R01, const TVector2d& R02, const TVector2d& R12, TVector2d& IntrsctPo, TLinesIntrsctCase& IntrsctCase)
{// This is used at subdivision
	const double t_Toler = 5.E-13;
	const double V_Toler = 5.E-13;

	TVector2d V2 = R12 - R02;

	double AbsV2x = fabs(V2.x), AbsV2y = fabs(V2.y);
	double MaxR = (AbsV2x > AbsV2y)? AbsV2x : AbsV2y;
	double D_Toler = MaxR*V_Toler;

	double V1yV2x = V1.y*V2.x;
	double V1xV2y = V1.x*V2.y;
	double D = V1xV2y - V1yV2x;

	if(fabs(D) > D_Toler)
	{
		IntrsctPo.x = -(-R01.y*V1.x*V2.x + R02.y*V1.x*V2.x + R01.x*V1yV2x - R02.x*V1xV2y)/D;
		IntrsctPo.y = -(R02.y*V1yV2x - R01.y*V1xV2y + R01.x*V1.y*V2.y - R02.x*V1.y*V2.y)/D;
		double t_Intrsct = (Abs(V2.x) > V_Toler)? (IntrsctPo.x - R02.x)/V2.x : (IntrsctPo.y - R02.y)/V2.y;
		IntrsctCase = ((t_Intrsct > t_Toler) && (t_Intrsct + t_Toler < 1.))? PointWithinBound : (((Abs(t_Intrsct) < t_Toler) || (Abs(t_Intrsct-1.) < t_Toler))? PointOnBoundEdge : PointOutsideBound);
	}
	else 
	{
		double LineCoinsToler = MaxR*MaxR*V_Toler;
		IntrsctCase = (fabs(V2.x*(R01.y-R02.y) - V2.y*(R01.x-R02.x)) < LineCoinsToler)? LineIsIntrsct : Zero; // To check !!!
		IntrsctPo = R02;
	}
}

//-------------------------------------------------------------------------

void radTPolygon::IntrsctOfTwoLines2(const TVector2d& R01, const TVector2d& R11, const TVector2d& R02, const TVector2d& R12, TVector2d& IntrsctPo, TLinesIntrsctCase& IntrsctCase)
{// This is used to determine self-intersection
	const double t_Toler = 5.E-12;

	TVector2d V1 = R11 - R01;
	TVector2d V2 = R12 - R02;

	double AbsV2x = fabs(V2.x), AbsV2y = fabs(V2.y);
	double V2Norm = (AbsV2x > AbsV2y)? AbsV2x : AbsV2y;
	double AbsV1x = fabs(V1.x), AbsV1y = fabs(V1.y);
	double V1Norm = (AbsV1x > AbsV1y)? AbsV1x : AbsV1y;
	double V_Toler = V1Norm*t_Toler;
	double D_Toler = V2Norm*V_Toler;

	double V1yV2x = V1.y*V2.x;
	double V1xV2y = V1.x*V2.y;
	double D = V1xV2y - V1yV2x;

	if(fabs(D) > D_Toler)
	{
		IntrsctPo.x = -(-R01.y*V1.x*V2.x + R02.y*V1.x*V2.x + R01.x*V1yV2x - R02.x*V1xV2y)/D;
		IntrsctPo.y = -(R02.y*V1yV2x - R01.y*V1xV2y + R01.x*V1.y*V2.y - R02.x*V1.y*V2.y)/D;

		double t_Intrsct2 = (Abs(V2.x) > Abs(V2.y))? (IntrsctPo.x - R02.x)/V2.x : (IntrsctPo.y - R02.y)/V2.y;
		IntrsctCase = ((t_Intrsct2 > t_Toler) && (t_Intrsct2 + t_Toler < 1.))? PointWithinBound : PointOutsideBound;
		if(IntrsctCase == PointWithinBound)
		{
			double t_Intrsct1 = (Abs(V1.x) > Abs(V1.y))? (IntrsctPo.x - R01.x)/V1.x : (IntrsctPo.y - R01.y)/V1.y;
			IntrsctCase = ((t_Intrsct1 > t_Toler) && (t_Intrsct1 + t_Toler < 1.))? PointWithinBound : PointOutsideBound;
		}

		if(IntrsctCase == PointWithinBound)
		{
			TVector2d IP_mi_R = IntrsctPo - R01;
			double AbsIP_mi_Rx = fabs(IP_mi_R.x), AbsIP_mi_Ry = fabs(IP_mi_R.y);
			double NormIP_mi_R = (AbsIP_mi_Rx > AbsIP_mi_Ry)? AbsIP_mi_Rx : AbsIP_mi_Ry;
			double AbsRx = fabs(R01.x), AbsRy = fabs(R01.x);
			double RNorm = ((AbsRx > AbsRy)? AbsRx : AbsRy);
			if(NormIP_mi_R < RNorm*t_Toler) { IntrsctCase = PointOnBoundEdge; return;}

			IP_mi_R = IntrsctPo - R11;
			AbsIP_mi_Rx = fabs(IP_mi_R.x); AbsIP_mi_Ry = fabs(IP_mi_R.y);
			NormIP_mi_R = (AbsIP_mi_Rx > AbsIP_mi_Ry)? AbsIP_mi_Rx : AbsIP_mi_Ry;
			AbsRx = fabs(R11.x); AbsRy = fabs(R11.x);
			RNorm = ((AbsRx > AbsRy)? AbsRx : AbsRy);
			if(NormIP_mi_R < RNorm*t_Toler) { IntrsctCase = PointOnBoundEdge; return;}

			IP_mi_R = IntrsctPo - R02;
			AbsIP_mi_Rx = fabs(IP_mi_R.x); AbsIP_mi_Ry = fabs(IP_mi_R.y);
			NormIP_mi_R = (AbsIP_mi_Rx > AbsIP_mi_Ry)? AbsIP_mi_Rx : AbsIP_mi_Ry;
			AbsRx = fabs(R02.x); AbsRy = fabs(R02.x);
			RNorm = ((AbsRx > AbsRy)? AbsRx : AbsRy);
			if(NormIP_mi_R < RNorm*t_Toler) { IntrsctCase = PointOnBoundEdge; return;}

			IP_mi_R = IntrsctPo - R12;
			AbsIP_mi_Rx = fabs(IP_mi_R.x), AbsIP_mi_Ry = fabs(IP_mi_R.y);
			NormIP_mi_R = (AbsIP_mi_Rx > AbsIP_mi_Ry)? AbsIP_mi_Rx : AbsIP_mi_Ry;
			AbsRx = fabs(R12.x); AbsRy = fabs(R12.x);
			RNorm = ((AbsRx > AbsRy)? AbsRx : AbsRy);
			if(NormIP_mi_R < RNorm*t_Toler) { IntrsctCase = PointOnBoundEdge; return;}
		}
	}
	else 
	{
		IntrsctPo = R02;

		double V3x = R01.x-R02.x, V3y = R01.y-R02.y;
		double AbsV3x = fabs(V3x), AbsV3y = fabs(V3y);
		double V3Norm = ((AbsV3x > AbsV3y)? AbsV3x : AbsV3y);

		double AbsR01x = fabs(R01.x), AbsR01y = fabs(R01.y);
		double R01Norm = ((AbsR01x > AbsR01y)? AbsR01x : AbsR01y);
		double AbsR02x = fabs(R02.x), AbsR02y = fabs(R02.y);
		double R02Norm = ((AbsR02x > AbsR02y)? AbsR02x : AbsR02y);
		double MaxNormR01R02 = (R01Norm > R02Norm)? R01Norm : R02Norm;

		if(V3Norm < MaxNormR01R02*t_Toler)
		{
			IntrsctCase = LineIsIntrsct; return;
		}

		double LineCoinsBufToler = V1Norm*t_Toler;
		double LineCoinsToler = LineCoinsBufToler*V3Norm;
		double CompareVal = fabs(V1.x*V3y - V1.y*V3x);
		IntrsctCase = (CompareVal < LineCoinsToler)? LineIsIntrsct : Zero; // To check !!!
	}
}

//-------------------------------------------------------------------------

void radTPolygon::FindStPointsForIntrsctLines(const TVector2d& V1, const TVector2d& V2, double k1, double q1, double k2, double q2, 
											  TVector2d* StPoForFirstSetOfIntrsctLines, TVector2d* StPoForSecondSetOfIntrsctLines)
{
	radTVect2dVect::iterator Iter = EdgePointsVector.begin();

	TVector2d* CurEdgePoPtr = &(*Iter);
	TVector2d* StEdgePoForV1Ptr;
	TVector2d* FiEdgePoForV1Ptr;
	TVector2d* StEdgePoForV2Ptr;
	TVector2d* FiEdgePoForV2Ptr;

	double MinProdV1 = 1.E+23;
	double MinProdV2 = 1.E+23;
	double MaxProdV1 = -1.E+23;
	double MaxProdV2 = -1.E+23;

	for(int i=1; i<=AmOfEdgePoints; i++)
	{
		double TestProdV1 = (*CurEdgePoPtr)*V1;
		double TestProdV2 = (*CurEdgePoPtr)*V2;

		if(TestProdV1 < MinProdV1)
		{
			MinProdV1 = TestProdV1; StEdgePoForV1Ptr = CurEdgePoPtr;
		}
		if(TestProdV1 > MaxProdV1)
		{
			MaxProdV1 = TestProdV1; FiEdgePoForV1Ptr = CurEdgePoPtr;
		}

		if(TestProdV2 < MinProdV2)
		{
			MinProdV2 = TestProdV2; StEdgePoForV2Ptr = CurEdgePoPtr;
		}
		if(TestProdV2 > MaxProdV2)
		{
			MaxProdV2 = TestProdV2; FiEdgePoForV2Ptr = CurEdgePoPtr;
		}

		if(i!= AmOfEdgePoints) CurEdgePoPtr = &(*(++Iter));
	}
	
	double GapForK2 = (V1.x!=0)? (FiEdgePoForV1Ptr->x - StEdgePoForV1Ptr->x)/V1.x : (FiEdgePoForV1Ptr->y - StEdgePoForV1Ptr->y)/V1.y;
	double GapForK1 = (V2.x!=0)? (FiEdgePoForV2Ptr->x - StEdgePoForV2Ptr->x)/V2.x : (FiEdgePoForV2Ptr->y - StEdgePoForV2Ptr->y)/V2.y;

	const double AbsZeroTol = 5.E-12;

	double q0_1 = (fabs(k1-1.)>AbsZeroTol)? pow(q1, 1./(k1-1.)) : q1;
	double Buf1 = q1*q0_1 - 1.;
	double a1_1 = (fabs(Buf1) > AbsZeroTol)? GapForK1*(q0_1 - 1.)/Buf1 : GapForK1/k1;

	double q0_2 = (fabs(k2-1.)>AbsZeroTol)? pow(q2, 1./(k2-1.)) : q2;
	double Buf2 = q2*q0_2 - 1.;
	double a1_2 = (fabs(Buf2) > AbsZeroTol)? GapForK2*(q0_2 - 1.)/Buf2 : GapForK2/k2;

	double InstSmallGap1 = a1_1, InstSmallGap2 = a1_2;

// The FirstSetOfIntrsctLines (k1+1 Lines incl. two edge lines) is parallel to V1
	StPoForFirstSetOfIntrsctLines[0] = *StEdgePoForV2Ptr;
	StPoForSecondSetOfIntrsctLines[0] = *StEdgePoForV1Ptr;

	for(int i1=1; i1<=int(k1); i1++)
	{
		if(i1 == int(k1)) StPoForFirstSetOfIntrsctLines[i1] = StPoForFirstSetOfIntrsctLines[0] + GapForK1*V2;
		else StPoForFirstSetOfIntrsctLines[i1] = StPoForFirstSetOfIntrsctLines[i1-1] + InstSmallGap1*V2;
		InstSmallGap1 *= q0_1;
	}
	for(int i2=1; i2<=int(k2); i2++)
	{
		if(i2 == int(k2)) StPoForSecondSetOfIntrsctLines[i2] = StPoForSecondSetOfIntrsctLines[0] + GapForK2*V1;
		else StPoForSecondSetOfIntrsctLines[i2] = StPoForSecondSetOfIntrsctLines[i2-1] + InstSmallGap2*V1;
		InstSmallGap2 *= q0_2;
	}
}

//-------------------------------------------------------------------------

void radTPolygon::FillInIntrsctInfoStruct(const TVector2d& V, const TVector2d* StPoForSetOfParLines, int k, 
										  radTPolyg2dIntrsctInfoVect* IntrsctInfoVectsForPolygBounds, 
										  radTPolyg2dIntrsctInfoVect* IntrsctInfoVectsForParLines, double& AbsPrecToler)
{
	const double RelPrecToler = 5.E-13; // To recognize coinsidence with edge points
	double MaxAbsPrecToler = 0.;

	TVector2d V_ParLin = V;
	TVector2d V_Bound, R0_ParLin, R0_Bound, R1_Bound, IntrsctPo;
	TLinesIntrsctCase IntrsctCase;

	int AmOfEdgePoints_m_1 = AmOfEdgePoints - 1;

	R0_Bound = EdgePointsVector[0];
	int i_Bound, i_ParLin;
	for(i_Bound = 0; i_Bound < AmOfEdgePoints; i_Bound++)
	{
		int NextPoNo = (i_Bound != AmOfEdgePoints_m_1)? (i_Bound + 1) : 0;
		R1_Bound = EdgePointsVector[NextPoNo];

		TVector2d BoundVect = R1_Bound - R0_Bound;
		double AbsX_BoundVect = fabs(BoundVect.x), AbsY_BoundVect = fabs(BoundVect.y);
		double NormBoundVect = (AbsX_BoundVect > AbsY_BoundVect)? AbsX_BoundVect : AbsY_BoundVect;
		AbsPrecToler = NormBoundVect*RelPrecToler;
		if(MaxAbsPrecToler < AbsPrecToler) MaxAbsPrecToler = AbsPrecToler;

		short StartEdgeAlreadyMarked = 0, EndEdgeAlreadyMarked = 0;
		for(i_ParLin = 0; i_ParLin <= k; i_ParLin++)
		{
			R0_ParLin = StPoForSetOfParLines[i_ParLin];

			IntrsctOfTwoLines(V_ParLin, R0_ParLin, R0_Bound, R1_Bound, IntrsctPo, IntrsctCase);

			radTPolyg2dIntrsctInfo IntrsctInfoForBound, IntrsctInfoForParLine;
			short CurIntrsctIsAtBoundEnd = 0;
			if((IntrsctCase == PointWithinBound) || (IntrsctCase == PointOnBoundEdge))
			{
				IntrsctInfoForBound.IntrsctPoint = IntrsctPo;
				IntrsctInfoForBound.NoOfIntrsctItem = i_ParLin;
				IntrsctInfoForBound.TypeOfIntrsctItem = ParLine;
				//NoOfIntrsctPoOnIntrsctItem to fill later

				IntrsctInfoForParLine.IntrsctPoint = IntrsctPo;
				IntrsctInfoForParLine.NoOfIntrsctItem = i_Bound;
				IntrsctInfoForParLine.TypeOfIntrsctItem = Bound;
				//NoOfIntrsctPoOnIntrsctItem to fill later

				if((Abs(IntrsctPo.x - R0_Bound.x) < AbsPrecToler) && (Abs(IntrsctPo.y - R0_Bound.y) < AbsPrecToler))
				{
					StartEdgeAlreadyMarked = 1;
					IntrsctInfoForBound.IntrsctPoint.x = R0_Bound.x;
					IntrsctInfoForBound.IntrsctPoint.y = R0_Bound.y;
				}
				if((Abs(IntrsctPo.x - R1_Bound.x) < AbsPrecToler) && (Abs(IntrsctPo.y - R1_Bound.y) < AbsPrecToler))
				{
					EndEdgeAlreadyMarked = 1; CurIntrsctIsAtBoundEnd = 1;
					IntrsctInfoForBound.IntrsctPoint.x = R1_Bound.x;
					IntrsctInfoForBound.IntrsctPoint.y = R1_Bound.y;
				}

				IntrsctInfoForParLine.IntrsctMultiplicity = (IntrsctCase == PointWithinBound)? 1 : 2;
				IntrsctInfoForBound.IntrsctMultiplicity = IntrsctInfoForParLine.IntrsctMultiplicity;
				IntrsctInfoForParLine.IntrsctIsLine = IntrsctInfoForBound.IntrsctIsLine = 0;
			}
			else if(IntrsctCase == LineIsIntrsct)
			{
				IntrsctInfoForBound.IntrsctPoint = R0_Bound;
				IntrsctInfoForBound.NoOfIntrsctItem = i_ParLin;
				IntrsctInfoForBound.TypeOfIntrsctItem = ParLine;
				//NoOfIntrsctPoOnIntrsctItem to fill later

				IntrsctInfoForParLine.IntrsctPoint = R0_Bound;
				IntrsctInfoForParLine.NoOfIntrsctItem = i_Bound;
				IntrsctInfoForParLine.TypeOfIntrsctItem = Bound;
				//NoOfIntrsctPoOnIntrsctItem to fill later

				StartEdgeAlreadyMarked = 1;

				IntrsctInfoForParLine.IntrsctMultiplicity = IntrsctInfoForBound.IntrsctMultiplicity = 2;
				IntrsctInfoForParLine.IntrsctIsLine = IntrsctInfoForBound.IntrsctIsLine = 1;
			}
			if((IntrsctCase == PointWithinBound) || (IntrsctCase == PointOnBoundEdge) || (IntrsctCase == LineIsIntrsct))
			{
				IntrsctInfoVectsForPolygBounds[i_Bound].push_back(IntrsctInfoForBound);
				if(!CurIntrsctIsAtBoundEnd) IntrsctInfoVectsForParLines[i_ParLin].push_back(IntrsctInfoForParLine);
			}
		}
		if(!StartEdgeAlreadyMarked)
		{
			radTPolyg2dIntrsctInfo IntrsctInfoForBound;
			IntrsctInfoForBound.IntrsctPoint = R0_Bound;
			IntrsctInfoForBound.NoOfIntrsctItem = (i_Bound != 0)? (i_Bound - 1) : AmOfEdgePoints_m_1;
			IntrsctInfoForBound.TypeOfIntrsctItem = Bound;
			//NoOfIntrsctPoOnIntrsctItem to fill later: the last one
			IntrsctInfoVectsForPolygBounds[i_Bound].push_back(IntrsctInfoForBound);
		}
		if(!EndEdgeAlreadyMarked)
		{
			radTPolyg2dIntrsctInfo IntrsctInfoForBound;
			IntrsctInfoForBound.IntrsctPoint = R1_Bound;
			IntrsctInfoForBound.NoOfIntrsctItem = (i_Bound != AmOfEdgePoints_m_1)? (i_Bound + 1) : 0;
			IntrsctInfoForBound.TypeOfIntrsctItem = Bound;
			IntrsctInfoForBound.NoOfIntrsctPoOnIntrsctItem = 0;
			IntrsctInfoVectsForPolygBounds[i_Bound].push_back(IntrsctInfoForBound);
		}
		R0_Bound = R1_Bound;
	}

// Sorting vectors of IntrsctInfo:
	for(i_Bound = 0; i_Bound < AmOfEdgePoints; i_Bound++)
	{
		radTPolyg2dIntrsctInfoVect&  IntrsctInfoVect = IntrsctInfoVectsForPolygBounds[i_Bound];
		TVector2d BoundDir = (i_Bound != AmOfEdgePoints_m_1)? (EdgePointsVector[i_Bound+1] - EdgePointsVector[i_Bound]) : (EdgePointsVector[0] - EdgePointsVector[i_Bound]);
		sort(IntrsctInfoVect.begin(), IntrsctInfoVect.end(), TCompareIntrsctInfo(BoundDir));
	}
	TVector2d dR0 = StPoForSetOfParLines[1] - StPoForSetOfParLines[0];
	for(i_ParLin = 0; i_ParLin <= k; i_ParLin++)
	{
		radTPolyg2dIntrsctInfoVect&  IntrsctInfoVect = IntrsctInfoVectsForParLines[i_ParLin];
		if(!IntrsctInfoVect.empty())
		{
			TVector2d mV_ParLin(-V_ParLin.x, -V_ParLin.y);
			TVector2d V_ForSort = (V_ParLin.x*dR0.y > dR0.x*V_ParLin.y)? V_ParLin : mV_ParLin;
			sort(IntrsctInfoVect.begin(), IntrsctInfoVect.end(), TCompareIntrsctInfo(V_ForSort));
		}
	}

// Filling-in NoOfIntrsctPoOnIntrsctItem Fields:
	radTPolyg2dIntrsctInfo* BufInfoPtr;
	radTPolyg2dIntrsctInfoVect* InfoVectPtr;
	for(i_Bound = 0; i_Bound < AmOfEdgePoints; i_Bound++)
	{
		radTPolyg2dIntrsctInfoVect&  IntrsctInfoVect = IntrsctInfoVectsForPolygBounds[i_Bound];
		for(radTPolyg2dIntrsctInfoVect::iterator iter = IntrsctInfoVect.begin(); iter != IntrsctInfoVect.end(); ++iter)
		{
			BufInfoPtr = &(*iter);
			TVector2d ThisPoint = BufInfoPtr->IntrsctPoint;
			InfoVectPtr = (BufInfoPtr->TypeOfIntrsctItem == Bound)? &(IntrsctInfoVectsForPolygBounds[BufInfoPtr->NoOfIntrsctItem]) : &(IntrsctInfoVectsForParLines[BufInfoPtr->NoOfIntrsctItem]);
			int ProbNo;
			for(ProbNo = 0; ProbNo < (int)(InfoVectPtr->size()); ProbNo++)
			{
				TVector2d& TestPoint = (*InfoVectPtr)[ProbNo].IntrsctPoint;
				if((Abs(ThisPoint.x - TestPoint.x) < AbsPrecToler) && (Abs(ThisPoint.y - TestPoint.y) < AbsPrecToler)) break;
			}
			BufInfoPtr->NoOfIntrsctPoOnIntrsctItem = ProbNo;
		}
	}
	for(i_ParLin = 0; i_ParLin <= k; i_ParLin++)
	{
		radTPolyg2dIntrsctInfoVect&  IntrsctInfoVect = IntrsctInfoVectsForParLines[i_ParLin];
		for(radTPolyg2dIntrsctInfoVect::iterator iter = IntrsctInfoVect.begin(); iter != IntrsctInfoVect.end(); ++iter)
		{
			BufInfoPtr = &(*iter);
			TVector2d ThisPoint = BufInfoPtr->IntrsctPoint;
			InfoVectPtr = &(IntrsctInfoVectsForPolygBounds[BufInfoPtr->NoOfIntrsctItem]);
			int ProbNo;
			for(ProbNo = 0; ProbNo < (int)(InfoVectPtr->size()); ProbNo++)
			{
				TVector2d& TestPoint = (*InfoVectPtr)[ProbNo].IntrsctPoint;
				if((Abs(ThisPoint.x - TestPoint.x) < AbsPrecToler) && (Abs(ThisPoint.y - TestPoint.y) < AbsPrecToler)) break;
			}
			BufInfoPtr->NoOfIntrsctPoOnIntrsctItem = ProbNo;
		}
	}
	AbsPrecToler = MaxAbsPrecToler;
}

//-------------------------------------------------------------------------

int radTPolygon::SubdivideBySetOfParallelLines(const TVector2d& V, TVector2d* StPoForSetOfLines, int k, radTvhg& VectOfHandlesToNewPolyg)
{
	const int MaxNewPointsInLoop = 200; // To abort too fine or erroneous subdivision
	radTSend Send;

	radTPolyg2dIntrsctInfoVect IntrsctInfoVectsForPolygBounds[100];
	radTPolyg2dIntrsctInfoVect IntrsctInfoVectsForParLines[100];

	double AbsPrecToler;
	FillInIntrsctInfoStruct(V, StPoForSetOfLines, k, IntrsctInfoVectsForPolygBounds, IntrsctInfoVectsForParLines, AbsPrecToler);

	radTVect2dVect BaseForBoundStartLoop, BaseForLineStartLoop;

#ifdef __GCC__
	vector<int> NosOfPointsPassedVect, NoOfPointOnCurLineVect;
#else
	vector<int, allocator<int> > NosOfPointsPassedVect, NoOfPointOnCurLineVect;
#endif

	radTPolyg2dIntrsctInfoVect* CurLinIntrsctInfoVectPtr = &(IntrsctInfoVectsForParLines[1]);

	short FirstTrial = 1;
	short SomeNewPolygonsRecognized = 0;

	for(int NoOfParLin = 1; NoOfParLin <= k; NoOfParLin++)
	{
		int AmOfIntrsctPointsOnTheLine = (int)(CurLinIntrsctInfoVectPtr->size());

		double dAmOfPolygToTest = 1.E-10;
		double Prev_dAmOfPolygToTest = dAmOfPolygToTest;

		radTPolyg2dIntrsctInfo* BufFirstInfoPtr = 0;
		radTPolyg2dIntrsctInfo* BufLastInfoPtr = 0;
		short IntrsctWasLine = 0;

		int AmOfIntrsctPointsOnTheLine_1 = AmOfIntrsctPointsOnTheLine - 1;
		for(int ic=0; ic<AmOfIntrsctPointsOnTheLine; ic++)
		{
			if(ic==0) NoOfPointOnCurLineVect.push_back(ic);
			radTPolyg2dIntrsctInfo* BufIntrsctInfoPtr = &((*CurLinIntrsctInfoVectPtr)[ic]);
			dAmOfPolygToTest += 0.5*(BufIntrsctInfoPtr->IntrsctMultiplicity - BufIntrsctInfoPtr->IntrsctIsLine);

			if(ic==0) BufFirstInfoPtr = BufIntrsctInfoPtr;
			if(ic==AmOfIntrsctPointsOnTheLine_1) BufLastInfoPtr = BufIntrsctInfoPtr;
			if(BufIntrsctInfoPtr->IntrsctIsLine) IntrsctWasLine = 1;

			if(((dAmOfPolygToTest - Prev_dAmOfPolygToTest) > 1.1) && (ic > 0))
			{
				NoOfPointOnCurLineVect.push_back(ic);
				Prev_dAmOfPolygToTest = dAmOfPolygToTest;
			}
		}
		if((BufFirstInfoPtr != 0) && (BufLastInfoPtr != 0))
			if((BufFirstInfoPtr != BufLastInfoPtr) && (!IntrsctWasLine) && (BufFirstInfoPtr->IntrsctMultiplicity==2) && (BufLastInfoPtr->IntrsctMultiplicity==2)) dAmOfPolygToTest -= 0.5;

		int AmOfPolygToTest = int(dAmOfPolygToTest);

		if((AmOfIntrsctPointsOnTheLine==1) && FirstTrial) goto BottomOfTheMainLoop;

		int i_PolOnLin;
		for(i_PolOnLin = 0; i_PolOnLin < AmOfPolygToTest; i_PolOnLin++)
		{
			int NoOfPointOnCurLine = NoOfPointOnCurLineVect[i_PolOnLin];

			radTPolyg2dIntrsctInfo* StartIntrsctInfoPtr = &((*CurLinIntrsctInfoVectPtr)[NoOfPointOnCurLine]);
			TVector2d CurFirstPoint = StartIntrsctInfoPtr->IntrsctPoint;

			radTPolyg2dIntrsctInfo* SecondIntrsctInfoPtr = 0;
			TVector2d SecondPointOnThisLine(1.e+23, 1.e+23);

			if((int)(CurLinIntrsctInfoVectPtr->size()) > (NoOfPointOnCurLine + 1)) //OC031107
			{
				SecondIntrsctInfoPtr = &((*CurLinIntrsctInfoVectPtr)[NoOfPointOnCurLine + 1]);
				SecondPointOnThisLine = SecondIntrsctInfoPtr->IntrsctPoint;
			}
			else //OC031107
			{
				SecondIntrsctInfoPtr = StartIntrsctInfoPtr;
				SecondPointOnThisLine = CurFirstPoint;
			}

// Loop starting from Bound:
			radTPolyg2dIntrsctInfo* IntrsctInfoPtr = StartIntrsctInfoPtr;
			radTPolyg2dIntrsctInfo* PrevIntrsctInfoPtr;

			short ThePointOnTheLineWasNotPassed = 1;
			for(radTVectInt::iterator Iter = NosOfPointsPassedVect.begin(); Iter != NosOfPointsPassedVect.end(); ++Iter)
			{
				if(*Iter == NoOfPointOnCurLine) ThePointOnTheLineWasNotPassed = 0; break;
			}
			if(ThePointOnTheLineWasNotPassed)
			{
				BaseForBoundStartLoop.push_back(CurFirstPoint);
				PrevIntrsctInfoPtr = IntrsctInfoPtr;

				IntrsctInfoPtr = &((IntrsctInfoVectsForPolygBounds[IntrsctInfoPtr->NoOfIntrsctItem])[IntrsctInfoPtr->NoOfIntrsctPoOnIntrsctItem + 1]);
				TVector2d CurPoint = IntrsctInfoPtr->IntrsctPoint;
				TVector2d PrevPoint;

				short SecondPointOnThisLineWasUsed = ((CurPoint.x==SecondPointOnThisLine.x) && (CurPoint.y==SecondPointOnThisLine.y))? 1 : 0;
				BaseForBoundStartLoop.push_back(CurPoint);

				int LocLoopPassCount = 0;
				for(;;)
				{
					radTPolyg2dIntrsctInfoVect* CurIntrsctInfoVectPtr = (IntrsctInfoPtr->TypeOfIntrsctItem == Bound)? &(IntrsctInfoVectsForPolygBounds[IntrsctInfoPtr->NoOfIntrsctItem]) : &(IntrsctInfoVectsForParLines[IntrsctInfoPtr->NoOfIntrsctItem]);

					if(CurIntrsctInfoVectPtr != CurLinIntrsctInfoVectPtr)
					{
						int BufIndx = IntrsctInfoPtr->NoOfIntrsctPoOnIntrsctItem;
						if(IntrsctInfoPtr->NoOfIntrsctPoOnIntrsctItem + 1 >= (int)(CurIntrsctInfoVectPtr->size()))
						{
							PrevIntrsctInfoPtr = IntrsctInfoPtr;

							IntrsctInfoPtr = &((*CurIntrsctInfoVectPtr)[BufIndx]);
							CurIntrsctInfoVectPtr = (IntrsctInfoPtr->TypeOfIntrsctItem == Bound)? &(IntrsctInfoVectsForPolygBounds[IntrsctInfoPtr->NoOfIntrsctItem]) : &(IntrsctInfoVectsForParLines[IntrsctInfoPtr->NoOfIntrsctItem]);
						}

						PrevIntrsctInfoPtr = IntrsctInfoPtr;
						IntrsctInfoPtr = &((*CurIntrsctInfoVectPtr)[BufIndx + 1]);
						
						if(IntrsctInfoPtr->SkipThisPoint)
						{
							PrevIntrsctInfoPtr = IntrsctInfoPtr;
							IntrsctInfoPtr = &((*CurIntrsctInfoVectPtr)[BufIndx + 2]);
						}
					}
					else
					{
						int NoOfPointPassed = IntrsctInfoPtr->NoOfIntrsctPoOnIntrsctItem - 1;
						radTPolyg2dIntrsctInfo* ProbIntrsctInfoPtr = &((*CurIntrsctInfoVectPtr)[NoOfPointPassed]);
						TVector2d ProbCurPoint = ProbIntrsctInfoPtr->IntrsctPoint;

						if((IntrsctInfoPtr->IntrsctMultiplicity == 2) && (!IntrsctInfoPtr->IntrsctIsLine))
						{
							radTPolyg2dIntrsctInfo* LocProbIntrsctInfoPtr = &((*((PrevIntrsctInfoPtr->NoOfIntrsctItem != (EdgePointsVector.size()-1))? &(IntrsctInfoVectsForPolygBounds[PrevIntrsctInfoPtr->NoOfIntrsctItem + 1]) : &(IntrsctInfoVectsForPolygBounds[0])))[1]);
							TVector2d ProbBoundCurPoint = LocProbIntrsctInfoPtr->IntrsctPoint;
							TVector2d ProbBoundCurPo_m_CurPo = ProbBoundCurPoint - CurPoint;
							double InvLen = 1./sqrt(ProbBoundCurPo_m_CurPo.x*ProbBoundCurPo_m_CurPo.x + ProbBoundCurPo_m_CurPo.y*ProbBoundCurPo_m_CurPo.y);
							ProbBoundCurPo_m_CurPo = InvLen*ProbBoundCurPo_m_CurPo;

							TVector2d ProbCurPo_m_CurPo = ProbCurPoint - CurPoint;
							InvLen = 1./sqrt(ProbCurPo_m_CurPo.x*ProbCurPo_m_CurPo.x + ProbCurPo_m_CurPo.y*ProbCurPo_m_CurPo.y);
							ProbCurPo_m_CurPo = InvLen*ProbCurPo_m_CurPo;

							TVector2d CurPo_m_PrevPo = CurPoint - PrevPoint;

							if((CurPo_m_PrevPo.x*ProbBoundCurPo_m_CurPo.x + CurPo_m_PrevPo.y*ProbBoundCurPo_m_CurPo.y) < (CurPo_m_PrevPo.x*ProbCurPo_m_CurPo.x + CurPo_m_PrevPo.y*ProbCurPo_m_CurPo.y))
							{
								NosOfPointsPassedVect.push_back(IntrsctInfoPtr->NoOfIntrsctPoOnIntrsctItem);

								PrevIntrsctInfoPtr = IntrsctInfoPtr;
								IntrsctInfoPtr = LocProbIntrsctInfoPtr;
								goto AssignNewCurPoint;
							}
						}

						int LocCheckCount=0;
CheckTheProbPoint:
						radTPolyg2dIntrsctInfoVect* ProbIntrsctInfoVectPtr = (ProbIntrsctInfoPtr->TypeOfIntrsctItem == Bound)? &(IntrsctInfoVectsForPolygBounds[ProbIntrsctInfoPtr->NoOfIntrsctItem]) : &(IntrsctInfoVectsForParLines[ProbIntrsctInfoPtr->NoOfIntrsctItem]);
						TVector2d NextProbPoint = ((*ProbIntrsctInfoVectPtr)[ProbIntrsctInfoPtr->NoOfIntrsctPoOnIntrsctItem + 1]).IntrsctPoint;
						TVector2d NextProbPo_m_ProbCurPo = NextProbPoint - ProbCurPoint;

						TVector2d R0cur_m_R0prev = StPoForSetOfLines[NoOfParLin] - StPoForSetOfLines[NoOfParLin-1];
						if(Sign(V.x*NextProbPo_m_ProbCurPo.y - V.y*NextProbPo_m_ProbCurPo.x) == Sign(V.x*R0cur_m_R0prev.y - V.y*R0cur_m_R0prev.x))
						//if((Sign(V.x*NextProbPo_m_ProbCurPo.y - V.y*NextProbPo_m_ProbCurPo.x) == Sign(V.x*R0cur_m_R0prev.y - V.y*R0cur_m_R0prev.x))
						//	&& (IntrsctInfoPtr->NoOfIntrsctPoOnIntrsctItem >= 2)) //OC03032017 (to make sure that [IntrsctInfoPtr->NoOfIntrsctPoOnIntrsctItem - 2] can be formally taken)
						{
							NosOfPointsPassedVect.push_back(NoOfPointPassed);

							NoOfPointPassed = IntrsctInfoPtr->NoOfIntrsctPoOnIntrsctItem - 2;
							ProbIntrsctInfoPtr = &((*CurIntrsctInfoVectPtr)[NoOfPointPassed]);
							ProbCurPoint = ProbIntrsctInfoPtr->IntrsctPoint;

							if((Abs(ProbCurPoint.x - CurFirstPoint.x) < AbsPrecToler) && (Abs(ProbCurPoint.y - CurFirstPoint.y) < AbsPrecToler)) break;
							else
							{
								if((++LocCheckCount) > MaxNewPointsInLoop)
								{
									Send.ErrorMessage("Radia::Error103"); return 0;
								}
								goto CheckTheProbPoint;
							}
						}
						else
						{
							PrevIntrsctInfoPtr = IntrsctInfoPtr;

							IntrsctInfoPtr = ProbIntrsctInfoPtr;
							NosOfPointsPassedVect.push_back(NoOfPointPassed);
						}
					}
AssignNewCurPoint:
					PrevPoint = CurPoint;
					CurPoint = IntrsctInfoPtr->IntrsctPoint;
					if((Abs(CurPoint.x - CurFirstPoint.x) < AbsPrecToler) && (Abs(CurPoint.y - CurFirstPoint.y) < AbsPrecToler)) break;

					if(SecondPointOnThisLineWasUsed)
						if((Abs(CurPoint.x - SecondPointOnThisLine.x) < AbsPrecToler) && (Abs(CurPoint.y - SecondPointOnThisLine.y) < AbsPrecToler))
						{
							SecondIntrsctInfoPtr->SkipThisPoint = 1;
							int LocSize = (int)(BaseForBoundStartLoop.size());
							for(int ii=1; ii<LocSize; ii++)	BaseForBoundStartLoop[ii-1] = BaseForBoundStartLoop[ii];
							BaseForBoundStartLoop.pop_back();
							break;
						}

					if((Abs(CurPoint.x - PrevPoint.x) > AbsPrecToler) || (Abs(CurPoint.y - PrevPoint.y) > AbsPrecToler)) 
						BaseForBoundStartLoop.push_back(CurPoint);

					if((++LocLoopPassCount) > MaxNewPointsInLoop)
					{
						Send.ErrorMessage("Radia::Error103"); return 0;
					}
				}

				RandomizeNonConvexEdgePoints(BaseForBoundStartLoop);

				radThg hg(new radTPolygon(BaseForBoundStartLoop));
				VectOfHandlesToNewPolyg.push_back(hg);
				SomeNewPolygonsRecognized = 1;
				FirstTrial = 0;
			}
			BaseForBoundStartLoop.erase(BaseForBoundStartLoop.begin(), BaseForBoundStartLoop.end());

// Loop starting from Line:
			if((AmOfIntrsctPointsOnTheLine > 1) && (NoOfParLin != k))
			{
				short PolygHasSegmentOfOtherParLine = 0;
				BaseForLineStartLoop.push_back(CurFirstPoint);

				PrevIntrsctInfoPtr = IntrsctInfoPtr;

				IntrsctInfoPtr = &((*CurLinIntrsctInfoVectPtr)[NoOfPointOnCurLine + 1]);
				TVector2d CurPoint = IntrsctInfoPtr->IntrsctPoint;

				if((Abs(CurPoint.x - CurFirstPoint.x) > AbsPrecToler) || (Abs(CurPoint.y - CurFirstPoint.y) > AbsPrecToler))
				{
					BaseForLineStartLoop.push_back(CurPoint);
					int LocLoopPassCount=0;
					for(;;)
					{
						radTPolyg2dIntrsctInfoVect* CurIntrsctInfoVectPtr = (IntrsctInfoPtr->TypeOfIntrsctItem == Bound)? &(IntrsctInfoVectsForPolygBounds[IntrsctInfoPtr->NoOfIntrsctItem]) : &(IntrsctInfoVectsForParLines[IntrsctInfoPtr->NoOfIntrsctItem]);

						PrevIntrsctInfoPtr = IntrsctInfoPtr;

						IntrsctInfoPtr = &((*CurIntrsctInfoVectPtr)[IntrsctInfoPtr->NoOfIntrsctPoOnIntrsctItem + 1]);
						CurPoint = IntrsctInfoPtr->IntrsctPoint;
						if(IntrsctInfoPtr->TypeOfIntrsctItem == ParLine)
						{
							if(&(IntrsctInfoVectsForParLines[IntrsctInfoPtr->NoOfIntrsctItem]) != CurLinIntrsctInfoVectPtr) { PolygHasSegmentOfOtherParLine = 1; break;}
							else if((Abs(CurPoint.x - CurFirstPoint.x) < AbsPrecToler) && (Abs(CurPoint.y - CurFirstPoint.y) < AbsPrecToler)) break;
						}
						BaseForLineStartLoop.push_back(CurPoint);

						if((++LocLoopPassCount) > MaxNewPointsInLoop)
						{
							Send.ErrorMessage("Radia::Error103"); return 0;
						}
					}
					if(!PolygHasSegmentOfOtherParLine)
					{
						RandomizeNonConvexEdgePoints(BaseForLineStartLoop);

						radThg hg(new radTPolygon(BaseForLineStartLoop));
						VectOfHandlesToNewPolyg.push_back(hg);
						SomeNewPolygonsRecognized = 1;
						FirstTrial = 0;
					}
				}
			}
			BaseForLineStartLoop.erase(BaseForLineStartLoop.begin(), BaseForLineStartLoop.end());

			NoOfPointOnCurLineVect.erase(NoOfPointOnCurLineVect.begin(), NoOfPointOnCurLineVect.end());
		}
BottomOfTheMainLoop:

		NosOfPointsPassedVect.erase(NosOfPointsPassedVect.begin(), NosOfPointsPassedVect.end());
		CurLinIntrsctInfoVectPtr++;
	}

	if(!SomeNewPolygonsRecognized)
	{
		radThg hg(new radTPolygon(EdgePointsVector));
		VectOfHandlesToNewPolyg.push_back(hg);
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTPolygon::SubdivideItself(double* SubdivArray, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions) 
{
	char SubdivideCoils = pSubdivOptions->SubdivideCoils;
	char PutNewStuffIntoGenCont = pSubdivOptions->PutNewStuffIntoGenCont;

	double k1 = *SubdivArray; // kx
	double q1 = *(SubdivArray+1);
	double k2 = *(SubdivArray+2); // ky
	double q2 = *(SubdivArray+3);

	if(pSubdivOptions->SubdivisionParamCode == 1)
	{
		double SizeX, SizeY;
		EstimateSize(SizeX, SizeY);

		k1 = (k1 < SizeX)? Round(SizeX/k1) : 1.;
		k2 = (k2 < SizeY)? Round(SizeY/k2) : 1.;
	}

	const double ZeroTol = 1.E-10;
	if((fabs(k1-1.)<ZeroTol) && (fabs(k2-1.)<ZeroTol)) return 1;

	TVector2d V1(0.,1.), V2(1.,0.); // V1 and V2 should be input in General Case

	radTGroup* GroupInPlaceOfThisPtr = new radTGroup;
	GroupInPlaceOfThisPtr->IsGroupMember = IsGroupMember;
	GroupInPlaceOfThisPtr->g3dListOfTransform = g3dListOfTransform;

	radThg NewHandle(GroupInPlaceOfThisPtr);
	TVector2d StPoForFirstSetOfIntrsctLines[200];
	TVector2d StPoForSecondSetOfIntrsctLines[200];

	radTvhg HandlersToPolygAfterFirstSubdiv, HandlersToFinalPolyg;

	FindStPointsForIntrsctLines(V1, V2, k1, q1, k2, q2, StPoForFirstSetOfIntrsctLines, StPoForSecondSetOfIntrsctLines);

	int LocSubdOK = SubdivideBySetOfParallelLines(V1, StPoForFirstSetOfIntrsctLines, int(k1), HandlersToPolygAfterFirstSubdiv);
	if(!LocSubdOK) return 0;

	int AmOfPolygAfterFirstSubdiv = (int)(HandlersToPolygAfterFirstSubdiv.size());

	int NewStuffCounter = 0;
	if(k2 > 1.+ZeroTol)
	{
		for(int i1=0; i1<AmOfPolygAfterFirstSubdiv; i1++)
		{
			int LocSubdOK = ((radTPolygon*)((HandlersToPolygAfterFirstSubdiv[i1]).rep))->SubdivideBySetOfParallelLines(V2, StPoForSecondSetOfIntrsctLines, int(k2), HandlersToFinalPolyg);
			if(!LocSubdOK) return 0;

			int AmOfFinalPolygInThePortion = (int)(HandlersToFinalPolyg.size());
			for(int i2=0; i2<AmOfFinalPolygInThePortion; i2++)
			{
				radThg hg = HandlersToFinalPolyg[i2];
				if(PutNewStuffIntoGenCont) GroupInPlaceOfThisPtr->AddElement(radPtr->AddElementToContainer(hg), hg);
				else GroupInPlaceOfThisPtr->AddElement(++NewStuffCounter, hg);
			}
			HandlersToFinalPolyg.erase(HandlersToFinalPolyg.begin(), HandlersToFinalPolyg.end());
		}
	}
	else
	{
		for(int i1=0; i1<AmOfPolygAfterFirstSubdiv; i1++)
		{
			radThg hg = HandlersToPolygAfterFirstSubdiv[i1];
			if(PutNewStuffIntoGenCont) GroupInPlaceOfThisPtr->AddElement(radPtr->AddElementToContainer(hg), hg);
			else GroupInPlaceOfThisPtr->AddElement(++NewStuffCounter, hg);
		}
	}
	HandlersToPolygAfterFirstSubdiv.erase(HandlersToPolygAfterFirstSubdiv.begin(), HandlersToPolygAfterFirstSubdiv.end());

	In_hg = NewHandle;
	return 1;
}

//-------------------------------------------------------------------------

void radTPolygon::ComputeCentrPoint(short& Out_InsideBlock)
{
	const double Max_k = 1.E+10;

	double AbsRandX = radCR.AbsRandMagnitude(CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(CentrPoint.y);
	if(AbsRandX == 0.) AbsRandX = 1.E-12;
	if(AbsRandY == 0.) AbsRandY = 1.E-12;

	int AmOfEdPoInBase = AmOfEdgePoints;
	int AmOfEdPoInBase_mi_1 = AmOfEdPoInBase - 1;

	radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
	TVector2d First2d = *BaseIter;
	
	double x1 = First2d.x;
	double y1 = First2d.y;
	double x2, y2;
	double x1e2 = x1*x1, x2e2;

	double SS=0., SXc=0., SYc=0.;

	int i;
	for(i=0; i<AmOfEdPoInBase; i++)
	{
		++BaseIter;
		if(i!=AmOfEdPoInBase_mi_1)
		{
			x2 = (*BaseIter).x; y2 = (*BaseIter).y;
		}
		else
		{
			x2 = First2d.x; y2 = First2d.y;
		}
		x2e2 = x2*x2;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);

		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			double k = (y2-y1)/(x2-x1), b = y1 - k*x1;
			double x1px2 = x1 + x2;

			SS += x2mx1*(2.*b + k*x1px2);
			SXc += -3.*b*x1e2 - 2.*k*x1e2*x1 + x2e2*(3.*b + 2.*k*x2);
			SYc += x2mx1*(3.*(b*b + b*k*x1px2) + k*k*(x1e2 + x1*x2 + x2e2));
		}
		x1 = x2; y1 = y2;
		x1e2 = x2e2;
	}
	double Square = 0.5*SS;
	double One_d_6Square = 1./Square/6.;
	CentrPoint.x = SXc*One_d_6Square;
	CentrPoint.y = SYc*One_d_6Square;

	int X_In_Count=0, X_Out_Count=0;
	int Y_In_Count=0, Y_Out_Count=0;
	short OutsideGateX=1, OutsideGateY=1;

	First2d = First2d - CentrPoint;
	x1 = First2d.x;
	y1 = First2d.y;

	// Artificial shift (to correctly determine if center point is inside the polygon):
	if(x1==0.) x1 = AbsRandX;
	if(y1==0.) y1 = AbsRandY;

	BaseIter = EdgePointsVector.begin();

	for(i=0; i<AmOfEdPoInBase; i++)
	{
		++BaseIter;
		if(i!=AmOfEdPoInBase_mi_1)
		{
			x2 = (*BaseIter).x - CentrPoint.x; 
			y2 = (*BaseIter).y - CentrPoint.y;
		}
		else
		{
			x2 = First2d.x; y2 = First2d.y;
		}

		// Artificial shift (to correctly determine if center point is inside the polygon):
		if(x2==0.) x2 = AbsRandX;
		if(y2==0.) y2 = AbsRandY;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);
		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			double k = (y2-y1)/(x2-x1), b = y1 - k*x1;
			if(x1*x2 <= 0.)
			{
				short LocInside = ((x2>x1)? (b<=0) : (b>=0));
				X_In_Count += LocInside? 1 : 0;
				X_Out_Count += (!LocInside)? 1 : 0;

				OutsideGateX = 0;
			}
			if(y1*y2 <= 0.)	OutsideGateY = 0;
		}
		else
		{
			if(y1*y2 <= 0.)
			{
				short LocInside = ((y2>y1)? (x1>=0) : (x1<=0));
				Y_In_Count += LocInside? 1 : 0;
				Y_Out_Count += (!LocInside)? 1 : 0;

				OutsideGateY = 0;
			}
		}
		x1 = x2; y1 = y2;
	}
	int X_In = X_In_Count - X_Out_Count;
	int Y_In = X_In_Count - X_Out_Count;
	short InsideBlock = (CheckIfPosEven(X_In) && CheckIfPosEven(Y_In) && (!OutsideGateX) && (!OutsideGateY));

	Out_InsideBlock = InsideBlock;
}

//-------------------------------------------------------------------------

//void radTPolygon::SimpleComputeCentrPoint(short& Out_InsideBlock)
void radTPolygon::SimpleComputeCentrPoint(short& Out_InsideBlock, char& Out_PointsAreDifferent) //OC 220902
{
	Out_PointsAreDifferent = 1; //OC 220902
	Out_InsideBlock = 1; //OC 080108

	const double Max_k = 1.E+10;

	double AbsRandX = radCR.AbsRandMagnitude(CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(CentrPoint.y);
	if(AbsRandX == 0.) AbsRandX = 1.E-12;
	if(AbsRandY == 0.) AbsRandY = 1.E-12;

	int AmOfEdPoInBase = AmOfEdgePoints;
	int AmOfEdPoInBase_mi_1 = AmOfEdPoInBase - 1;

	double SumX = 0., SumY = 0.;
	double xFirstP=0, yFirstP=0, xPrevP=0, yPrevP=0; //OC 220902

	for(int kk=0; kk<AmOfEdgePoints; kk++)
	{
		TVector2d& CurrentP = EdgePointsVector[kk];

		//SumX += CurrentP.x;
		//SumY += CurrentP.y;
		double& xCurrentP = CurrentP.x; //OC 220902
		double& yCurrentP = CurrentP.y; //OC 220902
		SumX += xCurrentP;
		SumY += yCurrentP;

		if(kk==0) //OC 220902
		{ 
			xFirstP = xCurrentP; yFirstP = yCurrentP;
		}
		else 
		{ 
			if(((Abs(xCurrentP - xFirstP) < AbsRandX) && (Abs(yCurrentP - yFirstP) < AbsRandY)) ||
			   ((Abs(xCurrentP - xPrevP) < AbsRandX) && (Abs(yCurrentP - yPrevP) < AbsRandY)))
			   Out_PointsAreDifferent = 0;
		}
		xPrevP = xCurrentP; yPrevP = yCurrentP;
	}
	double InvAmOfEdgePoints = 1./AmOfEdgePoints;
	CentrPoint.x = SumX*InvAmOfEdgePoints;
	CentrPoint.y = SumY*InvAmOfEdgePoints;

	int X_In_Count=0, X_Out_Count=0;
	int Y_In_Count=0, Y_Out_Count=0;
	short OutsideGateX=1, OutsideGateY=1;

	radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
	TVector2d First2d = *BaseIter;

	First2d = First2d - CentrPoint;
	double x1 = First2d.x;
	double y1 = First2d.y;
	double x2, y2;

	// Artificial shift (to correctly determine if center point is inside the polygon):
	if(x1==0.) x1 = AbsRandX;
	if(y1==0.) y1 = AbsRandY;

	BaseIter = EdgePointsVector.begin();

	for(int i=0; i<AmOfEdPoInBase; i++)
	{
		++BaseIter;
		if(i!=AmOfEdPoInBase_mi_1)
		{
			x2 = (*BaseIter).x - CentrPoint.x; 
			y2 = (*BaseIter).y - CentrPoint.y;
		}
		else
		{
			x2 = First2d.x; y2 = First2d.y;
		}

		// Artificial shift (to correctly determine if center point is inside the polygon):
		if(x2==0.) x2 = AbsRandX;
		if(y2==0.) y2 = AbsRandY;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);
		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			double k = (y2-y1)/(x2-x1), b = y1 - k*x1;
			if(x1*x2 <= 0.)
			{
				short LocInside = ((x2>x1)? (b<=0) : (b>=0));
				X_In_Count += LocInside? 1 : 0;
				X_Out_Count += (!LocInside)? 1 : 0;

				OutsideGateX = 0;
			}
			if(y1*y2 <= 0.)	OutsideGateY = 0;
		}
		else
		{
			if(y1*y2 <= 0.)
			{
				short LocInside = ((y2>y1)? (x1>=0) : (x1<=0));
				Y_In_Count += LocInside? 1 : 0;
				Y_Out_Count += (!LocInside)? 1 : 0;

				OutsideGateY = 0;
			}
		}
		x1 = x2; y1 = y2;
	}
	if(AmOfEdgePoints <= 3) //OC 080108
	{
		Out_InsideBlock = 1; //OC 080108
		return;
	}

	int X_In = X_In_Count - X_Out_Count;
	int Y_In = X_In_Count - X_Out_Count;
	short InsideBlock = (CheckIfPosEven(X_In) && CheckIfPosEven(Y_In) && (!OutsideGateX) && (!OutsideGateY));

	Out_InsideBlock = InsideBlock;
}

//-------------------------------------------------------------------------

int radTPolygon::CheckIfConvex()
{
	radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
	TVector2d p1 = *BaseIter;
	TVector2d p2 = *(++BaseIter);

	TVector2d v1 = p2 - p1;
	TVector2d v2;
	p1 = p2;

	short IsNonCovex = 0;

	int AmOfEdgePoints_m_1 = AmOfEdgePoints - 1;
	for(int i=1; i<=AmOfEdgePoints; i++)
	{
		if(i==AmOfEdgePoints_m_1) BaseIter = EdgePointsVector.begin();
		else ++BaseIter;

		p2 = *BaseIter;
		v2 = p2 - p1;

		double VectProd = v1.x*v2.y - v2.x*v1.y;
		if(VectProd < -1.E-13)
		{
			IsNonCovex = 1;
		}

		p1 = p2; v1 = v2;
	}
	return (IsNonCovex)? 0 : 1;
}

//-------------------------------------------------------------------------

int radTPolygon::RandomizeNonConvexEdgePoints(radTVect2dVect& LocEdgePointsVector)
{
	radTVect2dVect::iterator BaseIter = LocEdgePointsVector.begin();
	TVector2d p1 = *BaseIter;
	TVector2d p2 = *(++BaseIter);

	TVector2d v1 = p2 - p1;
	TVector2d v2;
	p1 = p2;

	short IsNonCovex = 0;

	int LocAmOfEdgePoints = (int)(LocEdgePointsVector.size());

	int LocAmOfEdgePoints_m_1 = LocAmOfEdgePoints - 1;
	for(int i=1; i<=LocAmOfEdgePoints; i++)
	{
		if(i==LocAmOfEdgePoints_m_1) BaseIter = LocEdgePointsVector.begin();
		else ++BaseIter;

		p2 = *BaseIter;
		v2 = p2 - p1;

		double VectProd = v1.x*v2.y - v2.x*v1.y;
		if(VectProd < -1.E-13)
		{
			IsNonCovex = 1;
		}

		p1 = p2; v1 = v2;
	}
	return (IsNonCovex)? 0 : 1;
}

//-------------------------------------------------------------------------

void radTPolygon::CheckAndRearrangeEdgePoints(TVector2d* InEdgePointsArray, int InAmOfEdgePoints)
{
	int CheckCount = 0;

	TVector2d p1 = InEdgePointsArray[0];
	TVector2d p2 = InEdgePointsArray[1];

	TVector2d v1 = p2 - p1;
	TVector2d v2;
	p1 = p2;

	int AmOfEdgePoints_p_1 = InAmOfEdgePoints + 1;
	for(int i=2; i<=AmOfEdgePoints_p_1; i++)
	{
		p2 = (i==InAmOfEdgePoints)? InEdgePointsArray[0] : ((i!=AmOfEdgePoints_p_1)? InEdgePointsArray[i] : InEdgePointsArray[1]);
		v2 = p2 - p1;

		CheckCount += (v1.x*v2.y - v2.x*v1.y < -1.E-13)? -1 : 1;

		p1 = p2; v1 = v2;
	}

	if(CheckCount<0)
	{
		TVector2d* BufArray = new TVector2d[InAmOfEdgePoints];
		*BufArray = *InEdgePointsArray;
		for(int i=1; i<InAmOfEdgePoints; i++) BufArray[i] = InEdgePointsArray[InAmOfEdgePoints-i];
		for(int k=0; k<InAmOfEdgePoints; k++) InEdgePointsArray[k] = BufArray[k];
		delete[] BufArray;
	}
}

//-------------------------------------------------------------------------

void radTPolygon::CheckAndRearrangeEdgePoints(radTVect2dVect& InEdgePointsArray)
{
	int InAmOfEdgePoints = (int)(InEdgePointsArray.size());
	int CheckCount = 0;

	TVector2d p1 = InEdgePointsArray[0];
	TVector2d p2 = InEdgePointsArray[1];

	TVector2d v1 = p2 - p1;
	TVector2d v2;
	p1 = p2;

	int AmOfEdgePoints_p_1 = InAmOfEdgePoints + 1;
	for(int i=2; i<=AmOfEdgePoints_p_1; i++)
	{
		p2 = (i==InAmOfEdgePoints)? InEdgePointsArray[0] : ((i!=AmOfEdgePoints_p_1)? InEdgePointsArray[i] : InEdgePointsArray[1]);
		v2 = p2 - p1;

		CheckCount += (v1.x*v2.y - v2.x*v1.y < -1.E-13)? -1 : 1;
		p1 = p2; v1 = v2;
	}

	if(CheckCount<0) // Crasy bad style! Change sometime.
	{
		radTVect2dVect BufArray;
		BufArray.push_back(InEdgePointsArray[0]);

		for(int i=1; i<InAmOfEdgePoints; i++) BufArray.push_back(InEdgePointsArray[InAmOfEdgePoints-i]);
		for(int k=0; k<InAmOfEdgePoints; k++) InEdgePointsArray[k] = BufArray[k];
	}
}

//-------------------------------------------------------------------------

int radTPolygon::CheckIfNotSelfIntersecting()
{
	if(AmOfEdgePoints <= 3) return 1; //OC291003

	TVector2d r0Gen = EdgePointsVector[0];
	TVector2d IntrsctPo;
	int AmOfEdgePoints_m_1 = AmOfEdgePoints-1;
	TLinesIntrsctCase IntrsctCase;

	radTSend Send;
	for(int i=0; i<AmOfEdgePoints-2; i++)
	{
		int i_Aux = i+1;
		TVector2d r1Gen = EdgePointsVector[i_Aux++];

		TVector2d r0Tmp = EdgePointsVector[i_Aux];
		for(int k=i_Aux; k<AmOfEdgePoints; k++)
		{
			TVector2d r1Tmp = EdgePointsVector[(k!=AmOfEdgePoints_m_1)? k+1 : 0];
			
			IntrsctOfTwoLines2(r0Gen, r1Gen, r0Tmp, r1Tmp, IntrsctPo, IntrsctCase);
			if(IntrsctCase == PointWithinBound) 
			{ 
				Send.ErrorMessage("Radia::Error104"); return 0;
			}

			r0Tmp = r1Tmp;
		}
		r0Gen = r1Gen;
	}
	return 1;
}

//-------------------------------------------------------------------------
