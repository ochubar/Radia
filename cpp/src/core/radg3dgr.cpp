/*-------------------------------------------------------------------------
*
* File name:      radg3dgr.cpp
*
* Project:        RADIA
*
* Description:    Graphical representations of 3D magnetic field sources
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radg3dgr.h"
#include "radappl.h"
#include "radcast.h"
#include "radexpgn.h"
#include "radvlpgn.h"
#include "radsbdep.h"
#include "radsbdvp.h"
#include "radsbdac.h"

#include <math.h>
#include <stdio.h>

//-------------------------------------------------------------------------

radTSend radTg3dGraphPresent::Send;

//#ifdef _WITH_QD3D
radTVectOfDrawAttr radTg3dGraphPresent::DrawAttrStack;
radRGB radTg3dGraphPresent::SbdLineColor = radRGB(0,0,0); //ensure BLACK lines
//#endif

//-------------------------------------------------------------------------

radTg3dGraphPresent::radTg3dGraphPresent(radTg3d* Ing3dPtr) 
{ 
	g3dPtr = Ing3dPtr; DrawAttrAreSet = 0;

	TVector3d St0(1.,0.,0.), St1(0.,1.,0.), St2(0.,0.,1.);
	TMatrix3d E(St0, St1, St2);
	TVector3d ZeroVect(0.,0.,0.);
	radTrans IdentTrans(E, E, ZeroVect, 1., 1.);
	GenTrans = IdentTrans;

	ShowInternalFacesAfterCut = true;
	GraphPresOptions.ShowSymmetryChilds = 1;
	GraphPresOptions.doDebug = 0;

	DrawFacilityInd = 0;
	ShowGridTicks = 1;
}

//-------------------------------------------------------------------------

void radTg3dGraphPresent::NestedFor_Draw(radTrans* BaseTransPtr, const radTlphg::iterator& Iter)
{
	radTrans* TransPtr = (radTrans*)(((*Iter).Handler_g).rep);
	radTlphg::iterator LocalNextIter = Iter;
	LocalNextIter++;

	radTrans LocTotTrans = *BaseTransPtr;

	if((*Iter).m == 1)
	{
		LocTotTrans = Product(LocTotTrans, *TransPtr);
		DrawOrNestedFor(&LocTotTrans, LocalNextIter);
	}
	else
	{
		DrawOrNestedFor(&LocTotTrans, LocalNextIter);
		int Mult = (*Iter).m;
		for(int km = 1; km < Mult; km++)
		{
			LocTotTrans = Product(LocTotTrans, *TransPtr);
			DrawOrNestedFor(&LocTotTrans, LocalNextIter);
		}
	}
}

//-------------------------------------------------------------------------

void radTg3dGraphPresent::GenDraw()
{
//#ifdef _WITH_QD3D
	//if(DrawAttrAreSet) DrawAttrStack.push_back(DrawAttr);
	if(DrawAttrAreSet && (!GraphPresOptions.doDebug)) DrawAttrStack.push_back(DrawAttr); //OC130709
//#endif

	//radTCast Cast;
	if(GraphPresOptions.ShowSymmetryChilds)
	{
		int NumberOfSymChilds_pl_Orig = 1;
		for(radTlphg::iterator iter = g3dPtr->g3dListOfTransform.begin();
			iter != g3dPtr->g3dListOfTransform.end(); ++iter) NumberOfSymChilds_pl_Orig *= (*iter).m;

		if(radTCast::FlmLinCurCast(g3dPtr) != NULL) // To separate Linear primitives (since DrawAttr apply differently to them)
			Send.InitDrawLinElem(DrawAttrAreSet, DrawAttr, NumberOfSymChilds_pl_Orig, DrawFacilityInd);
		else Send.InitDrawSurfElem(DrawAttrAreSet, DrawAttr, NumberOfSymChilds_pl_Orig, DrawFacilityInd);

		if(g3dPtr->g3dListOfTransform.empty()) Draw(&GenTrans);
		else NestedFor_Draw(&GenTrans, g3dPtr->g3dListOfTransform.begin());
	}
	else
	{
		if(radTCast::FlmLinCurCast(g3dPtr) != NULL) // To separate Linear primitives (since DrawAttr apply diferently to them)
			Send.InitDrawLinElem(DrawAttrAreSet, DrawAttr, 1, DrawFacilityInd);
		else Send.InitDrawSurfElem(DrawAttrAreSet, DrawAttr, 1, DrawFacilityInd);
		Draw(&GenTrans);
	}

//#ifdef _WITH_QD3D
	//if(DrawAttrAreSet) DrawAttrStack.pop_back();
	if(DrawAttrAreSet && (!GraphPresOptions.doDebug)) DrawAttrStack.pop_back(); //OC130709
//#endif
}

//-------------------------------------------------------------------------

int radTg3dGraphPresent::RetrieveDrawAttr(int ElemKey) 
{ 
	radTMapOfDrawAttr::iterator iter = MapOfDrawAttrPtr->find(ElemKey);
	if(!(iter == MapOfDrawAttrPtr->end()))
	{
		DrawAttr = (*iter).second; 
		DrawAttrAreSet = 1; 
		return 1;
	}
	else return 0;
}

//-------------------------------------------------------------------------

void radTg3dGraphPresent::SetupRotation(const TVector3d& PoiOnAxVect, const TVector3d& InAxVect, double Angle, radTrans& Rotation)
{
	double NormFact = 1./sqrt(InAxVect.x*InAxVect.x+InAxVect.y*InAxVect.y+InAxVect.z*InAxVect.z);
	TVector3d AxVect = NormFact*InAxVect;
	double VxVx, VyVy, VzVz;
	VxVx=AxVect.x*AxVect.x; VyVy=AxVect.y*AxVect.y; VzVz=AxVect.z*AxVect.z;

	double cosAng, sinAng, One_m_cos;
	cosAng = cos(Angle); sinAng = sin(Angle); One_m_cos = 1. - cosAng;
	double One_m_cosVxVy, One_m_cosVxVz, One_m_cosVyVz, sinVx, sinVy, sinVz;
	One_m_cosVxVy = One_m_cos*AxVect.x*AxVect.y;
	One_m_cosVxVz = One_m_cos*AxVect.x*AxVect.z;
	One_m_cosVyVz = One_m_cos*AxVect.y*AxVect.z;
	sinVx = sinAng*AxVect.x; sinVy = sinAng*AxVect.y; sinVz = sinAng*AxVect.z;

	TVector3d St0(VxVx+cosAng*(VyVy+VzVz), One_m_cosVxVy-sinVz, One_m_cosVxVz+sinVy);
	TVector3d St1(One_m_cosVxVy+sinVz, VyVy+cosAng*(VxVx+VzVz), One_m_cosVyVz-sinVx);
	TVector3d St2(One_m_cosVxVz-sinVy, One_m_cosVyVz+sinVx, VzVz+cosAng*(VxVx+VyVy));
	TMatrix3d M(St0, St1, St2);
	TVector3d St00(1.-St0.x, -St0.y, -St0.z);
	TVector3d St01(-St1.x, 1.-St1.y, -St1.z);
	TVector3d St02(-St2.x, -St2.y, 1.-St2.z);
	TMatrix3d M0(St00, St01, St02);

	Rotation = radTrans(M, M0*PoiOnAxVect, 1., 1.);
}

//-------------------------------------------------------------------------

void radTg3dGraphPresent::SetupPlaneSym(const TVector3d& PoiOnPlaneVect, const TVector3d& InN, radTrans& PlaneSym)
{
	double NormFact = 1./sqrt(InN.x*InN.x+InN.y*InN.y+InN.z*InN.z);
	TVector3d N = NormFact*InN;

	TVector3d St0(1.-2.*N.x*N.x, -2.*N.x*N.y, -2.*N.x*N.z);
	TVector3d St1(St0.y, 1.-2.*N.y*N.y, -2.*N.y*N.z);
	TVector3d St2(St0.z, St1.z, 1.-2.*N.z*N.z);
	TMatrix3d M(St0, St1, St2);
	TVector3d St00(1.-St0.x, -St0.y, -St0.z);
	TVector3d St01(-St1.x, 1.-St1.y, -St1.z);
	TVector3d St02(-St2.x, -St2.y, 1.-St2.z);
	TMatrix3d M0(St00, St01, St02);

	PlaneSym = radTrans(M, M0*PoiOnPlaneVect, -1., 1., 3); // ID_No = 3
}

//-------------------------------------------------------------------------

void radTg3dGraphPresent::SetupTranslation(const TVector3d& StPoi, const TVector3d& EndPoi, radTrans& Translation)
{
	TVector3d St0(1.,0.,0.), St1(0.,1.,0.), St2(0.,0.,1.);
	TMatrix3d E(St0, St1, St2);
	TVector3d TranslVect = EndPoi - StPoi;

	Translation = radTrans(E, E, TranslVect, 1., 1.);
}

//-------------------------------------------------------------------------

void radTg3dGraphPresent::DrawFrameLines()
{//to call after g3dGraphPresentPtr->GenDraw() was called
	char PrevShowLines = Send.ShowLines;
    Send.ShowLines = 1;

	const double RelFrameLinesOffset = 0.07;

	const double RelArrowHeight = 0.08;
	const double ArrowBottomRadToHeightRatio = 0.15;

	const double RelCharHeight = 0.05;
	const double CharWidthToHeigthRatio = 0.6;

	radRGB FrameLinesColor(1,1,1); //Frame lines are white
	radTDrawAttr FrameLinesDrawAttr;
	FrameLinesDrawAttr.RGB_col = FrameLinesColor;
	DrawAttrStack.push_back(FrameLinesDrawAttr);

	double OrigLimits3D[6];
	double *tOrigLimits3D = OrigLimits3D, *tLimits3D = Send.Limits3D;
	for(int k=0; k<6; k++) *(tOrigLimits3D++) = *(tLimits3D++);

	double Sx = OrigLimits3D[1] - OrigLimits3D[0];
	double Sy = OrigLimits3D[3] - OrigLimits3D[2];
	double Sz = OrigLimits3D[5] - OrigLimits3D[4];
	double Smax = (Sx > Sy)? ((Sx > Sz)? Sx : Sz) : ((Sy > Sz)? Sy : Sz);

	double Dx = Sx*RelFrameLinesOffset;
	double Dy = Sy*RelFrameLinesOffset;
	double Dz = Sz*RelFrameLinesOffset;

// Frame contour lines
	TVector3d Side[5], Segment[2];
	Side[0].x = OrigLimits3D[0] - Dx; Side[0].y = OrigLimits3D[2] - Dy; Side[0].z = OrigLimits3D[4] - Dz; 
	Side[1].x = Side[0].x; Side[1].y = Side[0].y; Side[1].z = OrigLimits3D[5] + Dz; 
	Side[2].x = Side[0].x; Side[2].y = OrigLimits3D[3] + Dy; Side[2].z = Side[1].z; 
	Side[3].x = Side[0].x; Side[3].y = Side[2].y; Side[3].z = Side[0].z;
	Side[4] = Side[0];
	Segment[0] = Side[0];
	Send.Line(Side, 5, DrawFacilityInd);
	Side[0].x = OrigLimits3D[1] + Dx;
	Side[4].x = Side[3].x = Side[2].x = Side[1].x = Side[0].x;
	Send.Line(Side, 5, DrawFacilityInd);

	Segment[1] = Side[0];
	Send.Line(Segment, 2, DrawFacilityInd);
	Segment[0].y = Side[1].y; Segment[0].z = Side[1].z;
	Segment[1] = Side[1];
	Send.Line(Segment, 2, DrawFacilityInd);
	Segment[0].y = Side[2].y; Segment[0].z = Side[2].z;
	Segment[1] = Side[2];
	Send.Line(Segment, 2, DrawFacilityInd);
	Segment[0].y = Side[3].y; Segment[0].z = Side[3].z;
	Segment[1] = Side[3];
	Send.Line(Segment, 2, DrawFacilityInd);

// Frame arrows
	double ArrowHeight = RelArrowHeight*Smax;
	double ArrowBottomRad = ArrowBottomRadToHeightRatio*ArrowHeight;

	TVector3d PyramidArrowInfo[4];
	TVector3d ArrowHeightVect(ArrowHeight,0.,0.);
	TVector3d ArrowBottomRad1(0., -ArrowBottomRad, 0.), ArrowBottomRad2(0., 0., -ArrowBottomRad);
	if(ArrowHeight < Sx)
	{
		PyramidArrowInfo[0].x = 0.5*(OrigLimits3D[0] + OrigLimits3D[1]) + 0.5*ArrowHeight;
		PyramidArrowInfo[0].y = OrigLimits3D[2] - Dy;
		PyramidArrowInfo[0].z = OrigLimits3D[4] - Dz;
		PyramidArrowInfo[1] = ArrowHeightVect;
		PyramidArrowInfo[2] = ArrowBottomRad1; PyramidArrowInfo[3] = ArrowBottomRad2;
		DrawPyramidArrow(PyramidArrowInfo, DrawFacilityInd);
	}
	if(ArrowHeight < Sy)
	{
		PyramidArrowInfo[0].x = OrigLimits3D[1] + Dx;
		PyramidArrowInfo[0].y = 0.5*(OrigLimits3D[2] + OrigLimits3D[3]) + 0.5*ArrowHeight;
		PyramidArrowInfo[0].z = OrigLimits3D[4] - Dz;
		PyramidArrowInfo[1].x = 0.; PyramidArrowInfo[1].y = ArrowHeight; PyramidArrowInfo[1].z = 0.;
		PyramidArrowInfo[2].x = ArrowBottomRad; PyramidArrowInfo[2].y = PyramidArrowInfo[2].z = 0.;
		PyramidArrowInfo[3].x = PyramidArrowInfo[3].y = 0.; PyramidArrowInfo[3].z = -ArrowBottomRad;
		DrawPyramidArrow(PyramidArrowInfo, DrawFacilityInd);
	}
	if(ArrowHeight < Sz)
	{
		PyramidArrowInfo[0].x = OrigLimits3D[1] + Dx;
		PyramidArrowInfo[0].y = OrigLimits3D[3] + Dy;
		PyramidArrowInfo[0].z = 0.5*(OrigLimits3D[4] + OrigLimits3D[5]) + 0.5*ArrowHeight;
		PyramidArrowInfo[1].x = PyramidArrowInfo[1].y = 0.; PyramidArrowInfo[1].z = ArrowHeight;
		PyramidArrowInfo[2].x = ArrowBottomRad; PyramidArrowInfo[2].y = PyramidArrowInfo[2].z = 0.;
		PyramidArrowInfo[3].x = 0.; PyramidArrowInfo[3].y = ArrowBottomRad; PyramidArrowInfo[3].z = 0.;
		DrawPyramidArrow(PyramidArrowInfo, DrawFacilityInd);
	}

// Frame characters
	double AbsCharHeigth = RelCharHeight*Smax;
	double AbsCharWidth = CharWidthToHeigthRatio*AbsCharHeigth;

	TVector3d CharInfo3D[3];
	CharInfo3D[0].x = 0.5*(OrigLimits3D[0] + OrigLimits3D[1]) - 0.5*AbsCharWidth;
	CharInfo3D[0].y = OrigLimits3D[2] - Dy;
	CharInfo3D[0].z = OrigLimits3D[4] - Dz - 1.5*AbsCharHeigth;
	CharInfo3D[1].x = 0.; CharInfo3D[1].y = -1.; CharInfo3D[1].z = 0.;
	CharInfo3D[2].x = 0.; CharInfo3D[2].y = 0.; CharInfo3D[2].z = AbsCharHeigth;
	DrawCharacter('X', CharWidthToHeigthRatio, CharInfo3D, DrawFacilityInd);

	CharInfo3D[0].x = OrigLimits3D[1] + Dx;
	CharInfo3D[0].y = 0.5*(OrigLimits3D[2] + OrigLimits3D[3]) - 0.5*AbsCharWidth;
	CharInfo3D[0].z = OrigLimits3D[4] - Dz - 1.5*AbsCharHeigth;
	CharInfo3D[1].x = 1.; CharInfo3D[1].y = 0.; CharInfo3D[1].z = 0.;
	CharInfo3D[2].x = 0.; CharInfo3D[2].y = 0.; CharInfo3D[2].z = AbsCharHeigth;
	DrawCharacter('Y', CharWidthToHeigthRatio, CharInfo3D, DrawFacilityInd);

	CharInfo3D[0].x = OrigLimits3D[1] + Dx;
	CharInfo3D[0].y = OrigLimits3D[3] + Dy + AbsCharWidth;
	CharInfo3D[0].z = 0.5*(OrigLimits3D[4] + OrigLimits3D[5]) - 0.5*AbsCharHeigth;
	CharInfo3D[1].x = 1.; CharInfo3D[1].y = 0.; CharInfo3D[1].z = 0.;
	CharInfo3D[2].x = 0.; CharInfo3D[2].y = 0.; CharInfo3D[2].z = AbsCharHeigth;
	DrawCharacter('Z', CharWidthToHeigthRatio, CharInfo3D, DrawFacilityInd);

// Frame ticks and labels
	if(ShowGridTicks)
	{
		TVector3d OffsetVect(Dx, Dy, Dz);
		DrawGridTicks(OrigLimits3D, OffsetVect, AbsCharHeigth, DrawFacilityInd);
	}

	DrawAttrStack.pop_back();
	Send.ShowLines = PrevShowLines;
}

//-------------------------------------------------------------------------

void radTg3dGraphPresent::DrawGridTicks(double* OrigLimits3D, TVector3d& OffsetVect, double AbsCharHeight, char LocDrawFacilityInd)
{
    TVector3d GridLimits[2];
    GridLimits->x = OrigLimits3D[0] - OffsetVect.x;
	GridLimits->y = OrigLimits3D[2] - OffsetVect.y;
	GridLimits->z = OrigLimits3D[4] - OffsetVect.z;
	GridLimits[1].x = OrigLimits3D[1] + OffsetVect.x;
	GridLimits[1].y = OrigLimits3D[3] + OffsetVect.y;
	GridLimits[1].z = OrigLimits3D[5] + OffsetVect.z;

	vector<double> MainTickPositions[3];
    FindGridTickPositions(GridLimits, MainTickPositions);
	double TickLength = FindGridTickLength(GridLimits);

	TVector3d Ex(1,0,0), Ey(0,1,0), Ez(0,0,1), mEx(-1,0,0), mEy(0,-1,0), mEz(0,0,-1);

    TVector3d AxesDirections[12], AxeTr1[12], AxeTr2[12];
	AxesDirections[0] = Ex; AxeTr1[0] = Ey; AxeTr2[0] = Ez;
	AxesDirections[1] = Ex; AxeTr1[1] = mEy; AxeTr2[1] = Ez;
	AxesDirections[2] = Ex; AxeTr1[2] = mEy; AxeTr2[2] = mEz;
	AxesDirections[3] = Ex; AxeTr1[3] = Ey; AxeTr2[3] = mEz;
	AxesDirections[4] = Ey; AxeTr1[4] = Ex; AxeTr2[4] = Ez;
	AxesDirections[5] = Ey; AxeTr1[5] = Ex; AxeTr2[5] = mEz;
	AxesDirections[6] = Ey; AxeTr1[6] = mEx; AxeTr2[6] = mEz;
	AxesDirections[7] = Ey; AxeTr1[7] = mEx; AxeTr2[7] = Ez;
	AxesDirections[8] = Ez; AxeTr1[8] = Ex; AxeTr2[8] = Ey;
	AxesDirections[9] = Ez; AxeTr1[9] = mEx; AxeTr2[9] = Ey;
	AxesDirections[10] = Ez; AxeTr1[10] = mEx; AxeTr2[10] = mEy;
	AxesDirections[11] = Ez; AxeTr1[11] = Ex; AxeTr2[11] = mEy;

	TVector3d P0(GridLimits[0]), P1(GridLimits[0].x, GridLimits[1].y, GridLimits[0].z), P2(GridLimits[0].x, GridLimits[1].y, GridLimits[1].z), P3(GridLimits[0].x, GridLimits[0].y, GridLimits[1].z);
	TVector3d P4(GridLimits[1].x, GridLimits[0].y, GridLimits[0].z), P5(GridLimits[1].x, GridLimits[1].y, GridLimits[0].z), P6(GridLimits[1]), P7(GridLimits[1].x, GridLimits[0].y, GridLimits[1].z);

    TVector3d AxesStartPoints[12];
	AxesStartPoints[0] = P0; AxesStartPoints[1] = P1; AxesStartPoints[2] = P2; AxesStartPoints[3] = P3;
	AxesStartPoints[4] = P0; AxesStartPoints[5] = P3; AxesStartPoints[6] = P7; AxesStartPoints[7] = P4;
	AxesStartPoints[8] = P0; AxesStartPoints[9] = P4; AxesStartPoints[10] = P5; AxesStartPoints[11] = P1;

	int IndsTickPos[] = {0,0,0,0, 1,1,1,1, 2,2,2,2};
	double FrameInitPositions[] = {GridLimits->x,GridLimits->x,GridLimits->x,GridLimits->x, GridLimits->y,GridLimits->y,GridLimits->y,GridLimits->y, GridLimits->z,GridLimits->z,GridLimits->z,GridLimits->z};

	double MinTickInterval = 1.e+23;

	for(int i=0; i<12; i++)
	{
		TVector3d& CurStartPt = AxesStartPoints[i];
		TVector3d& CurAxDir = AxesDirections[i];
		TVector3d CurAxTr1 = TickLength*(AxeTr1[i]);
		TVector3d CurAxTr2 = TickLength*(AxeTr2[i]);

		vector<double>& CurMainTickPositions = MainTickPositions[IndsTickPos[i]];
		int CurNumTicks = (int)(CurMainTickPositions.size());
		double CurFrameInitPos = FrameInitPositions[i];
		TVector3d TickLine[3];

		double CurTickInterv = 1.e+23;
		if(CurNumTicks > 1) CurTickInterv = CurMainTickPositions[1] - CurMainTickPositions[0];
		if(MinTickInterval > CurTickInterv) MinTickInterval = CurTickInterv;

		for(int k=0; k<CurNumTicks; k++)
		{
			double CurRelTickPos = CurMainTickPositions[k] - CurFrameInitPos;

			TVector3d CurTickPos = (CurStartPt + (CurRelTickPos*CurAxDir));
			TickLine[0] = CurTickPos + CurAxTr1;
			TickLine[1] = CurTickPos;
			TickLine[2] = CurTickPos + CurAxTr2;

			Send.Line(TickLine, 3, LocDrawFacilityInd);
		}
	}

	AbsCharHeight *= 0.7;

	double TickIntervalX = MinTickInterval, TickIntervalY = MinTickInterval, TickIntervalZ = MinTickInterval;
	if(MainTickPositions[0].size() > 1) TickIntervalX = (MainTickPositions[0])[1] - (MainTickPositions[0])[0];
	if(MainTickPositions[1].size() > 1) TickIntervalY = (MainTickPositions[1])[1] - (MainTickPositions[1])[0];
	if(MainTickPositions[2].size() > 1) TickIntervalZ = (MainTickPositions[2])[1] - (MainTickPositions[2])[0];

	double NumCharHeightX = AbsCharHeight, NumCharHeightY = AbsCharHeight, NumCharHeightZ = AbsCharHeight;
	double MaxNumCharHeightX = 0.7*TickIntervalX; if(NumCharHeightX > MaxNumCharHeightX) NumCharHeightX = MaxNumCharHeightX;
	double MaxNumCharHeightY = 0.7*TickIntervalY; if(NumCharHeightY > MaxNumCharHeightY) NumCharHeightY = MaxNumCharHeightY;
	double MaxNumCharHeightZ = 0.7*TickIntervalZ; if(NumCharHeightZ > MaxNumCharHeightZ) NumCharHeightZ = MaxNumCharHeightZ;

	DrawTickNumbers(P3, Ex, Ez, -1, MainTickPositions[0], GridLimits->x, NumCharHeightX, LocDrawFacilityInd);
	DrawTickNumbers(P7, Ey, Ez, -1, MainTickPositions[1], GridLimits->y, NumCharHeightY, LocDrawFacilityInd);
	DrawTickNumbers(P4, Ez, mEy, -1, MainTickPositions[2], GridLimits->z, NumCharHeightZ, LocDrawFacilityInd);
}

//-------------------------------------------------------------------------

void radTg3dGraphPresent::DrawTickNumbers(TVector3d& P0, TVector3d& Vpar, TVector3d& Vper, int PerDir, vector<double>& TickNumPositions, double TickOffsetPos, double AbsCharHeight, char LocDrawFacilityInd)
{
	const int AmOfNumChar = 6;
	const double AmOfCharOffset = 0.7;
	const double WidthToHeightRat = 0.6;
	const double IntervToWidthRat = 0.3;

	double AbsCharWidth = WidthToHeightRat*AbsCharHeight;
	double HalfAllCharWidth = 0.5*AmOfNumChar*(1 + IntervToWidthRat)*AbsCharWidth;
	double AbsOffsetFromAx = AmOfCharOffset*AbsCharWidth + HalfAllCharWidth;

	TVector3d UnitVpar = (1./sqrt(Vpar.x*Vpar.x + Vpar.y*Vpar.y + Vpar.z*Vpar.z))*Vpar;
	TVector3d UnitVper = (1./sqrt(Vper.x*Vper.x + Vper.y*Vper.y + Vper.z*Vper.z))*Vper;
	TVector3d StartCenP = P0 + (AbsOffsetFromAx*UnitVper); 

	TVector3d CharInfo3D[3];
	CharInfo3D[1] = UnitVpar^Vper;
	CharInfo3D[2] = AbsCharHeight*UnitVpar;

	TVector3d WidthVect = AbsCharWidth*UnitVper;
	TVector3d HalfHeightVect = (0.5*AbsCharHeight)*UnitVpar;

	TVector3d VectFromCenToLowerLeft = (HalfAllCharWidth*UnitVper) - HalfHeightVect;
	TVector3d CharTranslVect = ((-1 - IntervToWidthRat)*AbsCharWidth)*UnitVper;
	char FormStr[6];
	strcpy(FormStr, "%");
	if(PerDir > 0) 
	{
		VectFromCenToLowerLeft = (-1)*VectFromCenToLowerLeft;
		CharTranslVect = (-1)*CharTranslVect;
		strcat(FormStr, "-");
	}
	TVector3d StartLowerLeftP = StartCenP + VectFromCenToLowerLeft;

	char ExtraStrBuf[3];
	sprintf(ExtraStrBuf, "%d", AmOfNumChar);
	strcat(FormStr, ExtraStrBuf);
	strcat(FormStr, "g");

	char NumStrBuf[10];
	int AmOfTicks = (int)(TickNumPositions.size());
	for(int i=0; i<AmOfTicks; i++)
	{
		double CurTickNum = TickNumPositions[i];
		double CurTickOffset = CurTickNum - TickOffsetPos;
		TVector3d CurLowerLeftP = StartLowerLeftP + (CurTickOffset*UnitVpar);

		int AmOfChar = sprintf(NumStrBuf, FormStr, CurTickNum);
		char *tStrBuf = NumStrBuf;

		for(int k=0; k<AmOfChar; k++)
		{
			CharInfo3D[0] = CurLowerLeftP;
			DrawCharacter(*tStrBuf, WidthToHeightRat, CharInfo3D, LocDrawFacilityInd);

			tStrBuf++;
			CurLowerLeftP += CharTranslVect;
		}
	}
}

//-------------------------------------------------------------------------

void radTg3dGraphPresent::DrawPyramidArrow(TVector3d* PyramidArrowInfo, char LocDrawFacilityInd)
{// PyramidArrowInfo[0] - main vertex of pyramid
// PyramidArrowInfo[1] - height vector of pyramid
// PyramidArrowInfo[2] - bottom radius_1 vect
// PyramidArrowInfo[3] - bottom radius_2 vect

	TVector3d ArrowVerices[5], Face[4];
	*ArrowVerices = *PyramidArrowInfo;
	TVector3d ArrowOrigin = PyramidArrowInfo[0] - PyramidArrowInfo[1];
	TVector3d& ArrowBottomRad1 = PyramidArrowInfo[2];
	TVector3d& ArrowBottomRad2 = PyramidArrowInfo[3];
	ArrowVerices[1] = ArrowOrigin + ArrowBottomRad1;
	ArrowVerices[2] = ArrowOrigin + ArrowBottomRad2;
	ArrowVerices[3] = ArrowOrigin - ArrowBottomRad1;
	ArrowVerices[4] = ArrowOrigin - ArrowBottomRad2;
	Face[0] = ArrowVerices[0]; Face[1] = ArrowVerices[1]; Face[2] = ArrowVerices[2];
	Send.Polygon(Face, 3, LocDrawFacilityInd);
	Face[1] = ArrowVerices[2]; Face[2] = ArrowVerices[3];
	Send.Polygon(Face, 3, LocDrawFacilityInd);
	Face[1] = ArrowVerices[3]; Face[2] = ArrowVerices[4];
	Send.Polygon(Face, 3, LocDrawFacilityInd);
	Face[1] = ArrowVerices[4]; Face[2] = ArrowVerices[1];
	Send.Polygon(Face, 3, LocDrawFacilityInd);
	Face[0] = ArrowVerices[1]; Face[1] = ArrowVerices[4]; Face[2] = ArrowVerices[3]; Face[3] = ArrowVerices[2];
	Send.Polygon(Face, 4, LocDrawFacilityInd);
}

//-------------------------------------------------------------------------

void radTg3dGraphPresent::DrawCharacter(char Ch, double Ratio, TVector3d* Info3D, char LocDrawFacilityInd)
{//Don't use this ! Use radTg3dGraphPresent::DrawCharacter() instead !
// Info3D[0], Info3D[1] - point and normal vector of a plane to which the character belongs
// Info3D[2] - height vector of a character (the Height is the length of Info3D[2])
// Width = Ratio*Height

	TVector3d& Normal = Info3D[1];
	Normal = (1./sqrt(Normal.x*Normal.x + Normal.y*Normal.y + Normal.z*Normal.z))*Normal;

	TVector3d& HeightVect = Info3D[2];
	TVector3d LowerLeft = *Info3D;
	TVector3d WidthVect = Ratio*(HeightVect^Normal);
	TVector3d LowerRight = LowerLeft + WidthVect;
	TVector3d UpperLeft = LowerLeft + HeightVect;
	TVector3d UpperRight = UpperLeft + WidthVect;
	TVector3d Middle, SmallHorVect, SmallVertVect;

	TVector3d Segment[6];
	Ch = (char)toupper(Ch);

	switch(Ch) 
	{
        case 'E':
			Segment[0] = UpperRight; Segment[1] = UpperLeft; Segment[2] = LowerLeft; Segment[3] = LowerRight;
			Send.Line(Segment, 4, LocDrawFacilityInd);
			Segment[0] = LowerLeft + (0.5*HeightVect); Segment[1] = LowerRight + (0.5*HeightVect);
			Send.Line(Segment, 2, LocDrawFacilityInd);
			break;
		case 'X':
			Segment[0] = LowerLeft; Segment[1] = UpperRight;
			Send.Line(Segment, 2, LocDrawFacilityInd);
			Segment[0] = UpperLeft; Segment[1] = LowerRight;
			Send.Line(Segment, 2, LocDrawFacilityInd);
			break;
        case 'Y':
			Middle = 0.25*(LowerLeft + UpperLeft + LowerRight + UpperRight);
			Segment[0] = UpperLeft; Segment[1] = Middle; Segment[2] = UpperRight;
			Send.Line(Segment, 3, LocDrawFacilityInd);
			Segment[0] = 0.5*(LowerLeft + LowerRight); Segment[1] = Middle;
			Send.Line(Segment, 2, LocDrawFacilityInd);
			break;
        case 'Z':
			Segment[0] = UpperLeft; Segment[1] = UpperRight; Segment[2] = LowerLeft; Segment[3] = LowerRight;
			Send.Line(Segment, 4, LocDrawFacilityInd);
			break;
        case '0':
			Segment[0] = LowerLeft; Segment[1] = LowerRight; Segment[2] = UpperRight; Segment[3] = UpperLeft; Segment[4] = LowerLeft; 
			Send.Line(Segment, 5, LocDrawFacilityInd);
			break;
        case '1':
			Segment[0] = LowerRight; Segment[1] = UpperRight; Segment[2] = LowerLeft + (0.5*HeightVect);
			Send.Line(Segment, 3, LocDrawFacilityInd);
			break;
        case '2':
			Segment[0] = UpperLeft; Segment[1] = UpperRight; Segment[2] = LowerRight + (0.5*HeightVect); Segment[3] = LowerLeft; Segment[4] = LowerRight; 
			Send.Line(Segment, 5, LocDrawFacilityInd);
			break;
        case '3':
			Segment[0] = UpperLeft; Segment[1] = UpperRight; Segment[2] = LowerLeft + (0.5*HeightVect); Segment[3] = LowerRight + (0.5*HeightVect); Segment[4] = LowerLeft; 
			Send.Line(Segment, 5, LocDrawFacilityInd);
			break;
        case '4':
			Segment[0] = UpperLeft; Segment[1] = LowerLeft + (0.5*HeightVect); Segment[2] = LowerRight + (0.5*HeightVect); 
			Send.Line(Segment, 3, LocDrawFacilityInd);
			Segment[0] = LowerRight; Segment[1] = UpperRight; 
			Send.Line(Segment, 2, LocDrawFacilityInd);
			break;
        case '5':
			Segment[0] = LowerLeft; Segment[1] = LowerRight; Segment[2] = LowerRight + (0.5*HeightVect); Segment[3] = LowerLeft + (0.5*HeightVect); Segment[4] = UpperLeft; Segment[5] = UpperRight;
			Send.Line(Segment, 6, LocDrawFacilityInd);
			break;
        case '6':
			Segment[0] = UpperRight; Segment[1] = LowerLeft + (0.5*HeightVect); Segment[2] = LowerLeft; Segment[3] = LowerRight; Segment[4] = LowerRight + (0.5*HeightVect); Segment[5] = Segment[1];
			Send.Line(Segment, 6, LocDrawFacilityInd);
			break;
        case '7':
			Segment[0] = UpperLeft; Segment[1] = UpperRight; Segment[2] = LowerLeft + (0.5*HeightVect); Segment[3] = LowerLeft; 
			Send.Line(Segment, 4, LocDrawFacilityInd);
			break;
        case '8':
			Segment[0] = LowerLeft; Segment[1] = LowerRight; Segment[2] = UpperRight; Segment[3] = UpperLeft; Segment[4] = LowerLeft; 
			Send.Line(Segment, 5, LocDrawFacilityInd);
			Segment[0] = LowerLeft + (0.5*HeightVect); Segment[1] = LowerRight + (0.5*HeightVect);
			Send.Line(Segment, 2, LocDrawFacilityInd);
			break;
        case '9':
			Segment[0] = LowerRight + (0.5*HeightVect); Segment[1] = LowerLeft + (0.5*HeightVect); Segment[2] = UpperLeft; Segment[3] = UpperRight; Segment[4] = Segment[0]; Segment[5] = LowerLeft;
			Send.Line(Segment, 6, LocDrawFacilityInd);
			break;
        case '.':
			SmallHorVect = 0.05*WidthVect;
            SmallVertVect = (sqrt((SmallHorVect.x*SmallHorVect.x + SmallHorVect.y*SmallHorVect.y + SmallHorVect.z*SmallHorVect.z)/(HeightVect.x*HeightVect.x + HeightVect.y*HeightVect.y + HeightVect.z*HeightVect.z)))*HeightVect;
			Middle = 0.5*(LowerLeft + LowerRight);
			Segment[0] = Middle - (SmallHorVect + SmallVertVect); Segment[1] = (Middle - SmallVertVect) + SmallHorVect; Segment[2] = Middle + SmallVertVect + SmallHorVect;
            Segment[3] = (Middle + SmallVertVect) - SmallHorVect; Segment[4] = Segment[0]; 
			Send.Line(Segment, 5, LocDrawFacilityInd);
			break;
        case '-':
			Middle = 0.25*(LowerLeft + UpperLeft + LowerRight + UpperRight);
			SmallHorVect = 0.4*WidthVect;
			Segment[0] = Middle - SmallHorVect; Segment[1] = Middle + SmallHorVect;
			Send.Line(Segment, 2, LocDrawFacilityInd);
			break;
        case '+':
			Middle = 0.25*(LowerLeft + UpperLeft + LowerRight + UpperRight);
			SmallHorVect = 0.4*WidthVect;
            SmallVertVect = (sqrt((SmallHorVect.x*SmallHorVect.x + SmallHorVect.y*SmallHorVect.y + SmallHorVect.z*SmallHorVect.z)/(HeightVect.x*HeightVect.x + HeightVect.y*HeightVect.y + HeightVect.z*HeightVect.z)))*HeightVect;
			Segment[0] = Middle - SmallHorVect; Segment[1] = Middle + SmallHorVect;
			Send.Line(Segment, 2, LocDrawFacilityInd);
			Segment[0] = Middle - SmallVertVect; Segment[1] = Middle + SmallVertVect;
			Send.Line(Segment, 2, LocDrawFacilityInd);
			break;
	}
}

//-------------------------------------------------------------------------

void radTg3dGraphPresent::FindGridTickPositions(TVector3d* GridLimits, vector<double>* MainTickPositions)
{//Calculates absolute tick positions
	double Lx = GridLimits[1].x - GridLimits[0].x;
	double Ly = GridLimits[1].y - GridLimits[0].y;
	double Lz = GridLimits[1].z - GridLimits[0].z;

	double Dims[] = {Lx, Ly, Lz};
	double StartLim[] = { GridLimits[0].x,  GridLimits[0].y,  GridLimits[0].z};
	double EndLim[] = { GridLimits[1].x,  GridLimits[1].y,  GridLimits[1].z};

	for(int i=0; i<3; i++)
	{
        double TickDelta = ChooseGridTickInterval(Dims[i]);
		double CurStartLim = StartLim[i];
        double TickPos = FindGridTickFirstPosition(CurStartLim, TickDelta);
		double CurEndPos = EndLim[i];

		vector<double>& CurTickPosVect = MainTickPositions[i];
		for(int k=0; k<1000; k++)
		{
			CurTickPosVect.push_back(TickPos);
			TickPos += TickDelta;

			if(::fabs(TickPos) < 0.001*TickDelta) TickPos = 0;

			if(TickPos > CurEndPos) break;
		}
	}
}

//-------------------------------------------------------------------------

double radTg3dGraphPresent::FindGridTickFirstPosition(double Start, double TickDelta)
{
	double dInit = floor(Start/TickDelta) + 1;
	//if(Start < 0) dInit += 1.;
	return dInit*TickDelta;
}

//-------------------------------------------------------------------------

double radTg3dGraphPresent::ChooseGridTickInterval(double Range)
{
	double PossibleNormDeltas[] = {1, 2, 5};
	int LenPossibleDeltas = 3;
	const int AvgNumMainTicks = 6; //to tune

	double StrictDelta = Range/((double)(AvgNumMainTicks - 1));

	double MinFractOrder = 1.e+23;
	long iResOrder = 0, IndNormDelta = 0;
	for(int i=0; i<LenPossibleDeltas; i++)
	{
		double CurNormDelta = PossibleNormDeltas[i];
		double dOrderShifted = log10(StrictDelta/CurNormDelta) + 1000.;

		long iOrderShifted = (long)dOrderShifted;
		double FractOrder = dOrderShifted - iOrderShifted;
		if(FractOrder > 0.5) 
		{
			FractOrder = 1. - FractOrder;
			iOrderShifted++;
		}
		long iOrder = iOrderShifted - 1000;

		if(MinFractOrder > FractOrder)
		{
			MinFractOrder = FractOrder;
			iResOrder = iOrder;
			IndNormDelta = i;
		}
	}
	
	return PossibleNormDeltas[IndNormDelta]*pow(10., (double)iResOrder);
}

//-------------------------------------------------------------------------

double radTg3dGraphPresent::FindGridTickLength(TVector3d* GridLimits)
{
	const double TickLimRatio = 0.013;
	double Dims[] = {GridLimits[1].x - GridLimits[0].x, GridLimits[1].y - GridLimits[0].y, GridLimits[1].z - GridLimits[0].z};
	double MaxLength = 0;
	for(int i=0; i<3; i++)
	{
		double CurLen = TickLimRatio*Dims[i];
		if(MaxLength < CurLen) MaxLength = CurLen;
	}
	return MaxLength;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRecMagGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTRecMag* RecMagP = (radTRecMag*)g3dPtr;
	TVector3d CornPoi(RecMagP->CentrPoint.x - 0.5 * RecMagP->Dimensions.x,
					  RecMagP->CentrPoint.y - 0.5 * RecMagP->Dimensions.y,
					  RecMagP->CentrPoint.z - 0.5 * RecMagP->Dimensions.z);
	TVector3d CornVect_x(RecMagP->Dimensions.x, 0., 0.);
	TVector3d CornVect_y(0., RecMagP->Dimensions.y, 0.);
	TVector3d CornVect_z(0., 0., RecMagP->Dimensions.z);

	CornPoi = BaseTransPtr->TrPoint(CornPoi);
	CornVect_x = BaseTransPtr->TrBiPoint(CornVect_x);
	CornVect_y = BaseTransPtr->TrBiPoint(CornVect_y);
	CornVect_z = BaseTransPtr->TrBiPoint(CornVect_z);

	TVector3d CP_pl_CVx = CornPoi + CornVect_x;
	TVector3d CP_pl_CVy = CornPoi + CornVect_y;
	TVector3d CP_pl_CVz = CornPoi + CornVect_z;
	TVector3d CP_pl_CVx_pl_CVy = CP_pl_CVx + CornVect_y;
	TVector3d CP_pl_CVx_pl_CVz = CP_pl_CVx + CornVect_z;
	TVector3d CP_pl_CVy_pl_CVz = CP_pl_CVy + CornVect_z;
	TVector3d CP_pl_CVx_pl_CVy_pl_CVz = CP_pl_CVx_pl_CVy + CornVect_z;

	short AmOfNonInternalFaces = 6;
	//short InternalFaces[6];
	short InternalFaces[] = {0,0,0,0,0,0};
	if(!ShowInternalFacesAfterCut)
	{
		AmOfNonInternalFaces = 0;
		RecMagP->ListFacesInternalAfterCut(InternalFaces);
		for(int k=0; k<6; k++) if(!InternalFaces[k]) AmOfNonInternalFaces++;
	}

	//char SendExtraContourLines = (ShowEdgeLinesInQD3D && (DrawFacilityInd == 1));
	char SendExtraContourLines = (ShowEdgeLines && ((DrawFacilityInd == 1) || (DrawFacilityInd == 2)));

	if(AmOfNonInternalFaces > 0)
	{
		Send.InitOutList(AmOfNonInternalFaces, DrawFacilityInd);
		TVector3d Side[5];

		int LocParity = BaseTransPtr->ShowParity();

		if(!InternalFaces[4] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVy; Side[2]=CP_pl_CVx_pl_CVy; Side[3]=CP_pl_CVx;
			}
			else
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVx; Side[2]=CP_pl_CVx_pl_CVy; Side[3]=CP_pl_CVy;
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		if(!InternalFaces[2] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVx; Side[2]=CP_pl_CVx_pl_CVz; Side[3]=CP_pl_CVz; 
			}
			else
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVz; Side[2]=CP_pl_CVx_pl_CVz; Side[3]=CP_pl_CVx; 
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		if(!InternalFaces[0] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVz; Side[2]=CP_pl_CVy_pl_CVz; Side[3]=CP_pl_CVy; 
			}
			else
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVy; Side[2]=CP_pl_CVy_pl_CVz; Side[3]=CP_pl_CVz; 
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		if(!InternalFaces[5] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CP_pl_CVz; Side[1]=CP_pl_CVx_pl_CVz; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVy_pl_CVz;	
			}
			else
			{
				Side[0]=CP_pl_CVz; Side[1]=CP_pl_CVy_pl_CVz; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVx_pl_CVz;	
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		if(!InternalFaces[3] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CP_pl_CVy; Side[1]=CP_pl_CVy_pl_CVz; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVx_pl_CVy;
			}
			else
			{
				Side[0]=CP_pl_CVy; Side[1]=CP_pl_CVx_pl_CVy; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVy_pl_CVz;	
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		if(!InternalFaces[1] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CP_pl_CVx; Side[1]=CP_pl_CVx_pl_CVy; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVx_pl_CVz;	
			}
			else
			{
				Side[0]=CP_pl_CVx; Side[1]=CP_pl_CVx_pl_CVz; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVx_pl_CVy;	
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTExtrPolygonGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTExtrPolygon* ExtrPolygonP = (radTExtrPolygon*)g3dPtr;
	radTPolygon* BasePolygonPtr = (radTPolygon*)(ExtrPolygonP->BasePolygonHandle.rep);

	radTPolygonGraphPresent* PolygonGraphPresentPtr = (radTPolygonGraphPresent*)(BasePolygonPtr->CreateGraphPresent());
	PolygonGraphPresentPtr->DrawFacilityInd = DrawFacilityInd;

    PolygonGraphPresentPtr->ShowEdgeLines = ShowEdgeLines; //OC071002
	PolygonGraphPresentPtr->ShowFaces = ShowFaces; //OC071002
	
	radTrans TransForBasePgn, AuxTrans1, AuxTrans2, TransForNextBasePgn;
	TVector3d ZeroVect(0.,0.,0.), AxVectX(1.,0.,0.), AxVectY(0.,1.,0.), AxVectZ(0.,0.,1.);
	const double Pi = 3.141592653589793238;
	const double HalfPi = 0.5*Pi;
	TVector3d NextPoint = ExtrPolygonP->FirstPoint;
	if(ExtrPolygonP->AxOrnt==ParallelToX)
	{
		SetupRotation(ZeroVect, AxVectZ, HalfPi, AuxTrans1);
		SetupRotation(ZeroVect, AxVectY, HalfPi, AuxTrans2);
		TrProduct(&AuxTrans2, &AuxTrans1, TransForBasePgn);
		
		NextPoint.x += ExtrPolygonP->Thickness;
	}
	else if(ExtrPolygonP->AxOrnt==ParallelToY)
	{
		SetupRotation(ZeroVect, AxVectZ, Pi, AuxTrans1);
		SetupRotation(ZeroVect, AxVectX, -HalfPi, AuxTrans2);
		TrProduct(&AuxTrans2, &AuxTrans1, TransForBasePgn);

		NextPoint.y += ExtrPolygonP->Thickness;
	}
	else
	{
		radIdentTrans IdentTrans;
		TransForBasePgn = IdentTrans;

		NextPoint.z += ExtrPolygonP->Thickness;
	}

	Send.InitOutList(2 + BasePolygonPtr->AmOfEdgePoints, DrawFacilityInd);

	TVector2d AuxVect2d = BasePolygonPtr->EdgePointsVector[0];
	TVector3d FirstEdgePoint(AuxVect2d.x, AuxVect2d.y, 0.);

	SetupTranslation(TransForBasePgn.TrPoint(FirstEdgePoint), ExtrPolygonP->FirstPoint, AuxTrans1);
	TrProduct(&AuxTrans1, &TransForBasePgn, AuxTrans2);

	//OC 060902
	//adding mirror symmetry to TransForBasePgn, to ensure correct drawing
	radTrans AuxPlSym, AuxTransLoc, AuxTransForBasePgn;
	TVector3d* pFirstPoint = &(ExtrPolygonP->FirstPoint);
	TVector3d AuxN(NextPoint.x - pFirstPoint->x, NextPoint.y - pFirstPoint->y, NextPoint.z - pFirstPoint->z);
	SetupPlaneSym(ExtrPolygonP->FirstPoint, AuxN, AuxPlSym);
	TrProduct(&AuxPlSym, &AuxTrans2, AuxTransLoc);
	TrProduct(BaseTransPtr, &AuxTransLoc, AuxTransForBasePgn);
	//END OC

	TrProduct(BaseTransPtr, &AuxTrans2, TransForBasePgn);

	SetupTranslation(ExtrPolygonP->FirstPoint, NextPoint, AuxTrans1);
	TrProduct(&AuxTrans1, &AuxTrans2, TransForNextBasePgn);
	TrProduct(BaseTransPtr, &TransForNextBasePgn, TransForNextBasePgn);

	//PolygonGraphPresentPtr->Draw(&TransForBasePgn);
	PolygonGraphPresentPtr->Draw(&AuxTransForBasePgn); //OC 060902
	PolygonGraphPresentPtr->Draw(&TransForNextBasePgn);
	
	TVector3d Side[5], AuxSide[5];
	TVector3d ThcnVect = BaseTransPtr->TrBiPoint(NextPoint - ExtrPolygonP->FirstPoint);
	TVector3d AuxVect(0.,0.,0.);

	int LocParity = BaseTransPtr->ShowParity();
	//char SendExtraContourLines = (ShowEdgeLinesInQD3D && (DrawFacilityInd == 1));
	char SendExtraContourLines = (ShowEdgeLines && ((DrawFacilityInd == 1) || (DrawFacilityInd == 2)));

	FirstEdgePoint = TransForBasePgn.TrPoint(FirstEdgePoint);
	Side[0] = FirstEdgePoint;
	Side[3] = Side[0] + ThcnVect;
	for(int k=1; k<BasePolygonPtr->AmOfEdgePoints; k++)
	{
		AuxVect2d = BasePolygonPtr->EdgePointsVector[k];
		AuxVect.x = AuxVect2d.x; AuxVect.y = AuxVect2d.y;
		Side[1] = TransForBasePgn.TrPoint(AuxVect);
		Side[2] = Side[1] + ThcnVect;

		if(LocParity > 0)
		{
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		else
		{
			*AuxSide = *Side; AuxSide[1] = Side[3]; AuxSide[2] = Side[2]; AuxSide[3] = Side[1]; 
			Send.Polygon(AuxSide, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				AuxSide[4] = *AuxSide; Send.Line(AuxSide, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}

		Side[0] = Side[1]; Side[3] = Side[2];
	}
	Side[1] = FirstEdgePoint; Side[2] = Side[1] + ThcnVect;

	if(LocParity > 0)
	{
		Send.Polygon(Side, 4, DrawFacilityInd);
		if(SendExtraContourLines) 
		{ 
			SetCurrentColorInStack(SbdLineColor);
			Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
			RemoveCurrentColorFromStack();
		}
	}
	else
	{
		*AuxSide = *Side; AuxSide[1] = Side[3]; AuxSide[2] = Side[2]; AuxSide[3] = Side[1]; 
		Send.Polygon(AuxSide, 4, DrawFacilityInd);
		if(SendExtraContourLines) 
		{ 
			SetCurrentColorInStack(SbdLineColor);
			AuxSide[4] = *AuxSide; Send.Line(AuxSide, 5, DrawFacilityInd);
			RemoveCurrentColorFromStack();
		}
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTPolyhedronGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTPolyhedron* VolLimByPgnsP = (radTPolyhedron*)g3dPtr;
	radTVectHandlePgnAndTrans& VectHandlePgnAndTrans = VolLimByPgnsP->VectHandlePgnAndTrans;

	//TVector3d& CenPointVolLimByPgns = VolLimByPgnsP->CentrPoint; //OC090908

	int TotAmOfFaces = VolLimByPgnsP->AmOfFaces;
	int AmOfFacesToDraw = TotAmOfFaces;

	if(!ShowInternalFacesAfterCut)
	{
		CountNonInternalFaces();
		AmOfFacesToDraw = AmOfNonInternalFaces;
	}

			//test
			//char ErrorMesTitle[] = "SRW Debug";
			//char ErrorStr[100];
			//int j = sprintf(ErrorStr, "Polyhedron: AmOfFacesToDraw: %d", AmOfFacesToDraw);
			//UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
			//int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle); 
			//end test

	if(AmOfFacesToDraw > 0)
	{
		Send.InitOutList(AmOfFacesToDraw, DrawFacilityInd);

		for(int i=0; i<TotAmOfFaces; i++)
		{
			radTHandlePgnAndTrans HandlePgnAndTrans = VectHandlePgnAndTrans[i];
			if(ShowInternalFacesAfterCut || (!HandlePgnAndTrans.FaceIsInternalAfterCut))
			{
				radTrans* NativeFaceTrans = HandlePgnAndTrans.TransHndl.rep;
				//radTrans NativeFaceTrans(*(HandlePgnAndTrans.TransHndl.rep)); //OC090908
				//NativeFaceTrans.SetVector(CenPointVolLimByPgns); //OC090908

				radTrans TotFaceTrans;
				TrProduct(BaseTransPtr, NativeFaceTrans, TotFaceTrans);
				//TrProduct(BaseTransPtr, &NativeFaceTrans, TotFaceTrans); //OC090908

				radTg3dGraphPresent* g3dGraphPresentPtr = HandlePgnAndTrans.PgnHndl.rep->CreateGraphPresent();
				g3dGraphPresentPtr->DrawFacilityInd = DrawFacilityInd;

				g3dGraphPresentPtr->ShowEdgeLines = ShowEdgeLines; //OC071002
				g3dGraphPresentPtr->ShowFaces = ShowFaces; //OC071002

				g3dGraphPresentPtr->GenTrans = TotFaceTrans;
				g3dGraphPresentPtr->DrawAttrAreSet = DrawAttrAreSet;
				g3dGraphPresentPtr->DrawAttr = DrawAttr;
				g3dGraphPresentPtr->GenDraw();
				delete g3dGraphPresentPtr;
			}
		}
	}
}

//-------------------------------------------------------------------------

void radTPolyhedronGraphPresent::CountNonInternalFaces()
{
	radTPolyhedron* PolyhedronP = (radTPolyhedron*)g3dPtr;
	radTVectHandlePgnAndTrans& VectHandlePgnAndTrans = PolyhedronP->VectHandlePgnAndTrans;
	int AmOfFacesToDraw = 0;
	for(int i=0; i<PolyhedronP->AmOfFaces; i++)
		if(!VectHandlePgnAndTrans[i].FaceIsInternalAfterCut) AmOfFacesToDraw++;
	AmOfNonInternalFaces = AmOfFacesToDraw;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTrans* radTArcCurGraphPresent::CreateRotation(TVector3d& PoiOnAxVect, TVector3d& InAxVect, double Angle)
{
	// Later apply try-catch here
	double NormFact = 1./sqrt(InAxVect.x*InAxVect.x+InAxVect.y*InAxVect.y+InAxVect.z*InAxVect.z);
	TVector3d AxVect = NormFact*InAxVect;
	double VxVx, VyVy, VzVz;
	VxVx=AxVect.x*AxVect.x; VyVy=AxVect.y*AxVect.y; VzVz=AxVect.z*AxVect.z;

	double cosAng, sinAng, One_m_cos;
	cosAng = cos(Angle); sinAng = sin(Angle); One_m_cos = 1. - cosAng;
	double One_m_cosVxVy, One_m_cosVxVz, One_m_cosVyVz, sinVx, sinVy, sinVz;
	One_m_cosVxVy = One_m_cos*AxVect.x*AxVect.y;
	One_m_cosVxVz = One_m_cos*AxVect.x*AxVect.z;
	One_m_cosVyVz = One_m_cos*AxVect.y*AxVect.z;
	sinVx = sinAng*AxVect.x; sinVy = sinAng*AxVect.y; sinVz = sinAng*AxVect.z;

	TVector3d St0(VxVx+cosAng*(VyVy+VzVz), One_m_cosVxVy-sinVz, One_m_cosVxVz+sinVy);
	TVector3d St1(One_m_cosVxVy+sinVz, VyVy+cosAng*(VxVx+VzVz), One_m_cosVyVz-sinVx);
	TVector3d St2(One_m_cosVxVz-sinVy, One_m_cosVyVz+sinVx, VzVz+cosAng*(VxVx+VyVy));
	TMatrix3d M(St0, St1, St2);
	TVector3d St00(1.-St0.x, -St0.y, -St0.z);
	TVector3d St01(-St1.x, 1.-St1.y, -St1.z);
	TVector3d St02(-St2.x, -St2.y, 1.-St2.z);
	TMatrix3d M0(St00, St01, St02);

	return new radTrans(M, M0*PoiOnAxVect, 1., 1.);
}

//-------------------------------------------------------------------------

void radTArcCurGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTArcCur* ArcCurP = (radTArcCur*)g3dPtr;

	short AmOfNonInternalFaces = 6;
	//short InternalFaces[6];
	short InternalFaces[] = {0,0,0,0,0,0};
	if(!ShowInternalFacesAfterCut)
	{
		AmOfNonInternalFaces = 0;
		ArcCurP->ListFacesInternalAfterCut(InternalFaces);
		for(int k=0; k<6; k++) if(!InternalFaces[k]) AmOfNonInternalFaces++;
	}
	if(AmOfNonInternalFaces == 0) return;

	int LocParity = 1;

	char Face0 = (!InternalFaces[0]) || ShowInternalFacesAfterCut;
	char Face1 = (!InternalFaces[1]) || ShowInternalFacesAfterCut;
	char Face2 = (!InternalFaces[2]) || ShowInternalFacesAfterCut;
	char Face3 = (!InternalFaces[3]) || ShowInternalFacesAfterCut;
	char Face4 = (!InternalFaces[4]) || ShowInternalFacesAfterCut;
	char Face5 = (!InternalFaces[5]) || ShowInternalFacesAfterCut;

	if(!OutlineSegments) { Send.DrawEdgesSuppression(DrawFacilityInd); Send.InitOutList(2, DrawFacilityInd);}

	TVector3d CenPoi = ArcCurP->CircleCentrPoint;
	TVector3d AxVect(0., 0., ArcCurP->Height);

	radTrans* AuxRotPtr = CreateRotation(CenPoi, AxVect, ArcCurP->Phi_min);
	TVector3d SmallRadVect(ArcCurP->R_min, 0., 0.);
	SmallRadVect = AuxRotPtr->TrBiPoint(SmallRadVect);
	delete AuxRotPtr;

	CenPoi = BaseTransPtr->TrPoint(CenPoi);
	AxVect = BaseTransPtr->TrAxialVect(AxVect);
	SmallRadVect = BaseTransPtr->TrBiPoint(SmallRadVect);

	TVector3d SideBottom[4], SideTop[4], SideInt[4], SideExt[4], SideFrontEnd[4], AuxPgn[4];

	double RadRatio = ArcCurP->R_max/ArcCurP->R_min;

	AuxRotPtr = CreateRotation(CenPoi, AxVect, (ArcCurP->Phi_max - ArcCurP->Phi_min)/ArcCurP->NumberOfSectors);
	TVector3d SmallAuxTranslVect = AuxRotPtr->TrBiPoint(SmallRadVect) - SmallRadVect;

	SideBottom[0] = CenPoi + (-0.5)*AxVect + SmallRadVect;
	SideBottom[3] = CenPoi + (-0.5)*AxVect + RadRatio*SmallRadVect;
	SideTop[0] = SideBottom[0] + AxVect;
	SideTop[3] = SideBottom[3] + AxVect;
	SideInt[0] = SideBottom[0];
	SideInt[3] = SideTop[0];
	SideExt[0] = SideBottom[3];
	SideExt[3] = SideTop[3];

	int AmOfNonInternalLongFaces = (Face2? 1 : 0) + (Face3? 1 : 0) + (Face4? 1 : 0) + (Face5? 1 : 0);
	int AmOfNonInternalFrontEndFaces = (Face0? 1 : 0) + (Face1? 1 : 0);
	Send.InitOutList(AmOfNonInternalLongFaces*ArcCurP->NumberOfSectors + AmOfNonInternalFrontEndFaces, DrawFacilityInd);

	SideFrontEnd[0] = SideBottom[0];
	SideFrontEnd[1] = SideTop[0];
	SideFrontEnd[2] = SideTop[3];
	SideFrontEnd[3] = SideBottom[3];

	if(Face0)
	{
		if(LocParity > 0)
		{
			AuxPgn[0] = SideFrontEnd[0]; AuxPgn[1] = SideFrontEnd[3]; AuxPgn[2] = SideFrontEnd[2]; AuxPgn[3] = SideFrontEnd[1]; 
			Send.Polygon(AuxPgn, 4, DrawFacilityInd);
		}
		else
		{
			Send.Polygon(SideFrontEnd, 4, DrawFacilityInd);
		}
	}

	int lenEdgeLine = ArcCurP->NumberOfSectors + 1;
	TVector3d *LowerCloseEdgeLine=0, *LowerFarEdgeLine=0, *UpperCloseEdgeLine=0, *UpperFarEdgeLine=0;
	TVector3d FrontLine[] = { SideFrontEnd[0], SideFrontEnd[1], SideFrontEnd[2], SideFrontEnd[3], SideFrontEnd[0] };
	if(!OutlineSegments)
	{
		if(AmOfNonInternalLongFaces > 0)
		{
			LowerCloseEdgeLine = new TVector3d[lenEdgeLine];
			LowerFarEdgeLine = new TVector3d[lenEdgeLine];
			UpperCloseEdgeLine = new TVector3d[lenEdgeLine];
			UpperFarEdgeLine = new TVector3d[lenEdgeLine];

			LowerCloseEdgeLine[0] = SideBottom[0];
			LowerFarEdgeLine[0] = SideBottom[3];
			UpperCloseEdgeLine[0] = SideTop[0];
			UpperFarEdgeLine[0] = SideTop[3];
		}
	}

	for(int i = 1; i <= ArcCurP->NumberOfSectors; i++)
	{
		SideBottom[1] = SideBottom[0] + SmallAuxTranslVect;
		SideBottom[2] = SideBottom[3] + RadRatio*SmallAuxTranslVect;
		SideTop[1] = SideBottom[1] + AxVect;
		SideTop[2] = SideBottom[2] + AxVect;
		SideInt[1] = SideBottom[1];
		SideInt[2] = SideTop[1];
		SideExt[1] = SideBottom[2];
		SideExt[2] = SideTop[2];

		if(!OutlineSegments)
		{
			if(Face4 || Face2) LowerCloseEdgeLine[i] = SideBottom[1];
			if(Face4 || Face3) LowerFarEdgeLine[i] = SideBottom[2];
			if(Face5 || Face2) UpperCloseEdgeLine[i] = SideTop[1];
			if(Face5 || Face3) UpperFarEdgeLine[i] = SideTop[2];
		}

		if(Face4) 
		{
			if(LocParity > 0)
			{
				Send.Polygon(SideBottom, 4, DrawFacilityInd);
			}
			else
			{
				AuxPgn[0] = SideBottom[0]; AuxPgn[1] = SideBottom[3]; AuxPgn[2] = SideBottom[2]; AuxPgn[3] = SideBottom[1]; 
				Send.Polygon(AuxPgn, 4, DrawFacilityInd);
			}
		}
		if(Face5) 
		{
			if(LocParity > 0)
			{
				AuxPgn[0] = SideTop[0]; AuxPgn[1] = SideTop[3]; AuxPgn[2] = SideTop[2]; AuxPgn[3] = SideTop[1]; 
				Send.Polygon(AuxPgn, 4, DrawFacilityInd);
			}
			else
			{
				Send.Polygon(SideTop, 4, DrawFacilityInd);
			}
		}
		if(Face2) 
		{
			if(LocParity > 0)
			{
				AuxPgn[0] = SideInt[0]; AuxPgn[1] = SideInt[3]; AuxPgn[2] = SideInt[2]; AuxPgn[3] = SideInt[1]; 
				Send.Polygon(AuxPgn, 4, DrawFacilityInd);
			}
			else
			{
				Send.Polygon(SideInt, 4, DrawFacilityInd);
			}
		}
		if(Face3) 
		{
			if(LocParity > 0)
			{
				Send.Polygon(SideExt, 4, DrawFacilityInd);
			}
			else
			{
				AuxPgn[0] = SideExt[0]; AuxPgn[1] = SideExt[3]; AuxPgn[2] = SideExt[2]; AuxPgn[3] = SideExt[1]; 
				Send.Polygon(AuxPgn, 4, DrawFacilityInd);
			}
		}

		SideBottom[0] = SideBottom[1];
		SideBottom[3] = SideBottom[2];
		SideTop[0] = SideTop[1];
		SideTop[3] = SideTop[2];
		SideInt[0] = SideInt[1];
		SideInt[3] = SideInt[2];
		SideExt[0] = SideExt[1];
		SideExt[3] = SideExt[2];
		
		SmallAuxTranslVect = AuxRotPtr->TrBiPoint(SmallAuxTranslVect);
	}

	SideFrontEnd[0] = SideBottom[0];
	SideFrontEnd[1] = SideTop[0];
	SideFrontEnd[2] = SideTop[3];
	SideFrontEnd[3] = SideBottom[3];
	if(Face1)
	{
		if(LocParity > 0)
		{
			Send.Polygon(SideFrontEnd, 4, DrawFacilityInd);
		}
		else
		{
			AuxPgn[0] = SideFrontEnd[0]; AuxPgn[1] = SideFrontEnd[3]; AuxPgn[2] = SideFrontEnd[2]; AuxPgn[3] = SideFrontEnd[1]; 
			Send.Polygon(AuxPgn, 4, DrawFacilityInd);
		}
	}

	delete AuxRotPtr;

	if(!OutlineSegments)
	{
		TVector3d EndLine[] = { SideFrontEnd[0], SideFrontEnd[1], SideFrontEnd[2], SideFrontEnd[3], SideFrontEnd[0] };

		Send.InitDrawLineWithThickness(DrawAttrAreSet, DrawAttr, DrawFacilityInd);
		char BufChar = Face2 || Face3 || Face4 || Face5;
		char Face0Line = Face0 || BufChar;
		char Face1Line = Face1 || BufChar;

		int AmOfLinesToSend = (Face0Line? 1 : 0) + (Face1Line? 1 : 0) + ((Face2 || Face4)? 1 : 0) 
			+ ((Face3 || Face4)? 1 : 0) + ((Face2 || Face5)? 1 : 0) + ((Face3 || Face5)? 1 : 0);
		Send.InitOutList(AmOfLinesToSend, DrawFacilityInd);

		if(Face0Line) Send.Line(FrontLine, 5, DrawFacilityInd);
		if(Face4 || Face2) Send.Line(LowerCloseEdgeLine, lenEdgeLine, DrawFacilityInd);
		if(Face4 || Face3) Send.Line(LowerFarEdgeLine, lenEdgeLine, DrawFacilityInd);
		if(Face5 || Face2) Send.Line(UpperCloseEdgeLine, lenEdgeLine, DrawFacilityInd);
		if(Face5 || Face3) Send.Line(UpperFarEdgeLine, lenEdgeLine, DrawFacilityInd);
		if(Face1Line) Send.Line(EndLine, 5, DrawFacilityInd);

		if(LowerCloseEdgeLine != 0) delete[] LowerCloseEdgeLine;
		if(LowerFarEdgeLine != 0) delete[] LowerFarEdgeLine;
		if(UpperCloseEdgeLine != 0) delete[] UpperCloseEdgeLine;
		if(UpperFarEdgeLine != 0) delete[] UpperFarEdgeLine;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTGroupGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTGroup* GroupPtr = (radTGroup*)g3dPtr;

	int AmOfElemToDraw = (int)(GroupPtr->GroupMapOfHandlers.size());
	int *ItemIsNotFullyInternal, *tItemIsNotFullyInternal;
	ItemIsNotFullyInternal = tItemIsNotFullyInternal = 0;

	ItemIsNotFullyInternal = new int[AmOfElemToDraw]; //OC071002
	for(int k=0; k<AmOfElemToDraw; k++) ItemIsNotFullyInternal[k] = 1; //OC071002

	if(!ShowInternalFacesAfterCut)
	{
		//ItemIsNotFullyInternal = new int[AmOfElemToDraw]; //OC071002
		tItemIsNotFullyInternal = ItemIsNotFullyInternal;
		AmOfElemToDraw = 0;
		radTmhg& GrMapOfHndl = GroupPtr->GroupMapOfHandlers;
		radTmhg::const_iterator iter;
		for(iter = GrMapOfHndl.begin(); iter != GrMapOfHndl.end(); ++iter)
		{
			radTg3d* Loc_g3dPtr = (radTg3d*)(((*iter).second).rep);
			int CurrentItemIsNotInternal = Loc_g3dPtr->ItemIsNotFullyInternalAfterCut();
			*(tItemIsNotFullyInternal++) = CurrentItemIsNotInternal;
			AmOfElemToDraw += CurrentItemIsNotInternal;
		}
	}

		//test
		//char ErrorMesTitle[] = "SRW Debug";
		//char ErrorStr[100];
		//int j = sprintf(ErrorStr, "AmOfElemToDraw: %d", AmOfElemToDraw);
		//UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
		//int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle); 
		//end test

	if(AmOfElemToDraw > 0)
	{
		Send.InitOutList(AmOfElemToDraw, DrawFacilityInd);
		tItemIsNotFullyInternal = ItemIsNotFullyInternal;

		for(radTmhg::const_iterator iter = GroupPtr->GroupMapOfHandlers.begin();
			iter != GroupPtr->GroupMapOfHandlers.end(); ++iter)
		{
			if(ShowInternalFacesAfterCut || (*(tItemIsNotFullyInternal++)))
			{
				radTg3d* g3dPtrLoc = (radTg3d*)(((*iter).second).rep);
				radTg3dGraphPresent* g3dGraphPresentPtr = g3dPtrLoc->CreateGraphPresent();
				if(g3dGraphPresentPtr != NULL)
				{
					if(!ShowInternalFacesAfterCut) g3dGraphPresentPtr->ShowInternalFacesAfterCut = false;

					g3dGraphPresentPtr->DrawFacilityInd = DrawFacilityInd;

					g3dGraphPresentPtr->ShowEdgeLines = ShowEdgeLines;
					g3dGraphPresentPtr->ShowFaces = ShowFaces;

					g3dGraphPresentPtr->GraphPresOptions = GraphPresOptions;
					g3dGraphPresentPtr->MapOfDrawAttrPtr = MapOfDrawAttrPtr;

					if(!DrawAttrAreSet) g3dGraphPresentPtr->RetrieveDrawAttr((*iter).first);
					else
					{
						g3dGraphPresentPtr->DrawAttrAreSet = DrawAttrAreSet;
						g3dGraphPresentPtr->DrawAttr = DrawAttr;
					}
					g3dGraphPresentPtr->GenTrans = *BaseTransPtr;
					g3dGraphPresentPtr->GenDraw();

					delete g3dGraphPresentPtr;
				}
			}
		}
	}
	if(ItemIsNotFullyInternal != 0) delete[] ItemIsNotFullyInternal;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTSubdivRecMagGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTSubdividedRecMag* SubdividedRecMagP = (radTSubdividedRecMag*)((radTGroup*)g3dPtr);
	//TVector3d CornPoi(SubdividedRecMagP->CentrPoint.x - 0.5 * SubdividedRecMagP->Dimensions.x,
	//				  SubdividedRecMagP->CentrPoint.y - 0.5 * SubdividedRecMagP->Dimensions.y,
	//				  SubdividedRecMagP->CentrPoint.z - 0.5 * SubdividedRecMagP->Dimensions.z);

	TVector3d &LocCentrPoint = ((radTRecMag*)g3dPtr)->CentrPoint; //OC061008
	TVector3d CornPoi(LocCentrPoint.x - 0.5 * SubdividedRecMagP->Dimensions.x,
					  LocCentrPoint.y - 0.5 * SubdividedRecMagP->Dimensions.y,
					  LocCentrPoint.z - 0.5 * SubdividedRecMagP->Dimensions.z);

	TVector3d CornVect_x(SubdividedRecMagP->Dimensions.x, 0., 0.);
	TVector3d CornVect_y(0., SubdividedRecMagP->Dimensions.y, 0.);
	TVector3d CornVect_z(0., 0., SubdividedRecMagP->Dimensions.z);

	CornPoi = BaseTransPtr->TrPoint(CornPoi);
	CornVect_x = BaseTransPtr->TrBiPoint(CornVect_x);
	CornVect_y = BaseTransPtr->TrBiPoint(CornVect_y);
	CornVect_z = BaseTransPtr->TrBiPoint(CornVect_z);

	TVector3d CP_pl_CVx = CornPoi + CornVect_x;
	TVector3d CP_pl_CVy = CornPoi + CornVect_y;
	TVector3d CP_pl_CVz = CornPoi + CornVect_z;
	TVector3d CP_pl_CVx_pl_CVy = CP_pl_CVx + CornVect_y;
	TVector3d CP_pl_CVx_pl_CVz = CP_pl_CVx + CornVect_z;
	TVector3d CP_pl_CVy_pl_CVz = CP_pl_CVy + CornVect_z;
	TVector3d CP_pl_CVx_pl_CVy_pl_CVz = CP_pl_CVx_pl_CVy + CornVect_z;

	short AmOfNonInternalFaces = 6;
	short InternalFaces[6];
	if(!ShowInternalFacesAfterCut)
	{
		AmOfNonInternalFaces = 0;
		SubdividedRecMagP->ListFacesInternalAfterCut(InternalFaces);
		for(int k=0; k<6; k++) if(!InternalFaces[k]) AmOfNonInternalFaces++;
	}

	if(AmOfNonInternalFaces > 0)
	{
		Send.InitOutList(AmOfNonInternalFaces + 1, DrawFacilityInd);
		TVector3d Side[5];

		//char SendExtraContourLines = (ShowEdgeLinesInQD3D && (DrawFacilityInd == 1));
		char SendExtraContourLines = (ShowEdgeLines && ((DrawFacilityInd == 1) || (DrawFacilityInd == 2)));

		int LocParity = BaseTransPtr->ShowParity();
		if(!InternalFaces[4] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVy; Side[2]=CP_pl_CVx_pl_CVy; Side[3]=CP_pl_CVx;
			}
			else
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVx; Side[2]=CP_pl_CVx_pl_CVy; Side[3]=CP_pl_CVy;
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		if(!InternalFaces[2] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVx; Side[2]=CP_pl_CVx_pl_CVz; Side[3]=CP_pl_CVz; 
			}
			else
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVz; Side[2]=CP_pl_CVx_pl_CVz; Side[3]=CP_pl_CVx; 
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		if(!InternalFaces[0] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVz; Side[2]=CP_pl_CVy_pl_CVz; Side[3]=CP_pl_CVy; 
			}
			else
			{
				Side[0]=CornPoi; Side[1]=CP_pl_CVy; Side[2]=CP_pl_CVy_pl_CVz; Side[3]=CP_pl_CVz; 
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		if(!InternalFaces[5] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CP_pl_CVz; Side[1]=CP_pl_CVx_pl_CVz; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVy_pl_CVz;	
			}
			else
			{
				Side[0]=CP_pl_CVz; Side[1]=CP_pl_CVy_pl_CVz; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVx_pl_CVz;	
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		if(!InternalFaces[3] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CP_pl_CVy; Side[1]=CP_pl_CVy_pl_CVz; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVx_pl_CVy;
			}
			else
			{
				Side[0]=CP_pl_CVy; Side[1]=CP_pl_CVx_pl_CVy; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVy_pl_CVz;	
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		if(!InternalFaces[1] || ShowInternalFacesAfterCut)
		{
			if(LocParity > 0)
			{
				Side[0]=CP_pl_CVx; Side[1]=CP_pl_CVx_pl_CVy; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVx_pl_CVz;	
			}
			else
			{
				Side[0]=CP_pl_CVx; Side[1]=CP_pl_CVx_pl_CVz; Side[2]=CP_pl_CVx_pl_CVy_pl_CVz; Side[3]=CP_pl_CVx_pl_CVy;	
			}
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}

		TVector3d BasicVect3dArray[] = {CornPoi, CornVect_x, CornVect_y, CornVect_z};
		DrawSubdivisionLines(SubdividedRecMagP, BasicVect3dArray);
	}
}

//-------------------------------------------------------------------------

void radTSubdivRecMagGraphPresent::DrawSubdivisionLines(radTSubdividedRecMag* SubdividedRecMagP, TVector3d* BasicVect3dArray)
{
	SetCurrentColorInStack(SbdLineColor);

	//double SmallDist = 1.E-05;

	TVector3d CornPoi = BasicVect3dArray[0];
	TVector3d CornVect_x = BasicVect3dArray[1];
	TVector3d CornVect_y = BasicVect3dArray[2];
	TVector3d CornVect_z = BasicVect3dArray[3];

// To remove?
	//TVector3d UnitDiagVect = CornVect_x + CornVect_y + CornVect_z;
	//UnitDiagVect = (1./sqrt(UnitDiagVect.x*UnitDiagVect.x + UnitDiagVect.y*UnitDiagVect.y + UnitDiagVect.z*UnitDiagVect.z))*UnitDiagVect;
	//CornPoi = CornPoi - SmallDist*UnitDiagVect;
	//CornVect_x = (1. + 2.*SmallDist/SubdividedRecMagP->Dimensions.x)*CornVect_x;
	//CornVect_y = (1. + 2.*SmallDist/SubdividedRecMagP->Dimensions.y)*CornVect_y;
	//CornVect_z = (1. + 2.*SmallDist/SubdividedRecMagP->Dimensions.z)*CornVect_z;
// End To remove?

	radTDrawAttr LocDrawAttr;
	LocDrawAttr.RGB_col.Red = LocDrawAttr.RGB_col.Green = LocDrawAttr.RGB_col.Blue = 0.;
	LocDrawAttr.LineThickness = DrawAttr.LineThickness;

	int AmOfLocSubdRecMags = 0;

	//radTCast Cast;
	for(radTmhg::const_iterator Iter = SubdividedRecMagP->GroupMapOfHandlers.begin();
		Iter != SubdividedRecMagP->GroupMapOfHandlers.end(); ++Iter)
	{
		radTg3d* g3dPtr = (radTg3d*)((*Iter).second.rep);
		radTGroup* GroupPtr = radTCast::GroupCast(g3dPtr); // Because Subdivided ExtrPolygons (and Subd. RecMags) are placed to the general container through the Group branch
		radTSubdividedRecMag* SubdividedRecMagPtr = (GroupPtr!=0)? radTCast::SubdividedRecMagCast(GroupPtr) : 0;
		if(SubdividedRecMagPtr != 0) AmOfLocSubdRecMags++;
	}

	int AmOfExtraItems = int(SubdividedRecMagP->kx) + int(SubdividedRecMagP->ky) + int(SubdividedRecMagP->kz) - 3 + AmOfLocSubdRecMags;
	Send.InitDrawLinElem(DrawAttrAreSet, LocDrawAttr, AmOfExtraItems, DrawFacilityInd);

	TVector3d LineLoop[5];
	TVector3d TranslVect;

	const double AbsZeroTol = 5.E-12;
	if(SubdividedRecMagP->kx > 1)
	{
		double q0x = (fabs(SubdividedRecMagP->kx-1.)>AbsZeroTol)? pow(SubdividedRecMagP->qx, 1./(SubdividedRecMagP->kx-1.)) : SubdividedRecMagP->qx;
		double BufX = SubdividedRecMagP->qx*q0x - 1.;
		double BufRatX = (fabs(BufX) > AbsZeroTol)? (q0x - 1.)/BufX : 1./SubdividedRecMagP->kx;

		TranslVect = BufRatX*CornVect_x;

		LineLoop[0] = CornPoi + TranslVect;
		LineLoop[1] = LineLoop[0] + CornVect_y;
		LineLoop[2] = LineLoop[1] + CornVect_z;
		LineLoop[3] = LineLoop[0] + CornVect_z;
		LineLoop[4] = LineLoop[0];
		Send.Line(LineLoop, 5, DrawFacilityInd);
		for(int ix=2; ix<SubdividedRecMagP->kx; ix++)
		{
			BufRatX *= q0x;
			TranslVect = BufRatX*CornVect_x;

			for(int i=0; i<4; i++) 	LineLoop[i] += TranslVect;
			LineLoop[4] = LineLoop[0];
			Send.Line(LineLoop, 5, DrawFacilityInd);
		}
	}
	if(SubdividedRecMagP->ky > 1)
	{
		double q0y = (fabs(SubdividedRecMagP->ky-1.)>AbsZeroTol)? pow(SubdividedRecMagP->qy, 1./(SubdividedRecMagP->ky-1.)) : SubdividedRecMagP->qy;
		double BufY = SubdividedRecMagP->qy*q0y - 1.;
		double BufRatY = (fabs(BufY) > AbsZeroTol)? (q0y - 1.)/BufY : 1./SubdividedRecMagP->ky;

		TranslVect = BufRatY*CornVect_y;

		LineLoop[0] = CornPoi + TranslVect;
		LineLoop[1] = LineLoop[0] + CornVect_x;
		LineLoop[2] = LineLoop[1] + CornVect_z;
		LineLoop[3] = LineLoop[0] + CornVect_z;
		LineLoop[4] = LineLoop[0];
		Send.Line(LineLoop, 5, DrawFacilityInd);
		for(int iy=2; iy<SubdividedRecMagP->ky; iy++)
		{
			BufRatY *= q0y;
			TranslVect = BufRatY*CornVect_y;

			for(int i=0; i<4; i++) 	LineLoop[i] += TranslVect;
			LineLoop[4] = LineLoop[0];
			Send.Line(LineLoop, 5, DrawFacilityInd);
		}
	}
	if(SubdividedRecMagP->kz > 1)
	{
		double q0z = (fabs(SubdividedRecMagP->kz-1.)>AbsZeroTol)? pow(SubdividedRecMagP->qz, 1./(SubdividedRecMagP->kz-1.)) : SubdividedRecMagP->qz;
		double BufZ = SubdividedRecMagP->qz*q0z - 1.;
		double BufRatZ = (fabs(BufZ) > AbsZeroTol)? (q0z - 1.)/BufZ : 1./SubdividedRecMagP->kz;

		TranslVect = BufRatZ*CornVect_z;

		LineLoop[0] = CornPoi + TranslVect;
		LineLoop[1] = LineLoop[0] + CornVect_x;
		LineLoop[2] = LineLoop[1] + CornVect_y;
		LineLoop[3] = LineLoop[0] + CornVect_y;
		LineLoop[4] = LineLoop[0];
		Send.Line(LineLoop, 5, DrawFacilityInd);
		for(int iz=2; iz<SubdividedRecMagP->kz; iz++)
		{
			BufRatZ *= q0z;
			TranslVect = BufRatZ*CornVect_z;

			for(int i=0; i<4; i++) 	LineLoop[i] += TranslVect;
			LineLoop[4] = LineLoop[0];
			Send.Line(LineLoop, 5, DrawFacilityInd);
		}
	}
	for(radTmhg::const_iterator iter = SubdividedRecMagP->GroupMapOfHandlers.begin();
		iter != SubdividedRecMagP->GroupMapOfHandlers.end(); ++iter)
	{
		radTg3d* g3dPtr = (radTg3d*)((*iter).second.rep);
		radTGroup* GroupPtr = radTCast::GroupCast(g3dPtr); // Because Subdivided ExtrPolygons (and Subd. RecMags) are placed to the general container through the Group branch
		radTSubdividedRecMag* SubdividedRecMagPtr = (GroupPtr!=0)? radTCast::SubdividedRecMagCast(GroupPtr) : 0;
		if(SubdividedRecMagPtr != 0)
		{
			TVector3d LocCornVect_x = (SubdividedRecMagPtr->Dimensions.x/SubdividedRecMagP->Dimensions.x)*CornVect_x;
			TVector3d LocCornVect_y = (SubdividedRecMagPtr->Dimensions.y/SubdividedRecMagP->Dimensions.y)*CornVect_y;
			TVector3d LocCornVect_z = (SubdividedRecMagPtr->Dimensions.z/SubdividedRecMagP->Dimensions.z)*CornVect_z;

			//TVector3d OrigShift = (SubdividedRecMagPtr->CentrPoint - 0.5*SubdividedRecMagPtr->Dimensions) 
			//					- (SubdividedRecMagP->CentrPoint - 0.5*SubdividedRecMagP->Dimensions);
			TVector3d &LocCentrPoint = ((radTRecMag*)g3dPtr)->CentrPoint; //OC061008
			TVector3d OrigShift = (LocCentrPoint - 0.5*SubdividedRecMagPtr->Dimensions) 
								- (LocCentrPoint - 0.5*SubdividedRecMagP->Dimensions);

			TVector3d Shift = (OrigShift.x/SubdividedRecMagP->Dimensions.x) * CornVect_x
							+ (OrigShift.y/SubdividedRecMagP->Dimensions.y) * CornVect_y
							+ (OrigShift.z/SubdividedRecMagP->Dimensions.z) * CornVect_z;
			TVector3d LocCornPoi = CornPoi + Shift;
			TVector3d LocBasicVect3dArray[] = {LocCornPoi, LocCornVect_x, LocCornVect_y, LocCornVect_z};

			DrawSubdivisionLines(SubdividedRecMagPtr, LocBasicVect3dArray);
		}
	}
	RemoveCurrentColorFromStack();
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTSubdivExtrPolygGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTSubdividedExtrPolygon* SubdExtrPolygonP = (radTSubdividedExtrPolygon*)((radTGroup*)g3dPtr);
	radTPolygon* BasePolygonPtr = (radTPolygon*)(SubdExtrPolygonP->BasePolygonHandle.rep);
	radTPolygonGraphPresent* PolygonGraphPresentPtr = (radTPolygonGraphPresent*)(BasePolygonPtr->CreateGraphPresent());
	PolygonGraphPresentPtr->DrawFacilityInd = DrawFacilityInd;
	
	radTrans TransForBasePgn, AuxTrans1, AuxTrans2, TransForNextBasePgn;
	TVector3d ZeroVect(0.,0.,0.), AxVectX(1.,0.,0.), AxVectY(0.,1.,0.), AxVectZ(0.,0.,1.);

	const double Pi = 3.141592653589793238;
	const double HalfPi = 0.5*Pi;

	TVector3d NextPoint = SubdExtrPolygonP->FirstPoint;
	if(SubdExtrPolygonP->AxOrnt==ParallelToX)
	{
		SetupRotation(ZeroVect, AxVectZ, HalfPi, AuxTrans1);
		SetupRotation(ZeroVect, AxVectY, HalfPi, AuxTrans2);
		TrProduct(&AuxTrans2, &AuxTrans1, TransForBasePgn);
		
		NextPoint.x += SubdExtrPolygonP->Thickness;
	}
	else if(SubdExtrPolygonP->AxOrnt==ParallelToY)
	{
		SetupRotation(ZeroVect, AxVectZ, Pi, AuxTrans1);
		SetupRotation(ZeroVect, AxVectX, -HalfPi, AuxTrans2);
		TrProduct(&AuxTrans2, &AuxTrans1, TransForBasePgn);

		NextPoint.y += SubdExtrPolygonP->Thickness;
	}
	else
	{
		radIdentTrans IdentTrans;
		TransForBasePgn = IdentTrans;

		NextPoint.z += SubdExtrPolygonP->Thickness;
	}

	Send.InitOutList(2 + BasePolygonPtr->AmOfEdgePoints + 1, DrawFacilityInd);

	TVector2d AuxVect2d = BasePolygonPtr->EdgePointsVector[0];
	TVector3d FirstEdgePoint(AuxVect2d.x, AuxVect2d.y, 0.);

	SetupTranslation(TransForBasePgn.TrPoint(FirstEdgePoint), SubdExtrPolygonP->FirstPoint, AuxTrans1);
	TrProduct(&AuxTrans1, &TransForBasePgn, AuxTrans2);

	//OC 200902
	//adding mirror symmetry to TransForBasePgn, to ensure correct drawing
	radTrans AuxPlSym, AuxTransLoc, AuxTransForBasePgn;
	TVector3d* pFirstPoint = &(SubdExtrPolygonP->FirstPoint);
	TVector3d AuxN(NextPoint.x - pFirstPoint->x, NextPoint.y - pFirstPoint->y, NextPoint.z - pFirstPoint->z);
	SetupPlaneSym(SubdExtrPolygonP->FirstPoint, AuxN, AuxPlSym);
	TrProduct(&AuxPlSym, &AuxTrans2, AuxTransLoc);
	TrProduct(BaseTransPtr, &AuxTransLoc, AuxTransForBasePgn);
	//END OC

	TrProduct(BaseTransPtr, &AuxTrans2, TransForBasePgn);

	SetupTranslation(SubdExtrPolygonP->FirstPoint, NextPoint, AuxTrans1);
	TrProduct(&AuxTrans1, &AuxTrans2, TransForNextBasePgn);
	TrProduct(BaseTransPtr, &TransForNextBasePgn, TransForNextBasePgn);

	//PolygonGraphPresentPtr->Draw(&TransForBasePgn);
	PolygonGraphPresentPtr->Draw(&AuxTransForBasePgn); //OC 200902
	PolygonGraphPresentPtr->Draw(&TransForNextBasePgn);

	TVector3d Side[5], AuxSide[5];
	TVector3d ThcnVect = BaseTransPtr->TrBiPoint(NextPoint - SubdExtrPolygonP->FirstPoint);

	TVector3d AuxVect(0.,0.,0.);

	int LocParity = BaseTransPtr->ShowParity();
	//char SendExtraContourLines = (ShowEdgeLinesInQD3D && (DrawFacilityInd == 1));
	char SendExtraContourLines = (ShowEdgeLines && ((DrawFacilityInd == 1) || (DrawFacilityInd == 2)));

	FirstEdgePoint = TransForBasePgn.TrPoint(FirstEdgePoint);
	Side[0] = FirstEdgePoint;
	Side[3] = Side[0] + ThcnVect;
	for(int k=1; k<BasePolygonPtr->AmOfEdgePoints; k++)
	{
		AuxVect2d = BasePolygonPtr->EdgePointsVector[k];
		AuxVect.x = AuxVect2d.x; AuxVect.y = AuxVect2d.y;
		Side[1] = TransForBasePgn.TrPoint(AuxVect);
		Side[2] = Side[1] + ThcnVect;

		if(LocParity > 0)
		{
			Send.Polygon(Side, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{ 
				SetCurrentColorInStack(SbdLineColor);
				Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		else
		{
			*AuxSide = *Side; AuxSide[1] = Side[3]; AuxSide[2] = Side[2]; AuxSide[3] = Side[1]; 
			Send.Polygon(AuxSide, 4, DrawFacilityInd);
			if(SendExtraContourLines) 
			{
				SetCurrentColorInStack(SbdLineColor);
				AuxSide[4] = *AuxSide; Send.Line(AuxSide, 5, DrawFacilityInd);
				RemoveCurrentColorFromStack();
			}
		}
		Side[0] = Side[1]; Side[3] = Side[2];
	}
	Side[1] = FirstEdgePoint; Side[2] = Side[1] + ThcnVect;

	if(LocParity > 0)
	{
		Send.Polygon(Side, 4, DrawFacilityInd);
		if(SendExtraContourLines) 
		{
			SetCurrentColorInStack(SbdLineColor);
			Side[4] = *Side; Send.Line(Side, 5, DrawFacilityInd);
			RemoveCurrentColorFromStack();
		}
	}
	else
	{
		*AuxSide = *Side; AuxSide[1] = Side[3]; AuxSide[2] = Side[2]; AuxSide[3] = Side[1]; 
		Send.Polygon(AuxSide, 4, DrawFacilityInd);
		if(SendExtraContourLines) 
		{
			SetCurrentColorInStack(SbdLineColor);
			AuxSide[4] = *AuxSide; Send.Line(AuxSide, 5, DrawFacilityInd);
			RemoveCurrentColorFromStack();
		}
	}

	DrawSubdivisionLines(SubdExtrPolygonP, &TransForBasePgn, ThcnVect);
}

//-------------------------------------------------------------------------

void radTSubdivExtrPolygGraphPresent::DrawSubdivisionLines(radTSubdividedExtrPolygon* SubdExtrPolygonP, radTrans* TransForBasePgnPtr, TVector3d& ThcnVect)
{
	const double PointCoinsToler = 1.E-07;

	SetCurrentColorInStack(SbdLineColor);

	int kxLoc = (SubdExtrPolygonP->AxOrnt==ParallelToX)? int(SubdExtrPolygonP->kx) : ((SubdExtrPolygonP->AxOrnt==ParallelToY)? int(SubdExtrPolygonP->ky) : int(SubdExtrPolygonP->kz));
	int AmOfExtraItemsOnMantle = kxLoc - 1;
	int AmOfExtraItemsOnBase = int((double(SubdExtrPolygonP->AmOfSubElem) + 1.E-10)/double(kxLoc));
	if(AmOfExtraItemsOnBase==1) AmOfExtraItemsOnBase = 0;

	int AmOfSubdExtrPolygInTheGroup = 0;
	int AmOfExtrusionLines = 0;

	//radTCast Cast;
	for(radTmhg::const_iterator Iter = SubdExtrPolygonP->GroupMapOfHandlers.begin(); Iter != SubdExtrPolygonP->GroupMapOfHandlers.end(); ++Iter)
	{
		radTPolygon* PolygPtr;
		short SubPolygonsAreOnTheBase=0;

		radTg3d* g3dPtr = (radTg3d*)((*Iter).second.rep);
		radTGroup* GroupPtr = radTCast::GroupCast(g3dPtr); // Because Subdivided ExtrPolygons (and Subd. RecMags) are placed to the general container through the Group branch
		radTSubdividedExtrPolygon* SubdividedExtrPolygPtr = (GroupPtr!=0)? radTCast::SubdExtrPolygonCastFromGroup(GroupPtr) : 0;
		radTExtrPolygon* ExtrPolygPtr = (GroupPtr==0)? radTCast::ExtrPolygonCast((radTg3dRelax*)g3dPtr) : 0;

		if(ExtrPolygPtr != 0)
		{
			PolygPtr = (radTPolygon*)(ExtrPolygPtr->BasePolygonHandle.rep);
			if(SubdExtrPolygonP->AxOrnt==ParallelToX) SubPolygonsAreOnTheBase = (Abs(ExtrPolygPtr->FirstPoint.x - SubdExtrPolygonP->FirstPoint.x) < PointCoinsToler)? 1 : 0;
			else if(SubdExtrPolygonP->AxOrnt==ParallelToY) SubPolygonsAreOnTheBase = (Abs(ExtrPolygPtr->FirstPoint.y - SubdExtrPolygonP->FirstPoint.y) < PointCoinsToler)? 1 : 0;
			else SubPolygonsAreOnTheBase = (Abs(ExtrPolygPtr->FirstPoint.z - SubdExtrPolygonP->FirstPoint.z) < PointCoinsToler)? 1 : 0;
		}
		else if(SubdividedExtrPolygPtr != 0)
		{
			AmOfSubdExtrPolygInTheGroup++;

			PolygPtr = (radTPolygon*)(SubdividedExtrPolygPtr->BasePolygonHandle.rep);
			if(SubdExtrPolygonP->AxOrnt==ParallelToX) SubPolygonsAreOnTheBase = (Abs(SubdividedExtrPolygPtr->FirstPoint.x - SubdExtrPolygonP->FirstPoint.x) < PointCoinsToler)? 1 : 0;
			else if(SubdExtrPolygonP->AxOrnt==ParallelToY) SubPolygonsAreOnTheBase = (Abs(SubdividedExtrPolygPtr->FirstPoint.y - SubdExtrPolygonP->FirstPoint.y) < PointCoinsToler)? 1 : 0;
			else SubPolygonsAreOnTheBase = (Abs(SubdividedExtrPolygPtr->FirstPoint.z - SubdExtrPolygonP->FirstPoint.z) < PointCoinsToler)? 1 : 0;
		}
		if((AmOfExtraItemsOnBase>0) && SubPolygonsAreOnTheBase)	for(int i_Polyg = 0; i_Polyg < PolygPtr->AmOfEdgePoints; i_Polyg++) AmOfExtrusionLines++;
	}

	int AmOfExtraItems = AmOfExtraItemsOnMantle + 2*AmOfExtraItemsOnBase + AmOfExtrusionLines + AmOfSubdExtrPolygInTheGroup;

	radTDrawAttr LocDrawAttr;
	LocDrawAttr.RGB_col.Red = LocDrawAttr.RGB_col.Green = LocDrawAttr.RGB_col.Blue = 0.;
	LocDrawAttr.LineThickness = DrawAttr.LineThickness;

	Send.InitDrawLinElem(DrawAttrAreSet, LocDrawAttr, AmOfExtraItems, DrawFacilityInd);

// Preparing and Sending Some Mantle Lines
	TVector3d ExtraBaseArray[200]; // Hope it's sufficient
	TVector3d MantleArray[200];

	double qxLoc = (SubdExtrPolygonP->AxOrnt==ParallelToX)? SubdExtrPolygonP->qx : ((SubdExtrPolygonP->AxOrnt==ParallelToY)? SubdExtrPolygonP->qy : SubdExtrPolygonP->qz);

	const double AbsZeroTol = 5.E-12;
	double q0x = (fabs(kxLoc-1.)>AbsZeroTol)? pow(qxLoc, 1./(kxLoc-1.)) : qxLoc;
	double BufX = qxLoc*q0x - 1.;
	double a1x_R = (fabs(BufX) > AbsZeroTol)? (q0x - 1.)/BufX : 1./kxLoc;

	double BufRat = a1x_R;
	TVector3d SmallTrslVect = BufRat*ThcnVect;

	TVector2d* AuxVect2dPtr;
	TVector3d AuxVect(0.,0.,0.), ZeroVect(0.,0.,0.);
	radTPolygon* BasePolygonPtr = (radTPolygon*)(SubdExtrPolygonP->BasePolygonHandle.rep);

	for(int ii=0; ii<BasePolygonPtr->AmOfEdgePoints; ii++)
	{
		AuxVect2dPtr = &(BasePolygonPtr->EdgePointsVector[ii]);
		AuxVect.x = AuxVect2dPtr->x; AuxVect.y = AuxVect2dPtr->y;
		MantleArray[ii] = TransForBasePgnPtr->TrPoint(AuxVect);
	}
	for(int i_Mantle = 0; i_Mantle < kxLoc-1; i_Mantle++)
	{
		for(int i=0; i<BasePolygonPtr->AmOfEdgePoints; i++) MantleArray[i] += SmallTrslVect;
		MantleArray[BasePolygonPtr->AmOfEdgePoints] = MantleArray[0];
		Send.Line(MantleArray, BasePolygonPtr->AmOfEdgePoints + 1, DrawFacilityInd);

		BufRat *= q0x;
		SmallTrslVect = BufRat*ThcnVect;
	}

// Preparing and Sending Base Lines
	if((AmOfExtraItemsOnBase>0) || (AmOfSubdExtrPolygInTheGroup>0))
	{
		AuxVect = ZeroVect;
		for(radTmhg::const_iterator iter = SubdExtrPolygonP->GroupMapOfHandlers.begin();
			iter != SubdExtrPolygonP->GroupMapOfHandlers.end(); ++iter)
		{
			radTPolygon* PolygPtr;
			short SubPolygonsAreOnTheBase=0;

			radTg3d* g3dPtr = (radTg3d*)((*iter).second.rep);
			radTGroup* GroupPtr = radTCast::GroupCast(g3dPtr); // Because Subdivided ExtrPolygons (and Subd. RecMags) are placed to the general container through the Group branch
			radTSubdividedExtrPolygon* SubdividedExtrPolygPtr = (GroupPtr!=0)? radTCast::SubdExtrPolygonCastFromGroup(GroupPtr) : 0;
			radTExtrPolygon* ExtrPolygPtr = (GroupPtr==0)? radTCast::ExtrPolygonCast((radTg3dRelax*)g3dPtr) : 0;

			if(ExtrPolygPtr != 0)
			{
				PolygPtr = (radTPolygon*)(ExtrPolygPtr->BasePolygonHandle.rep);
				if(SubdExtrPolygonP->AxOrnt==ParallelToX) SubPolygonsAreOnTheBase = (Abs(ExtrPolygPtr->FirstPoint.x - SubdExtrPolygonP->FirstPoint.x) < PointCoinsToler)? 1 : 0;
				else if(SubdExtrPolygonP->AxOrnt==ParallelToY) SubPolygonsAreOnTheBase = (Abs(ExtrPolygPtr->FirstPoint.y - SubdExtrPolygonP->FirstPoint.y) < PointCoinsToler)? 1 : 0;
				else SubPolygonsAreOnTheBase = (Abs(ExtrPolygPtr->FirstPoint.z - SubdExtrPolygonP->FirstPoint.z) < PointCoinsToler)? 1 : 0;
			}
			else if(SubdividedExtrPolygPtr != 0)
			{
				PolygPtr = (radTPolygon*)(SubdividedExtrPolygPtr->BasePolygonHandle.rep);
				if(SubdExtrPolygonP->AxOrnt==ParallelToX) SubPolygonsAreOnTheBase = (Abs(SubdividedExtrPolygPtr->FirstPoint.x - SubdExtrPolygonP->FirstPoint.x) < PointCoinsToler)? 1 : 0;
				else if(SubdExtrPolygonP->AxOrnt==ParallelToY) SubPolygonsAreOnTheBase = (Abs(SubdividedExtrPolygPtr->FirstPoint.y - SubdExtrPolygonP->FirstPoint.y) < PointCoinsToler)? 1 : 0;
				else SubPolygonsAreOnTheBase = (Abs(SubdividedExtrPolygPtr->FirstPoint.z - SubdExtrPolygonP->FirstPoint.z) < PointCoinsToler)? 1 : 0;
				
				double Shift;
				if(SubdExtrPolygonP->AxOrnt==ParallelToX) Shift = SubdExtrPolygonP->FirstPoint.x - SubdividedExtrPolygPtr->FirstPoint.x;
				else if(SubdExtrPolygonP->AxOrnt==ParallelToY) Shift = SubdExtrPolygonP->FirstPoint.y - SubdividedExtrPolygPtr->FirstPoint.y;
				else if(SubdExtrPolygonP->AxOrnt==ParallelToZ) Shift = SubdExtrPolygonP->FirstPoint.z - SubdividedExtrPolygPtr->FirstPoint.z;

				TVector3d ExtraTranslVect = (-Shift/SubdExtrPolygonP->Thickness)*ThcnVect;

				radTrans ExtraTranslation, LocTransForBasePgn;
				SetupTranslation(ZeroVect, ExtraTranslVect, ExtraTranslation);
				TrProduct(&ExtraTranslation, TransForBasePgnPtr, LocTransForBasePgn);

				TVector3d LocThcnVect = (SubdividedExtrPolygPtr->Thickness/SubdExtrPolygonP->Thickness)*ThcnVect;

				DrawSubdivisionLines(SubdividedExtrPolygPtr, &LocTransForBasePgn, LocThcnVect);
			}
			if(SubPolygonsAreOnTheBase)
			{
				TVector3d ExtrusionLinesArray[2];
				int i_Polyg;
				for(i_Polyg = 0; i_Polyg < PolygPtr->AmOfEdgePoints; i_Polyg++)
				{
					AuxVect2dPtr = &(PolygPtr->EdgePointsVector[i_Polyg]);
					AuxVect.x = AuxVect2dPtr->x; AuxVect.y = AuxVect2dPtr->y;
					ExtraBaseArray[i_Polyg] = TransForBasePgnPtr->TrPoint(AuxVect);
					TVector3d& BasePo = ExtraBaseArray[i_Polyg];

					ExtrusionLinesArray[0] = BasePo;
					ExtrusionLinesArray[1] = BasePo + ThcnVect;
					Send.Line(ExtrusionLinesArray, 2, DrawFacilityInd);
				}
				ExtraBaseArray[PolygPtr->AmOfEdgePoints] = ExtraBaseArray[0];
				Send.Line(ExtraBaseArray, PolygPtr->AmOfEdgePoints + 1, DrawFacilityInd);

				for(i_Polyg = 0; i_Polyg < PolygPtr->AmOfEdgePoints; i_Polyg++) 
					ExtraBaseArray[i_Polyg] = ExtraBaseArray[i_Polyg] + ThcnVect;
				ExtraBaseArray[PolygPtr->AmOfEdgePoints] = ExtraBaseArray[0];
				Send.Line(ExtraBaseArray, PolygPtr->AmOfEdgePoints + 1, DrawFacilityInd);
			}
		}
	}
	RemoveCurrentColorFromStack();
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTSubdivPolyhedronGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTSubdividedPolyhedron* SubdPolyhedronP = (radTSubdividedPolyhedron*)((radTGroup*)g3dPtr);

	int AmOfElemToDraw = 0, ElemCount = 0;
	int* ItemIsNotFullyInternal = new int[SubdPolyhedronP->GroupMapOfHandlers.size()];

	radTmhg& GrMapOfHndl = SubdPolyhedronP->GroupMapOfHandlers;
	radTmhg::const_iterator iter;
	for(iter = GrMapOfHndl.begin(); iter != GrMapOfHndl.end(); ++iter)
	{
		radTg3d* Loc_g3dPtr = (radTg3d*)(((*iter).second).rep);
		int CurrentItemIsNotInternal = Loc_g3dPtr->ItemIsNotFullyInternalAfterCut();
		ItemIsNotFullyInternal[ElemCount++] = CurrentItemIsNotInternal;
		AmOfElemToDraw += CurrentItemIsNotInternal;
	}

			//test
			//char ErrorMesTitle[] = "SRW Debug";
			//char ErrorStr[100];
			//int j = sprintf(ErrorStr, "SubdivPolyhedron: AmOfElemToDraw: %d", AmOfElemToDraw);
			//UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
			//int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle); 
			//end test

	if(AmOfElemToDraw > 0)
	{
		Send.InitOutList(AmOfElemToDraw, DrawFacilityInd);
		ElemCount = 0;
		for(iter = GrMapOfHndl.begin(); iter != GrMapOfHndl.end(); ++iter)
		{
			if(ItemIsNotFullyInternal[ElemCount++])
			{
				radTg3d* Loc_g3dPtr = (radTg3d*)(((*iter).second).rep);
				radTg3dGraphPresent* g3dGraphPresentPtr = Loc_g3dPtr->CreateGraphPresent();
				if(g3dGraphPresentPtr != NULL)
				{
					g3dGraphPresentPtr->DrawFacilityInd = DrawFacilityInd;
					g3dGraphPresentPtr->ShowInternalFacesAfterCut = false;
				
					g3dGraphPresentPtr->GraphPresOptions = GraphPresOptions;
					g3dGraphPresentPtr->MapOfDrawAttrPtr = MapOfDrawAttrPtr;

					g3dGraphPresentPtr->ShowEdgeLines = ShowEdgeLines; //OC071002
					g3dGraphPresentPtr->ShowFaces = ShowFaces; //OC071002

					if(!DrawAttrAreSet) g3dGraphPresentPtr->RetrieveDrawAttr((*iter).first);
					else
					{
						g3dGraphPresentPtr->DrawAttrAreSet = DrawAttrAreSet;
						g3dGraphPresentPtr->DrawAttr = DrawAttr;
					}
					g3dGraphPresentPtr->GenTrans = *BaseTransPtr;
					g3dGraphPresentPtr->GenDraw();

					delete g3dGraphPresentPtr;
				}
			}
		}
	}
	if(ItemIsNotFullyInternal != 0) delete[] ItemIsNotFullyInternal;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTSubdivArcCurGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTSubdividedArcCur* SubdArcCurP = (radTSubdividedArcCur*)((radTGroup*)g3dPtr);

	int AmOfElemToDraw = 0, ElemCount = 0;
	int* ItemIsNotFullyInternal = new int[SubdArcCurP->GroupMapOfHandlers.size()];

	radTmhg& GrMapOfHndl = SubdArcCurP->GroupMapOfHandlers;
	radTmhg::const_iterator iter;
	for(iter = GrMapOfHndl.begin(); iter != GrMapOfHndl.end(); ++iter)
	{
		radTg3d* Loc_g3dPtr = (radTg3d*)(((*iter).second).rep);
		int CurrentItemIsNotInternal = Loc_g3dPtr->ItemIsNotFullyInternalAfterCut();
		ItemIsNotFullyInternal[ElemCount++] = CurrentItemIsNotInternal;
		AmOfElemToDraw += CurrentItemIsNotInternal;
	}

	if(AmOfElemToDraw > 0)
	{
		Send.InitOutList(AmOfElemToDraw, DrawFacilityInd);
		ElemCount = 0;
		for(iter = GrMapOfHndl.begin(); iter != GrMapOfHndl.end(); ++iter)
		{
			if(ItemIsNotFullyInternal[ElemCount++])
			{
				radTg3d* Loc_g3dPtr = (radTg3d*)(((*iter).second).rep);
				radTg3dGraphPresent* g3dGraphPresentPtr = Loc_g3dPtr->CreateGraphPresent();
				if(g3dGraphPresentPtr != NULL)
				{
					g3dGraphPresentPtr->DrawFacilityInd = DrawFacilityInd;
					g3dGraphPresentPtr->ShowInternalFacesAfterCut = false;
				
					g3dGraphPresentPtr->GraphPresOptions = GraphPresOptions;
					g3dGraphPresentPtr->MapOfDrawAttrPtr = MapOfDrawAttrPtr;

					if(!DrawAttrAreSet) g3dGraphPresentPtr->RetrieveDrawAttr((*iter).first);
					else
					{
						g3dGraphPresentPtr->DrawAttrAreSet = DrawAttrAreSet;
						g3dGraphPresentPtr->DrawAttr = DrawAttr;
					}
					g3dGraphPresentPtr->GenTrans = *BaseTransPtr;
					g3dGraphPresentPtr->GenDraw();
					delete g3dGraphPresentPtr;
				}
			}
		}
	}
	if(ItemIsNotFullyInternal != 0) delete[] ItemIsNotFullyInternal;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRectangleGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTRectangle* RectangleP = (radTRectangle*)g3dPtr;

	TVector3d CornPoi(RectangleP->CentrPoint.x - 0.5 * RectangleP->Dimensions.x,
					  RectangleP->CentrPoint.y - 0.5 * RectangleP->Dimensions.y,
					  RectangleP->CentrPoint.z);
	TVector3d CornVect_x(RectangleP->Dimensions.x, 0., 0.);
	TVector3d CornVect_y(0., RectangleP->Dimensions.y, 0.);

	CornPoi = BaseTransPtr->TrPoint(CornPoi);
	CornVect_x = BaseTransPtr->TrBiPoint(CornVect_x);
	CornVect_y = BaseTransPtr->TrBiPoint(CornVect_y);

	TVector3d CP_pl_CVx = CornPoi + CornVect_x;
	TVector3d CP_pl_CVy = CornPoi + CornVect_y;

	TVector3d Side[4];
	Side[0]=CornPoi; Side[1]=CP_pl_CVx; Side[2]=CP_pl_CVx + CornVect_y; Side[3]=CP_pl_CVy; Send.Polygon(Side, 4, DrawFacilityInd);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTPolygonGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTPolygon* PolygonP = (radTPolygon*)g3dPtr;
	TVector3d* PolygonArray = NULL;
	PolygonArray = new TVector3d[PolygonP->AmOfEdgePoints + 1];

	TVector2d AuxVect2d;
	TVector3d AuxVect(0.,0., PolygonP->CoordZ);

	int LocParity = BaseTransPtr->ShowParity();
	if(LocParity > 0)
	{
		for(int i=0; i<PolygonP->AmOfEdgePoints; i++)
		{
			AuxVect2d = PolygonP->EdgePointsVector[i];
			AuxVect.x = AuxVect2d.x; AuxVect.y = AuxVect2d.y;
			PolygonArray[i] = BaseTransPtr->TrPoint(AuxVect);
		}
	}
	else
	{
		int AmPo_mi_1 = PolygonP->AmOfEdgePoints - 1;
		for(int i=0; i<PolygonP->AmOfEdgePoints; i++)
		{
			AuxVect2d = PolygonP->EdgePointsVector[AmPo_mi_1 - i];
			AuxVect.x = AuxVect2d.x; AuxVect.y = AuxVect2d.y;
			PolygonArray[i] = BaseTransPtr->TrPoint(AuxVect);
		}
	}

	Send.Polygon(PolygonArray, PolygonP->AmOfEdgePoints, DrawFacilityInd);

	//char SendExtraContourLines = (ShowEdgeLinesInQD3D && (DrawFacilityInd == 1));
	char SendExtraContourLines = (ShowEdgeLines && ((DrawFacilityInd == 1) || (DrawFacilityInd == 2)));

	if(SendExtraContourLines) 
	{ 
		SetCurrentColorInStack(SbdLineColor);
		PolygonArray[PolygonP->AmOfEdgePoints] = *PolygonArray; 
		Send.Line(PolygonArray, PolygonP->AmOfEdgePoints + 1, DrawFacilityInd);
		RemoveCurrentColorFromStack();
	}

	delete[] PolygonArray;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTFlmLinCurGraphPresent::Draw(radTrans* BaseTransPtr)
{
	radTFlmLinCur* FlmLinCurP = (radTFlmLinCur*)g3dPtr;
	TVector3d aLine[2];
	aLine[0] = BaseTransPtr->TrPoint(FlmLinCurP->StartPoint);
	aLine[1] = BaseTransPtr->TrPoint(FlmLinCurP->EndPoint);

#ifdef _WITH_QD3D
	//char OldShowLinesInQD3D = Send.ShowLinesInQD3D;
	//Send.ShowLinesInQD3D = 1;
	char OldShowLines = Send.ShowLines;
	Send.ShowLines = 1;
#endif

	Send.Line(aLine, 2, DrawFacilityInd);

#ifdef _WITH_QD3D
	//Send.ShowLinesInQD3D = OldShowLinesInQD3D;
	Send.ShowLines = OldShowLines;
#endif
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
