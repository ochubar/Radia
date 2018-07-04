/*-------------------------------------------------------------------------
*
* File name:      radexpgn.h
*
* Project:        RADIA
*
* Description:    Magnetic field source: extruded polygon (prism)
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADEXPGN_H
#define __RADEXPGN_H

#include "radg3d.h"
#include "radplnr.h"
#include "radtrans.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTExtrPolygon : public radTg3dRelax {
public:
	TVector3d FirstPoint;
	TAxisOrient AxOrnt;
	radThg BasePolygonHandle;
	double Thickness;

	radTExtrPolygon(const TVector3d& InFirstPoint, const TAxisOrient& InAxOrnt, double InThickness,
					const radThg& InBasePolygonHandle, const TVector3d& InMagn, 
					const radThg& InMaterHandle) : radTg3dRelax(InFirstPoint, InMagn, InMaterHandle)
	{
		FirstPoint = InFirstPoint; AxOrnt = InAxOrnt; Thickness = InThickness; Magn = InMagn;
		BasePolygonHandle = InBasePolygonHandle;
		DefineCentrPoint();
	}
	radTExtrPolygon(const TVector3d& InFirstPoint, const TAxisOrient& InAxOrnt, double InThickness,
					TVector2d* ArrayOfPoints2d, int lenArrayOfPoints2d, const TVector3d& InMagn) 
					: radTg3dRelax(InMagn)
	{
		FirstPoint = InFirstPoint; AxOrnt = InAxOrnt; Thickness = InThickness; //Magn = InMagn;
		radThg hg(new radTPolygon(ArrayOfPoints2d, lenArrayOfPoints2d));
		BasePolygonHandle = hg;
		DefineCentrPoint();
	}
	radTExtrPolygon(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{
		DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);
		DumpBinParse_g3dRelax(inStr, mKeysOldNew, gMapOfHandlers);
		DumpBinParse_ExtrPolygon(inStr);
	}
	radTExtrPolygon() : radTg3dRelax() {}

	int Type_g3dRelax() { return 2;}

	void B_comp(radTField*);
	void B_intComp(radTField*);
	void B_intCompSpecCases(radTField*, const TSpecCaseID&);

	void Dump(std::ostream&, int ShortSign =0);
	void DumpPureObjInfo(std::ostream&, int);
	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey);
	void DumpBin_ExtrPolygon(CAuxBinStrVect& oStr);
	void DumpBinParse_ExtrPolygon(CAuxBinStrVect& inStr);

	radTg3dGraphPresent* CreateGraphPresent();

	double Volume() { return Thickness*(((radTPolygon*)(BasePolygonHandle.rep))->Area());}
	void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique);

	int SubdivideItself(double*, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideItselfByOneSetOfParPlanes(TVector3d& InPlanesNormal, TVector3d* InPointsOnCuttingPlanes, int AmOfPieces_mi_1, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions, radTvhg* pVectOfHgChanged)
	{
		radThg& NewHandle = In_hg;
		radThg OldHandle = In_hg;
		if(!((radTg3d*)(OldHandle.rep))->ConvertToPolyhedron(NewHandle, radPtr, pSubdivOptions->PutNewStuffIntoGenCont)) return 0;
		OldHandle = NewHandle;
		if(!((radTg3d*)(OldHandle.rep))->SubdivideItselfByOneSetOfParPlanes(InPlanesNormal, InPointsOnCuttingPlanes, AmOfPieces_mi_1, NewHandle, radPtr, pSubdivOptions, pVectOfHgChanged)) return 0;
		return 1;
	}
	int SubdivideItselfByParPlanes(double* SubdivArray, int AmOfDir, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{
		radThg& NewHandle = In_hg;
		radThg OldHandle = In_hg;
		if(!((radTg3d*)(OldHandle.rep))->ConvertToPolyhedron(NewHandle, radPtr, pSubdivOptions->PutNewStuffIntoGenCont)) return 0;
		OldHandle = NewHandle;
		if(!((radTg3d*)(OldHandle.rep))->SubdivideItselfByParPlanes(SubdivArray, AmOfDir, NewHandle, radPtr, pSubdivOptions)) return 0;
		return 1;
	}
	int CutItself(TVector3d* InCuttingPlane, radThg& In_hg, radTPair_int_hg& LowerNewPair_int_hg, radTPair_int_hg& UpperNewPair_int_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{
		radThg& NewHandle = In_hg;
		radThg OldHandle = In_hg;
		if(!((radTg3d*)(OldHandle.rep))->ConvertToPolyhedron(NewHandle, radPtr, pSubdivOptions->PutNewStuffIntoGenCont)) return 0;
		OldHandle = NewHandle;
		if(!((radTg3d*)(OldHandle.rep))->CutItself(InCuttingPlane, NewHandle, LowerNewPair_int_hg, UpperNewPair_int_hg, radPtr, pSubdivOptions)) return 0;
		return 1;
	}
	int SubdivideItselfByEllipticCylinder(double* SubdivArray, radTCylindricSubdivSpec* pSubdivSpec, radThg& In_hg, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{ 
		radThg& NewHandle = In_hg;
		radThg OldHandle = In_hg;
		if(!((radTg3d*)(OldHandle.rep))->ConvertToPolyhedron(NewHandle, radPtr, pSubdivOptions->PutNewStuffIntoGenCont)) return 0;
		OldHandle = NewHandle;
		if(!((radTg3d*)(OldHandle.rep))->SubdivideItselfByEllipticCylinder(SubdivArray, pSubdivSpec, NewHandle, radPtr, pSubdivOptions)) return 0;
		return 1;
	}

	int FindLowestAndUppestVertices(TVector3d&, radTSubdivOptions*, TVector3d&, TVector3d&, radTrans&, char&, char&);
	int ConvertToPolyhedron(radThg&, radTApplication*, char);

	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{
		return FinishDuplication(new radTExtrPolygon(*this), hg);
	}
	void DefineCentrPoint()
	{
		TVector2d& CentrPo2d = ((radTPolygon*)(BasePolygonHandle.rep))->CentrPoint;
		TVector2d& BaseFirstPo2d = ((radTPolygon*)(BasePolygonHandle.rep))->EdgePointsVector[0];
		CentrPoint.x = FirstPoint.x + 0.5*Thickness;
		CentrPoint.y = FirstPoint.y + CentrPo2d.x - BaseFirstPo2d.x;
		CentrPoint.z = FirstPoint.z + CentrPo2d.y - BaseFirstPo2d.y;
	}
	double ArcTanTwo(double x, double y)
	{
		return atan(y/x) + (x>0. ? 0. : (y<0 ? -3.14159265358979 : 3.14159265358979));
	}
	double PhCorrForArcTanTwo(double x, double y)
	{
		return x>0. ? 0. : (y<0 ? -3.14159265358979 : 3.14159265358979);
	}
	double SumOfFourAtansTwo(double X1, double Y1, double X2, double Y2, double X3, double Y3, double X4, double Y4)
	{
		double PiMult1=0., PiMult2=0., PiMult3=0.;
		return atan(TransAtans(TransAtans(Y1/X1, Y2/X2, PiMult1), TransAtans(Y3/X3, Y4/X4, PiMult2), PiMult3))
			   + PhCorrForArcTanTwo(X1, Y1) + PhCorrForArcTanTwo(X2, Y2) + PhCorrForArcTanTwo(X3, Y3) + PhCorrForArcTanTwo(X4, Y4)
			   + (PiMult1+PiMult2+PiMult3)*3.14159265358979;
	}
	int NumberOfDegOfFreedom() { return 3;}
	int SizeOfThis()
	{
		int GenSize = sizeof(radTExtrPolygon);
		GenSize += sizeof(radTPolygon);
		GenSize += ((radTPolygon*)(BasePolygonHandle.rep))->AmOfEdgePoints * sizeof(TVector2d);
		return GenSize;
	}
};

//-------------------------------------------------------------------------

#endif
