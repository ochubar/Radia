/*-------------------------------------------------------------------------
*
* File name:      radplnr.h
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

#ifndef __RADPLNR_H
#define __RADPLNR_H

#include "radsend.h"
#include "radg3d.h"

#include <algorithm>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radTRectangleSurfIntData {
	TVector3d PointOnSurface;
	int IntegrandLen;
	void (radTg3d::*IntegrandFunPtr)(radTField*);
	radTField Field;

	double* InnerAbsPrecAndLimitsArray;
	short* InnerElemCompNotFinished;
	TVector3d** InnerIntegVal;
};

//-------------------------------------------------------------------------

class radTRectangle : public radTg3d {
	radTRectangleSurfIntData* SurfIntDataPtr;

public:
	//TVector3d CentrPoint; //moved to radTg3d
	TVector2d Dimensions;
	// Something concerning material

	radTRectangle(const TVector3d& InCentrPoint, const TVector2d& InDimensions)
	{
		CentrPoint = InCentrPoint; Dimensions = InDimensions;
	}
	radTRectangle() {}

	int Type_g3d() { return 5;}

	void IntOverShape(radTField* FieldPtr) 
	{
		if(FieldPtr->ShapeIntDataPtr->IntOverSurf_) IntOverSurf(FieldPtr);
		else if(FieldPtr->ShapeIntDataPtr->IntOverLine_) IntOverLine(FieldPtr);
	}
	void IntOverSurf(radTField*);
	inline void FunForOuterIntAtSurfInt(double, TVector3d*);
	inline void FunForInnerIntAtSurfInt(double, TVector3d*);
	void IntOverLine(radTField*) {}

	int NumberOfDegOfFreedom() { return 0;}

	void Dump(std::ostream& o, int ShortSign) // Porting
	{
		radTg3d::Dump(o, ShortSign);
		o << "Rectangle";

		if(ShortSign) return;

		o << endl;
		o << "   {x,y,z}= {" << CentrPoint.x << ',' << CentrPoint.y << ',' << CentrPoint.z << "}" << endl;
		o << "   {wx,wy}= {" << Dimensions.x << ',' << Dimensions.y << "}";

		DumpTransApplied(o);

		o << endl;
		o << "   Memory occupied: " << SizeOfThis() << " bytes";
	}
	radTg3dGraphPresent* CreateGraphPresent();

	int SizeOfThis() { return sizeof(radTRectangle);}
};
//-------------------------------------------------------------------------

inline void radTRectangle::FunForOuterIntAtSurfInt(double Arg, TVector3d* VectArray)
{
	const double PrecEnhFact = 1.; // Don't make it >1 : it's dangerous for convergence of outer itegral !
	double* OuterIntPrecArray = SurfIntDataPtr->Field.ShapeIntDataPtr->AbsPrecArray;

	double SmallPositive = 1.E-10;

	for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++)
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[i] = PrecEnhFact*OuterIntPrecArray[i]/Dimensions.y;
	(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen] = CentrPoint.x - 0.5*Dimensions.x + SmallPositive;
	(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen + 1] = CentrPoint.x + 0.5*Dimensions.x;
	SurfIntDataPtr->PointOnSurface.y = Arg;

	FormalOneFoldInteg(this, &radTRectangle::FunForInnerIntAtSurfInt, SurfIntDataPtr->IntegrandLen, 
					   SurfIntDataPtr->InnerAbsPrecAndLimitsArray, 
					   SurfIntDataPtr->InnerElemCompNotFinished, SurfIntDataPtr->InnerIntegVal);

	for(int ii=0; ii<SurfIntDataPtr->IntegrandLen; ii++) VectArray[ii] = ((SurfIntDataPtr->InnerIntegVal)[0])[ii];
}

//-------------------------------------------------------------------------

inline void radTRectangle::FunForInnerIntAtSurfInt(double Arg, TVector3d* VectArray)
{
	SurfIntDataPtr->PointOnSurface.x = Arg;

	SurfIntDataPtr->Field.P = SurfIntDataPtr->PointOnSurface;
	(((radTg3d*)this)->*(SurfIntDataPtr->IntegrandFunPtr))(&(SurfIntDataPtr->Field));

	for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++) 
		VectArray[i] = (SurfIntDataPtr->Field.ShapeIntDataPtr->VectArray)[i];
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

enum TBoundOrLine { Bound, ParLine };
enum TLinesIntrsctCase { PointWithinBound, PointOnBoundEdge, PointOutsideBound, LineIsIntrsct, Zero };

//-------------------------------------------------------------------------

struct radTPolyg2dIntrsctInfo {
	TVector2d IntrsctPoint;
	int NoOfIntrsctItem;
	int NoOfIntrsctPoOnIntrsctItem;
	TBoundOrLine TypeOfIntrsctItem;
	short SkipThisPoint;

	int IntrsctMultiplicity;
	short IntrsctIsLine;

	radTPolyg2dIntrsctInfo() { SkipThisPoint = 0;}

	inline friend int operator <(const radTPolyg2dIntrsctInfo&, const radTPolyg2dIntrsctInfo&);
	inline friend int operator ==(const radTPolyg2dIntrsctInfo&, const radTPolyg2dIntrsctInfo&);
};

//-------------------------------------------------------------------------

inline int operator <(const radTPolyg2dIntrsctInfo&, const radTPolyg2dIntrsctInfo&) { return 1;}

//-------------------------------------------------------------------------

inline int operator ==(const radTPolyg2dIntrsctInfo& inf1, const radTPolyg2dIntrsctInfo& inf2) 
{ 
	return (inf1.IntrsctPoint == inf2.IntrsctPoint) && (inf1.NoOfIntrsctItem == inf2.NoOfIntrsctItem)
		&& (inf1.NoOfIntrsctPoOnIntrsctItem == inf2.NoOfIntrsctPoOnIntrsctItem) 
		&& (inf1.TypeOfIntrsctItem == inf2.TypeOfIntrsctItem)
		&& (inf1.SkipThisPoint == inf2.SkipThisPoint) 
		&& (inf1.IntrsctMultiplicity == inf2.IntrsctMultiplicity)
		&& (inf1.IntrsctIsLine == inf2.IntrsctIsLine);
}

//-------------------------------------------------------------------------

#ifdef __MWERKS__
/*
null_template
struct iterator_traits <radTPolyg2dIntrsctInfo*> {
     typedef ptrdiff_t difference_type;
     typedef radTPolyg2dIntrsctInfo value_type;
     typedef radTPolyg2dIntrsctInfo* pointer;
     typedef radTPolyg2dIntrsctInfo& reference;
     typedef random_access_iterator_tag iterator_category;
};
*/
#endif

//-------------------------------------------------------------------------

class TCompareIntrsctInfo {
public:
	TVector2d V;
	double V_Toler; // To recognize Zero V.x | V.y

	TCompareIntrsctInfo(const TVector2d& InV) { V = InV; V_Toler = 5.E-14;}

	bool operator()(const radTPolyg2dIntrsctInfo& Info1, const radTPolyg2dIntrsctInfo& Info2)
	{ // To check !
		double vAbsE2 = V.x*V.x + V.y*V.y; //OC04032017
		double V_TolerE2 = V_Toler*V_Toler;
		if(vAbsE2 < V_TolerE2) return (Info1.NoOfIntrsctItem > Info2.NoOfIntrsctItem); //OC04032017

		if((Info1.IntrsctPoint.x == Info2.IntrsctPoint.x) && (Info1.IntrsctPoint.y == Info2.IntrsctPoint.y))
		{
			//return (V.y > V_Toler)? (Info1.NoOfIntrsctItem < Info2.NoOfIntrsctItem) : ((V.y < -V_Toler)? (Info1.NoOfIntrsctItem > Info2.NoOfIntrsctItem) : ((V.x > V_Toler)? (Info1.NoOfIntrsctItem > Info2.NoOfIntrsctItem) : (Info1.NoOfIntrsctItem < Info2.NoOfIntrsctItem)));
			return (Info1.NoOfIntrsctItem < Info2.NoOfIntrsctItem); //OC04032017
		}
		else
		{
		//	if(V.x > V_Toler) return Info1.IntrsctPoint.x < Info2.IntrsctPoint.x;
		//	else if(V.x < -V_Toler) return Info1.IntrsctPoint.x > Info2.IntrsctPoint.x;
		//	else return (V.y > V_Toler)? (Info1.IntrsctPoint.y < Info2.IntrsctPoint.y) : (Info1.IntrsctPoint.y > Info2.IntrsctPoint.y);
			//OC04032017
			//Compare projections of all points on V
			return ((Info1.IntrsctPoint*V) < (Info2.IntrsctPoint*V));
		}
	}
};

//-------------------------------------------------------------------------

#ifdef __GCC__
typedef vector<radTPolyg2dIntrsctInfo> radTPolyg2dIntrsctInfoVect;
typedef vector<TVector2d> radTVect2dVect;
typedef vector<int> radTVectInt;
#else
typedef vector<radTPolyg2dIntrsctInfo, allocator<radTPolyg2dIntrsctInfo> > radTPolyg2dIntrsctInfoVect;
typedef vector<TVector2d, allocator<TVector2d> > radTVect2dVect;
typedef vector<int, allocator<int> > radTVectInt;
#endif

#ifdef __MWERKS__
/*
null_template
struct iterator_traits <TVector2d*> {
     typedef ptrdiff_t difference_type;
     typedef TVector2d value_type;
     typedef TVector2d* pointer;
     typedef TVector2d& reference;
     typedef random_access_iterator_tag iterator_category;
};
*/
#endif

//-------------------------------------------------------------------------

enum TSpecCaseID { ZeroVxVy, ZeroVxVz, ZeroVyVz, NoSpecCase };

//-------------------------------------------------------------------------

class radTPolygon : public radTg3d {
public:
	radTVect2dVect EdgePointsVector;
	int AmOfEdgePoints;
	TVector2d CentrPoint;

	double CoordZ;
	TVector3d Magn;

	short SomethingIsWrong;
	char IsConvex;

	radTPolygon(TVector2d*, int);
	radTPolygon(const radTVect2dVect&);
	radTPolygon(double, TVector2d*, int, const TVector3d&);
	radTPolygon(radTVect2dVect&, double, const TVector3d&);
	radTPolygon(CAuxBinStrVect& inStr);
	radTPolygon() { CoordZ =0.; SomethingIsWrong = 0;}

	virtual int Type_g3d() { return 6;}

	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{
		return FinishDuplication(new radTPolygon(*this), hg);
	}

	int SubdivideItself(double*, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideBySetOfParallelLines(const TVector2d&, TVector2d*, int, radTvhg&);
	void FindStPointsForIntrsctLines(const TVector2d&, const TVector2d&, double, double, double, double, TVector2d*, TVector2d*);
	void FillInIntrsctInfoStruct(const TVector2d&, const TVector2d*, int, radTPolyg2dIntrsctInfoVect*, radTPolyg2dIntrsctInfoVect*, double&);
	void IntrsctOfTwoLines(const TVector2d&, const TVector2d&, const TVector2d&, const TVector2d&, TVector2d&, TLinesIntrsctCase&);
	void IntrsctOfTwoLines2(const TVector2d&, const TVector2d&, const TVector2d&, const TVector2d&, TVector2d&, TLinesIntrsctCase&);
	void ComputeCentrPoint(short&);
	void SimpleComputeCentrPoint(short&, char& Out_PointsAreDifferent);
	int CheckIfConvex();
	int RandomizeNonConvexEdgePoints(radTVect2dVect&);

	int CheckIfNotSelfIntersecting();
	void CheckAndRearrangeEdgePoints(TVector2d*, int);
	void CheckAndRearrangeEdgePoints(radTVect2dVect&);

	void B_comp(radTField*);
	void B_intComp(radTField*);
	void B_intCompSpecCases(radTField*, const TSpecCaseID&);

	//void B_comp_frJ(radTField*);

	radTg3dGraphPresent* CreateGraphPresent();
	void Dump(std::ostream& o, int ShortSign) // Porting
	{
		radTg3d::Dump(o, ShortSign);
		o << "Polygon";

		o << endl;
		o << "   Memory occupied: " << SizeOfThis() << " bytes";
	}

	void DumpBin_Polygon(CAuxBinStrVect& oStr)
	{
		//int AmOfEdgePoints;
		oStr << AmOfEdgePoints;

		//radTVect2dVect EdgePointsVector;
		for(int i=0; i<AmOfEdgePoints; i++)
		{
			oStr << EdgePointsVector[i];
		}

		//TVector2d CentrPoint;
		oStr << radTPolygon::CentrPoint;

		//double CoordZ;
		oStr << CoordZ;

		//TVector3d Magn;
		oStr << Magn;

		//short SomethingIsWrong;
		oStr << SomethingIsWrong;

		//char IsConvex;
		oStr << IsConvex;
	}

	int NumberOfDegOfFreedom() { return 0;}
	int SizeOfThis() 
	{
		int GenSize = sizeof(radTPolygon);
		GenSize += AmOfEdgePoints*sizeof(TVector2d);
		return GenSize;
	}

	double Area()
	{
		const double Max_k = 1.E+10;
		radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
		int AmOfEdPoInBase_mi_1 = AmOfEdgePoints - 1;
		TVector2d First2d = *BaseIter;
		double x1 = First2d.x, y1 = First2d.y, x2, y2;
		double x1e2 = x1*x1, x2e2, SS = 0.;
		for(int i=0; i<AmOfEdgePoints; i++)
		{
			if(i!=AmOfEdPoInBase_mi_1) { x2 = (*(++BaseIter)).x; y2 = (*BaseIter).y;}
			else { x2 = First2d.x; y2 = First2d.y;}
			x2e2 = x2*x2;
			double x2mx1 = x2-x1, y2my1 = y2-y1;
			double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);
			if(abs_x2mx1*Max_k > abs_y2my1)
			{
				double k = (y2-y1)/(x2-x1), b = y1 - k*x1, x1px2 = x1 + x2; 
				SS += x2mx1*(2.*b + k*x1px2);
			}
			x1 = x2; y1 = y2; x1e2 = x2e2;
		}
		return Abs(0.5*SS);
	}
	
	int IndexOfGoodThirdPoint() // i.e. not belonging to one line with first two points
	{
		const double RelTol = 1.E-09;
		TVector2d p0 = EdgePointsVector[0], p1 = EdgePointsVector[1], p2;
		TVector2d v0 = p1 - p0, v1;
		v0 = (1./sqrt(v0.x*v0.x + v0.y*v0.y))*v0;
		for(int k=2; k<=AmOfEdgePoints; k++)
		{
			v1 = EdgePointsVector[(k==AmOfEdgePoints)? 0 : k] - p1;
			v1 = (1./sqrt(v1.x*v1.x + v1.y*v1.y))*v1;
			if(fabs(1. - fabs(v0*v1)) > RelTol) return k;
		}
		return 0;
	}

	void EstimateSize(double& SizeX, double& SizeY)
	{
		double MinX = 1.E+23, MaxX = -1.E+23, MinY = 1.E+23, MaxY = -1.E+23;
		for(int k=0; k<AmOfEdgePoints; k++)
		{
			TVector2d& P = EdgePointsVector[k];
			if(P.x < MinX) MinX = P.x;
			if(P.x > MaxX) MaxX = P.x;
			if(P.y < MinY) MinY = P.y;
			if(P.y > MaxY) MaxY = P.y;
		}
		SizeX = MaxX - MinX; SizeY = MaxY - MinY;
	}

	double EstimateTypSize()
	{
		double SizeX, SizeY;
		EstimateSize(SizeX, SizeY);
		return (SizeX > SizeY)? SizeX : SizeY;
	}

	//void ApplyHomothety(double coefHom)
	void ApplyHomothety(double kxH, double kyH, double phi)
	{//assumes that CentrPoint was defined
		//if(coefHom == 0) return;
		if((kxH == 0) || (kyH == 0)) return;

		double cosPhi = cos(phi), sinPhi = sin(phi);
		TVector2d vHomLocX(cosPhi, sinPhi), vHomLocY(-sinPhi, cosPhi);

		for(int i=0; i<AmOfEdgePoints; i++)
		{
			TVector2d vRelVertex = EdgePointsVector[i] - CentrPoint;

			//EdgePointsVector[i] = CentrPoint + (coefHom*vRelVertex);
			EdgePointsVector[i] = CentrPoint + ((kxH*(vHomLocX*vRelVertex))*vHomLocX) + ((kyH*(vHomLocY*vRelVertex))*vHomLocY);
		}
	}

	int GetNormalSign()
	{//to distinguish between "left" and "right" polygons
	 //assumes that Center Point was defined
		TVector2d v1 = EdgePointsVector[0] - CentrPoint;
		TVector2d v2 = EdgePointsVector[1] - CentrPoint;
		double vectProd = v1.x*v2.y - v2.x*v1.y;
		return (vectProd >= 0)? 1 : -1; //to modify
	}
};

//-------------------------------------------------------------------------

#endif
