/*-------------------------------------------------------------------------
*
* File name:      radg3dgr.h
*
* Project:        RADIA
*
* Description:    Graphical representations of 3D magnetic field sources
*                 (front-end / interface dependent)
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADG3DGR_H
#define __RADG3DGR_H

#include "radsend.h"
#include "gmvect.h"
#include "radrec.h"
#include "radarccu.h"
#include "radgroup.h"
#include "radflm.h"
//#include "radtrans.h"
#include "gmtrans.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radRGB { 
	double Red, Green, Blue; 
	radRGB(double InRed =0, double InGreen =0, double InBlue =0) 
	{ 
		Red = InRed; Green = InGreen; Blue = InBlue; 
	}
};

//-------------------------------------------------------------------------

struct radTDrawAttr { 
	radRGB RGB_col;
	double LineThickness;
};

//-------------------------------------------------------------------------

struct radGraphPresOptions {
	char ShowSymmetryChilds;
	char doDebug;
	radGraphPresOptions(char InShowSymmetryChilds =1, char inDoDebug =0)
	{
		ShowSymmetryChilds = InShowSymmetryChilds;
		doDebug = inDoDebug;
	}
};

//-------------------------------------------------------------------------

#ifdef __GCC__
typedef	map <int, radTDrawAttr, less<int> >  radTMapOfDrawAttr;

//#ifdef _WITH_QD3D
typedef	vector<radTDrawAttr> radTVectOfDrawAttr;
//#endif

#else
typedef	map <int, radTDrawAttr, less<int> >  radTMapOfDrawAttr;

//#ifdef _WITH_QD3D
typedef	vector<radTDrawAttr, allocator<radTDrawAttr> > radTVectOfDrawAttr;
//#endif
#endif

#ifdef __MWERKS__
#ifdef _WITH_QD3D
/*
null_template
struct iterator_traits <radTDrawAttr*> {
     typedef ptrdiff_t difference_type;
     typedef radTDrawAttr value_type;
     typedef radTDrawAttr* pointer;
     typedef radTDrawAttr& reference;
     typedef random_access_iterator_tag iterator_category;
};
*/
#endif
#endif

//-------------------------------------------------------------------------

class radTg3dGraphPresent {
public:

//#ifdef _WITH_QD3D
	static radTVectOfDrawAttr DrawAttrStack;
//#endif

	static radRGB SbdLineColor;

	radTg3d* g3dPtr;
	radTMapOfDrawAttr* MapOfDrawAttrPtr;

	radGraphPresOptions GraphPresOptions;
	radTDrawAttr DrawAttr;
	static radTSend Send;

	bool ShowInternalFacesAfterCut;
	int DrawAttrAreSet;
	char DrawFacilityInd; // 0- Mathematica, 1- QuickDraw3D, 2- OpenGL
	char ShowEdgeLines; //ShowEdgeLinesInQD3D;
	char ShowFaces; //ShowFacesInQD3D;
	char ShowGridTicks;

	radTrans GenTrans;

	radTg3dGraphPresent(radTg3d*); 
	radTg3dGraphPresent() 
	{ 
		ShowInternalFacesAfterCut = true; DrawFacilityInd = 0;
		ShowEdgeLines = 0; //ShowEdgeLinesInQD3D = 0;
		ShowFaces = 1; //ShowFacesInQD3D = 1;
		ShowGridTicks = 1;
	}

	int RetrieveDrawAttr(int ElemKey);
	void SetGraphPresOptions(const radGraphPresOptions& InGraphPresOptions)
	{
		GraphPresOptions = InGraphPresOptions;
	}
	void SetGraphPresOptionsExt(const radGraphPresOptions& InGraphPresOptions, char DoShowLines, char DoShowFaces)
	{
		GraphPresOptions = InGraphPresOptions;

        Send.ShowLines = ShowEdgeLines = DoShowLines;
        Send.ShowFaces = ShowFaces = DoShowFaces;

		Send.GeomPolygons.reserve(100);
        Send.GeomLines.reserve(100);
        Send.InitLimits3D();
	}

	virtual void Draw(radTrans*) {}
	void NestedFor_Draw(radTrans* BaseTransPtr, const radTlphg::iterator& Iter);
	void DrawOrNestedFor(radTrans* BaseTransPtr, const radTlphg::iterator& Iter)
	{
		if(Iter == g3dPtr->g3dListOfTransform.end()) Draw(BaseTransPtr);
		else NestedFor_Draw(BaseTransPtr, Iter);
	}
	void GenDraw();

	void SetCurrentColorInStack(double r, double g, double b)
	{
		radRGB LinesColor(r, g, b);
		radTDrawAttr LinesDrawAttr;
		LinesDrawAttr.RGB_col = LinesColor;
		DrawAttrStack.insert(DrawAttrStack.begin(), LinesDrawAttr);
	}
	void SetCurrentColorInStack(radRGB& aColor)
	{
		radTDrawAttr LinesDrawAttr;
		LinesDrawAttr.RGB_col = aColor;
		DrawAttrStack.insert(DrawAttrStack.begin(), LinesDrawAttr);
	}
	void RemoveCurrentColorFromStack()
	{
		DrawAttrStack.erase(DrawAttrStack.begin());
	}

	void SetupRotation(const TVector3d&, const TVector3d&, double, radTrans&);
	void SetupTranslation(const TVector3d&, const TVector3d&, radTrans&);
	void SetupPlaneSym(const TVector3d&, const TVector3d&, radTrans&);

	void DrawFrameLines();
	static void DrawPyramidArrow(TVector3d* PyramidArrowInfo, char LocDrawFacilityInd);
    static void DrawCharacter(char Ch, double Ratio, TVector3d* Info3D, char LocDrawFacilityInd);

	static void DrawGridTicks(double* OrigLimits3D, TVector3d& OffsetVect, double AbsCharHeight, char LocDrawFacilityInd);
    static void DrawTickNumbers(TVector3d& P0, TVector3d& Vpar, TVector3d& Vper, int PerDir, vector<double>& TickNumPositions, double TickOffsetPos, double AbsCharHeight, char LocDrawFacilityInd);

	static void FindGridTickPositions(TVector3d* GridLimits, vector<double>* MainTickPositions);
	static double ChooseGridTickInterval(double Lx);
	static double FindGridTickFirstPosition(double Lx, double TickDeltaX);
	static double FindGridTickLength(TVector3d* GridLimits);

};

//-------------------------------------------------------------------------
	
class radTRecMagGraphPresent : public radTg3dGraphPresent {
public:
	radTRecMagGraphPresent(radTg3d* InRecMagPtr) : radTg3dGraphPresent(InRecMagPtr) {}
	radTRecMagGraphPresent() {}

	void Draw(radTrans*);
};

//-------------------------------------------------------------------------

class radTExtrPolygonGraphPresent : public radTg3dGraphPresent {
public:
	radTExtrPolygonGraphPresent(radTg3d* InExtrPolygonPtr) : radTg3dGraphPresent(InExtrPolygonPtr) {}
	radTExtrPolygonGraphPresent() {}

	void Draw(radTrans*);
};

//-------------------------------------------------------------------------

class radTPolyhedronGraphPresent : public radTg3dGraphPresent {
public:
	int AmOfNonInternalFaces;

	radTPolyhedronGraphPresent(radTg3d* InPolyhedronPtr) : radTg3dGraphPresent(InPolyhedronPtr) { ShowInternalFacesAfterCut = true;}
	radTPolyhedronGraphPresent() { ShowInternalFacesAfterCut = true;}

	void Draw(radTrans*);
	void CountNonInternalFaces();
};

//-------------------------------------------------------------------------

class radTArcCurGraphPresent : public radTg3dGraphPresent {
public:
	char OutlineSegments;

	radTArcCurGraphPresent(radTg3d* InArcCurPtr) : radTg3dGraphPresent(InArcCurPtr) { OutlineSegments = 0;}
	radTArcCurGraphPresent() { OutlineSegments = 0;}
	
	void Draw(radTrans*);
	radTrans* CreateRotation(TVector3d&, TVector3d&, double);
};

//-------------------------------------------------------------------------

class radTGroupGraphPresent : public radTg3dGraphPresent {
public:
	radTGroupGraphPresent(radTg3d* InGroupPtr) : radTg3dGraphPresent(InGroupPtr) {}
	radTGroupGraphPresent() {}

	void Draw(radTrans*);
};

//-------------------------------------------------------------------------

class radTSubdividedRecMag;

//-------------------------------------------------------------------------

class radTSubdivRecMagGraphPresent : public radTg3dGraphPresent {
public:
	radTSubdivRecMagGraphPresent(radTg3d* InSubdivRecMagPtr) : radTg3dGraphPresent(InSubdivRecMagPtr) {}
	radTSubdivRecMagGraphPresent() {}

	void Draw(radTrans*);
	void DrawSubdivisionLines(radTSubdividedRecMag*, TVector3d*);
};

//-------------------------------------------------------------------------

class radTSubdividedExtrPolygon;

//-------------------------------------------------------------------------

class radTSubdivExtrPolygGraphPresent : public radTg3dGraphPresent {
public:
	radTSubdivExtrPolygGraphPresent(radTg3d* InSubdivExtrPolygPtr) : radTg3dGraphPresent(InSubdivExtrPolygPtr) {}
	radTSubdivExtrPolygGraphPresent() {}

	void Draw(radTrans*);
	void DrawSubdivisionLines(radTSubdividedExtrPolygon*, radTrans*, TVector3d&);
};

//-------------------------------------------------------------------------

class radTSubdividedPolyhedron;

//-------------------------------------------------------------------------

class radTSubdivPolyhedronGraphPresent : public radTg3dGraphPresent {
public:
	radTSubdivPolyhedronGraphPresent(radTg3d* InSubdivPolyhedronPtr) : radTg3dGraphPresent(InSubdivPolyhedronPtr) {}
	radTSubdivPolyhedronGraphPresent() {}

	void Draw(radTrans*);
};

//-------------------------------------------------------------------------

class radTSubdividedArcCur;

//-------------------------------------------------------------------------

class radTSubdivArcCurGraphPresent : public radTg3dGraphPresent {
public:
	radTSubdivArcCurGraphPresent(radTg3d* InSubdivArcCurPtr) : radTg3dGraphPresent(InSubdivArcCurPtr) {}
	radTSubdivArcCurGraphPresent() {}

	void Draw(radTrans*);
};

//-------------------------------------------------------------------------

class radTRectangleGraphPresent : public radTg3dGraphPresent {
public:
	radTRectangleGraphPresent(radTg3d* InRectanglePtr) : radTg3dGraphPresent(InRectanglePtr) {}
	radTRectangleGraphPresent() {}

	void Draw(radTrans*);
};

//-------------------------------------------------------------------------

class radTPolygonGraphPresent : public radTg3dGraphPresent {
public:
	radTPolygonGraphPresent(radTg3d* InPolygonPtr) : radTg3dGraphPresent(InPolygonPtr) {}
	radTPolygonGraphPresent() {}

	void Draw(radTrans*);
};

//-------------------------------------------------------------------------

class radTFlmLinCurGraphPresent : public radTg3dGraphPresent {
public:
	radTFlmLinCurGraphPresent(radTg3d* InFlmLinCurPtr) : radTg3dGraphPresent(InFlmLinCurPtr) {}
	radTFlmLinCurGraphPresent() {}

	void Draw(radTrans*);
};

//-------------------------------------------------------------------------

class radTBackgroundFldSrcGraphPresent : public radTg3dGraphPresent {
public:
	radTBackgroundFldSrcGraphPresent(radTg3d* InBackgroundFldSrcPtr) : radTg3dGraphPresent(InBackgroundFldSrcPtr) {}
	radTBackgroundFldSrcGraphPresent() {}

	void Draw(radTrans* BaseTransPtr)
	{
		Send.InitOutList(0, DrawFacilityInd);
	}
};

//-------------------------------------------------------------------------

#endif
