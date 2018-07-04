/*-------------------------------------------------------------------------
*
* File name:      radsend.h
*
* Project:        RADIA
*
* Description:    Interface functions (data input / output)
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __RADSEND_H
#define __RADSEND_H

#ifdef _WITH_QD3D
#include "radq3ld.h"
#endif

#include "radauxst.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radRGB;
struct radTDrawAttr;
class radTrans;
class radTField;

//-------------------------------------------------------------------------

class radTSend {
public:

	double Limits3D[6];

	char ShowLines;
	char ShowFaces;

	radTVectGeomPolygon GeomPolygons;
	radTVectGeomPolygon GeomLines;

#ifdef _WITH_QD3D
	TQ3GroupObject QD3D_GroupToDraw;
	//char ShowLinesInQD3D;
	//char ShowFacesInQD3D;

	radTrans* pSmallRotForQD3D;
#endif

	radTSend() 
	{
		ShowLines = 0; ShowFaces = 1;
		InitLimits3D();

#ifdef _WITH_QD3D
		QD3D_GroupToDraw = 0; 
		pSmallRotForQD3D = 0;
#endif
	}

	static void ErrorMessage(const char*);
	void OrdinaryMessage(const char*);
	static void WarningMessage(const char*);

	void String(const char*);
	void ByteString(const unsigned char* MessageString, long len);
	void Vector3d(const TVector3d*);
	void Vector3d(const TVector3df*);
	void ArrayOfVector3d(const TVector3d*, int);
	void Matrix3d(const TMatrix3d*);
	void Matrix3d(const TMatrix3df*);
	void MatrixOfMatrix3d(TMatrix3d**, int, int);
	void MatrixOfMatrix3d(TMatrix3df**, int, int);
	void Long(long);
	void Int(int);
	void IntList(int*, int);
	void Double(double);
	void DoubleList(double*, int);
	void ArbNestedArrays(double*, int*, int);
	void SubArbNestedArrays(double*, int*, int, int&);

	void MultiDimArrayOfDouble(double*, int*, int);

	void ArrayOfPairOfVect3d(radTVectPairOfVect3d* pVectPairOfVect3d);
	void OutFieldForceOrTorqueThroughEnergyCompRes(char* ForceComponID, TVector3d& Vect, char ID);
	void OutFieldIntCompRes(char* FieldIntChar, radTField* FieldPtr, double* ArgArray = 0, int Np = 1);
	void OutFieldCompRes(char* FieldChar, radTField* FieldArray, double* ArgArray, int Np);
	void OutRelaxResultsInfo(double* RelaxStatusParamArray, int lenRelaxStatusParamArray, int ActualIterNum);
	void OutMagnetizCompRes(char* MagnChar, TVector3d& M_vect);

	void MyMLPutDouble(double);

	void GenInitDraw(char =0);
	void InitDrawSurfElem(int, radTDrawAttr&, int, char =0);
	void InitDrawLinElem(int, radTDrawAttr&, int, char =0);
	void InitOutList(int, char =0);
	void DrawEdgesSuppression(char =0);
	void InitDrawLineWithThickness(int, radTDrawAttr&, char =0);

	void Color(const radRGB&, char =0);
	void Polygon(const TVector3d*, int, char =0);
	void Line(const TVector3d*, int, char =0);

	//void FrameLines(char);
	void InitLimits3D()
	{
		Limits3D[0] = 1.E+23; Limits3D[1] = -1.E+23; 
		Limits3D[2] = 1.E+23; Limits3D[3] = -1.E+23; 
		Limits3D[4] = 1.E+23; Limits3D[5] = -1.E+23; 
	}
	int InitSmallRotForQD3D();
	void DelSmallRotForQD3D();

	//void DrawCharacter(char, double, TVector3d*, char);
	//void DrawPyramidArrow(TVector3d*, char);

	int GetInteger(int&);
	int GetDouble(double&);
	int GetString(const char*&);
	//int GetString(char*&);
	void DisownString(char* Str);
	int GetArbitraryListOfVector3d(radTVectorOfVector3d&, radTVectInputCell&);

	int GetVector3d(TVector3d& vect3d);
	int GetVector2d(TVector2d& vect2d);

	int GetArrayOfVector3d(TVector3d*&, int&);
	int GetArrayOfVector2d(TVector2d*&, int&);
	int GetArrayOfVector2dVersion2(TVector2d*&, int&);
	int GetArrayOfArrayOfVector3d(TVector3d**&, int*&, int&);
	int GetArrayOfArrayOfInt(int**&, int*&, int&);

	int GetArrayOfDouble(double*&, long&);

	void AddGeomPolygon(const TVector3d* Side, int lenSide, radTVectGeomPolygon& VectGeomPolygons);
	void DeallocateGeomPolygonData();
};

//-------------------------------------------------------------------------

#endif

