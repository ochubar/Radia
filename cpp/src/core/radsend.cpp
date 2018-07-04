/*-------------------------------------------------------------------------
*
* File name:      radsend.cpp
*
* Project:        RADIA
*
* Description:    Interface functions
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "radsend.h"
#include <stdio.h>
#include <string.h>

//#ifdef __JAVA__
//#ifndef __SEND2JAVA_H
//#include "Send2Java.h"
//#endif
//extern CSendToJava gSendToJava;
//#endif

//#ifdef __DLLVBA__
//#ifndef __SEND2VBA_H
//#include "Send2VBA.h"
//#endif
//extern radTSendToVBA gSendToVBA;
//#endif

#include "radiobuf.h"
extern radTIOBuffer ioBuffer;
#include "radg3dgr.h"

#ifdef __MATHEMATICA__
extern "C" {
//#include <mathlink.h>
#include "mathlink_wrap.h" //OC091015
}
#endif

#ifdef _WITH_QD3D
#ifdef _WINDOWS
extern "C" {
/*extern*/ radTFQ3Object_Dispose radQ3Object_Dispose;
/*extern*/ radTFQ3Group_AddObject radQ3Group_AddObject;
/*extern*/ radTFQ3Polygon_New radQ3Polygon_New;
/*extern*/ radTFQ3PolyLine_New radQ3PolyLine_New;
/*extern*/ radTFQ3AttributeSet_New radQ3AttributeSet_New;
/*extern*/ radTFQ3AttributeSet_Add radQ3AttributeSet_Add;

/*extern*/ radTFQ3Cone_New radQ3Cone_New;
}
#endif
#endif

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTSend::ErrorMessage(const char* MessageString)
{
#ifdef __MATHEMATICA__

	char ErrDescr[2048];
	ioBuffer.PrepErrWarnMesageForMathematica(ErrDescr, MessageString, 'e');

	char err_msg[300];
	sprintf(err_msg, "%s%s%s", "Message[",MessageString,"]");

	MLClearError(stdlink);
	MLNewPacket(stdlink);
	MLEvaluate(stdlink, ErrDescr);
	MLNextPacket(stdlink);
	MLNewPacket(stdlink);
	MLEvaluate(stdlink, err_msg);
	MLNextPacket(stdlink);
	MLNewPacket(stdlink);
	MLPutSymbol(stdlink, "$Failed");

	//char err_msg[300];
	//sprintf(err_msg, "%s%s%s", "Message[",MessageString,"]");
	//MLClearError(stdlink);
	//MLNewPacket(stdlink);
	//MLEvaluate(stdlink, err_msg);
	//MLNextPacket(stdlink);
	//MLNewPacket(stdlink);
	//MLPutSymbol(stdlink, "$Failed");
	return;
#endif
//#ifdef __JAVA__
//	gSendToJava.SendErrorMessage(MessageString);
//#endif
//#ifdef __DLLVBA__
//	gSendToVBA.SendErrorMessage(MessageString);
//#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreErrorMessage(MessageString);
#endif
}

//-------------------------------------------------------------------------

void radTSend::WarningMessage(const char* MessageString)
{
#ifdef __MATHEMATICA__
	
	char ErrDescr[2048];
	ioBuffer.PrepErrWarnMesageForMathematica(ErrDescr, MessageString, 'w');

	char err_msg[300];
	sprintf(err_msg, "%s%s%s", "Message[",MessageString,"]");

	MLNewPacket(stdlink);
	MLEvaluate(stdlink, ErrDescr);
	MLNextPacket(stdlink);
	MLNewPacket(stdlink);
	MLEvaluate(stdlink, err_msg);
	MLNextPacket(stdlink);
	MLNewPacket(stdlink);

	//char err_msg[300];
	//sprintf(err_msg, "%s%s%s", "Message[",MessageString,"]");
	//MLNewPacket(stdlink);
	//MLEvaluate(stdlink, err_msg);
	//MLNextPacket(stdlink);
	//MLNewPacket(stdlink);
	//return;
#endif
//#ifdef __JAVA__
//	gSendToJava.SendWarningMessage(MessageString);
//#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreWarningMessage(MessageString);
#endif
}

//-------------------------------------------------------------------------

void radTSend::OrdinaryMessage(const char* MessageString)
{
#ifdef __MATHEMATICA__
	char InfoMessage[500];
	sprintf(InfoMessage, "%s%s%s", "Message[",MessageString,"]");
	MLEvaluate(stdlink, InfoMessage);
	MLNextPacket(stdlink);
	MLNewPacket(stdlink);
	MLPutSymbol(stdlink, "Null");
#endif
}

//-------------------------------------------------------------------------

void radTSend::String(const char* MessageString)
{
#ifdef __MATHEMATICA__
	MLPutString(stdlink, MessageString);
#endif
#ifdef __JAVA__
	gSendToJava.SendString(MessageString);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreString(MessageString);
#endif
}

//-------------------------------------------------------------------------

void radTSend::ByteString(const unsigned char* MessageString, long len)
{
#ifdef __MATHEMATICA__
	MLPutByteString(stdlink, MessageString, len);
#endif
//#ifdef __JAVA__
//	gSendToJava.SendString(MessageString);
//#endif
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreByteString((const char*)MessageString, len);
#endif
}

//-------------------------------------------------------------------------

void radTSend::Double(double d)
{
#ifdef __MATHEMATICA__
	if(d==0.) d=1.E-17; 
	MLPutDouble(stdlink, d);
#endif
#ifdef __JAVA__
	gSendToJava.SendDouble(d);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreDouble(d);
#endif
}

//-------------------------------------------------------------------------

void radTSend::MyMLPutDouble(double d)
{
#ifdef __MATHEMATICA__
	if(d==0.) d=1.E-17; 
	MLPutDouble(stdlink, d);
#endif
}

//-------------------------------------------------------------------------

void radTSend::DoubleList(double* ArrayOfDouble, int lenArrayOfDouble)
{
#ifdef __MATHEMATICA__
	InitOutList(lenArrayOfDouble);
	for(int i=0; i<lenArrayOfDouble; i++) MyMLPutDouble(ArrayOfDouble[i]);
#endif
//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__
	int Dims[] = {lenArrayOfDouble};
	MultiDimArrayOfDouble(ArrayOfDouble, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::Long(long LongIntValue)
{
#ifdef __MATHEMATICA__
	MLPutLongInteger(stdlink, LongIntValue);
#endif
#ifdef __JAVA__
	gSendToJava.SendLong(LongIntValue);
#endif
}

//-------------------------------------------------------------------------

void radTSend::Int(int IntValue)
{
#ifdef __MATHEMATICA__
	MLPutInteger(stdlink, IntValue);
#endif
#ifdef __JAVA__
	gSendToJava.SendInt(IntValue);
#endif
#ifdef __DLLVBA__
	gSendToVBA.SendInt(IntValue);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreInt(IntValue);
#endif
}

//-------------------------------------------------------------------------

void radTSend::IntList(int* ArrayOfInt, int lenArrayOfInt)
{
#ifdef __MATHEMATICA__
	InitOutList(lenArrayOfInt);
	for(int i=0; i<lenArrayOfInt; i++) Int(ArrayOfInt[i]);
#endif
#ifdef __JAVA__
	int Dims[] = { lenArrayOfInt};
	gSendToJava.SendMultiDimArrayOfInt(ArrayOfInt, Dims, 1);
#endif
#ifdef __DLLVBA__
	int Dims[] = { lenArrayOfInt};
	gSendToVBA.SendMultiDimArrayOfInt(ArrayOfInt, Dims, 1);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	int Dims[] = { lenArrayOfInt};
	ioBuffer.StoreMultiDimArrayOfInt(ArrayOfInt, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::InitOutList(int NumberOfElem, char DrawFacilityInd)
{
#ifdef __MATHEMATICA__
	if(DrawFacilityInd == 0)
	{
		MLPutFunction(stdlink, "List", NumberOfElem);
	}
#endif
}

//-------------------------------------------------------------------------

void radTSend::Vector3d(const TVector3d* VectorPtr)
{
#ifdef __MATHEMATICA__
	MLPutFunction(stdlink, "List", 3);
		MyMLPutDouble(VectorPtr->x);
		MyMLPutDouble(VectorPtr->y);
		MyMLPutDouble(VectorPtr->z);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	double TotOutArray[] = {VectorPtr->x, VectorPtr->y, VectorPtr->z};
	int Dims[] = {3};
	MultiDimArrayOfDouble(TotOutArray, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::Vector3d(const TVector3df* VectorPtr)
{
#ifdef __MATHEMATICA__
	MLPutFunction(stdlink, "List", 3);
		MyMLPutDouble(VectorPtr->x);
		MyMLPutDouble(VectorPtr->y);
		MyMLPutDouble(VectorPtr->z);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	double TotOutArray[] = {VectorPtr->x, VectorPtr->y, VectorPtr->z};
	int Dims[] = {3};
	MultiDimArrayOfDouble(TotOutArray, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::ArrayOfVector3d(const TVector3d* ArrayOfVector3d, int lenArray)
{
#ifdef __MATHEMATICA__
	MLPutFunction(stdlink, "List", lenArray);
	for(int i=0; i<lenArray; i++) Vector3d(&(ArrayOfVector3d[i]));
#endif
}

//-------------------------------------------------------------------------

void radTSend::Matrix3d(const TMatrix3d* MatrixPtr)
{
#ifdef __MATHEMATICA__
	MLPutFunction(stdlink, "List", 3);
		Vector3d(&(MatrixPtr->Str0));
		Vector3d(&(MatrixPtr->Str1));
		Vector3d(&(MatrixPtr->Str2));
#endif
}

//-------------------------------------------------------------------------

void radTSend::Matrix3d(const TMatrix3df* MatrixPtr)
{
#ifdef __MATHEMATICA__
	MLPutFunction(stdlink, "List", 3);
		Vector3d(&(MatrixPtr->Str0));
		Vector3d(&(MatrixPtr->Str1));
		Vector3d(&(MatrixPtr->Str2));
#endif
}

//-------------------------------------------------------------------------

void radTSend::MatrixOfMatrix3d(TMatrix3d** MatrixOfMatrix3d, int AmOfStr, int AmOfCol)
{
#ifdef __MATHEMATICA__
	MLPutFunction(stdlink, "List", AmOfStr);
	for(int StrNo=0; StrNo<AmOfStr; StrNo++)
	{
		MLPutFunction(stdlink, "List", AmOfCol);
		for(int ColNo=0; ColNo<AmOfCol; ColNo++) Matrix3d(&(MatrixOfMatrix3d[StrNo][ColNo]));
	}
#endif
}

//-------------------------------------------------------------------------

void radTSend::MatrixOfMatrix3d(TMatrix3df** MatrixOfMatrix3d, int AmOfStr, int AmOfCol)
{
#ifdef __MATHEMATICA__
	MLPutFunction(stdlink, "List", AmOfStr);
	for(int StrNo=0; StrNo<AmOfStr; StrNo++)
	{
		MLPutFunction(stdlink, "List", AmOfCol);
		for(int ColNo=0; ColNo<AmOfCol; ColNo++) Matrix3d(&(MatrixOfMatrix3d[StrNo][ColNo]));
	}
#endif
}

//-------------------------------------------------------------------------

void radTSend::SubArbNestedArrays(double* Data, int* Dims, int Depth, int& CntData)
{
#ifdef __MATHEMATICA__
	for(int i=0; i<Dims[Depth-1]; i++)
	{
		if(Depth>1)
		{
			MLPutFunction(stdlink, "List", Dims[Depth-2]);
			SubArbNestedArrays(Data, Dims, Depth-1, CntData);
		}
		else
		{
			double Buf = Data[CntData];
			MyMLPutDouble(Data[CntData++]);
		}
	}
#endif
}

//-------------------------------------------------------------------------

void radTSend::ArbNestedArrays(double* Data, int* Dims, int Depth)
{
#ifdef __MATHEMATICA__
	MLPutFunction(stdlink, "List", Dims[Depth-1]);
	int CntData =0;
	SubArbNestedArrays(Data, Dims, Depth, CntData);
#endif
//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__
	MultiDimArrayOfDouble(Data, Dims, Depth);
#endif
}

//-------------------------------------------------------------------------

void radTSend::Color(const radRGB& RGB_color, char DrawFacilityInd)
{// This is not for Lines !!!
#ifdef __MATHEMATICA__
	MLPutFunction(stdlink, "SurfaceColor", 1);
		MLPutFunction(stdlink, "RGBColor", 3);
			MyMLPutDouble(RGB_color.Red);
			MyMLPutDouble(RGB_color.Green);
			MyMLPutDouble(RGB_color.Blue);
#endif
}

//-------------------------------------------------------------------------

void radTSend::GenInitDraw(char DrawFacilityInd)
{
	if(DrawFacilityInd == 0) 
	{
#ifdef __MATHEMATICA__
		MLPutFunction(stdlink, "List", 1);
#endif
	}
	else if(DrawFacilityInd == 1) 
	{
#ifdef _WITH_QD3D

#endif
	}
}

//-------------------------------------------------------------------------

void radTSend::InitDrawSurfElem(int DrawAttrAreSet, radTDrawAttr& DrawAttr, int NumberOfSymChilds_pl_Orig, char DrawFacilityInd)
{
	if(DrawFacilityInd == 0)
	{
#ifdef __MATHEMATICA__
		if(DrawAttrAreSet) 
		{
			MLPutFunction(stdlink, "List", 3);
				Color(DrawAttr.RGB_col);
				MLPutFunction(stdlink, "EdgeForm", 1);
					MLPutFunction(stdlink, "Thickness", 1);
						MyMLPutDouble(DrawAttr.LineThickness);
		}
		if(NumberOfSymChilds_pl_Orig > 1) MLPutFunction(stdlink, "List", NumberOfSymChilds_pl_Orig);
#endif
	}
	else if(DrawFacilityInd == 1) 
	{
#ifdef _WITH_QD3D

#endif
	}
}

//-------------------------------------------------------------------------

void radTSend::InitDrawLinElem(int DrawAttrAreSet, radTDrawAttr& DrawAttr, int NumberOfSymChilds_pl_Orig, char DrawFacilityInd)
{
	if(DrawFacilityInd == 0)
	{
#ifdef __MATHEMATICA__
		if(DrawAttrAreSet) 
		{
			MLPutFunction(stdlink, "List", 2);
				//MLPutFunction(stdlink, "LineForm", 2);
				MLPutFunction(stdlink, "List", 2); //OC240907
					MLPutFunction(stdlink, "Thickness", 1);
						MyMLPutDouble(DrawAttr.LineThickness);
					MLPutFunction(stdlink, "RGBColor", 3);
						MyMLPutDouble(DrawAttr.RGB_col.Red);
						MyMLPutDouble(DrawAttr.RGB_col.Green);
						MyMLPutDouble(DrawAttr.RGB_col.Blue);
		}
		if(NumberOfSymChilds_pl_Orig > 1) MLPutFunction(stdlink, "List", NumberOfSymChilds_pl_Orig);
#endif
	}
	else if(DrawFacilityInd == 1) 
	{
#ifdef _WITH_QD3D

#endif
	}
}

//-------------------------------------------------------------------------

void radTSend::InitDrawLineWithThickness(int DrawAttrAreSet, radTDrawAttr& DrawAttr, char DrawFacilityInd)
{
	if(DrawFacilityInd == 0)
	{
#ifdef __MATHEMATICA__
		if(DrawAttrAreSet) 
		{
			MLPutFunction(stdlink, "List", 2);
				//MLPutFunction(stdlink, "LineForm", 1); //OC240907
					MLPutFunction(stdlink, "Thickness", 1);
						MyMLPutDouble(DrawAttr.LineThickness);
		}
#endif
	}
	else if(DrawFacilityInd == 1)
	{
#ifdef _WITH_QD3D

#endif
	}
}

//-------------------------------------------------------------------------

void radTSend::DrawEdgesSuppression(char DrawFacilityInd)
{
	if(DrawFacilityInd == 0)
	{
#ifdef __MATHEMATICA__
		MLPutFunction(stdlink, "List", 2);
		MLPutFunction(stdlink, "EdgeForm", 0);
#endif
	}
	else if(DrawFacilityInd == 1)
	{
#ifdef _WITH_QD3D

#endif
	}
}

//-------------------------------------------------------------------------

int radTSend::InitSmallRotForQD3D()
{
#ifdef _WITH_QD3D

	pSmallRotForQD3D = 0;
	pSmallRotForQD3D = new radTrans();
	if(pSmallRotForQD3D == 0) { ErrorMessage("Radia::Error900"); return 0;}

	TVector3d PoiOnAxVect(0.,0.,0.), AxVect(1.,2.,3.);
	const double Angle = 0.0005;
	pSmallRotForQD3D->SetupRotation(PoiOnAxVect, AxVect, Angle);
	return 1;

#endif
#ifndef _WITH_QD3D

	return 1;

#endif
}

//-------------------------------------------------------------------------

void radTSend::DelSmallRotForQD3D()
{
#ifdef _WITH_QD3D

	if(pSmallRotForQD3D != 0) delete pSmallRotForQD3D;
	pSmallRotForQD3D = 0;

#endif
}

//-------------------------------------------------------------------------

void radTSend::Polygon(const TVector3d* Side, int lenSide, char DrawFacilityInd)
{
	if(DrawFacilityInd == 0)
	{
#ifdef __MATHEMATICA__
		MLPutFunction(stdlink, "Polygon", 1);
			MLPutFunction(stdlink, "List", lenSide);
				TVector3d* SideTravers = (TVector3d*)Side;
				for(int i = 0; i < lenSide; i++) Vector3d(SideTravers++);
#endif
	}
	else if(DrawFacilityInd == 1)
	{
#ifdef _WITH_QD3D

		//if(!ShowFacesInQD3D) return;
		if(!ShowFaces) return;

		TQ3Vertex3D* Vertices = new TQ3Vertex3D[lenSide];
		TQ3Point3D* Points = new TQ3Point3D[lenSide];

		TVector3d* tSide = (TVector3d*)Side;
		TQ3Point3D* tPoints = Points;
		TQ3Vertex3D* tVertices = Vertices;
		for(int i = 0; i < lenSide; i++)
		{
			TVector3d P = pSmallRotForQD3D->TrPoint(*tSide);

			if(P.x < Limits3D[0]) Limits3D[0] = P.x;
			else if(P.x > Limits3D[1]) Limits3D[1] = P.x;
			if(P.y < Limits3D[2]) Limits3D[2] = P.y;
			else if(P.y > Limits3D[3]) Limits3D[3] = P.y;
			if(P.z < Limits3D[4]) Limits3D[4] = P.z;
			else if(P.z > Limits3D[5]) Limits3D[5] = P.z;

			tPoints->x = float(P.y); tPoints->y = float(P.z); tPoints->z = float(P.x);
			tVertices->point = *tPoints; 
			tVertices->attributeSet = NULL;

			tSide++; tPoints++; tVertices++;
		}

		TQ3AttributeSet PgnAttrSet = NULL;

		TQ3ColorRGB ColorRGB;
		if(!(radTg3dGraphPresent::DrawAttrStack).empty())
		{
			radRGB& RGB = ((radTg3dGraphPresent::DrawAttrStack).begin())->RGB_col;
			ColorRGB.r = float(RGB.Red); ColorRGB.g = float(RGB.Green); ColorRGB.b = float(RGB.Blue);

			PgnAttrSet = radQ3AttributeSet_New();
			radQ3AttributeSet_Add(PgnAttrSet, kQ3AttributeTypeDiffuseColor, &ColorRGB);
		}

		TQ3PolygonData PolygonData;
		PolygonData.numVertices = lenSide;
		PolygonData.vertices = Vertices;
		PolygonData.polygonAttributeSet = PgnAttrSet;

		TQ3GeometryObject MyPolygon = radQ3Polygon_New(&PolygonData);
		TQ3GroupPosition aPosition = radQ3Group_AddObject(QD3D_GroupToDraw, MyPolygon);
		
		if(PgnAttrSet != NULL) radQ3Object_Dispose(PgnAttrSet);
		if(MyPolygon) radQ3Object_Dispose(MyPolygon);

		if(Vertices!=0) delete[] Vertices;
		if(Points!=0) delete[] Points;
#endif
	}
	else if(DrawFacilityInd == 2) // OpenGL
	{
		if(!ShowFaces) return;
		AddGeomPolygon(Side, lenSide, GeomPolygons);
	}
}

//-------------------------------------------------------------------------

void radTSend::AddGeomPolygon(const TVector3d* Side, int lenSide, radTVectGeomPolygon& VectGeomPolygons)
{
	if((Side == 0) || (lenSide == 0)) return;

	TVector3d* tSide = (TVector3d*)Side;
	double* NewVertCoords = new double[lenSide*3];
	double* tVertCoords = NewVertCoords;

	for(int i = 0; i < lenSide; i++)
	{
		*(tVertCoords++) = tSide->x;
		*(tVertCoords++) = tSide->y;
		*(tVertCoords++) = tSide->z;

		if(tSide->x < Limits3D[0]) Limits3D[0] = tSide->x;
		else if(tSide->x > Limits3D[1]) Limits3D[1] = tSide->x;
		if(tSide->y < Limits3D[2]) Limits3D[2] = tSide->y;
		else if(tSide->y > Limits3D[3]) Limits3D[3] = tSide->y;
		if(tSide->z < Limits3D[4]) Limits3D[4] = tSide->z;
		else if(tSide->z > Limits3D[5]) Limits3D[5] = tSide->z;

		tSide++;
	}

	radTGeomPolygon aNewGeomPolygon;
	aNewGeomPolygon.VertCoords = NewVertCoords;
	aNewGeomPolygon.Nv = lenSide;

	if(!(radTg3dGraphPresent::DrawAttrStack).empty())
	{
		radRGB& RGB = ((radTg3dGraphPresent::DrawAttrStack).begin())->RGB_col;
		aNewGeomPolygon.ColRGB[0] = float(RGB.Red); 
		aNewGeomPolygon.ColRGB[1] = float(RGB.Green); 
		aNewGeomPolygon.ColRGB[2] = float(RGB.Blue);
	}

	VectGeomPolygons.push_back(aNewGeomPolygon);
}

//-------------------------------------------------------------------------

void radTSend::Line(const TVector3d* EdgePoints, int lenEdgePoints, char DrawFacilityInd)
{
	if(DrawFacilityInd == 0)
	{
#ifdef __MATHEMATICA__
		MLPutFunction(stdlink, "Line", 1);
		MLPutFunction(stdlink, "List", lenEdgePoints);
		TVector3d* Travers = (TVector3d*)EdgePoints;
		for(int i = 0; i < lenEdgePoints; i++) Vector3d(Travers++);
#endif
	}
	else if(DrawFacilityInd == 1)
	{
#ifdef _WITH_QD3D

		//if(!ShowLinesInQD3D) return;
		if(!ShowLines) 
		{
			return;
		}
		
		TQ3Vertex3D* Vertices = new TQ3Vertex3D[lenEdgePoints];
		TQ3Point3D* Points = new TQ3Point3D[lenEdgePoints];

		TVector3d* tSide = (TVector3d*)EdgePoints;
		TQ3Point3D* tPoints = Points;
		TQ3Vertex3D* tVertices = Vertices;
		for(int i = 0; i < lenEdgePoints; i++)
		{
			TVector3d P = pSmallRotForQD3D->TrPoint(*tSide);

			if(P.x < Limits3D[0]) Limits3D[0] = P.x;
			else if(P.x > Limits3D[1]) Limits3D[1] = P.x;
			if(P.y < Limits3D[2]) Limits3D[2] = P.y;
			else if(P.y > Limits3D[3]) Limits3D[3] = P.y;
			if(P.z < Limits3D[4]) Limits3D[4] = P.z;
			else if(P.z > Limits3D[5]) Limits3D[5] = P.z;

			tPoints->x = float(P.y); tPoints->y = float(P.z); tPoints->z = float(P.x);

			tVertices->point = *tPoints; tVertices->attributeSet = NULL;
			tSide++; tPoints++; tVertices++;
		}

		TQ3ColorRGB LineColor = {float(0.), float(0.), float(0.)};
		TQ3AttributeSet LineAttrSet = radQ3AttributeSet_New();
		radQ3AttributeSet_Add(LineAttrSet, kQ3AttributeTypeDiffuseColor, &LineColor);

		TQ3PolyLineData PolyLineData;
		PolyLineData.numVertices = lenEdgePoints;
		PolyLineData.vertices = Vertices;
		PolyLineData.segmentAttributeSet = NULL;
		PolyLineData.polyLineAttributeSet = LineAttrSet;

		TQ3GeometryObject PolyLine = radQ3PolyLine_New(&PolyLineData);
		TQ3GroupPosition aPosition = radQ3Group_AddObject(QD3D_GroupToDraw, PolyLine);

		if(LineAttrSet) radQ3Object_Dispose(LineAttrSet);
		if(PolyLine) radQ3Object_Dispose(PolyLine);

		if(Vertices!=0) delete[] Vertices;
		if(Points!=0) delete[] Points;
#endif
	}
	else if(DrawFacilityInd == 2)
	{
		if(!ShowLines) 
		{
			return;
		}

		AddGeomPolygon(EdgePoints, lenEdgePoints, GeomLines);
	}
}

//-------------------------------------------------------------------------
/**
void radTSend::FrameLines(char DrawFacilityInd)
{//Don't use this ! Use radTg3dGraphPresent::DrawFrameLines() instead !
	const double RelFrameLinesOffset = 0.07;

	const double RelArrowHeight = 0.08;
	const double ArrowBottomRadToHeightRatio = 0.15;

	const double RelCharHeight = 0.05;
	const double CharWidthToHeigthRatio = 0.6;

	radRGB FrameLinesColor(1,1,1); //Frame lines are white
	radTDrawAttr FrameLinesDrawAttr;
	FrameLinesDrawAttr.RGB_col = FrameLinesColor;
	radTg3dGraphPresent::DrawAttrStack.push_back(FrameLinesDrawAttr);

	double OrigLimits3D[6];
	double *tOrigLimits3D = OrigLimits3D, *tLimits3D = Limits3D;
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
	Line(Side, 5, DrawFacilityInd);
	Side[0].x = OrigLimits3D[1] + Dx;
	Side[4].x = Side[3].x = Side[2].x = Side[1].x = Side[0].x;
	Line(Side, 5, DrawFacilityInd);

	Segment[1] = Side[0];
	Line(Segment, 2, DrawFacilityInd);
	Segment[0].y = Side[1].y; Segment[0].z = Side[1].z;
	Segment[1] = Side[1];
	Line(Segment, 2, DrawFacilityInd);
	Segment[0].y = Side[2].y; Segment[0].z = Side[2].z;
	Segment[1] = Side[2];
	Line(Segment, 2, DrawFacilityInd);
	Segment[0].y = Side[3].y; Segment[0].z = Side[3].z;
	Segment[1] = Side[3];
	Line(Segment, 2, DrawFacilityInd);

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

	radTg3dGraphPresent::DrawAttrStack.pop_back();
}
**/
//-------------------------------------------------------------------------
/**
void radTSend::DrawPyramidArrow(TVector3d* PyramidArrowInfo, char DrawFacilityInd)
{//Don't use this ! Use radTg3dGraphPresent::DrawPyramidArrow() instead !
// PyramidArrowInfo[0] - main vertex of pyramid
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
	Polygon(Face, 3, DrawFacilityInd);
	Face[1] = ArrowVerices[2]; Face[2] = ArrowVerices[3];
	Polygon(Face, 3, DrawFacilityInd);
	Face[1] = ArrowVerices[3]; Face[2] = ArrowVerices[4];
	Polygon(Face, 3, DrawFacilityInd);
	Face[1] = ArrowVerices[4]; Face[2] = ArrowVerices[1];
	Polygon(Face, 3, DrawFacilityInd);
	Face[0] = ArrowVerices[1]; Face[1] = ArrowVerices[4]; Face[2] = ArrowVerices[3]; Face[3] = ArrowVerices[2];
	Polygon(Face, 4, DrawFacilityInd);
}
**/
//-------------------------------------------------------------------------
/**
void radTSend::DrawCharacter(char Ch, double Ratio, TVector3d* Info3D, char DrawFacilityInd)
{//Don't use this ! Use radTg3dGraphPresent::DrawCharacter() instead !
// Info3D[0], Info3D[1] - point and normal vector of a plane in which the character lies
// Info3D[2] - height vector of a character (the Height is the length of Info3D[2])
// Width = Ratio*Height

	TVector3d& Normal = Info3D[1];
	Normal = (1./sqrt(Normal.x*Normal.x + Normal.y*Normal.y + Normal.z*Normal.z))*Normal;

	TVector3d LowerLeft = *Info3D;
	TVector3d WidthVect = Ratio*(Info3D[2]^Normal);
	TVector3d LowerRight = LowerLeft + WidthVect;
	TVector3d UpperLeft = LowerLeft + Info3D[2];
	TVector3d UpperRight = UpperLeft + WidthVect;

	TVector3d Segment[4];

	if((Ch == 'X') || (Ch == 'x'))
	{
		Segment[0] = LowerLeft; Segment[1] = UpperRight;
		Line(Segment, 2, DrawFacilityInd);
		Segment[0] = UpperLeft; Segment[1] = LowerRight;
		Line(Segment, 2, DrawFacilityInd);
	}
	else if((Ch == 'Y') || (Ch == 'y'))
	{
		TVector3d Middle = 0.25*(LowerLeft + UpperLeft + LowerRight + UpperRight);
		Segment[0] = UpperLeft; Segment[1] = Middle; Segment[2] = UpperRight;
		Line(Segment, 3, DrawFacilityInd);
		Segment[0] = 0.5*(LowerLeft + LowerRight); Segment[1] = Middle;
		Line(Segment, 2, DrawFacilityInd);
	}
	else if((Ch == 'Z') || (Ch == 'z'))
	{
		Segment[0] = UpperLeft; Segment[1] = UpperRight; Segment[2] = LowerLeft; Segment[3] = LowerRight;
		Line(Segment, 4, DrawFacilityInd);
	}
}
**/
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int radTSend::GetArrayOfDouble(double*& Data, long& lenData)
{
#ifdef __MATHEMATICA__

	long* Dims;
	char** Heads;
	long Depth;
	double* DataPtr;

	int ArrayReadOK = MLGetDoubleArray(stdlink, &DataPtr, &Dims, &Heads, &Depth);
	if(!ArrayReadOK || Depth!=1 || Dims[0]<=0) { ErrorMessage("Radia::Error000"); return 0;}

	lenData = Dims[0];
    Data = new double[lenData];
	double *tData = Data, *tDataPtr = DataPtr;
	for(long i=0; i<lenData; i++) *(tData++) = *(tDataPtr++);

	//MLDisownDoubleArray(stdlink, DataPtr, Dims, Heads, Depth);
	MLWrapDeleteDoubleArray(stdlink, DataPtr, Dims, Heads, Depth); //OC091015

#endif

	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetArrayOfVector3d(TVector3d*& ArrayOfVector3d, int& lenArrayOfVector3d)
{
#ifdef __MATHEMATICA__

	long* Dims;
	char** Heads;
	long Depth;
	double* DataPtr;

	int ArrayReadOK = MLGetDoubleArray(stdlink, &DataPtr, &Dims, &Heads, &Depth);
	if(!ArrayReadOK || Dims[1]!=3 || Depth!=2) { ErrorMessage("Radia::Error000"); return 0;}

	lenArrayOfVector3d = Dims[0];
	//try
	//{
		ArrayOfVector3d = new TVector3d[lenArrayOfVector3d];
	//}
	//catch (radTException* radExceptionPtr)
	//{
	//	ErrorMessage(radExceptionPtr->what()); return 0;
	//}
	//catch (...)
	//{
	//	ErrorMessage("Radia::Error999"); return 0;
	//}

	// Make some other check for memory allocation !!!

	for(int i=0; i<lenArrayOfVector3d; i++)
	{
		int i3 = i*3;
		ArrayOfVector3d[i] = TVector3d(DataPtr[i3], DataPtr[i3+1], DataPtr[i3+2]);
	}

	//MLDisownDoubleArray(stdlink, DataPtr, Dims, Heads, Depth);
	MLWrapDeleteDoubleArray(stdlink, DataPtr, Dims, Heads, Depth); //OC091015

#endif

	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetVector3d(TVector3d& vect3d)
{
#ifdef __MATHEMATICA__

	long* Dims;
	char** Heads;
	long Depth;
	double* DataPtr;

	int ArrayReadOK = MLGetDoubleArray(stdlink, &DataPtr, &Dims, &Heads, &Depth);
	if(!ArrayReadOK || Dims[0]!=3 || Depth!=1) { ErrorMessage("Radia::Error000"); return 0;}

	vect3d.x = *DataPtr; vect3d.y = *(DataPtr + 1); vect3d.z = *(DataPtr + 2);

	//MLDisownDoubleArray(stdlink, DataPtr, Dims, Heads, Depth);
	MLWrapDeleteDoubleArray(stdlink, DataPtr, Dims, Heads, Depth); //OC091015

#endif

	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetVector2d(TVector2d& vect2d)
{
#ifdef __MATHEMATICA__

	long* Dims;
	char** Heads;
	long Depth;
	double* DataPtr;

	int ArrayReadOK = MLGetDoubleArray(stdlink, &DataPtr, &Dims, &Heads, &Depth);
	if(!ArrayReadOK || Dims[0]!=2 || Depth!=1) { ErrorMessage("Radia::Error000"); return 0;}

	vect2d.x = *DataPtr; vect2d.y = *(DataPtr + 1);

	//MLDisownDoubleArray(stdlink, DataPtr, Dims, Heads, Depth);
	MLWrapDeleteDoubleArray(stdlink, DataPtr, Dims, Heads, Depth); //OC091015

#endif

	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetArrayOfVector2d(TVector2d*& ArrayOfVector2d, int& lenArrayOfVector2d)
{
#ifdef __MATHEMATICA__

	long* Dims;
	char** Heads;
	long Depth;
	double* DataPtr;

	int ArrayReadOK = MLGetDoubleArray(stdlink, &DataPtr, &Dims, &Heads, &Depth);
	if(!ArrayReadOK || Dims[1]!=2 || Depth!=2) { ErrorMessage("Radia::Error000"); return 0;}

	lenArrayOfVector2d = Dims[0];
	ArrayOfVector2d = new TVector2d[lenArrayOfVector2d];
	if(ArrayOfVector2d == 0) { ErrorMessage("Radia::Error900"); return 0;}

	for(int i=0; i<lenArrayOfVector2d; i++)
	{
		int i2 = i*2;
		ArrayOfVector2d[i] = TVector2d(DataPtr[i2], DataPtr[i2+1]);
	}

	//MLDisownDoubleArray(stdlink, DataPtr, Dims, Heads, Depth);
	MLWrapDeleteDoubleArray(stdlink, DataPtr, Dims, Heads, Depth); //OC091015

#endif

	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetArrayOfVector2dVersion2(TVector2d*& ArrayOfVector2d, int& lenArrayOfVector2d)
{
#ifdef __MATHEMATICA__

	const char *FunName;
	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;

	int ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	if((!ReadOK) || strcmp(FunName, "List") || (ArgNum <= 0)) { ErrorMessage("Radia::Error000"); return 0;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015

	lenArrayOfVector2d = ArgNum;
	ArrayOfVector2d = new TVector2d[lenArrayOfVector2d];
	if(ArrayOfVector2d == 0) { ErrorMessage("Radia::Error000"); return 0;}

	TVector2d* tArrayOfVector2d = ArrayOfVector2d;
	for(int i=0; i<lenArrayOfVector2d; i++)
	{
		ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		if((!ReadOK) || strcmp(FunName, "List") || (ArgNum != 2)) { ErrorMessage("Radia::Error000"); return 0;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015

		if(!MLGetDouble(stdlink, &(tArrayOfVector2d->x))) { ErrorMessage("Radia::Error000"); return 0;}
		if(!MLGetDouble(stdlink, &(tArrayOfVector2d->y))) { ErrorMessage("Radia::Error000"); return 0;}
		tArrayOfVector2d++;
	}

#endif

	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetArrayOfArrayOfVector3d(TVector3d**& ArrayOfArrayOfVector3d, int*& ArrayOfLengths, int& lenArrayOfArrayOfVector3d)
{
#ifdef __MATHEMATICA__

#if(MLVERSION >= 3)
	const char *FunName;
#else
	char *FunName;
#endif

	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;

	int Next = MLGetNext(stdlink);
	int ReadOK = 0;
	if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	else { ErrorMessage("Radia::Error000"); return 0;}
	if((!ReadOK) || strcmp(FunName, "List")) { ErrorMessage("Radia::Error000"); return 0;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015
	
	lenArrayOfArrayOfVector3d = (int)ArgNum;
	ArrayOfArrayOfVector3d = new TVector3d*[lenArrayOfArrayOfVector3d];
	ArrayOfLengths = new int[lenArrayOfArrayOfVector3d];
	if(ArrayOfArrayOfVector3d==0) { ErrorMessage("Radia::Error900"); return 0;}
	if(ArrayOfLengths==0) { ErrorMessage("Radia::Error900"); return 0;}

	for(int kk=0; kk<lenArrayOfArrayOfVector3d; kk++)
	{
		Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		else { ErrorMessage("Radia::Error000"); return 0;}
		if((!ReadOK) || strcmp(FunName, "List")) { ErrorMessage("Radia::Error000"); return 0;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015

		int CurrentLength = (int)ArgNum;
		ArrayOfLengths[kk] = CurrentLength;
		TVector3d*& CurrentArrayOfVector3d = ArrayOfArrayOfVector3d[kk];
		CurrentArrayOfVector3d = new TVector3d[CurrentLength];
		if(CurrentArrayOfVector3d==0) { ErrorMessage("Radia::Error900"); return 0;}

		for(int ii=0; ii<CurrentLength; ii++)
		{
			Next = MLGetNext(stdlink);
			if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
			else { ErrorMessage("Radia::Error000"); return 0;}
			if((!ReadOK) || strcmp(FunName, "List") || ArgNum!=3) { ErrorMessage("Radia::Error000"); return 0;}
			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC091015

			double Array3d[3];
			for(int jj=0; jj<3; jj++)
			{
				Next = MLGetNext(stdlink);
				ReadOK = 0;
				if((Next==MLTKINT) || (Next==MLTKREAL)) ReadOK = MLGetDouble(stdlink, &(Array3d[jj]));
				if(!ReadOK) { ErrorMessage("Radia::Error000"); return 0;}
			}
			TVector3d aVect3d(Array3d[0], Array3d[1], Array3d[2]);
			CurrentArrayOfVector3d[ii] = aVect3d;
		}
	}

#endif

	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetArrayOfArrayOfInt(int**& ArrayOfArrayOfInt, int*& ArrayOfLengths, int& lenArrayOfArrayOfInt)
{
#ifdef __MATHEMATICA__

#if(MLVERSION >= 3)
	const char *FunName;
#else
	char *FunName;
#endif

	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;

	int Next = MLGetNext(stdlink);
	int ReadOK = 0;
	if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	else { ErrorMessage("Radia::Error000"); return 0;}
	if((!ReadOK) || strcmp(FunName, "List")) { ErrorMessage("Radia::Error000"); return 0;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015
	
	lenArrayOfArrayOfInt = (int)ArgNum;
	ArrayOfArrayOfInt = new int*[lenArrayOfArrayOfInt];

	ArrayOfLengths = new int[lenArrayOfArrayOfInt];
	if(lenArrayOfArrayOfInt==0) { ErrorMessage("Radia::Error900"); return 0;}
	if(ArrayOfLengths==0) { ErrorMessage("Radia::Error900"); return 0;}

	for(int kk=0; kk<lenArrayOfArrayOfInt; kk++)
	{
		Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		else { ErrorMessage("Radia::Error000"); return 0;}
		if((!ReadOK) || strcmp(FunName, "List")) { ErrorMessage("Radia::Error000"); return 0;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015

		int CurrentLength = (int)ArgNum;
		ArrayOfLengths[kk] = CurrentLength;
		int*& CurrentArrayOfInt = ArrayOfArrayOfInt[kk];
		CurrentArrayOfInt = new int[CurrentLength];
		if(CurrentArrayOfInt==0) { ErrorMessage("Radia::Error900"); return 0;}

		for(int ii=0; ii<CurrentLength; ii++)
		{
			Next = MLGetNext(stdlink);
			ReadOK = 0;
			if(Next==MLTKINT) ReadOK = MLGetInteger(stdlink, &(CurrentArrayOfInt[ii]));
			if(!ReadOK) { ErrorMessage("Radia::Error000"); return 0;}
		}
	}

#endif

	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetInteger(int& Value)
{
#ifdef __MATHEMATICA__
	return MLGetInteger(stdlink, &Value);
#endif

	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetDouble(double& Value)
{
#ifdef __MATHEMATICA__
	return MLGetDouble(stdlink, &Value);
#endif

	return 1;
}

//-------------------------------------------------------------------------

int radTSend::GetString(const char*& Str)
//int radTSend::GetString(char*& Str)
{
#ifdef __MATHEMATICA__
	//const char** StrBuf=0;
	return MLGetString(stdlink, &Str);
	//int res = MLGetString(stdlink, StrBuf);
	//Str = StrBuf;
	//return res;
#endif

	return 1;
}

//-------------------------------------------------------------------------

void radTSend::DisownString(char* Str)
{
#ifdef __MATHEMATICA__
	//MLDisownString(stdlink, Str);
	MLWrapDeleteString(stdlink, Str); //OC091015
#endif

}

//-------------------------------------------------------------------------

int radTSend::GetArbitraryListOfVector3d(radTVectorOfVector3d& VectorOfVector3d, radTVectInputCell& VectInputCell)
{
	bool symbWasObtainedInsteadOfNumber = false; //OC28092010

#ifdef __MATHEMATICA__

	const char *FunName;
	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;
	int Next;

	if(!MLGetFunction(stdlink, &FunName, &ArgNum)) { ErrorMessage("Radia::Error000"); return 0;}
	if(strcmp(FunName, "List")) { ErrorMessage("Radia::Error000"); return 0;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015

	Next = MLGetNext(stdlink);
	//if((ArgNum == 3) && ((Next == MLTKINT) || (Next == MLTKREAL)))
	if((ArgNum == 3) && ((Next == MLTKINT) || (Next == MLTKREAL) || (Next == MLTKSYM)))
	{
		VectInputCell.push_back(radTInputCell('P', 3));

		TVector3d Vector3d;
		//if(!MLGetDouble(stdlink, &(Vector3d.x))) { ErrorMessage("Radia::Error000"); return 0;}
		//if(!MLGetDouble(stdlink, &(Vector3d.y))) { ErrorMessage("Radia::Error000"); return 0;}
		//if(!MLGetDouble(stdlink, &(Vector3d.z))) { ErrorMessage("Radia::Error000"); return 0;}

		//OC240907 //port to math6
		const char* AuxNameStr = 0;
		//int Next = MLGetNext(stdlink);
		if((Next == MLTKINT) || (Next == MLTKREAL))
		{
			if(!MLGetDouble(stdlink, &(Vector3d.x))) { ErrorMessage("Radia::Error000"); return 0;}
		}
		else//the following code was added to walk around probable math6 Plot[] bug, i.e. sending eventually symbolic "y" in cases like:
			//Plot[radFld[Grp,"Bz",{0,y,0}],{y,0,300}]
		{
			if(!MLGetSymbol(stdlink, &AuxNameStr)) { ErrorMessage("Radia::Error000"); return 0;}
			
			//MLDisownSymbol(stdlink, AuxNameStr);
			MLWrapDeleteSymbol(stdlink, AuxNameStr); //OC091015
			//WarningMessage("Radia::Warning017");
			
			symbWasObtainedInsteadOfNumber = true; //OC28092010
		}

		Next = MLGetNext(stdlink);
		if((Next == MLTKINT) || (Next == MLTKREAL))
		{
			if(!MLGetDouble(stdlink, &(Vector3d.y))) { ErrorMessage("Radia::Error000"); return 0;}
		}
		else
		{
			if(!MLGetSymbol(stdlink, &AuxNameStr)) { ErrorMessage("Radia::Error000"); return 0;}
			
			//MLDisownSymbol(stdlink, AuxNameStr);
			MLWrapDeleteSymbol(stdlink, AuxNameStr); //OC091015
			//WarningMessage("Radia::Warning017");
			
			symbWasObtainedInsteadOfNumber = true;
		}

		Next = MLGetNext(stdlink);
		if((Next == MLTKINT) || (Next == MLTKREAL))
		{
			if(!MLGetDouble(stdlink, &(Vector3d.z))) { ErrorMessage("Radia::Error000"); return 0;}
		}
		else
		{
			if(!MLGetSymbol(stdlink, &AuxNameStr)) { ErrorMessage("Radia::Error000"); return 0;}
			
			//MLDisownSymbol(stdlink, AuxNameStr);
			MLWrapDeleteSymbol(stdlink, AuxNameStr); //OC091015
			//WarningMessage("Radia::Warning017");
			
			symbWasObtainedInsteadOfNumber = true;
		}
		VectorOfVector3d.push_back(Vector3d);
	}
	else if(Next == MLTKFUNC)
	{
		VectInputCell.push_back(radTInputCell('L', ArgNum));
		
		for(long i=0; i<ArgNum; i++)
			if(!GetArbitraryListOfVector3d(VectorOfVector3d, VectInputCell)) return 0;
	}
	else { ErrorMessage("Radia::Error000"); return 0;}

#endif

	if(symbWasObtainedInsteadOfNumber) return 2;
	else return 1;
}

//-------------------------------------------------------------------------

void radTSend::MultiDimArrayOfDouble(double* Array, int* Dims, int NumDims)
{
#ifdef __JAVA__
	gSendToJava.SendMultiDimArrayOfDouble(Array, Dims, NumDims);
#endif
#ifdef __DLLVBA__
	gSendToVBA.SendMultiDimArrayOfDouble(Array, Dims, NumDims);
#endif
//#ifdef ALPHA__DLL__
#if defined ALPHA__DLL__ || defined ALPHA__LIB__
	ioBuffer.StoreMultiDimArrayOfDouble(Array, Dims, NumDims);
#endif
}

//-------------------------------------------------------------------------

void radTSend::ArrayOfPairOfVect3d(radTVectPairOfVect3d* pVectPairOfVect3d)
{
#ifdef __MATHEMATICA__
	int AmOfPoints = (int)pVectPairOfVect3d->size();
	InitOutList(AmOfPoints, 0);
	for(int i=0; i<AmOfPoints; i++)
	{
		InitOutList(2, 0);
		radTPairOfVect3d& Pair = (*pVectPairOfVect3d)[i];
		Vector3d(&(Pair.V1));
		Vector3d(&(Pair.V2));
	}
#endif
//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__

	int AmOfPoints = (int)pVectPairOfVect3d->size();
	int NumDims = 3;
	int Dims[] = {3,2,AmOfPoints};

	long TotLen = Dims[0]*Dims[1]*Dims[2];
	double *TotArray = new double[TotLen];
	if(TotArray == 0) { ErrorMessage("Radia::Error900"); return;}
	double *tTotArray = TotArray;
	for(int k=0; k<AmOfPoints; k++)
	{
		radTPairOfVect3d& aPair = (*pVectPairOfVect3d)[k];
		TVector3d &V1 = aPair.V1, &V2 = aPair.V2;
		*(tTotArray++) = V1.x; *(tTotArray++) = V1.y; *(tTotArray++) = V1.z;
		*(tTotArray++) = V2.x; *(tTotArray++) = V2.y; *(tTotArray++) = V2.z;
	}
	MultiDimArrayOfDouble(TotArray, Dims, NumDims);

#endif
}

//-------------------------------------------------------------------------

void radTSend::OutFieldForceOrTorqueThroughEnergyCompRes(char* ForceComponID, TVector3d& Vect, char ID)
{// This is only for Force and Torque!
	char* BufChar = ForceComponID;
	//char* EqEmptyStr = (ID=='f')? "FxFyFz" : "TxTyTz";
	//char EqEmptyStr[6];
	char EqEmptyStr[10]; //OC150505
	strcpy(EqEmptyStr, "TxTyTz");
	if(ID=='f') strcpy(EqEmptyStr, "FxFyFz");

	char SmallID = ID;
	char CapitalID = (SmallID=='f')? 'F' : 'T';

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			char* BufChar_pl_1 = BufChar+1;
			if((((*BufChar==CapitalID) || (*BufChar==SmallID)) && 
			   (*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) ||
			   (*BufChar == 'X') || (*BufChar == 'x') ||
			   (*BufChar == 'Y') || (*BufChar == 'y') ||
			   (*BufChar == 'Z') || (*BufChar == 'z')) ItemCount++;
			BufChar++;
		}
		BufChar = ForceComponID;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 3;
	}

#ifdef __MATHEMATICA__
	if(ItemCount > 1) InitOutList(ItemCount);

	while(*BufChar != '\0') 
	{
		if((*(BufChar)==CapitalID) || (*(BufChar)==SmallID))
		{
			char* BufChar_pl_1 = BufChar+1;
			if((*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) Vector3d(&Vect);
		}
		else if((*(BufChar)=='X') || (*(BufChar)=='x')) Double(Vect.x);
		else if((*(BufChar)=='Y') || (*(BufChar)=='y')) Double(Vect.y);
		else if((*(BufChar)=='Z') || (*(BufChar)=='z')) Double(Vect.z);
		BufChar++;
	}
#endif
//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__

	double TotOutArray[10];
	double *t = TotOutArray;
	int nv = 0;

	while(*BufChar != '\0') 
	{
		if((*(BufChar)==CapitalID) || (*(BufChar)==SmallID))
		{
			char* BufChar_pl_1 = BufChar+1;
			if((*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z'))
			{ *(t++) = Vect.x; *(t++) =Vect.y; *(t++) = Vect.z; nv += 3;}
		}
		else if((*(BufChar)=='X') || (*(BufChar)=='x')) { *(t++) = Vect.x; nv++;}
		else if((*(BufChar)=='Y') || (*(BufChar)=='y')) { *(t++) = Vect.y; nv++;}
		else if((*(BufChar)=='Z') || (*(BufChar)=='z')) { *(t++) = Vect.z; nv++;}
		BufChar++;
	}
	int Dims[] = { nv};
	MultiDimArrayOfDouble(TotOutArray, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::OutFieldCompRes(char* FieldChar, radTField* FieldArray, double* ArgArray, int Np)
{
	char* BufChar = FieldChar;
	//char* EqEmptyStr = "BHAM";
	char EqEmptyStr[] = "BHAM"; //OC01052013

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while(*BufChar != '\0') 
		{
			if((*BufChar == 'B') || (*BufChar == 'b') || 
			   (*BufChar == 'H') || (*BufChar == 'h') ||
			   (*BufChar == 'A') || (*BufChar == 'a') ||
			   (*BufChar == 'M') || (*BufChar == 'm') ||
			   (*BufChar == 'J') || (*BufChar == 'j') ||
			   (*BufChar == 'P') || (*BufChar == 'p')) ItemCount++;
			BufChar++;
		}
		BufChar = FieldChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 4;
	}
	char* ActualInitCharPtr = BufChar;

#ifdef __MATHEMATICA__
	if(Np > 1) InitOutList(Np);

	radTField* FieldPtr = FieldArray;
	for(int i=0; i<Np; i++)
	{
		if(ArgArray != NULL) // Argument Needed
		{
			InitOutList(2);
			Double(ArgArray[i]);
		}

		if(ItemCount > 1) InitOutList(ItemCount);
		while(*BufChar != '\0') 
		{
			char* BufChar_p_1 = BufChar+1;
			if(*(BufChar)=='B' || *(BufChar)=='b')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Double(FieldPtr->B.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Double(FieldPtr->B.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Double(FieldPtr->B.z);
				else Vector3d(&(FieldPtr->B));
			}
			else if(*(BufChar)=='H' || *(BufChar)=='h')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Double(FieldPtr->H.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Double(FieldPtr->H.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Double(FieldPtr->H.z);
				else Vector3d(&(FieldPtr->H));
			}
			else if(*(BufChar)=='A' || *(BufChar)=='a')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Double(FieldPtr->A.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Double(FieldPtr->A.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Double(FieldPtr->A.z);
				else Vector3d(&(FieldPtr->A));
			}
			else if(*(BufChar)=='M' || *(BufChar)=='m')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Double(FieldPtr->M.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Double(FieldPtr->M.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Double(FieldPtr->M.z);
				else Vector3d(&(FieldPtr->M));
			}
			else if(*(BufChar)=='J' || *(BufChar)=='j')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Double(FieldPtr->J.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Double(FieldPtr->J.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Double(FieldPtr->J.z);
				else Vector3d(&(FieldPtr->J));
			}
			else if(*(BufChar)=='P' || *(BufChar)=='p')	Double(FieldPtr->Phi);
			BufChar++;
		}
		FieldPtr++;
		BufChar = ActualInitCharPtr;
	}
#endif
//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__

	double *TotOutArray = new double[14*Np];
	if(TotOutArray == 0) { ErrorMessage("Radia::Error900"); return;}
	double *t = TotOutArray;
	int nv = 0;

	radTField* FieldPtr = FieldArray;
	for(int i=0; i<Np; i++)
	{
		nv = 0;
		if(ArgArray != NULL) // Argument Needed
		{
			*(t++) = ArgArray[i]; nv++;
		}

		while(*BufChar != '\0') 
		{
			char* BufChar_p_1 = BufChar+1;
			if(*(BufChar)=='B' || *(BufChar)=='b')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->B.x; nv++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->B.y; nv++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->B.z; nv++;}
				else { *(t++) = FieldPtr->B.x; *(t++) = FieldPtr->B.y; *(t++) = FieldPtr->B.z; nv += 3;}
			}
			else if(*(BufChar)=='H' || *(BufChar)=='h')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->H.x; nv++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->H.y; nv++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->H.z; nv++;}
				else { *(t++) = FieldPtr->H.x; *(t++) = FieldPtr->H.y; *(t++) = FieldPtr->H.z; nv += 3;}
			}
			else if(*(BufChar)=='A' || *(BufChar)=='a')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->A.x; nv++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->A.y; nv++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->A.z; nv++;}
				else { *(t++) = FieldPtr->A.x; *(t++) = FieldPtr->A.y; *(t++) = FieldPtr->A.z; nv += 3;}
			}
			else if(*(BufChar)=='M' || *(BufChar)=='m')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->M.x; nv++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->M.y; nv++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->M.z; nv++;}
				else { *(t++) = FieldPtr->M.x; *(t++) = FieldPtr->M.y; *(t++) = FieldPtr->M.z; nv += 3;}
			}
			else if(*(BufChar)=='J' || *(BufChar)=='j')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->J.x; nv++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->J.y; nv++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->J.z; nv++;}
				else { *(t++) = FieldPtr->J.x; *(t++) = FieldPtr->J.y; *(t++) = FieldPtr->J.z; nv += 3;}
			}
			else if(*(BufChar)=='P' || *(BufChar)=='p')	{ *(t++) = FieldPtr->Phi; nv++;}
			BufChar++;
		}
		FieldPtr++;
		BufChar = ActualInitCharPtr;
	}
	int Dims[] = { nv, Np};
	MultiDimArrayOfDouble(TotOutArray, Dims, 2);

	if(TotOutArray != 0) delete[] TotOutArray;
#endif
}

//-------------------------------------------------------------------------

//void radTSend::OutFieldIntCompRes(char* FieldIntChar, radTField* FieldPtr, double* ArgArray, int Np)
void radTSend::OutFieldIntCompRes(char* FieldIntChar, radTField* FieldArray, double* ArgArray, int Np)
{
	char* BufChar = FieldIntChar;
	char* BufCharPrev = NULL;
	//char* EqEmptyStr = "Ib";
	char EqEmptyStr[] = "Ib"; //OC01052013

	short I_used = 0;
	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if(((*BufChar == 'B') || (*BufChar == 'b') || 
			    (*BufChar == 'H') || (*BufChar == 'h')) ||
			   (((*BufChar == 'X') || (*BufChar == 'x') ||
			     (*BufChar == 'Y') || (*BufChar == 'y') ||
				 (*BufChar == 'Z') || (*BufChar == 'z')) &&
				(*BufCharPrev != 'B') && (*BufCharPrev != 'b') &&
				(*BufCharPrev != 'H') && (*BufCharPrev != 'h'))) ItemCount++;

			if((*BufChar == 'I') || (*BufChar == 'i')) I_used = 1;
			BufCharPrev = BufChar;
			BufChar++;
		}
		BufChar = FieldIntChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 1;
	}
	if(I_used && (ItemCount == 0))
	{
		BufChar = EqEmptyStr;
		ItemCount = 1;
	}

	char* ActualInitCharPtr = BufChar;

#ifdef __MATHEMATICA__

	if(Np > 1) InitOutList(Np);

	radTField* FieldPtr = FieldArray;
	for(int i=0; i<Np; i++)
	{
		if(ArgArray != NULL) // Argument Needed
		{
			InitOutList(2);
			Double(ArgArray[i]);
		}

		if(ItemCount > 1) InitOutList(ItemCount);
		while(*BufChar != '\0') 
		{
			char* BufChar_pl_1 = BufChar+1;
			char* BufChar_mi_1 = BufChar-1;

			if((*BufChar =='I') || (*BufChar == 'i'))
			{
				if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) Double(FieldPtr->Ib.x);
				else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) Double(FieldPtr->Ib.y);
				else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) Double(FieldPtr->Ib.z);
				else if((*BufChar_pl_1 != 'B') && (*BufChar_pl_1 != 'b') &&
					(*BufChar_pl_1 != 'H') && (*BufChar_pl_1 != 'h') &&
					(*BufChar_pl_1 != 'X') && (*BufChar_pl_1 != 'x') &&
					(*BufChar_pl_1 != 'Y') && (*BufChar_pl_1 != 'y') &&
					(*BufChar_pl_1 != 'Z') && (*BufChar_pl_1 != 'z')) { Vector3d(&(FieldPtr->Ib)); break;}
			}
			else if((*BufChar == 'B') || (*BufChar == 'b'))
			{
				if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) Double(FieldPtr->Ib.x);
				else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) Double(FieldPtr->Ib.y);
				else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) Double(FieldPtr->Ib.z);
				else Vector3d(&(FieldPtr->Ib));
			}
			else if((*BufChar == 'H') || (*BufChar == 'h'))
			{
				if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) Double(FieldPtr->Ih.x);
				else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) Double(FieldPtr->Ih.y);
				else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) Double(FieldPtr->Ih.z);
				else Vector3d(&(FieldPtr->Ih));
			}
			else if(((*BufChar == 'X') || (*BufChar == 'x')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) Double(FieldPtr->Ib.x);
			else if(((*BufChar == 'Y') || (*BufChar == 'y')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) Double(FieldPtr->Ib.y);
			else if(((*BufChar == 'Z') || (*BufChar == 'z')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) Double(FieldPtr->Ib.z);
			BufChar++;
		}
		FieldPtr++;
		BufChar = ActualInitCharPtr;
	}

#endif
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__

	//double TotOutArray[10];
	//double *t = TotOutArray;
	//int nv = 0;

	double *TotOutArray = new double[10*Np];
	if(TotOutArray == 0) { ErrorMessage("Radia::Error900"); return;}
	double *t = TotOutArray;
	int nv = 0;

	radTField* FieldPtr = FieldArray;
	for(int i=0; i<Np; i++)
	{
		nv = 0;
		if(ArgArray != NULL) // Argument Needed
		{
			*(t++) = ArgArray[i]; nv++;
		}

		while(*BufChar != '\0') 
		{
			char* BufChar_pl_1 = BufChar+1;
			char* BufChar_mi_1 = BufChar-1;

			if((*BufChar =='I') || (*BufChar == 'i'))
			{
				if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) { *(t++) = FieldPtr->Ib.x; nv++;}
				else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) { *(t++) = FieldPtr->Ib.y; nv++;}
				else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) { *(t++) = FieldPtr->Ib.z; nv++;}
				else if((*BufChar_pl_1 != 'B') && (*BufChar_pl_1 != 'b') &&
					(*BufChar_pl_1 != 'H') && (*BufChar_pl_1 != 'h') &&
					(*BufChar_pl_1 != 'X') && (*BufChar_pl_1 != 'x') &&
					(*BufChar_pl_1 != 'Y') && (*BufChar_pl_1 != 'y') &&
					(*BufChar_pl_1 != 'Z') && (*BufChar_pl_1 != 'z')) 
				{ *(t++) = FieldPtr->Ib.x; *(t++) = FieldPtr->Ib.y; *(t++) = FieldPtr->Ib.z; nv += 3; break;}
			}
			else if((*BufChar == 'B') || (*BufChar == 'b'))
			{
				if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) { *(t++) = FieldPtr->Ib.x; nv++;}
				else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) { *(t++) = FieldPtr->Ib.y; nv++;}
				else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) { *(t++) = FieldPtr->Ib.z; nv++;}
				else { *(t++) = FieldPtr->Ib.x; *(t++) = FieldPtr->Ib.y; *(t++) = FieldPtr->Ib.z; nv += 3;}
			}
			else if((*BufChar == 'H') || (*BufChar == 'h'))
			{
				if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) { *(t++) = FieldPtr->Ih.x; nv++;}
				else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) { *(t++) = FieldPtr->Ih.y; nv++;}
				else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) { *(t++) = FieldPtr->Ih.z; nv++;}
				else { *(t++) = FieldPtr->Ih.x; *(t++) = FieldPtr->Ih.y; *(t++) = FieldPtr->Ih.z; nv += 3;}
			}
			else if(((*BufChar == 'X') || (*BufChar == 'x')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) { *(t++) = FieldPtr->Ib.x; nv++;}
			else if(((*BufChar == 'Y') || (*BufChar == 'y')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) { *(t++) = FieldPtr->Ib.y; nv++;}
			else if(((*BufChar == 'Z') || (*BufChar == 'z')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) { *(t++) = FieldPtr->Ib.z; nv++;}
				BufChar++;
		}
		FieldPtr++;
		BufChar = ActualInitCharPtr;
	}

	//int Dims[] = { nv};
	//MultiDimArrayOfDouble(TotOutArray, Dims, 1);

	int Dims[] = { nv, Np};
	MultiDimArrayOfDouble(TotOutArray, Dims, 2);

	if(TotOutArray != 0) delete[] TotOutArray;
#endif
}

//-------------------------------------------------------------------------

void radTSend::OutRelaxResultsInfo(double* RelaxStatusParamArray, int lenRelaxStatusParamArray, int ActualIterNum)
{
#ifdef __MATHEMATICA__
	InitOutList(lenRelaxStatusParamArray + 1);

	double Buf;
	for(int i=0; i<lenRelaxStatusParamArray; i++) 
	{
		Buf = RelaxStatusParamArray[i];
		Double(Buf);
	}
	Int(ActualIterNum);
#endif
//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__
	int TotOutElem = lenRelaxStatusParamArray + 1;
	double *TotOutArray = new double[TotOutElem];
	if(TotOutArray == 0) { ErrorMessage("Radia::Error900"); return;}
	double *t = TotOutArray;
	double *tRelaxStatusParamArray = RelaxStatusParamArray;
	for(int i=0; i<lenRelaxStatusParamArray; i++) *(t++) = *(tRelaxStatusParamArray++);
	*t = ActualIterNum;

	int Dims[] = { TotOutElem};
	MultiDimArrayOfDouble(TotOutArray, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::OutMagnetizCompRes(char* MagnChar, TVector3d& M_vect)
{
	char* BufChar = MagnChar;
	//char* EqEmptyStr = "MxMyMz";
	char EqEmptyStr[] = "MxMyMz";

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			char* BufChar_pl_1 = BufChar+1;
			if((((*BufChar == 'M') || (*BufChar == 'm')) && 
			   (*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) ||
			   (*BufChar == 'X') || (*BufChar == 'x') ||
			   (*BufChar == 'Y') || (*BufChar == 'y') ||
			   (*BufChar == 'Z') || (*BufChar == 'z')) ItemCount++;
			BufChar++;
		}
		BufChar = MagnChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 3;
	}

#ifdef __MATHEMATICA__
	if(ItemCount > 1) InitOutList(ItemCount);

	while (*BufChar != '\0') 
	{
		if((*(BufChar)=='M') || (*(BufChar)=='m'))
		{
			char* BufChar_pl_1 = BufChar+1;
			if((*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) Vector3d(&M_vect);
		}
		else if((*(BufChar)=='X') || (*(BufChar)=='x')) Double(M_vect.x);
		else if((*(BufChar)=='Y') || (*(BufChar)=='y')) Double(M_vect.y);
		else if((*(BufChar)=='Z') || (*(BufChar)=='z')) Double(M_vect.z);

		BufChar++;
	}
#endif
//#ifdef __JAVA__
#if defined __JAVA__ || defined ALPHA__DLL__ || defined ALPHA__LIB__

	double TotOutArray[10];
	double *t = TotOutArray;
	int nv = 0;

	while (*BufChar != '\0') 
	{
		if((*(BufChar)=='M') || (*(BufChar)=='m'))
		{
			char* BufChar_pl_1 = BufChar+1;
			if((*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z'))
			{ *(t++) = M_vect.x; *(t++) =M_vect.y; *(t++) = M_vect.z; nv += 3;}
		}
		else if((*(BufChar)=='X') || (*(BufChar)=='x')) { *(t++) = M_vect.x; nv++;}
		else if((*(BufChar)=='Y') || (*(BufChar)=='y')) { *(t++) = M_vect.y; nv++;}
		else if((*(BufChar)=='Z') || (*(BufChar)=='z')) { *(t++) = M_vect.z; nv++;}
		BufChar++;
	}
	int Dims[] = { nv};
	MultiDimArrayOfDouble(TotOutArray, Dims, 1);
#endif
}

//-------------------------------------------------------------------------

void radTSend::DeallocateGeomPolygonData()
{
	int AmOfGeomPolygons = (int)GeomPolygons.size();
	int AmOfGeomLines = (int)GeomLines.size();

	if((AmOfGeomPolygons == 0) && (AmOfGeomLines == 0)) return;

	int k;
	for(k=0; k<AmOfGeomPolygons; k++)
	{
		double* pCoords = GeomPolygons[k].VertCoords;
		if(pCoords != 0) delete[] pCoords;
	}
	GeomPolygons.clear();

	for(k=0; k<AmOfGeomLines; k++)
	{
		double* pCoords = GeomLines[k].VertCoords;
		if(pCoords != 0) delete[] pCoords;
	}
	GeomLines.clear();
}

//-------------------------------------------------------------------------
