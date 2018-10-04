/*-------------------------------------------------------------------------
*
* File name:      radinter.cpp
*
* Project:        RADIA
*
* Description:    C-style interface, used by Mathematica (Mathlink)
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifdef __MATHEMATICA__
extern "C" {
//#include <mathlink.h>
#include "mathlink_wrap.h" //OC091015
}
#endif

#ifdef __MWERKS__
#ifndef WIN32
#include "radexcep.h"
#include <profiler.h>
#endif
#endif

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <string.h>

#include "radappl.h"
#include "gmvect.h"
#include "radiobuf.h"

//-------------------------------------------------------------------------

extern "C" {

void RecMag( double,double,double, double,double,double, double,double,double );
void ExtrudedPolygon();
void ExtrudedPolygon2();
void ExtrudedPolygonOpt( double, double, double**, int, double* );
void ExtrudedPolygonDLL( double, double, double*, int, char, double* );
void PlanarPolygon();
void Polyhedron1();
void PolyhedronOpt( double**, int, int**, int*, int, double* );
void PolyhedronDLL( double*, int, int*, int*, int, double*, double*, double*, double* );
void Polyhedron2();
void MultGenExtrPolygon();
void MultGenExtrPolygonOpt( double**, double*, int*, int, double* );
void MultGenExtrPolygonDLL( double*, int*, double*, int, double* );
void MultGenExtrPolygonCur();
void MultGenExtrPolygonMag();
void MultGenExtrRectangle();
void MultGenExtrRectangleOpt( double**, int, double* );
void MultGenExtrRectangleDLL( double*, double*, int, double* );
void MultGenExtrTriangle();
//void MultGenExtrTriangleDLL( double, double, double*, double*, int, char, double*, const char*,const char*,const char* );
void MultGenExtrTriangleDLL( double, double, double*, double*, int, char, double*, const char*,const char*,const char*,const char* ); //OC30072018

//void ArcMag( double,double,double, double,double, double,double, double, int, double,double,double, char* );
void ArcMag( double,double,double, double,double, double,double, double, int, char*, double,double,double );
void ArcPolygon();
void ArcPolygonDLL( double,double, char, double*, int, double,double, int, char, double,double,double );
//void CylMag( double,double,double, double, double, int, double,double,double, char* );
void CylMag( double,double,double, double, double, int, char*, double,double,double );
void RecCur( double,double,double, double,double,double, double,double,double );
void ArcCur( double,double,double, double,double, double,double, double, int, double, char*, char* );
void RaceTrack( double,double,double, double,double, double,double, double, int, double, char*, char* );
void FlmCur();
void FlmCurOpt( double**, long, double );
void FlmCurDLL( double*, int, double );
void Rectngl( double,double,double, double,double );
void Group( int*, long );
void AddToGroup(int, int*, long);
void OutGroupSize( int );
void OutGroupSubObjectKeys( int );

void DuplicateElementG3D();
void DuplicateElementG3DOpt( int, const char* );
void CreateFromG3DObjectWithSymmetries( int );
void SubdivideElementG3D();
void SubdivideElementG3DOpt( int, double*, char, double*, int, const char*, const char*, const char* );
void CutElementG3D();
void CutElementG3DOpt1( int, double,double,double, double,double,double, const char* );
void CutElementG3DOpt0( int, double,double,double, double,double,double );
void CutElementG3DOpt( int, double,double,double, double,double,double, const char* );
void SubdivideElementG3DByParPlanes();
void GeometricalVolume( int );
void GeometricalLimits( int );
void FldCmpMetForSubdRecMag( int, int, int );
void SetLocMgnInSbdRecMag();
void RecMagsAsExtrPolygons( char* );
void RecMagsAsPolyhedrons( char* );
void RecognizeRecMags( char* );
void ExtPgnsAsPolyhedrons( char* );
void NumberOfDegOfFreedom( int );
void MagnOfObj( int );
void ObjField( int, char* );
void SetObjMagn( int, double,double,double );
void BackgroundFieldSource( double,double,double );
void ScaleCurInObj( int,double );

void Translation( double,double,double );
void Rotation( double,double,double, double,double,double, double );
void PlaneSym( double,double,double, double,double,double );
void FieldInversion();
void CombineTransformLeft( int, int );
void CombineTransformRight( int, int );
void TransformObject( int, int );
void ApplySymmetry( int, int, int );

void LinearMaterial( double,double, double,double,double );
void LinearMaterial2( double,double, double );
void MaterialStd( char*, double );
void NonlinearIsotropMaterial( double,double,double, double,double,double );
void NonlinearIsotropMaterial2( double,double, double,double, double,double );
void NonlinearIsotropMaterial3();
void NonlinearIsotropMaterial3Opt( double**, long );
void NonlinearLaminatedMaterialML();
void NonlinearLaminatedMaterialFrm( double*,double*,double*, double, double* );
void NonlinearLaminatedMaterialTab( double*, int, double, double* );
void NonlinearAnisotropMaterial();
void NonlinearAnisotropMaterialOpt0( double*, int, double*, int );
void NonlinearAnisotropMaterialOpt1( double**, double** );
void NonlinearAnisotropMaterialOpt2( double**, double );
void NonlinearAnisotropMaterialOpt3( double, double** );
void ApplyMaterial( int, int );
void MvsH( int, char*, double,double,double );

void Field( int, char*, double,double,double, double,double,double, int, char*, double );
//void FieldArbitraryPointsStruct();
void FieldArbitraryPointsStruct( int, char* ); //OCTEST18042016
void FieldArbitraryPointsArray( long, const char*, double**, long );
void FieldInt( int, char*, char*, double,double,double, double,double,double );
void FieldForce( int, int );
void FieldEnergy( int, int, int,int,int );
void FieldForceThroughEnergy( int, int, char*, int,int,int );
void FieldTorqueThroughEnergy( int, int, char*, double,double,double, int,int,int );
void CompCriterium( double, double, double, double, double,double );
void CompPrecision();
void CompPrecisionOpt( const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char* );
void MultipoleThresholds(double, double, double, double); // Maybe to be removed later
void PreRelax( int, int );
void ShowInteractMatrix(int);
void ShowInteractVector(int, char*);
void ManualRelax( int, int, int, double );
//void AutoRelax( int, double, int, int );
void AutoRelax();
void AutoRelaxOpt( int, double, int, int, const char* );
void UpdateSourcesForRelax( int );
void SolveGen( int, double, int, int );
void ParticleTrajectory( int, double, double,double,double,double, double,double, int );
void FocusingPotential( int, double,double,double, double,double,double, int );
//void FocusingKickPer( int, double,double,double, double,double,double, double,int, double,double,double, double,int,double,int, const char*, int,int,double,double, const char*, double );
void FocusingKickPer( int, double,double,double, double,double,double, double,double, double,double,double, double,int,double,int, const char*, int,int,double,double, const char*, double, const char* );
void FocusingKickPerFormStrRep( double*,double*,double*,double*,double*, int,int, double, int, const char* );
//void FocusingKick( int, double,double,double, double,double,double, double*,long,int, double,double,double, double,int,double,int, const char*, double,double );
void FocusingKickML();
void ShimSignature( int, char*, double,double,double, double,double,double, double,double,double, int, double,double,double );

void TolForConvergence( double, double, double );
void RandomizationOnOrOff( char* );
void PhysicalUnits();

void GraphicsForElemWithoutSymChilds(int);
//void GraphicsForElemWithSymChilds(int);
void GraphicsForElemWithSymChildsExt();
void GraphicsForAllWithoutSymChilds();
void GraphicsForAllWithSymChilds();
void ApplyDrawAttrToElem( int, double,double,double, double );
void ApplyColorToElem( int, double,double,double );
void RemoveDrawAttrFromElem(int);
void QuickDraw3D_Viewer();
void QuickDraw3D_ViewerOpt( int, const char*, const char*, const char* );
void QuickDraw3D_ViewerOpt3( int, const char*, const char*, const char* );
void QuickDraw3D_ViewerOpt2( int, const char*, const char* );
void QuickDraw3D_ViewerOpt1( int, const char* );
void QuickDraw3D_ViewerOpt0( int );
void OpenGL_3D_Viewer();
void OpenGL_3D_ViewerOpt( int, const char*, const char*, const char* );

void DeleteElement( int );
void DeleteAllElements1();
void InterruptTime( double );
void RadiaVersion();
//void DumpElem( int );
void DumpElem();
void DumpElemOpt( int*, int, const char* );
void DumpElemParse();
void DumpElemParseOpt( const unsigned char*, int );
void GenDump();
void DeleteAllElements2();

void StartProf(int, int, int);
void StopProf();
void OutCommandsInfo();
void ReturnInput(double, int);
void MemAllocMethForIntrctMatr(char*);

//int AuxSetOptionNameAndValue(const char* OptTot, char* OptName, char** OptValue);
int AuxSetOptionNameAndValue(const char* OptTot, char* OptName, const char** OptValue);

}

//-------------------------------------------------------------------------

radTApplication rad;
radTIOBuffer ioBuffer;
radTYield radYield;
radTConvergRepair& radCR = rad.CnRep;

//-------------------------------------------------------------------------

//This is a "missing" function in VC++8 SDK;
//To allow linking DEBUG configuration with Libs compiled in release mode
//copied from Microsoft's "invarg.c"
//OC101015: Commented-out; no longer required for VC++14 (?)
//#ifdef WIN32
//#ifdef _DEBUG
//extern "C" _CRTIMP void __cdecl _invalid_parameter_noinfo(void)
//{
//	_invalid_parameter(NULL, NULL, NULL, 0, 0);
//}
//#endif
//#endif

//-------------------------------------------------------------------------

int AuxSetOptionNameAndValue(const char* OptTot, char* OptName, const char** OptValue)
{
	if(*OptTot != '\0') 
	{
		strcpy(OptName, OptTot);
		char *pEndOptName = strrchr(OptName, '-');
		if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return 0;}
		*pEndOptName = '\0';
		*OptValue = strrchr(OptTot, '>') + 1;
	}
	return 1;
}

//-------------------------------------------------------------------------

//void AuxParseOptionNamesAndValues(int Nopt, const char** NonParsedOpts, const char** OptionNames, const char** OptionValues, int& OptionCount)
void AuxParseOptionNamesAndValues(const char** NonParsedOpts, const char** OptionNames, const char** OptionValues, int& OptionCount)
{
	int Nopt = OptionCount;
	OptionCount = 0;
	for(int i=0; i<Nopt; i++)
	{
		const char* Opt = NonParsedOpts[i];
		if(Opt != 0)
		{
			if(*Opt != '\0') 
			{
				//if(!AuxSetOptionNameAndValue(Opt, (char*)(OptionNames[OptionCount]), (char**)(OptionValues + OptionCount))) return;
				if(!AuxSetOptionNameAndValue(Opt, (char*)(OptionNames[OptionCount]), OptionValues + OptionCount)) return;
				OptionCount++;
			}
		}
	}
}

//-------------------------------------------------------------------------

#ifdef __MATHEMATICA__
int AuxGetMathematicaOptions(int& Next, const char** arOptionNames, const char** arOptionValues, int& numOptions)
{
	//if((arOptionNames == 0) || (arOptionValues == 0) || (numOptions <= 0)) return 0;
	if((arOptionNames == 0) || (arOptionValues == 0) || (numOptions <= 0)) return 1; //OC240612
	const char *FunName;
	int ArgNum;

	int OptionCount = 0;
	for(int ii=0; ii<numOptions; ii++)
	{
		if(Next <= 0) Next = MLGetNext(stdlink);
		if(Next == MLTKFUNC)
		{
			int ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
			if(!ReadOK) { rad.Send.ErrorMessage("Radia::Error062"); return 0;}
			if((strcmp(FunName, "Rule") != 0) || (ArgNum != 2)) { rad.Send.ErrorMessage("Radia::Error062"); return 0;}

			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC091015
			if(!MLGetSymbol(stdlink, arOptionNames + OptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); return 0;}
			//if(!MLGetSymbol(stdlink, arOptionValues + (OptionCount++))) { rad.Send.ErrorMessage("Radia::Error062"); return 0;}
			if(!MLGetString(stdlink, arOptionValues + (OptionCount++))) { rad.Send.ErrorMessage("Radia::Error062"); return 0;}

			//Next = MLGetNext(stdlink);
			//if(Next == MLTKSYM) 
			//	if(!MLGetSymbol(stdlink, arOptionValues + (OptionCount++))) { rad.Send.ErrorMessage("Radia::Error062"); return 0;}
			//else if(Next == MLTKSTR) 
			//	if(!MLGetString(stdlink, arOptionValues + (OptionCount++))) { rad.Send.ErrorMessage("Radia::Error062"); return 0;}
			//else if((Next == MLTKREAL) || (Next == MLTKINT))
			//{
			//	double d;
			//	if(!MLGetDouble(stdlink, &d)) { rad.Send.ErrorMessage("Radia::Error062"); return 0;}

			//	std::ostringstream o;
			//	o << d;
			//	const char *sd = o.str().c_str();
			//}

			Next = 0;
		}
		//else { MLNewPacket(stdlink); break;}
		else if((Next == MLTKREAL) || (Next == MLTKINT))
		{//OC250713
			double dummy;
			if(!MLGetDouble(stdlink, &dummy)) { rad.Send.ErrorMessage("Radia::Error062"); return 0;}
			//since optional variables by default are set to 0 in .tm file

			Next = 0;
		}
		else break;
	}
	numOptions = OptionCount;
	//return (numOptions == 0)? 0 : 1;
	return 1; //OC240612
}
#endif

//-------------------------------------------------------------------------

void RecMag(double xc, double yc, double zc, 
			double Lx, double Ly, double Lz, 
			double Mx, double My, double Mz)
{
	double J[] = {0.,0.,0.};
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	//double CPoi[] = {xc, yc, zc};

	short OldActOnDoubles = radCR.ActOnDoubles;
	if(rad.TreatRecMagsAsExtrPolygons) radCR.ActOnDoubles = 0;

	double Dims[] = {radCR.Double(Lx), radCR.Double(Ly), radCR.Double(Lz)};
	//double Dims[] = {Lx, Ly, Lz};
	double Magn[] = {Mx, My, Mz};
	rad.SetRecMag(CPoi, 3, Dims, 3, Magn, 3, J, 3, 0);

	if(rad.TreatRecMagsAsExtrPolygons) radCR.ActOnDoubles = OldActOnDoubles;
}

//-------------------------------------------------------------------------

void RecCur(double xc, double yc, double zc, 
			double Lx, double Ly, double Lz, 
			double Jx, double Jy, double Jz)
{
	double Magn[] = {0.,0.,0.};
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Dims[] = {radCR.Double(Lx), radCR.Double(Ly), radCR.Double(Lz)};
	double J[] = {Jx, Jy, Jz};
	rad.SetRecMag(CPoi, 3, Dims, 3, Magn, 3, J, 3, 1);
}

//-------------------------------------------------------------------------

void SetExtrPolygFirstPoint(double xc, double Lx, TVector2d& FirstPoint2d, char a, double* FirstPoi) //OC040306
{
	if((a == 'x') || (a == 'X'))
	{
		FirstPoi[0] = xc - 0.5*Lx;
		FirstPoi[1] = FirstPoint2d.x;
		FirstPoi[2] = FirstPoint2d.y;
	}
	else if((a == 'y') || (a == 'Y'))
	{
		FirstPoi[0] = FirstPoint2d.y;
		FirstPoi[1] = xc - 0.5*Lx;
		FirstPoi[2] = FirstPoint2d.x;
	}
	else
	{
		FirstPoi[0] = FirstPoint2d.x;
		FirstPoi[1] = FirstPoint2d.y;
		FirstPoi[2] = xc - 0.5*Lx;
	}
}

//-------------------------------------------------------------------------

void ExtrudedPolygon()
{
#ifdef __MATHEMATICA__
	double xc;
	if(!MLGetReal(stdlink, &xc)) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	double Lx;
	if(!MLGetReal(stdlink, &Lx)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	int lenArrayOfPoints2d;
	TVector2d* ArrayOfPoints2d = NULL;
	if(!rad.Send.GetArrayOfVector2d(ArrayOfPoints2d, lenArrayOfPoints2d)) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	double Magn[3];
	if(!MLGetReal(stdlink, &(Magn[0]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if(!MLGetReal(stdlink, &(Magn[1]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if(!MLGetReal(stdlink, &(Magn[2]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	const char* OrientStr;
	if(!MLGetString(stdlink, &OrientStr)) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	//double FirstPoi[] = {xc - 0.5*Lx, ArrayOfPoints2d[0].x, ArrayOfPoints2d[0].y};
	double FirstPoi[3];
	SetExtrPolygFirstPoint(xc, Lx, ArrayOfPoints2d[0], OrientStr[0], FirstPoi); //OC040306

	for(int i1=0; i1<3; i1++) FirstPoi[i1] = radCR.Double(FirstPoi[i1]);

	rad.SetExtrudedPolygon(FirstPoi, 3, radCR.Double(Lx), ArrayOfPoints2d, lenArrayOfPoints2d, Magn, 3, OrientStr);

	if(ArrayOfPoints2d != NULL) delete[] ArrayOfPoints2d;
	//MLDisownString(stdlink, OrientStr); //OC251108
	MLWrapDeleteString(stdlink, OrientStr); //OC091015 //OC251108
#endif
}

//-------------------------------------------------------------------------

void ExtrudedPolygon2()
{
#ifdef __MATHEMATICA__
	double xc;
	if(!MLGetReal(stdlink, &xc)) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	double Lx;
	if(!MLGetReal(stdlink, &Lx)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	int lenArrayOfPoints2d;
	TVector2d* ArrayOfPoints2d = NULL;
	if(!rad.Send.GetArrayOfVector2d(ArrayOfPoints2d, lenArrayOfPoints2d)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	const char* OrientStr;
	if(!MLGetString(stdlink, &OrientStr)) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	double Magn[3];
	if(!MLGetReal(stdlink, &(Magn[0]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if(!MLGetReal(stdlink, &(Magn[1]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if(!MLGetReal(stdlink, &(Magn[2]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	//double FirstPoi[] = {xc - 0.5*Lx, ArrayOfPoints2d[0].x, ArrayOfPoints2d[0].y};
	double FirstPoi[3];
	SetExtrPolygFirstPoint(xc, Lx, ArrayOfPoints2d[0], OrientStr[0], FirstPoi); //OC040306

	for(int i1=0; i1<3; i1++) FirstPoi[i1] = radCR.Double(FirstPoi[i1]);

	rad.SetExtrudedPolygon(FirstPoi, 3, radCR.Double(Lx), ArrayOfPoints2d, lenArrayOfPoints2d, Magn, 3, OrientStr);

	if(ArrayOfPoints2d != NULL) delete[] ArrayOfPoints2d;
	//MLDisownString(stdlink, OrientStr); //OC251108
	MLWrapDeleteString(stdlink, OrientStr); //OC091015 //OC251108
#endif
}

//-------------------------------------------------------------------------

void MultGenExtrTriangle()
{
#ifdef __MATHEMATICA__
	double xc;
	if(!MLGetReal(stdlink, &xc)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	double Lx;
	if(!MLGetReal(stdlink, &Lx)) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	int lenArrayOfPoints2d;
	TVector2d* ArrayOfPoints2d = NULL;
	if(!rad.Send.GetArrayOfVector2d(ArrayOfPoints2d, lenArrayOfPoints2d)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if(lenArrayOfPoints2d <= 0) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	long* SubdDims;
	char** SubdHeads;
	long SubdDepth;
	double *SubdDataPtr=0, *AuxSubdDataPtr=0;
	int ArrayReadOK = MLGetDoubleArray(stdlink, &SubdDataPtr, &SubdDims, &SubdHeads, &SubdDepth);
	if(!ArrayReadOK) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if((SubdDepth < 1) || (SubdDepth > 2)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	int numSubdParam = SubdDims[0];
	if(numSubdParam <= 0) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if(numSubdParam != lenArrayOfPoints2d) { rad.Send.ErrorMessage("Radia::Error097"); return;}

//Parse subdivision parameters: {{k1,q1},{k2,q2},...} or {k1,k2,...}
	if(SubdDepth == 2) AuxSubdDataPtr = SubdDataPtr;
	else
	{
		AuxSubdDataPtr = new double[numSubdParam << 1];
		double *tAuxSubdDataPtr = AuxSubdDataPtr, *tSubdDataPtr = SubdDataPtr;
		for(int i=0; i<numSubdParam; i++)
		{
			*(tAuxSubdDataPtr++) = *(tSubdDataPtr++);
			*(tAuxSubdDataPtr++) = 1.;
		}
	}

	const char* OrientStr;
	if(!MLGetString(stdlink, &OrientStr)) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	double FirstPoi[3];
	SetExtrPolygFirstPoint(xc, Lx, ArrayOfPoints2d[0], OrientStr[0], FirstPoi); //OC040306
	for(int i1=0; i1<3; i1++) FirstPoi[i1] = radCR.Double(FirstPoi[i1]);

	double Magn[3];
	if(!MLGetReal(stdlink, &(Magn[0]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if(!MLGetReal(stdlink, &(Magn[1]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if(!MLGetReal(stdlink, &(Magn[2]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}

//Options
	const char *FunName;
	int ArgNum = 0;
	const int AmOfOptions = 4;
	const char* OptionNames[] = {0,0,0,0};
	const char* OptionValues[] = {0,0,0,0};
	int OptionCount = 0;
	for(int ii=0; ii<AmOfOptions; ii++)
	{
		int Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC)
		{
			int ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
			if((!ReadOK) || strcmp(FunName, "Rule") || (ArgNum!=2)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC091015

			if(!MLGetSymbol(stdlink, OptionNames + OptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			//if(!MLGetSymbol(stdlink, OptionValues + OptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			if(!MLGetString(stdlink, OptionValues + OptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			OptionCount++;
		}
		else { MLNewPacket(stdlink); break;}
	}

	rad.SetMultGenExtrTriangle(FirstPoi, 3, radCR.Double(Lx), ArrayOfPoints2d, lenArrayOfPoints2d, AuxSubdDataPtr, Magn, 3, OrientStr, OptionNames, OptionValues, OptionCount);

	if(ArrayOfPoints2d != NULL) delete[] ArrayOfPoints2d;
	if(SubdDepth != 2) delete[] AuxSubdDataPtr;
	//MLDisownDoubleArray(stdlink, SubdDataPtr, SubdDims, SubdHeads, SubdDepth);
	MLWrapDeleteDoubleArray(stdlink, SubdDataPtr, SubdDims, SubdHeads, SubdDepth); //OC091015

	for(int k=0; k<AmOfOptions; k++) 
	{
		//if(OptionNames[k] != 0) MLDisownSymbol(stdlink, OptionNames[k]);
		if(OptionNames[k] != 0) MLWrapDeleteSymbol(stdlink, OptionNames[k]); //OC091015
		//if(OptionValues[k] != 0) MLDisownString(stdlink, OptionValues[k]);
		if(OptionValues[k] != 0) MLWrapDeleteString(stdlink, OptionValues[k]); //OC091015
	}
	//MLDisownString(stdlink, OrientStr);
	MLWrapDeleteString(stdlink, OrientStr); //OC091015
#endif
}

//-------------------------------------------------------------------------

void ExtrudedPolygonOpt(double xc, double Lx, double** Polygon, int AmOfVertices, double* M)
{
	TVector2d* ArrayOfPoints2d = new TVector2d[AmOfVertices];
	if(ArrayOfPoints2d == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector2d* tArrayOfPoints2d = ArrayOfPoints2d;
	double** tPolygon = Polygon;

	for(int i=0; i<AmOfVertices; i++)
	{
		double* tPoint = *(tPolygon++);
		tArrayOfPoints2d->x = *(tPoint++);
		(tArrayOfPoints2d++)->y = *tPoint;
	}
	double FirstPoi[] = {xc - 0.5*Lx, ArrayOfPoints2d->x, ArrayOfPoints2d->y};

	rad.SetExtrudedPolygon(FirstPoi, 3, radCR.Double(Lx), ArrayOfPoints2d, AmOfVertices, M, 3, "x");

	if(ArrayOfPoints2d != 0) delete[] ArrayOfPoints2d;
}

//-------------------------------------------------------------------------

void ExtrudedPolygonDLL(double xc, double Lx, double* Polygon, int AmOfVertices, char a, double* M)
{
	TVector2d* ArrayOfPoints2d = new TVector2d[AmOfVertices];
	if(ArrayOfPoints2d == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector2d* tArrayOfPoints2d = ArrayOfPoints2d;
	double* tPolygon = Polygon;

	for(int i=0; i<AmOfVertices; i++)
	{
		tArrayOfPoints2d->x = *(tPolygon++);
		(tArrayOfPoints2d++)->y = *(tPolygon++);
	}

	//double FirstPoi[] = {xc - 0.5*Lx, ArrayOfPoints2d->x, ArrayOfPoints2d->y};
	double FirstPoi[3];
	SetExtrPolygFirstPoint(xc, Lx, ArrayOfPoints2d[0], a, FirstPoi); //OC040306

	for(int i1=0; i1<3; i1++) FirstPoi[i1] = radCR.Double(FirstPoi[i1]);

	rad.SetExtrudedPolygon(FirstPoi, 3, radCR.Double(Lx), ArrayOfPoints2d, AmOfVertices, M, 3, &a);

	if(ArrayOfPoints2d != 0) delete[] ArrayOfPoints2d;
}

//-------------------------------------------------------------------------

void MultGenExtrTriangleDLL(double xc, double lx, double* pFlatVert, double* pFlatSubd, int nv, char a, double* pM, const char* sOpt1, const char* sOpt2, const char* sOpt3, const char* sOpt4) //OC30072018
//void MultGenExtrTriangleDLL(double xc, double lx, double* pFlatVert, double* pFlatSubd, int nv, char a, double* pM, const char* sOpt1, const char* sOpt2, const char* sOpt3)
{
	TVector2d* ArrayOfPoints2d = new TVector2d[nv];
	if(ArrayOfPoints2d == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}

	TVector2d* tArrayOfPoints2d = ArrayOfPoints2d;
	double* tPolygon = pFlatVert;

	for(int i=0; i<nv; i++)
	{
		tArrayOfPoints2d->x = *(tPolygon++);
		(tArrayOfPoints2d++)->y = *(tPolygon++);
	}

	double FirstPoi[3];
	SetExtrPolygFirstPoint(xc, lx, ArrayOfPoints2d[0], a, FirstPoi); //OC040306
	for(int i1=0; i1<3; i1++) FirstPoi[i1] = radCR.Double(FirstPoi[i1]);

	//char CharBuf1[200], CharBuf2[200], CharBuf3[200];
	char CharBuf1[200], CharBuf2[200], CharBuf3[200], CharBuf4[200]; //OC30072018
	//const char* OptionNames[] = {CharBuf1, CharBuf2, CharBuf3};
	const char* OptionNames[] = {CharBuf1, CharBuf2, CharBuf3, CharBuf4}; //OC30072018
	//const char* OptionValues[] = {0,0,0};
	const char* OptionValues[] = {0,0,0,0}; //OC30072018
	//const char* NonParsedOpts[] = {sOpt1, sOpt2, sOpt3};
	const char* NonParsedOpts[] = {sOpt1, sOpt2, sOpt3, sOpt4}; //OC30072018
	//int OptionCount = 3;
	int OptionCount = 4; //OC30072018
	//AuxParseOptionNamesAndValues(3, NonParsedOpts, OptionNames, OptionValues, OptionCount);
	AuxParseOptionNamesAndValues(NonParsedOpts, OptionNames, OptionValues, OptionCount);
	
	rad.SetMultGenExtrTriangle(FirstPoi, 3, radCR.Double(lx), ArrayOfPoints2d, nv, pFlatSubd, pM, 3, &a, OptionNames, OptionValues, OptionCount);

	if(ArrayOfPoints2d != NULL) delete[] ArrayOfPoints2d;
}

//-------------------------------------------------------------------------

void PlanarPolygon()
{
#ifdef __MATHEMATICA__
	double zc;
	if(!MLGetReal(stdlink, &zc)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	int lenArrayOfPoints2d;
	TVector2d* ArrayOfPoints2d = NULL;
	if(!rad.Send.GetArrayOfVector2d(ArrayOfPoints2d, lenArrayOfPoints2d)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	double Magn[3];
	if(!MLGetReal(stdlink, &(Magn[0]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	if(!MLGetReal(stdlink, &(Magn[1]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	if(!MLGetReal(stdlink, &(Magn[2]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	rad.SetPlanarPolygon(radCR.Double(zc), ArrayOfPoints2d, lenArrayOfPoints2d, Magn, 3);

	if(ArrayOfPoints2d != NULL) delete[] ArrayOfPoints2d;
#endif
}

//-------------------------------------------------------------------------

void Polyhedron1()
{
#ifdef __MATHEMATICA__
	int lenArrayOfPoints;
	TVector3d* ArrayOfPoints = NULL;
	if(!rad.Send.GetArrayOfVector3d(ArrayOfPoints, lenArrayOfPoints)) return;

	int AmOfFaces;
	int** ArrayOfFaces = NULL;
	int* ArrayOfNumOfPoInFaces = NULL;
	if(!rad.Send.GetArrayOfArrayOfInt(ArrayOfFaces, ArrayOfNumOfPoInFaces, AmOfFaces)) return;

	double Magn[] = {0,0,0}, J[] = {0,0,0};
	double M_LinCoef[] = {0,0,0, 0,0,0, 0,0,0}, J_LinCoef[] = {0,0,0, 0,0,0, 0,0,0};

	//if(!MLGetReal(stdlink, &(Magn[0]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	//if(!MLGetReal(stdlink, &(Magn[1]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	//if(!MLGetReal(stdlink, &(Magn[2]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	bool MisDefined=false, JisDefined=false;
	bool M_LinCoefDefined=false, J_LinCoefDefined=false;

	double *pJ = 0;
	double *pJ_LinCoef = 0;
	double *pM_LinCoef = 0;
	double *pMorJ = 0, *pLinMorJ = 0; 
	double *tLinMorJ = 0;
	bool LocMisDefined = false;
	bool LocM_LinCoefDefined = false;
	bool LocJisDefined = false;
	bool LocJ_LinCoefDefined = false;

	const char *FunName;
	int ArgNum;
	int Next, ReadOK;

	//const int AmOfOptions = 2;
	//const char* OptionNames[] = {0, 0};
	const int AmOfOptions = 3;
	const char* OptionNames[] = {0,0,0};
	int OptionCount = 0;
	const char* RealOptionNames[] = {0,0,0};
	const char* RealOptionValues[] = {0,0,0};
	int RealOptionCount = 0;

	for(int ii=0; ii<AmOfOptions; ii++)
	{
		//int Next = MLGetNext(stdlink);
		Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC)
		{
			//int ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
			ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
			if(!ReadOK) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
			if((ii==0) && (strcmp(FunName, "List")==0) && (ArgNum==3)) //{mx,my,mz}
			{
				//MLDisownSymbol(stdlink, FunName);
				MLWrapDeleteSymbol(stdlink, FunName); //OC091015
				for(int jj=0; jj<3; jj++)
				{
					if(!MLGetDouble(stdlink, Magn + jj)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
				}
				MisDefined = true;
				OptionCount++;
			}
			else if((strcmp(FunName, "Rule")==0) && (ArgNum==2)) //M->{...} or J->{...}
			{
				//MLDisownSymbol(stdlink, FunName);
				MLWrapDeleteSymbol(stdlink, FunName); //OC091015
				ReadOK = MLGetSymbol(stdlink, OptionNames + OptionCount);
				//if((!ReadOK) || (!((strcmp(OptionNames[OptionCount], "M")==0) || (strcmp(OptionNames[OptionCount], "J")==0)))) 
				//{ rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
				if(!ReadOK) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}

				//double *pMorJ = 0, *pLinMorJ = 0; 
				pMorJ = 0; pLinMorJ = 0; 
				if(strcmp(OptionNames[OptionCount], "M")==0) 
				{
					if(MisDefined) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
					pMorJ = Magn; pLinMorJ = M_LinCoef;
					MisDefined = true;
				}
				if(strcmp(OptionNames[OptionCount], "J")==0) 
				{
					if(JisDefined) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
					pMorJ = J; pLinMorJ = J_LinCoef;
					JisDefined = true;
				}

				if((strcmp(OptionNames[OptionCount], "J") != 0) && (strcmp(OptionNames[OptionCount], "M") != 0)) 
				{
					RealOptionNames[RealOptionCount] = OptionNames[OptionCount];
					//if(!MLGetString(stdlink, RealOptionValues + RealOptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
					if(!MLGetSymbol(stdlink, RealOptionValues + RealOptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
					RealOptionCount++;
					OptionCount++;
					continue;
				}

				ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
				if((!ReadOK) || strcmp(FunName, "List") || ((ArgNum != 2) && (ArgNum != 3))) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
				//MLDisownSymbol(stdlink, FunName);
				MLWrapDeleteSymbol(stdlink, FunName); //OC091015

				if(ArgNum == 3) //M->{mx,my,mz} or J->{jx,jy,jz}
				{
					for(int jj=0; jj<3; jj++)
					{
						if(!MLGetDouble(stdlink, pMorJ + jj)) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
					}
				}
				else if(ArgNum == 2) //M->{{mx,my,mz},{{mlxx,mlxy,mlxz},{mlyx,mlyy,mlyz},{mlzx,mlzy,mlzz}}}
				{ //or J->{{jx,jy,jz},{{jlxx,jlxy,jlxz},{jlyx,jlyy,jlyz},{jlzx,jlzy,jlzz}}}
					ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
					if((!ReadOK) || strcmp(FunName, "List") || (ArgNum != 3)) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
					//MLDisownSymbol(stdlink, FunName);
					MLWrapDeleteSymbol(stdlink, FunName); //OC091015
					for(int jj=0; jj<3; jj++)
					{
						if(!MLGetDouble(stdlink, pMorJ + jj)) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
					}

					ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
					if((!ReadOK) || strcmp(FunName, "List") || (ArgNum != 3)) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
					//MLDisownSymbol(stdlink, FunName);
					MLWrapDeleteSymbol(stdlink, FunName); //OC091015

					//double *tLinMorJ = pLinMorJ;
					tLinMorJ = pLinMorJ;
					for(int j=0; j<3; j++)
					{
						ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
						if((!ReadOK) || strcmp(FunName, "List") || (ArgNum != 3)) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
						//MLDisownSymbol(stdlink, FunName);
						MLWrapDeleteSymbol(stdlink, FunName); //OC091015

						for(int jj=0; jj<3; jj++)
						{
							if(!MLGetDouble(stdlink, tLinMorJ++)) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
						}
					}
					M_LinCoefDefined = MisDefined;
					J_LinCoefDefined = JisDefined;
				}
				OptionCount++;
			}
			else { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
		}
		else { MLNewPacket(stdlink); break;}
	}

	if(MisDefined)
	{
		//bool LocMisDefined = false;
		LocMisDefined = false;
		for(int kk=0; kk<3; kk++) { if(Magn[kk] != 0.) { LocMisDefined = true; break;}}
		MisDefined = LocMisDefined;

		if(M_LinCoefDefined)
		{
			//bool LocM_LinCoefDefined = false;
			LocM_LinCoefDefined = false;
			for(int kk=0; kk<9; kk++) { if(M_LinCoef[kk] != 0.) { LocM_LinCoefDefined = true; break;}}
			M_LinCoefDefined = LocM_LinCoefDefined;
		}
	}
	if(JisDefined)
	{
		//bool LocJisDefined = false;
		LocJisDefined = false;
		for(int kk=0; kk<3; kk++) { if(J[kk] != 0.) { LocJisDefined = true; break;}}
		JisDefined = LocJisDefined;

		if(J_LinCoefDefined)
		{
			//bool LocJ_LinCoefDefined = false;
			LocJ_LinCoefDefined = false;
			for(int kk=0; kk<9; kk++) { if(J_LinCoef[kk] != 0.) { LocJ_LinCoefDefined = true; break;}}
			J_LinCoefDefined = LocJ_LinCoefDefined;
		}
	}

	if(MisDefined && JisDefined) { rad.Send.ErrorMessage("Radia::Error120"); goto Finish;}
	if(M_LinCoefDefined) { rad.Send.ErrorMessage("Radia::Error121"); goto Finish;}

	//double *pJ = JisDefined? J : 0;
	//double *pJ_LinCoef = J_LinCoefDefined? J_LinCoef : 0;
	//double *pM_LinCoef = M_LinCoefDefined? M_LinCoef : 0;
	if(JisDefined) pJ = J; //OC260712 (fix of jump error on GCC)
	if(J_LinCoefDefined) pJ_LinCoef = J_LinCoef;
	if(M_LinCoefDefined) pM_LinCoef = M_LinCoef;

	rad.SetPolyhedron1(ArrayOfPoints, lenArrayOfPoints, ArrayOfFaces, ArrayOfNumOfPoInFaces, AmOfFaces, Magn, 0, pJ, pJ_LinCoef, RealOptionNames, RealOptionValues, RealOptionCount);
	//rad.SetPolyhedron1(ArrayOfPoints, lenArrayOfPoints, ArrayOfFaces, ArrayOfNumOfPoInFaces, AmOfFaces, Magn, pM_LinCoef, pJ, pJ_LinCoef);

Finish:
	if(ArrayOfPoints != 0) delete[] ArrayOfPoints;
	if(ArrayOfFaces != 0)
	{
		for(int k=0; k<AmOfFaces; k++) delete[] (ArrayOfFaces[k]);
		delete[] ArrayOfFaces;
	}
	if(ArrayOfNumOfPoInFaces != 0) delete[] ArrayOfNumOfPoInFaces;

	for(int k=0; k<AmOfOptions; k++) 
	{
		//if(OptionNames[k] != 0) MLDisownSymbol(stdlink, OptionNames[k]);
		if(OptionNames[k] != 0) MLWrapDeleteSymbol(stdlink, OptionNames[k]); //OC091015
		
		//if(RealOptionNames[k] != 0) MLDisownSymbol(stdlink, RealOptionNames[k]);
		//if(RealOptionValues[k] != 0) MLDisownString(stdlink, RealOptionValues[k]);
		
		//if(RealOptionValues[k] != 0) MLDisownSymbol(stdlink, RealOptionValues[k]);
		if(RealOptionValues[k] != 0) MLWrapDeleteSymbol(stdlink, RealOptionValues[k]); //OC091015
	}

#endif
}

//-------------------------------------------------------------------------

void PolyhedronOpt(double** Vertices, int AmOfVertices, int** Faces, int* AmOfPoInFaces, int AmOfFaces, double* M)
{
	TVector3d* ArrayOfPoints = new TVector3d[AmOfVertices];
	if(ArrayOfPoints == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector3d* tArrayOfPoints = ArrayOfPoints;
	double** tVertices = Vertices;
	for(int i=0; i<AmOfVertices; i++)
	{
		double* tPoint = *(tVertices++);
		tArrayOfPoints->x = *(tPoint++);
		tArrayOfPoints->y = *(tPoint++);
		(tArrayOfPoints++)->z = *(tPoint++);
	}

	rad.SetPolyhedron1(ArrayOfPoints, AmOfVertices, Faces, AmOfPoInFaces, AmOfFaces, M, 0, 0, 0);
	//rad.SetPolyhedron1(ArrayOfPoints, AmOfVertices, Faces, AmOfPoInFaces, AmOfFaces, M, 3);

	if(ArrayOfPoints != 0) delete[] ArrayOfPoints;
}

//-------------------------------------------------------------------------

void PolyhedronDLL(double* Vertices, int AmOfVertices, int* InFaces, int* AmOfPoInFaces, int AmOfFaces, double* M, double* M_LinCoef, double* J, double* J_LinCoef)
{
	bool MisDefined = false, M_LinCoefDefined = false; //declare at te top because of use of goto
	bool JisDefined = false, J_LinCoefDefined = false;
	bool LocMisDefined = false, LocM_LinCoefDefined = false;
	bool LocJisDefined = false, LocJ_LinCoefDefined = false;

	double *tVertices=0;
	TVector3d *tArrayOfPoints=0;
	int **Faces=0, **tFaces=0, *tInFaces=0;
	double *pJ=0, *pJ_LinCoef=0;

	TVector3d* ArrayOfPoints = new TVector3d[AmOfVertices];
	if(ArrayOfPoints == 0) { rad.Send.ErrorMessage("Radia::Error900"); goto Finish;}
	//TVector3d* tArrayOfPoints = ArrayOfPoints;
	tArrayOfPoints = ArrayOfPoints;
	//double* tVertices = Vertices;
	tVertices = Vertices;
	for(int i=0; i<AmOfVertices; i++)
	{
		tArrayOfPoints->x = *(tVertices++);
		tArrayOfPoints->y = *(tVertices++);
		(tArrayOfPoints++)->z = *(tVertices++);
	}

	//int** Faces = new int*[AmOfFaces];
	Faces = new int*[AmOfFaces];
	if(Faces == 0) { rad.Send.ErrorMessage("Radia::Error900"); goto Finish;}
	//int** tFaces = Faces;
	tFaces = Faces;
	//int* tInFaces = InFaces;
	tInFaces = InFaces;

	for(int k=0; k<AmOfFaces; k++)
	{
		*(tFaces++) = tInFaces;
		tInFaces += AmOfPoInFaces[k];
	}

	if(M != 0)
	{
		//bool LocMisDefined = false;
		LocMisDefined = false;
		for(int kk=0; kk<3; kk++) { if(M[kk] != 0.) { LocMisDefined = true; break;}}
		MisDefined = LocMisDefined;
	}
	if(M_LinCoef != 0)
	{
		//bool LocM_LinCoefDefined = false;
		LocM_LinCoefDefined = false;
		for(int kk=0; kk<9; kk++) { if(M_LinCoef[kk] != 0.) { LocM_LinCoefDefined = true; break;}}
		M_LinCoefDefined = LocM_LinCoefDefined;
	}

	if(J != 0)
	{
		//bool LocJisDefined = false;
		LocJisDefined = false;
		for(int kk=0; kk<3; kk++) { if(J[kk] != 0.) { LocJisDefined = true; break;}}
		JisDefined = LocJisDefined;
	}
	if(J_LinCoef != 0)
	{
		//bool LocJ_LinCoefDefined = false;
		LocJ_LinCoefDefined = false;
		for(int kk=0; kk<9; kk++) { if(J_LinCoef[kk] != 0.) { LocJ_LinCoefDefined = true; break;}}
		J_LinCoefDefined = LocJ_LinCoefDefined;
	}
	
	if(MisDefined && JisDefined) { rad.Send.ErrorMessage("Radia::Error120"); goto Finish;}
	if(M_LinCoefDefined) { rad.Send.ErrorMessage("Radia::Error121"); goto Finish;}

	//double *pJ = JisDefined? J : 0;
	//double *pJ_LinCoef = J_LinCoefDefined? J_LinCoef : 0;
	if(JisDefined) pJ = J;
	if(J_LinCoefDefined) pJ_LinCoef = J_LinCoef;

	rad.SetPolyhedron1(ArrayOfPoints, AmOfVertices, Faces, AmOfPoInFaces, AmOfFaces, M, 0, pJ, pJ_LinCoef);

Finish:
	if(ArrayOfPoints != 0) delete[] ArrayOfPoints;
	if(Faces != 0) delete[] Faces;
}

//-------------------------------------------------------------------------

void Polyhedron2()
{
#ifdef __MATHEMATICA__
	int lenArrayOfFaces;
	TVector3d** ArrayOfFaces = NULL;
	int* ArrayOfNumOfPoInFaces = NULL;
	if(!rad.Send.GetArrayOfArrayOfVector3d(ArrayOfFaces, ArrayOfNumOfPoInFaces, lenArrayOfFaces)) return;

	double Magn[3];
	if(!MLGetReal(stdlink, &(Magn[0]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	if(!MLGetReal(stdlink, &(Magn[1]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	if(!MLGetReal(stdlink, &(Magn[2]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	rad.SetPolyhedron2(ArrayOfFaces, ArrayOfNumOfPoInFaces, lenArrayOfFaces, Magn, 3);

	if(ArrayOfFaces != 0)
	{
		for(int k=0; k<lenArrayOfFaces; k++) delete[] (ArrayOfFaces[k]);
		delete[] ArrayOfFaces;
	}
	if(ArrayOfNumOfPoInFaces != 0) delete[] ArrayOfNumOfPoInFaces;
#endif
}

//-------------------------------------------------------------------------

void RecMagsAsExtrPolygons(char* OnOrOff)
{
	rad.RecMagsAsExtrPolygons(OnOrOff);
}

//-------------------------------------------------------------------------

void RecMagsAsPolyhedrons(char* OnOrOff)
{
	rad.RecMagsAsPolyhedrons(OnOrOff);
}

//-------------------------------------------------------------------------

void RecognizeRecMags(char* OnOrOff)
{
	rad.RecognizeRecMags(OnOrOff);
}

//-------------------------------------------------------------------------

void ExtPgnsAsPolyhedrons(char* OnOrOff)
{
	rad.ExtPgnsAsPolyhedrons(OnOrOff);
}

//-------------------------------------------------------------------------

void MultGenExtrPolygon()
{
#ifdef __MATHEMATICA__
	const char *FunName;
	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;

	int ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	if((!ReadOK) || strcmp(FunName, "List") || (ArgNum <= 0)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015

	int AmOfLayerPolygons = ArgNum;
	TVector2d** LayerPolygons = new TVector2d*[AmOfLayerPolygons];
	if(LayerPolygons == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}

	int* PtsNumbersInLayerPgns = new int[AmOfLayerPolygons];
	if(PtsNumbersInLayerPgns == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}

	double* CoordsZ = new double[AmOfLayerPolygons];
	if(CoordsZ == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}

	for(int i=0; i<AmOfLayerPolygons; i++)
	{
		ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		if((!ReadOK) || strcmp(FunName, "List") || (ArgNum != 2)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015

		//if(!rad.Send.GetArrayOfVector2d(LayerPolygons[i], PtsNumbersInLayerPgns[i])) return;
		if(!rad.Send.GetArrayOfVector2dVersion2(LayerPolygons[i], PtsNumbersInLayerPgns[i])) return;
		if(!MLGetDouble(stdlink, &(CoordsZ[i]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}

		CoordsZ[i] = radCR.Double(CoordsZ[i]);
	}

	double Magn[3];
	if(!MLGetReal(stdlink, &(Magn[0]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	if(!MLGetReal(stdlink, &(Magn[1]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	if(!MLGetReal(stdlink, &(Magn[2]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	rad.SetMultGenExtrPolygon(LayerPolygons, PtsNumbersInLayerPgns, CoordsZ, AmOfLayerPolygons, Magn, 3);

	if(LayerPolygons != 0)
	{
		for(int k=0; k<AmOfLayerPolygons; k++) delete[] (LayerPolygons[k]);
		delete[] LayerPolygons;
	}
	if(PtsNumbersInLayerPgns != 0) delete[] PtsNumbersInLayerPgns;
	if(CoordsZ != 0) delete[] CoordsZ;
#endif
}

//-------------------------------------------------------------------------

void MultGenExtrPolygonOpt(double** Layers, double* Heights, int* AmOfPointsInLayers, int AmOfLayers, double* M)
{
	TVector2d** LayerPolygons = new TVector2d*[AmOfLayers];
	if(LayerPolygons == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector2d** tLayerPolygons = LayerPolygons;
	double** tLayers = Layers;

	for(int i=0; i<AmOfLayers; i++)
	{
		int AmOfPointsInTheLayer = AmOfPointsInLayers[i];
		*tLayerPolygons = new TVector2d[AmOfPointsInTheLayer];
		if(*tLayerPolygons == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}

		TVector2d* tPoints = *(tLayerPolygons++);
		double* tLayer = *(tLayers++);
		for(int k=0; k<AmOfPointsInTheLayer; k++)
		{
			tPoints->x = *(tLayer++); (tPoints++)->y = *(tLayer++); 
		}
	}
	rad.SetMultGenExtrPolygon(LayerPolygons, AmOfPointsInLayers, Heights, AmOfLayers, M, 3);

	if(LayerPolygons != 0)
	{
		for(int k=0; k<AmOfLayers; k++) delete[] (LayerPolygons[k]);
		delete[] LayerPolygons;
	}
}

//-------------------------------------------------------------------------

void MultGenExtrPolygonDLL(double* Layers, int* AmOfPointsInLayers, double* Heights, int AmOfLayers, double* M)
{
	TVector2d** LayerPolygons = new TVector2d*[AmOfLayers];
	if(LayerPolygons == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector2d** tLayerPolygons = LayerPolygons;

	double* tLayer = Layers;
	for(int i=0; i<AmOfLayers; i++)
	{
		int AmOfPointsInTheLayer = AmOfPointsInLayers[i];
		*tLayerPolygons = new TVector2d[AmOfPointsInTheLayer];
		if(*tLayerPolygons == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}

		TVector2d* tPoints = *(tLayerPolygons++);
		for(int k=0; k<AmOfPointsInTheLayer; k++)
		{
			tPoints->x = *(tLayer++); (tPoints++)->y = *(tLayer++); 
		}
		Heights[i] = radCR.Double(Heights[i]);
	}

	rad.SetMultGenExtrPolygon(LayerPolygons, AmOfPointsInLayers, Heights, AmOfLayers, M, 3);
	if(LayerPolygons != 0)
	{
		for(int k=0; k<AmOfLayers; k++) delete[] (LayerPolygons[k]);
		delete[] LayerPolygons;
	}
}

//-------------------------------------------------------------------------

#ifdef __MATHEMATICA__
int ReadRotParamsFromML(double* arRotParam)
{//Tries to read {x,y,z},{vx,vy,vz},phi from ML, without the external wrapping "{}"

	if(arRotParam == 0) return 1;
	const char *FunName;
	int NumPar = 0;
	
	int ReadOK = MLGetFunction(stdlink, &FunName, &NumPar);
	if((!ReadOK) || (strcmp(FunName, "List") != 0) || (NumPar != 3)) { rad.Send.ErrorMessage("Radia::Error000"); return 0;}
	for(int i=0; i<3; i++)
	{
		if(!MLGetReal(stdlink, arRotParam++)) { rad.Send.ErrorMessage("Radia::Error000"); return -1;}
	}

	ReadOK = MLGetFunction(stdlink, &FunName, &NumPar);
	if((!ReadOK) || (strcmp(FunName, "List") != 0) || (NumPar != 3)) { rad.Send.ErrorMessage("Radia::Error000"); return 0;}
	for(int i=0; i<3; i++)
	{
		if(!MLGetReal(stdlink, arRotParam++)) { rad.Send.ErrorMessage("Radia::Error000"); return 0;}
	}

	if(!MLGetReal(stdlink, arRotParam)) { rad.Send.ErrorMessage("Radia::Error000"); return 0;}

	return 1;
}
#endif

//-------------------------------------------------------------------------

#ifdef __MATHEMATICA__
int ReadParamsOfSeveralTrfsFromML(double **&pTrfParInExtrStep, int& numTrf, char *&strTrfOrderInExtrStep, bool readDeepNested = true, int next = 0)
{
	const char *funName;
	int numPar01=0, numPar02=0, readOK=0;

	double *pCurTrfPar=0;
	double arBuf1[3], arBuf2[3];
	double phi=0;
	double *pCurTrfPar_0=0, *pCurTrfPar_1=0;
	
	int dummyNumTrf=0;
	double **pTrfParInExtrStepAux=0;
	char *strTrfOrderAux=0;
	bool noOtherTrfExpected = true;

	//if((numTrf > 0) && (pTrfParInExtrStep != 0))
	if(numTrf > 0)
	{
		numPar01 = numTrf;
	}
	else
	{
		readOK = MLGetFunction(stdlink, &funName, &numPar01);
		if((!readOK) || (strcmp(funName, "List") != 0) || (numPar01 <= 0)) return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
		//MLDisownSymbol(stdlink, funName);
		MLWrapDeleteSymbol(stdlink, funName); //OC091015

		numTrf = 0;
		//next = MLGetNext(stdlink);
	}
	//bool readFinished = true;
	//for(int i=0; i<numPar01; i++)
	//{
	if(next == 0) next = MLGetNext(stdlink);

	if(((next==MLTKINT) || (next==MLTKREAL)) && (numPar01 == 3))
	{//assuming one translation: {x,y,z}
		if(pTrfParInExtrStep == 0) pTrfParInExtrStep = new double*[1]; 
		*pTrfParInExtrStep = new double[3];
		numTrf = 1;
		if(strTrfOrderInExtrStep == 0) strTrfOrderInExtrStep = new char[2];
		strTrfOrderInExtrStep[0] = 't';
		strTrfOrderInExtrStep[1] = '\0';

		//double *pCurTrfPar = *pTrfParInExtrStep;
		pCurTrfPar = *pTrfParInExtrStep;
		for(int ii=0; ii<3; ii++)
		{
			if(!MLGetReal(stdlink, pCurTrfPar + ii)) return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
		}
		//readFinished = true;
		//break;
	}
	else if(((next==MLTKINT) || (next==MLTKREAL)) && (numPar01 == 2))
	{//assuming one homothety: {kxH,kyH}
		if(pTrfParInExtrStep == 0) pTrfParInExtrStep = new double*[1]; 
		*pTrfParInExtrStep = new double[2];
		numTrf = 1;
		if(strTrfOrderInExtrStep == 0) strTrfOrderInExtrStep = new char[2];
		strTrfOrderInExtrStep[0] = 'h';
		strTrfOrderInExtrStep[1] = '\0';

		//double *pCurTrfPar = *pTrfParInExtrStep;
		pCurTrfPar = *pTrfParInExtrStep;
		for(int ii=0; ii<2; ii++)
		{
			if(!MLGetReal(stdlink, pCurTrfPar + ii)) return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
		}
		//readFinished = true;
		//break;
	}
	else if(next == MLTKFUNC)
	{//assuming Rotations or Translation or Homothety
		readOK = MLGetFunction(stdlink, &funName, &numPar02);
		if((!readOK) || (strcmp(funName, "List") != 0) || (numPar02 < 0)) return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
		//MLDisownSymbol(stdlink, funName);
		MLWrapDeleteSymbol(stdlink, funName); //OC091015

		//assuming {{...},...}
		next = MLGetNext(stdlink);
		if(((next==MLTKINT) || (next==MLTKREAL)) && (numPar02 == 3) && (numPar01 == 3))
		{//assuming one rotation {{x,y,z},{vx,vy,vz},phi} or one translation and two other trfs {{x,y,z},trf2,trf3}
			//double arBuf1[3], arBuf2[3];
			for(int ii=0; ii<3; ii++)
			{
				if(!MLGetReal(stdlink, arBuf1 + ii)) return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
			}

			readOK = MLGetFunction(stdlink, &funName, &numPar02);
			if((!readOK) || (strcmp(funName, "List") != 0) || (numPar02 <= 0)) return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
			//MLDisownSymbol(stdlink, funName);
			MLWrapDeleteSymbol(stdlink, funName); //OC091015

			next = MLGetNext(stdlink);
			if(((next==MLTKINT) || (next==MLTKREAL)) && (numPar02 == 3))
			{
				for(int ii=0; ii<3; ii++)
				{
					if(!MLGetReal(stdlink, arBuf2 + ii)) return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
				}
				next = MLGetNext(stdlink);
				if((next==MLTKINT) || (next==MLTKREAL))
				{//clearly one rotation
					//double phi=0;
					phi=0;
					if(!MLGetReal(stdlink, &phi)) return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}

					if(pTrfParInExtrStep == 0) pTrfParInExtrStep = new double*[1]; 
					*pTrfParInExtrStep = new double[7];
					numTrf = 1;
					if(strTrfOrderInExtrStep == 0) strTrfOrderInExtrStep = new char[2];
					strTrfOrderInExtrStep[0] = 'r';
					strTrfOrderInExtrStep[1] = '\0';

					//double *pCurTrfPar = *pTrfParInExtrStep;
					pCurTrfPar = *pTrfParInExtrStep;
					for(int ii=0; ii<3; ii++)
					{
						pCurTrfPar[ii] = arBuf1[ii];
						pCurTrfPar[ii + 3] = arBuf2[ii];
					}
					pCurTrfPar[6] = phi;
					//readFinished = true;
					//break;
				}
				else if((next == MLTKFUNC) && readDeepNested)
				{//treat already read data as two translations, and try to read another transformation

					if(pTrfParInExtrStep == 0) pTrfParInExtrStep = new double*[3]; 
					numTrf = 3;
					pTrfParInExtrStep[0] = new double[3];
					pTrfParInExtrStep[1] = new double[3];
					pTrfParInExtrStep[2] = 0;
					if(strTrfOrderInExtrStep == 0) strTrfOrderInExtrStep = new char[4];
					strTrfOrderInExtrStep[0] = 't';
					strTrfOrderInExtrStep[1] = 't';
					strTrfOrderInExtrStep[2] = '\0'; //still to define
					strTrfOrderInExtrStep[3] = '\0';

					//double *pCurTrfPar_0 = pTrfParInExtrStep[0];
					//double *pCurTrfPar_1 = pTrfParInExtrStep[1];
					pCurTrfPar_0 = pTrfParInExtrStep[0];
					pCurTrfPar_1 = pTrfParInExtrStep[1];

					for(int ii=0; ii<3; ii++)
					{
						pCurTrfPar_0[ii] = arBuf1[ii];
						pCurTrfPar_1[ii] = arBuf2[ii];
					}

					//int dummyNumTrf = 0;
					//double **pTrfParInExtrStepAux = pTrfParInExtrStep + 2;
					//char *strTrfOrderAux = strTrfOrderInExtrStep + 2;
					dummyNumTrf = 0;
					pTrfParInExtrStepAux = pTrfParInExtrStep + 2;
					strTrfOrderAux = strTrfOrderInExtrStep + 2;
					//try to read-in one more transformation parameters (using a separate function):
					if(!ReadParamsOfSeveralTrfsFromML(pTrfParInExtrStepAux, dummyNumTrf, strTrfOrderAux, false)) return 0; //nested call, but not deep
					//{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
					//readFinished = true;
					//break;
				}
				else return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
			}
			else if((next == MLTKFUNC) && readDeepNested)
			{//fill-in one translation and try to read two other transformations

				if(pTrfParInExtrStep == 0) pTrfParInExtrStep = new double*[3]; 
				numTrf = 3;
				pTrfParInExtrStep[0] = new double[3];
				pTrfParInExtrStep[1] = 0;
				pTrfParInExtrStep[2] = 0;
				if(strTrfOrderInExtrStep == 0) strTrfOrderInExtrStep = new char[4];
				strTrfOrderInExtrStep[0] = 't';
				strTrfOrderInExtrStep[1] = '\0'; //still to define
				strTrfOrderInExtrStep[2] = '\0'; //still to define
				strTrfOrderInExtrStep[3] = '\0';

				//double *pCurTrfPar_0 = pTrfParInExtrStep[0];
				pCurTrfPar_0 = pTrfParInExtrStep[0];
				for(int ii=0; ii<3; ii++)
				{
					pCurTrfPar_0[ii] = arBuf1[ii];
				}

				for(int jj=0; jj<2; jj++)
				{
					//double **pTrfParInExtrStepAux = pTrfParInExtrStep + 1 + jj;
					//char *strTrfOrderAux = strTrfOrderInExtrStep + 1 + jj;
					//int dummyNumTrf = 0;
					pTrfParInExtrStepAux = pTrfParInExtrStep + 1 + jj;
					strTrfOrderAux = strTrfOrderInExtrStep + 1 + jj;
					dummyNumTrf = 0;

					//try to read-in one more transformation parameters (using a separate function):
					if(!ReadParamsOfSeveralTrfsFromML(pTrfParInExtrStep, dummyNumTrf, strTrfOrderAux, false)) return 0; //nested call, but not deep
					//{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
				}
				//readFinished = true;
				//break;
			}
			else return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
		}
		else if(((next==MLTKINT) || (next==MLTKREAL)) && (numPar02 == 2) && (numPar01 <= 2))
		{//try to read one homothety {{kxH,kyH},phi}

			//double arBuf1[] = {0,0,0};
			arBuf1[0] = 0; arBuf1[1] = 0; arBuf1[2] = 0;
			if(!MLGetReal(stdlink, arBuf1)) return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
			if(!MLGetReal(stdlink, arBuf1 + 1)) return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}

			//bool noOtherTrfExpected = true;
			noOtherTrfExpected = true;
			if(numPar01 == 2)
			{
				next = MLGetNext(stdlink);
				if((next==MLTKINT) || (next==MLTKREAL))
				{
					if(!MLGetReal(stdlink, arBuf1 + 2)) return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
				}
				else if(next == MLTKFUNC)
				{
					noOtherTrfExpected = false;
				}
			}

			if(noOtherTrfExpected)
			{
				if(pTrfParInExtrStep == 0) pTrfParInExtrStep = new double*[1]; 
				pTrfParInExtrStep[0] = 0;
				numTrf = 1;
				if(strTrfOrderInExtrStep == 0) strTrfOrderInExtrStep = new char[2];
			}
			else
			{
				if(pTrfParInExtrStep == 0) pTrfParInExtrStep = new double*[2]; 
				pTrfParInExtrStep[0] = 0;
				pTrfParInExtrStep[1] = 0;
				numTrf = 2;
				if(strTrfOrderInExtrStep == 0) strTrfOrderInExtrStep = new char[3];
				strTrfOrderInExtrStep[2] = '\0';
			}

			*pTrfParInExtrStep = new double[3];
			strTrfOrderInExtrStep[0] = 'h';
			strTrfOrderInExtrStep[1] = '\0';
			//double *pCurTrfPar = pTrfParInExtrStep[0];
			pCurTrfPar = pTrfParInExtrStep[0];
			for(int ii=0; ii<3; ii++)
			{
				pCurTrfPar[ii] = arBuf1[ii];
			}

			if(!noOtherTrfExpected)
			{
				//int dummyNumTrf = 0;
				//double **pTrfParInExtrStepAux = pTrfParInExtrStep + 1;
				//char *strTrfOrderAux = strTrfOrderInExtrStep + 1;
				dummyNumTrf = 0;
				pTrfParInExtrStepAux = pTrfParInExtrStep + 1;
				strTrfOrderAux = strTrfOrderInExtrStep + 1;
				//try to read-in one more transformation parameters (using a separate function):
				if(!ReadParamsOfSeveralTrfsFromML(pTrfParInExtrStepAux, dummyNumTrf, strTrfOrderAux, false)) return 0; //nested call, but not deep
				//{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
			}
			//readFinished = true;
			//break;
		}
		else if(readDeepNested)
		{
			if(pTrfParInExtrStep == 0) pTrfParInExtrStep = new double*[numPar01]; 
			numTrf = numPar01;
			if(strTrfOrderInExtrStep == 0) strTrfOrderInExtrStep = new char[numPar01 + 1];
			for(int jj=0; jj<=numPar01; jj++) 
			{
				if(jj != numPar01) pTrfParInExtrStep[jj] = 0;
				strTrfOrderInExtrStep[jj] = '\0';
			}

			//int dummyNumTrf = numPar02; //to signal that the first list was already read and the number of elements is numPar02
			dummyNumTrf = numPar02; //to signal that the first list was already read and the number of elements is numPar02
			if(!ReadParamsOfSeveralTrfsFromML(pTrfParInExtrStep, dummyNumTrf, strTrfOrderInExtrStep, false, next)) return 0; //nested call, but not deep
			//{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}

			for(int jj=1; jj<numPar01; jj++)
			{
				//double **pTrfParInExtrStepAux = pTrfParInExtrStep + jj;
				//char *strTrfOrderAux = strTrfOrderInExtrStep + jj;
				pTrfParInExtrStepAux = pTrfParInExtrStep + jj;
				strTrfOrderAux = strTrfOrderInExtrStep + jj;
				dummyNumTrf = 0;
				if(!ReadParamsOfSeveralTrfsFromML(pTrfParInExtrStepAux, dummyNumTrf, strTrfOrderAux, false)) return 0; //nested call, but not deep
				//{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
			}
		}
		else return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
	}
	else return 0; //{ rad.Send.ErrorMessage("Radia::Error000"); return 0;}
	//}

	return 1;
}
#endif

//-------------------------------------------------------------------------

void MultGenExtrPolygonCurOrMag(char cur_or_mag)
//void MultGenExtrPolygonCur()
{//radObjMltExtPgnCur[z:0,a:"z",({{x1,y1},{x2,y2},...},{{R1,T1,H1},{R2,T2,H2},...}},I,Frame->Loc|Lab]
//or radObjMltExtPgnMag[z:0,a:"z",({{x1,y1},{x2,y2},...},{{R1,T1,H1},{R2,T2,H2},...}},{{mx1,my1,mz1},{mx2,my2,mz2},...},Frame->Loc|Lab]
//Ri: {{x,y,z},{vx,vy,vz},phi}
//Ti: {vx,vy,vz}
//Hi: k or {kx,ky} or {{kx,ky},phi}
#ifdef __MATHEMATICA__

	double avgCur = 0;
	double zc = 0; //default value
	int lenArrayOfPoints2d = 0;
	TVector2d *ArrayOfPoints2d = NULL;
	TVector2d *tArrayOfPoints2d = NULL;
	char defOrientStr[] = "z\0";
	const char* OrientStr = defOrientStr;
	bool origHeightDefined = false, origOrientStrDefined = false, origOrientSymDefined = false;

	double ***arPtrTrfParInExtrSteps = 0;
	int *arNumTrfInExtrSteps = 0; //there can be several transformations in each extrusion step
	char **arStrTrfOrderInExtrSteps = 0;
	double *arMagnCompInSteps=0;

	int NumSteps = 0, NumTrfInStep = 0;

	int numOptions = 0;
	const char* arOptionNames[] = {0,0,0};
	const char* arOptionValues[] = {0,0,0};

	const char *FunName;
	bool basePgnAndExtrusionDefinedSeparately = true;
	int lenList = 0;
	int lenListPolygSubdExtr = 0;

	long *SubdDims=0;
	char **SubdHeads=0;
	long SubdDepth=0;
	double *SubdDataPtr=0, *AuxSubdDataPtr=0, *tAuxSubdDataPtr=0, *tSubdDataPtr=0;
	int numSubdParam=0;

	int ReadOK=0, ArrayReadOK=0;
	int OptionCount = 0;

	int Next = MLGetNext(stdlink);
	if((Next == MLTKINT) || (Next == MLTKREAL))
	{//Assuming zc
		if(!MLGetReal(stdlink, &zc)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		origHeightDefined = true;
		Next = MLGetNext(stdlink);
	}

	if(Next == MLTKSTR)
	{//Assuming a
		if(!MLGetString(stdlink, &OrientStr)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
		origOrientStrDefined = true;
		Next = MLGetNext(stdlink);
	}
	else if(Next == MLTKSYM)
	{
		if(!MLGetSymbol(stdlink, &OrientStr)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
		origOrientSymDefined = true;
		Next = MLGetNext(stdlink);
	}

	if(Next != MLTKFUNC) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
	ReadOK = MLGetFunction(stdlink, &FunName, &lenList);
	if((!ReadOK) || (strcmp(FunName, "List") != 0) || (lenList < 2)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015

	lenListPolygSubdExtr = lenList;
	//if(lenList == 2)
	if(lenListPolygSubdExtr == 2)
	{//reading list of Vector2d manually, assuming {{x1,y1},{x2,y2},...},{{R1,T1,H1},{R2,T2,H2},...} form
		ReadOK = MLGetFunction(stdlink, &FunName, &lenList);
		if((!ReadOK) || (strcmp(FunName, "List") != 0) || (lenList <= 0)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015

		basePgnAndExtrusionDefinedSeparately = false; //OC260713
	}
	else if(lenListPolygSubdExtr == 3)
	{
		ReadOK = MLGetFunction(stdlink, &FunName, &lenList);
		if((!ReadOK) || (strcmp(FunName, "List") != 0) || (lenList <= 0)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015

		if(lenList == 2)
		{
			lenArrayOfPoints2d = 3;
			ArrayOfPoints2d = new TVector2d[3]; 
			if(!MLGetReal(stdlink, &(ArrayOfPoints2d->x))) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
			if(!MLGetReal(stdlink, &(ArrayOfPoints2d->y))) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
			//TVector2d *tArrayOfPoints2d = ArrayOfPoints2d + 1;
			tArrayOfPoints2d = ArrayOfPoints2d + 1;
			for(int ii=0; ii<2; ii++)
			{
				if(!rad.Send.GetVector2d(*(tArrayOfPoints2d++))) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
			}
			lenListPolygSubdExtr = 0;
		}
	}

	if(ArrayOfPoints2d == NULL)
	{
		lenArrayOfPoints2d = lenList;
		ArrayOfPoints2d = new TVector2d[lenList]; 
		//TVector2d *tArrayOfPoints2d = ArrayOfPoints2d;
		tArrayOfPoints2d = ArrayOfPoints2d;
		for(int ii=0; ii<lenList; ii++)
		{
			if(!rad.Send.GetVector2d(*(tArrayOfPoints2d++))) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
		}
	}

	if(lenListPolygSubdExtr == 3)
	{//reading triangulation parameters, assuming {{x1,y1},{x2,y2},...},{{k1,q1},{k2,q2},...},{{R1,T1,H1},{R2,T2,H2},...} form

		//int ArrayReadOK = MLGetDoubleArray(stdlink, &SubdDataPtr, &SubdDims, &SubdHeads, &SubdDepth);
		ArrayReadOK = MLGetDoubleArray(stdlink, &SubdDataPtr, &SubdDims, &SubdHeads, &SubdDepth);
		if(!ArrayReadOK) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		if((SubdDepth < 1) || (SubdDepth > 2)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		//int numSubdParam = SubdDims[0];
		numSubdParam = SubdDims[0];
		if(numSubdParam <= 0) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		if(numSubdParam != lenArrayOfPoints2d) { rad.Send.ErrorMessage("Radia::Error097"); return;}

		//Parse subdivision parameters: {{k1,q1},{k2,q2},...} or {k1,k2,...}
		if(SubdDepth == 2) AuxSubdDataPtr = SubdDataPtr;
		else
		{
			AuxSubdDataPtr = new double[numSubdParam << 1];
			//double *tAuxSubdDataPtr = AuxSubdDataPtr, *tSubdDataPtr = SubdDataPtr;
			tAuxSubdDataPtr = AuxSubdDataPtr; tSubdDataPtr = SubdDataPtr;
			for(int i=0; i<numSubdParam; i++)
			{
				*(tAuxSubdDataPtr++) = *(tSubdDataPtr++);
				*(tAuxSubdDataPtr++) = 1.;
			}
		}
		basePgnAndExtrusionDefinedSeparately = false; //OC260713
	}

	ReadOK = MLGetFunction(stdlink, &FunName, &NumSteps);
	if((!ReadOK) || (strcmp(FunName, "List") != 0) || (NumSteps <= 0)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015

	arPtrTrfParInExtrSteps = new double**[NumSteps];
	arNumTrfInExtrSteps = new int[NumSteps]; //there can be several rotations in each extrusion step
	arStrTrfOrderInExtrSteps = new char*[NumSteps];

	for(int j=0; j<NumSteps; j++) 
	{//to allow for easy delete
		arPtrTrfParInExtrSteps[j] = 0; 
		arNumTrfInExtrSteps[j] = 0;
		arStrTrfOrderInExtrSteps[j] = 0;
	}

	for(int k=0; k<NumSteps; k++)
	{
		ReadOK = MLGetFunction(stdlink, &FunName, &NumTrfInStep);
		if((!ReadOK) || (strcmp(FunName, "List") != 0) || (NumTrfInStep < 0)) 
		{ rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015

		double **&arPtrTrfParInExtrSteps_k = arPtrTrfParInExtrSteps[k];
		char *&arStrTrfOrderInExtrSteps_k = arStrTrfOrderInExtrSteps[k];

		arPtrTrfParInExtrSteps_k = 0;
		arStrTrfOrderInExtrSteps_k = 0;
		//arNumTrfInExtrSteps[k] = 0;
		arNumTrfInExtrSteps[k] = NumTrfInStep;

		if(!ReadParamsOfSeveralTrfsFromML(arPtrTrfParInExtrSteps_k, arNumTrfInExtrSteps[k], arStrTrfOrderInExtrSteps_k)) 
		{ rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
	}

	if(cur_or_mag ==  'c')
	{//Assuming one number I
		if(!MLGetReal(stdlink, &avgCur)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
		numOptions = 1; //OC250713
	}
	else if(cur_or_mag ==  'm')
	{//Assuming {{mx1,my1,mz1},{mx2,my2,mz2},...} or {mx1,my1,mz1} or 0
		int totNumComp = NumSteps*3;
		arMagnCompInSteps = new double[totNumComp];
		double *t_arMagnCompInSteps = arMagnCompInSteps;
		for(int i=0; i<totNumComp; i++) *(t_arMagnCompInSteps++) = 0;

		int numElem = 0;
		ReadOK = MLGetFunction(stdlink, &FunName, &numElem);

		//OC250713 (to allow for absence of M in arguments)
		//if((!ReadOK) || (strcmp(FunName, "List") != 0) || (numElem <= 0)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
		if(!ReadOK) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}

		if((strcmp(FunName, "List") == 0) && (numElem > 0))
		{
			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC091015

			Next = MLGetNext(stdlink);
			if(((Next == MLTKINT) || (Next == MLTKREAL)) && (numElem == 3))
			{//Assuming {mx1,my1,mz1}
				t_arMagnCompInSteps = arMagnCompInSteps;
				for(int i=0; i<3; i++)
				{
					//if(!MLGetReal(stdlink, t_arMagnCompInSteps++)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
					if(!MLGetReal(stdlink, t_arMagnCompInSteps++)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
				}
			}
			else if(Next == MLTKFUNC)
			{//Assuming {{mx1,my1,mz1},{mx2,my2,mz2},...}
				if(numElem != NumSteps) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}

				t_arMagnCompInSteps = arMagnCompInSteps;
				for(int i=0; i<NumSteps; i++)
				{
					ReadOK = MLGetFunction(stdlink, &FunName, &numElem);
					if((!ReadOK) || (strcmp(FunName, "List") != 0) || (numElem != 3)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
					//MLDisownSymbol(stdlink, FunName);
					MLWrapDeleteSymbol(stdlink, FunName); //OC091015
					for(int i=0; i<3; i++)
					{
						//if(!MLGetReal(stdlink, t_arMagnCompInSteps++)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
						if(!MLGetReal(stdlink, t_arMagnCompInSteps++)) { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
					}
				}
			}
			//else { rad.Send.ErrorMessage("Radia::Error000"); return;}
			else { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}
		}
		else if((strcmp(FunName, "Rule") == 0) && (numElem == 2))
		{//OC250713
			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC091015
			if(!MLGetSymbol(stdlink, arOptionNames + OptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
			if(!MLGetString(stdlink, arOptionValues + (OptionCount++))) { rad.Send.ErrorMessage("Radia::Error062"); goto Finish;}
		}
		else { rad.Send.ErrorMessage("Radia::Error000"); goto Finish;}

		numOptions = 5; //OC250713
	}

	//reading options "Frame->Loc|Lab" and eventually the Triangulation options
	//numOptions = 1;
	//numOptions = 5; //OC250713 (commented-out (moved up))
	if(!origHeightDefined) numOptions++;
	if(!(origOrientStrDefined || origOrientSymDefined)) numOptions++;
	if(!basePgnAndExtrusionDefinedSeparately) numOptions++;
	Next = 0;
	//if(!AuxGetMathematicaOptions(Next, arOptionNames, arOptionValues, numOptions)) goto Finish;
	if(!AuxGetMathematicaOptions(Next, arOptionNames + OptionCount, arOptionValues + OptionCount, numOptions)) goto Finish; //OC250713

	//TEST
	//Next = MLGetNext(stdlink);
	//END TEST

	//rad.SetMultGenExtrPolygonCur(zc, OrientStr, ArrayOfPoints2d, lenArrayOfPoints2d, arPtrTrfParInExtrSteps, arStrTrfOrderInExtrSteps, arNumTrfInExtrSteps, NumSteps, avgCur, arOptionNames, arOptionValues, numOptions);
	//rad.SetMultGenExtrPolygonCur(zc, OrientStr, ArrayOfPoints2d, lenArrayOfPoints2d, arPtrTrfParInExtrSteps, arStrTrfOrderInExtrSteps, arNumTrfInExtrSteps, NumSteps, avgCur, arMagnCompInSteps, arOptionNames, arOptionValues, numOptions);
	rad.SetMultGenExtrPolygonCur(zc, OrientStr, ArrayOfPoints2d, lenArrayOfPoints2d, AuxSubdDataPtr, arPtrTrfParInExtrSteps, arStrTrfOrderInExtrSteps, arNumTrfInExtrSteps, NumSteps, avgCur, arMagnCompInSteps, arOptionNames, arOptionValues, numOptions);

Finish:

	if(ArrayOfPoints2d != 0) delete[] ArrayOfPoints2d;

	//if(origOrientStrDefined) MLDisownString(stdlink, OrientStr);
	if(origOrientStrDefined) MLWrapDeleteString(stdlink, OrientStr); //OC091015
	//else if(origOrientSymDefined) MLDisownSymbol(stdlink, OrientStr);
	else if(origOrientSymDefined) MLWrapDeleteSymbol(stdlink, OrientStr); //OC091015

	if(SubdDepth != 2) delete[] AuxSubdDataPtr;
	//MLDisownDoubleArray(stdlink, SubdDataPtr, SubdDims, SubdHeads, SubdDepth);
	MLWrapDeleteDoubleArray(stdlink, SubdDataPtr, SubdDims, SubdHeads, SubdDepth); //OC091015

	for(int i=0; i<NumSteps; i++)
	{
		int curNumTrf = 0;
		if(arNumTrfInExtrSteps != 0) curNumTrf = arNumTrfInExtrSteps[i];

		if((arPtrTrfParInExtrSteps != 0) && (curNumTrf > 0))
		{
			double **pCurPtrTrfParInExtrSteps = arPtrTrfParInExtrSteps[i];
			if(pCurPtrTrfParInExtrSteps != 0)
			{
				for(int j=0; j<curNumTrf; j++) 
				{
					if(pCurPtrTrfParInExtrSteps[j] != 0) 
					{
						delete[] pCurPtrTrfParInExtrSteps[j];
						pCurPtrTrfParInExtrSteps[j] = 0;
					}
				}
				delete[] pCurPtrTrfParInExtrSteps;
				pCurPtrTrfParInExtrSteps = 0;
			}
		}
		if(arStrTrfOrderInExtrSteps != 0)
		{
			if(arStrTrfOrderInExtrSteps[i] != 0)
			{
				delete[] arStrTrfOrderInExtrSteps[i];
				arStrTrfOrderInExtrSteps[i] = 0;
			}
		}
	}

	if(arPtrTrfParInExtrSteps != 0) delete[] arPtrTrfParInExtrSteps;
	if(arStrTrfOrderInExtrSteps != 0) delete[] arStrTrfOrderInExtrSteps;
	if(arNumTrfInExtrSteps != 0) delete[] arNumTrfInExtrSteps;

	if(arMagnCompInSteps != 0) delete[] arMagnCompInSteps;

	for(int k=0; k<numOptions; k++) 
	{
		//if(arOptionNames[k] != 0) MLDisownSymbol(stdlink, arOptionNames[k]);
		if(arOptionNames[k] != 0) MLWrapDeleteSymbol(stdlink, arOptionNames[k]); //OC091015
		//if(arOptionValues[k] != 0) MLDisownString(stdlink, arOptionValues[k]);
		if(arOptionValues[k] != 0) MLWrapDeleteString(stdlink, arOptionValues[k]); //OC091015
	}
#endif
}

//-------------------------------------------------------------------------

void MultGenExtrPolygonCur()
{
	MultGenExtrPolygonCurOrMag('c');
}
void MultGenExtrPolygonMag()
{
	MultGenExtrPolygonCurOrMag('m');
}

//-------------------------------------------------------------------------

void MultGenExtrRectangle()
{
#ifdef __MATHEMATICA__
	const char *FunName;
	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;

	int ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	if((!ReadOK) || strcmp(FunName, "List") || (ArgNum <= 0)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015

	int AmOfLayerRect = ArgNum;

	TVector3d* RectCenPoints = new TVector3d[AmOfLayerRect];
	if(RectCenPoints == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector2d* RectDims = new TVector2d[AmOfLayerRect];
	if(RectDims == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}

	TVector3d* tRectCenPoints = RectCenPoints;
	TVector2d* tRectDims = RectDims;

	for(int i=0; i<AmOfLayerRect; i++)
	{
		ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		if((!ReadOK) || strcmp(FunName, "List") || (ArgNum != 2)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015

		ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		if((!ReadOK) || strcmp(FunName, "List") || (ArgNum != 3)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015
		if(!MLGetDouble(stdlink, &(tRectCenPoints->x))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		if(!MLGetDouble(stdlink, &(tRectCenPoints->y))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		if(!MLGetDouble(stdlink, &(tRectCenPoints->z))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		tRectCenPoints++;

		ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		if((!ReadOK) || strcmp(FunName, "List") || (ArgNum != 2)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015
		if(!MLGetDouble(stdlink, &(tRectDims->x))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		if(!MLGetDouble(stdlink, &(tRectDims->y))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		tRectDims++;
	}

	double Magn[3];
	if(!MLGetReal(stdlink, &(Magn[0]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	if(!MLGetReal(stdlink, &(Magn[1]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	if(!MLGetReal(stdlink, &(Magn[2]))) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	rad.SetMultGenExtrRectangle(RectCenPoints, RectDims, AmOfLayerRect, Magn, 3);

	if(RectDims != 0) delete[] RectDims;
	if(RectCenPoints != 0) delete[] RectCenPoints;
#endif
}

//-------------------------------------------------------------------------

void MultGenExtrRectangleOpt(double** Layers, int AmOfLayers, double* M)
{
	TVector3d* RectCenPoints = new TVector3d[AmOfLayers];
	if(RectCenPoints == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector2d* RectDims = new TVector2d[AmOfLayers];
	if(RectDims == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector3d* tRectCenPoints = RectCenPoints;
	TVector2d* tRectDims = RectDims;
	double** tLayers = Layers;

	for(int i=0; i<AmOfLayers; i++)
	{
		double* tCoord = *(tLayers++);
		tRectCenPoints->x = *(tCoord++); tRectCenPoints->y = *(tCoord++); (tRectCenPoints++)->z = *(tCoord++); 
		tRectDims->x = *(tCoord++); (tRectDims++)->y = *(tCoord++);
	}
	rad.SetMultGenExtrRectangle(RectCenPoints, RectDims, AmOfLayers, M, 3);

	if(RectDims != 0) delete[] RectDims;
	if(RectCenPoints != 0) delete[] RectCenPoints;
}

//-------------------------------------------------------------------------

void MultGenExtrRectangleDLL(double* pFlatCenPts, double* pFlatRtgSizes, int AmOfLayers, double* pM)
{
	TVector3d* RectCenPoints = new TVector3d[AmOfLayers];
	if(RectCenPoints == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector2d* RectDims = new TVector2d[AmOfLayers];
	if(RectDims == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}

	TVector3d* tRectCenPoints = RectCenPoints;
	TVector2d* tRectDims = RectDims;

	double* tCenPts = pFlatCenPts;
	double* tRtgSizes = pFlatRtgSizes;
	for(int i=0; i<AmOfLayers; i++)
	{
		tRectCenPoints->x = *(tCenPts++); tRectCenPoints->y = *(tCenPts++); (tRectCenPoints++)->z = *(tCenPts++); 
		tRectDims->x = *(tRtgSizes++); (tRectDims++)->y = *(tRtgSizes++);
	}
	rad.SetMultGenExtrRectangle(RectCenPoints, RectDims, AmOfLayers, pM, 3);

	if(RectDims != 0) delete[] RectDims;
	if(RectCenPoints != 0) delete[] RectCenPoints;
}

//-------------------------------------------------------------------------

void ArcMag(double xc, double yc, double zc, double rmin, double rmax, double phimin, double phimax, double Height, int NumberOfSegm, char* Orient, double mx, double my, double mz)
{
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Radii[] = {fabs(radCR.Double(rmin)), fabs(radCR.Double(rmax))};
	double Angles[] = {phimin, phimax}; // Consider shrinking this
	double Magn[] = {mx, my, mz};
	rad.SetArcMag(CPoi, 3, Radii, 2, Angles, 2, Height, NumberOfSegm, Magn, 3, Orient);
}

//-------------------------------------------------------------------------

void ArcPolygon()
{
#ifdef __MATHEMATICA__

	double xc, yc; //, zc;
	if(!MLGetReal(stdlink, &xc)) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	if(!MLGetReal(stdlink, &yc)) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	//if(!MLGetReal(stdlink, &zc)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	double CenP[] = {radCR.Double(xc), radCR.Double(yc)}; //, radCR.Double(zc)};

	const char* OrientStr = 0;
	if(!MLGetString(stdlink, &OrientStr)) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	int lenArrayOfPoints2d;
	TVector2d* ArrayOfPoints2d = NULL;
	if(!rad.Send.GetArrayOfVector2d(ArrayOfPoints2d, lenArrayOfPoints2d)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	//TVector2d vSmallArbShift(radCR.Double(0), radCR.Double(0));
	//for(int i=0; i<lenArrayOfPoints2d; i++)
	//{
    //       ArrayOfPoints2d[i] += vSmallArbShift;
	//}

	double Angles[] = {0, 3.1415926}; 
	if(!MLGetReal(stdlink, Angles)) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	if(!MLGetReal(stdlink, Angles + 1)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	int NumberOfSegm = 1;
	if(!MLGetInteger(stdlink, &NumberOfSegm)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	const char* SymOrNoSymStr = 0;
	if(!MLGetString(stdlink, &SymOrNoSymStr)) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	double Magn[3];
	if(!MLGetReal(stdlink, &(Magn[0]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if(!MLGetReal(stdlink, &(Magn[1]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if(!MLGetReal(stdlink, &(Magn[2]))) { rad.Send.ErrorMessage("Radia::Error000"); return;}

    rad.SetArcPolygon(CenP, OrientStr, ArrayOfPoints2d, lenArrayOfPoints2d, Angles, NumberOfSegm, SymOrNoSymStr, Magn);

	if(ArrayOfPoints2d != NULL) delete[] ArrayOfPoints2d;
	//if(SymOrNoSymStr != NULL) MLDisownString(stdlink, SymOrNoSymStr);
	if(SymOrNoSymStr != NULL) MLWrapDeleteString(stdlink, SymOrNoSymStr); //OC091015
	//if(OrientStr != NULL) MLDisownString(stdlink, OrientStr);
	if(OrientStr != NULL) MLWrapDeleteString(stdlink, OrientStr); //OC091015

#endif
}

//-------------------------------------------------------------------------

void ArcPolygonDLL(double xc, double yc, char a, double* pFlatVert, int nv, double PhiMin, double PhiMax, int nseg, char sym_no, double mx, double my, double mz)
{
	double CenP[] = {radCR.Double(xc), radCR.Double(yc)};

	TVector2d* ArrayOfPoints2d = new TVector2d[nv];
	double *tFlatVert = pFlatVert;
	TVector2d *tArrayOfPoints2d = ArrayOfPoints2d;
	for(int i=0; i<nv; i++)
	{
		tArrayOfPoints2d->x = *(tFlatVert++);
		tArrayOfPoints2d->y = *(tFlatVert++);
		tArrayOfPoints2d++;
	}

    double Angles[] = {radCR.DoublePlus(PhiMin), radCR.DoubleMinus(PhiMax)}; 

	char SymOrNoSymStr[20];
	if((sym_no == 's') || (sym_no == 'S')) strcpy(SymOrNoSymStr, "sym");
	else strcpy(SymOrNoSymStr, "nosym");

	double Magn[] = {mx, my, mz};

    rad.SetArcPolygon(CenP, &a, ArrayOfPoints2d, nv, Angles, nseg, SymOrNoSymStr, Magn);

	if(ArrayOfPoints2d != NULL) delete[] ArrayOfPoints2d;
}

//-------------------------------------------------------------------------

//void CylMag(double xc, double yc, double zc, double r, double h, int NumberOfSegm, double mx, double my, double mz, char* Orient)
//{
//	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
//	double Magn[] = {mx, my, mz};
//	rad.SetCylMag(CPoi, 3, r, h, NumberOfSegm, Magn, 3, Orient);
//}

//-------------------------------------------------------------------------

void CylMag(double xc, double yc, double zc, double r, double h, int NumberOfSegm, char* Orient, double mx, double my, double mz)
{
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Magn[] = {mx, my, mz};
	rad.SetCylMag(CPoi, 3, r, h, NumberOfSegm, Magn, 3, Orient);
}

//-------------------------------------------------------------------------

void ArcCur(double xc, double yc, double zc, double rmin, double rmax, double phimin, double phimax, double Height, int NumberOfSegm, double J_azim, char* ManOrAuto, char* Orient)
{
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Radii[] = {fabs(radCR.Double(rmin)), fabs(radCR.Double(rmax))};
	double Angles[] = {phimin, phimax}; // Consider shrinking this
	rad.SetArcCur(CPoi, 3, Radii, 2, Angles, 2, Height, J_azim, NumberOfSegm, ManOrAuto, Orient);
}

//-------------------------------------------------------------------------

void RaceTrack(double xc, double yc, double zc, double rmin, double rmax, double Lx, double Ly, double Height, int NumberOfSegm, double J_azim, char* ManOrAuto, char* Orient)
{
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Radii[] = {fabs(radCR.Double(rmin)), fabs(radCR.Double(rmax))};
	double StrPartDims[] = {(Lx==0.)? Lx : radCR.Double(Lx), (Ly==0.)? Ly : radCR.Double(Ly)};
	rad.SetRaceTrack(CPoi, 3, Radii, 2, StrPartDims, 2, Height, J_azim, NumberOfSegm, ManOrAuto, Orient);
}

//-------------------------------------------------------------------------

void FlmCur()
{
#ifdef __MATHEMATICA__
	TVector3d* ArrayOfPoints = NULL;
	int lenArrayOfPoints;
	if(rad.Send.GetArrayOfVector3d(ArrayOfPoints, lenArrayOfPoints))
	{
		for(int i1=0; i1<lenArrayOfPoints; i1++)
		{
			TVector3d& Vect3d = ArrayOfPoints[i1];
			Vect3d.x = radCR.Double(Vect3d.x);
			Vect3d.y = radCR.Double(Vect3d.y);
			Vect3d.z = radCR.Double(Vect3d.z);
		}
		double I;
		MLGetReal(stdlink, &I);

		rad.SetFlmCur(I, ArrayOfPoints, lenArrayOfPoints);
		delete[] ArrayOfPoints;
	}
#endif
}

//-------------------------------------------------------------------------

void FlmCurOpt(double** Points, long LenPoints, double Cur)
{
	TVector3d* ArrayOfPoints = new TVector3d[LenPoints];
	if(ArrayOfPoints == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector3d* tArrayOfPoints = ArrayOfPoints;
	double** tPoints = Points;

	for(long k=0; k<LenPoints; k++)
	{
		double *aPoint = *(tPoints++);
		*(tArrayOfPoints++) = TVector3d(aPoint[0], aPoint[1], aPoint[2]);
	}

	rad.SetFlmCur(Cur, ArrayOfPoints, (int)LenPoints);
	delete[] ArrayOfPoints;
}

//-------------------------------------------------------------------------

void FlmCurDLL(double* Points, int LenPoints, double Cur)
{
	TVector3d* ArrayOfPoints = new TVector3d[LenPoints];
	if(ArrayOfPoints == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector3d* tArrayOfPoints = ArrayOfPoints;
	double* tPoints = Points;

	for(long k=0; k<LenPoints; k++)
	{
		double x = *(tPoints++);
		double y = *(tPoints++);
		double z = *(tPoints++);
		*(tArrayOfPoints++) = TVector3d(radCR.Double(x), radCR.Double(y), radCR.Double(z));
	}
	rad.SetFlmCur(Cur, ArrayOfPoints, (int)LenPoints);
	delete[] ArrayOfPoints;
}

//-------------------------------------------------------------------------

void Rectngl(double xc, double yc, double zc, double Lx, double Ly)
{
	double CPoi[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	double Dims[] = {radCR.Double(Lx), radCR.Double(Ly)};
	rad.SetRectangle(CPoi, 3, Dims, 2);
}

//-------------------------------------------------------------------------

void BackgroundFieldSource(double Bx, double By, double Bz)
{
	double B[] = {Bx, By, Bz};
	rad.SetBackgroundFieldSource(B, 3);
}

//-------------------------------------------------------------------------

void Group(int* ArrayOfKeys, long lenArrayOfKeys)
{
	rad.SetGroup(ArrayOfKeys, lenArrayOfKeys);
}

//-------------------------------------------------------------------------

void AddToGroup(int GroupKey, int* ArrayOfKeys, long lenArrayOfKeys)
{
	rad.AddToGroup(GroupKey, ArrayOfKeys, lenArrayOfKeys);
}

//-------------------------------------------------------------------------

void OutGroupSize(int ElemKey)
//void OutGroupSize(int ElemKey, char deep=0)
{
	rad.OutGroupSize(ElemKey);
	//rad.OutGroupSize(ElemKey, deep);
}

//-------------------------------------------------------------------------

void OutGroupSubObjectKeys(int ElemKey)
{
	rad.OutGroupSubObjectKeys(ElemKey);
}

//-------------------------------------------------------------------------

void DuplicateElementG3D()
{
#ifdef __MATHEMATICA__
	int ElemKey;
	MLGetInteger(stdlink, &ElemKey);

	const char *FunName;
	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;

	const int AmOfOptions = 1;
	const char* OptionNames[] = {0};
	const char* OptionValues[] = {0};
	int OptionCount = 0;
	for(int ii=0; ii<AmOfOptions; ii++)
	{
		int Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC)
		{
			int ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
			if((!ReadOK) || strcmp(FunName, "Rule") || (ArgNum!=2)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC091015

			if(!MLGetSymbol(stdlink, OptionNames + OptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			if(!MLGetSymbol(stdlink, OptionValues + (OptionCount++))) { rad.Send.ErrorMessage("Radia::Error062"); return;}
		}
		else { MLNewPacket(stdlink); break;}
	}

	rad.DuplicateElement_g3d(ElemKey, OptionNames, OptionValues, OptionCount);

	for(int k=0; k<AmOfOptions; k++) 
	{
		//if(OptionNames[k] != 0) MLDisownSymbol(stdlink, OptionNames[k]);
		if(OptionNames[k] != 0) MLWrapDeleteSymbol(stdlink, OptionNames[k]); //OC091015
		//if(OptionValues[k] != 0) MLDisownSymbol(stdlink, OptionValues[k]);
		if(OptionValues[k] != 0) MLWrapDeleteSymbol(stdlink, OptionValues[k]); //OC091015
	}
#endif
}

//-------------------------------------------------------------------------

void DuplicateElementG3DOpt(int ElemKey, const char* Opt)
{
	const char* OptionNames[] = {0};
	const char* OptionValues[] = {0};
	int OptionCount = 0;

	char CharBuf[200];
	if((Opt != 0) && (*Opt != '\0')) 
	{
		strcpy(CharBuf, Opt);
		char *pEndOptName = strrchr(CharBuf, '-');
		if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
		*pEndOptName = '\0';
		OptionNames[0] = CharBuf;
		OptionValues[0] = strrchr(Opt, '>') + 1;
		OptionCount++;
	}

	rad.DuplicateElement_g3d(ElemKey, OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

void CreateFromG3DObjectWithSymmetries(int ElemKey)
{
	rad.CreateFromObj_g3dWithSym(ElemKey);
}

//-------------------------------------------------------------------------

void NumberOfDegOfFreedom(int ElemKey)
{
	rad.ComputeNumberOfDegOfFreedom(ElemKey);
}

//-------------------------------------------------------------------------

void MagnOfObj(int ElemKey)
{
	//rad.ComputeMagnInCenter(ElemKey);
	rad.ComputeMagnOrJ_InCenter(ElemKey, 'M');
}

//-------------------------------------------------------------------------

void ObjField(int ElemKey, char* fldType)
{
	rad.ComputeMagnOrJ_InCenter(ElemKey, *fldType);
}

//-------------------------------------------------------------------------

void ScaleCurInObj(int ElemKey, double scaleCoef)
{
	rad.ScaleCurrent(ElemKey, scaleCoef);
}

//-------------------------------------------------------------------------

void SetObjMagn(int ElemKey, double Mx, double My, double Mz)
{
	rad.SetObjMagn(ElemKey, Mx, My, Mz);
}

//-------------------------------------------------------------------------

void SubdivideElementG3D()
{
#ifdef __MATHEMATICA__
	int ElemKey;
	MLGetInteger(stdlink, &ElemKey);

	double SubdivArray[6];

	const char *FunName;
	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;
	double DoubleValue;
	int ReadOK = 0;

	int Next = MLGetNext(stdlink);
	if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	else { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if((!ReadOK) || strcmp(FunName, "List") || ArgNum!=3) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015

	for(int kk=0; kk<3; kk++)
	{
		Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC)
		{
			ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
			if((!ReadOK) || strcmp(FunName, "List") || ArgNum!=2) { rad.Send.ErrorMessage("Radia::Error000"); return;}
			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC091015
			Next = MLGetNext(stdlink);
			if((Next==MLTKINT) || (Next==MLTKREAL)) { if(MLGetDouble(stdlink, &DoubleValue)) SubdivArray[kk*2] = DoubleValue;}
			else { rad.Send.ErrorMessage("Radia::Error000"); return;}
			Next = MLGetNext(stdlink);
			if((Next==MLTKINT) || (Next==MLTKREAL)) { if(MLGetDouble(stdlink, &DoubleValue)) SubdivArray[kk*2 + 1] = DoubleValue;}
			else { rad.Send.ErrorMessage("Radia::Error000"); return;}
		}
		else if((Next==MLTKINT) || (Next==MLTKREAL))
		{
			if(MLGetDouble(stdlink, &DoubleValue)) SubdivArray[kk*2] = DoubleValue;
			SubdivArray[kk*2 + 1] = 1.;
		}
		else { rad.Send.ErrorMessage("Radia::Error000"); return;}
	}

// Optional parameters
// Additional Specifications (for cylindrical subdivision or subdivision by parallel planes)
	char StartReadingAdditionalSpecifications = 0;

	char StartReadingGeneralOptions = 0;
	Next = MLGetNext(stdlink);
	if((Next==MLTKFUNC) && MLGetFunction(stdlink, &FunName, &ArgNum))
	{
		if(!strcmp(FunName, "List")) 
		{
			const char* ExtraSpecNameStr;
			Next = MLGetNext(stdlink);
			if(Next==MLTKSTR)
			{
				if(!MLGetString(stdlink, &ExtraSpecNameStr)) { rad.Send.ErrorMessage("Radia::Error065"); return;}

				if((!strcmp(ExtraSpecNameStr, "cyl")) || (!strcmp(ExtraSpecNameStr, "CYL")) || (!strcmp(ExtraSpecNameStr, "Cyl"))) StartReadingAdditionalSpecifications = 1;
				else if((!strcmp(ExtraSpecNameStr, "pln")) || (!strcmp(ExtraSpecNameStr, "PLN")) || (!strcmp(ExtraSpecNameStr, "Pln"))) StartReadingAdditionalSpecifications = 2;
				else { rad.Send.ErrorMessage("Radia::Error065"); return;}

				//MLDisownString(stdlink, ExtraSpecNameStr);
				MLWrapDeleteString(stdlink, ExtraSpecNameStr); //OC091015
			}
			else if(Next==MLTKSYM)
			{
				if(!MLGetSymbol(stdlink, &ExtraSpecNameStr)) { rad.Send.ErrorMessage("Radia::Error065"); return;}

				if((!strcmp(ExtraSpecNameStr, "cyl")) || (!strcmp(ExtraSpecNameStr, "CYL")) || (!strcmp(ExtraSpecNameStr, "Cyl"))) StartReadingAdditionalSpecifications = 1;
				else if((!strcmp(ExtraSpecNameStr, "pln")) || (!strcmp(ExtraSpecNameStr, "PLN")) || (!strcmp(ExtraSpecNameStr, "Pln"))) StartReadingAdditionalSpecifications = 2;
				else { rad.Send.ErrorMessage("Radia::Error065"); return;}

				//MLDisownSymbol(stdlink, ExtraSpecNameStr);
				MLWrapDeleteSymbol(stdlink, ExtraSpecNameStr); //OC091015
			}
			ArgNum--;
		}
		else if((!strcmp(FunName, "Rule")) && (ArgNum==2)) StartReadingGeneralOptions = 1;
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015
	}

	char TypeExtraSpec = char(StartReadingAdditionalSpecifications);
	long gArgNum = ArgNum;

	double ExtraSpec[10];
	int LenExtraSpec = 0;
	for(int ss=0; ss<10; ss++) ExtraSpec[ss] = 0.;
	double* tExtraSpec = ExtraSpec;

	if(StartReadingAdditionalSpecifications == 1)
	{
		if(ArgNum==3) // Assuming {{{x,y,z},{vx,vy,vz}},{px,py,pz},rat}
		{
			LenExtraSpec = 10;
			//long ArgNum; //OC240907 (port to math6)
			int ArgNum;
			Next = MLGetNext(stdlink);
			if((Next==MLTKFUNC) && MLGetFunction(stdlink, &FunName, &ArgNum))
			{
				//MLDisownSymbol(stdlink, FunName);
				if((!strcmp(FunName, "List")) && (ArgNum==2))
				{
					//MLDisownSymbol(stdlink, FunName); //OC081007_BNL
					MLWrapDeleteSymbol(stdlink, FunName); //OC091015 //OC081007_BNL
					for(int j=0; j<2; j++)
					{
						//long ArgNum; //OC240907 (port to math6)
						int ArgNum;
						Next = MLGetNext(stdlink);
						if((Next==MLTKFUNC) && MLGetFunction(stdlink, &FunName, &ArgNum))
						{
							//MLDisownSymbol(stdlink, FunName);
							if((!strcmp(FunName, "List")) && (ArgNum==3))
							{
								//MLDisownSymbol(stdlink, FunName); //OC081007_BNL
								MLWrapDeleteSymbol(stdlink, FunName); //OC091015 //OC081007_BNL
								for(int i=0; i<3; i++)
								{
									Next = MLGetNext(stdlink);
									if(!(((Next==MLTKINT) || (Next==MLTKREAL)) && MLGetDouble(stdlink, tExtraSpec++))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
								}
							}
							else 
							{
								//MLDisownSymbol(stdlink, FunName);
								MLWrapDeleteSymbol(stdlink, FunName); //OC091015
								rad.Send.ErrorMessage("Radia::Error065"); return;
							}
						}
						else { rad.Send.ErrorMessage("Radia::Error065"); return;}
					}
				}
				else 
				{ 
					//MLDisownSymbol(stdlink, FunName); //OC081007_BNL
					MLWrapDeleteSymbol(stdlink, FunName); //OC091015 //OC081007_BNL
					rad.Send.ErrorMessage("Radia::Error065"); return;
				}
			}
			else { rad.Send.ErrorMessage("Radia::Error065"); return;}

			Next = MLGetNext(stdlink);
			if((Next==MLTKFUNC) && MLGetFunction(stdlink, &FunName, &ArgNum))
			{
				//MLDisownSymbol(stdlink, FunName);
				if((!strcmp(FunName, "List")) && (ArgNum==3))
				{
					//MLDisownSymbol(stdlink, FunName); //OC081007_BNL
					MLWrapDeleteSymbol(stdlink, FunName); //OC091015 //OC081007_BNL
					for(int i=0; i<3; i++)
					{
						Next = MLGetNext(stdlink);
						if(!(((Next==MLTKINT) || (Next==MLTKREAL)) && MLGetDouble(stdlink, tExtraSpec++))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
					}
				}
			}
			else 
			{ 
				//MLDisownSymbol(stdlink, FunName);
				MLWrapDeleteSymbol(stdlink, FunName); //OC091015
				rad.Send.ErrorMessage("Radia::Error065"); return;
			}

			Next = MLGetNext(stdlink);
			if(!(((Next==MLTKINT) || (Next==MLTKREAL)) && MLGetDouble(stdlink, tExtraSpec++))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		}
		else if(ArgNum==2) // Assuming {{x,y,z},{vx,vy,vz}}
		{
			LenExtraSpec = 6;
			for(int j=0; j<2; j++)
			{
				//long ArgNum; //OC240907 (port to math6)
				int ArgNum;
				Next = MLGetNext(stdlink);
				if((Next==MLTKFUNC) && MLGetFunction(stdlink, &FunName, &ArgNum))
				{
					//MLDisownSymbol(stdlink, FunName);
					if((!strcmp(FunName, "List")) && (ArgNum==3))
					{
						//MLDisownSymbol(stdlink, FunName); //OC081007_BNL
						MLWrapDeleteSymbol(stdlink, FunName); //OC091015 //OC081007_BNL
						for(int i=0; i<3; i++)
						{
							Next = MLGetNext(stdlink);
							if(!(((Next==MLTKINT) || (Next==MLTKREAL)) && MLGetDouble(stdlink, tExtraSpec++))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
						}
					}
					else { rad.Send.ErrorMessage("Radia::Error065"); return;}
				}
				else 
				{ 
					//MLDisownSymbol(stdlink, FunName); //OC081007_BNL
					MLWrapDeleteSymbol(stdlink, FunName); //OC091015 //OC081007_BNL
					rad.Send.ErrorMessage("Radia::Error065"); return;
				}
			}
		}
		else if(ArgNum==1) // Assuming {{{x,y,z},{vx,vy,vz}}}
		{
			LenExtraSpec = 6;
			//long ArgNum; //OC240907 (port to math6)
			int ArgNum;
			Next = MLGetNext(stdlink);
			if((Next==MLTKFUNC) && MLGetFunction(stdlink, &FunName, &ArgNum))
			{
				//MLDisownSymbol(stdlink, FunName);
				if((!strcmp(FunName, "List")) && (ArgNum==2))
				{
					//MLDisownSymbol(stdlink, FunName); //OC081007_BNL
					MLWrapDeleteSymbol(stdlink, FunName); //OC091015 //OC081007_BNL
					for(int j=0; j<2; j++)
					{
						//long ArgNum; //OC240907 (port to math6)
						int ArgNum;
						Next = MLGetNext(stdlink);
						if((Next==MLTKFUNC) && MLGetFunction(stdlink, &FunName, &ArgNum))
						{
							//MLDisownSymbol(stdlink, FunName);
							if((!strcmp(FunName, "List")) && (ArgNum==3))
							{
								//MLDisownSymbol(stdlink, FunName); //OC081007_BNL
								MLWrapDeleteSymbol(stdlink, FunName); //OC091015 //OC081007_BNL
								for(int i=0; i<3; i++)
								{
									Next = MLGetNext(stdlink);
									if(!(((Next==MLTKINT) || (Next==MLTKREAL)) && MLGetDouble(stdlink, tExtraSpec++))) { rad.Send.ErrorMessage("Radia::Error065"); return;}
								}
							}
							else 
							{ 
								//MLDisownSymbol(stdlink, FunName); //OC081007_BNL
								MLWrapDeleteSymbol(stdlink, FunName); //OC091015 //OC081007_BNL
								rad.Send.ErrorMessage("Radia::Error065"); return;
							}
						}
						else { rad.Send.ErrorMessage("Radia::Error065"); return;}
					}
				}
				else 
				{ 
					//MLDisownSymbol(stdlink, FunName); //OC081007_BNL
					MLWrapDeleteSymbol(stdlink, FunName); //OC091015 //OC081007_BNL
					rad.Send.ErrorMessage("Radia::Error065"); return;
				}
			}
		}
		else { rad.Send.ErrorMessage("Radia::Error065"); return;}

		StartReadingGeneralOptions = 0;
	}
	else if(StartReadingAdditionalSpecifications == 2)
	{
		//long ArgNum; //OC240907 (port to math6)
		int ArgNum;
		LenExtraSpec = gArgNum*3;
		for(int k=0; k<gArgNum; k++)
		{
			Next = MLGetNext(stdlink);
			if((Next==MLTKFUNC) && MLGetFunction(stdlink, &FunName, &ArgNum))
			{
				//MLDisownSymbol(stdlink, FunName);
				if((!strcmp(FunName, "List")) && (ArgNum==3))
				{
					//MLDisownSymbol(stdlink, FunName);
					MLWrapDeleteSymbol(stdlink, FunName); //OC091015
					for(int i=0; i<3; i++)
					{
						Next = MLGetNext(stdlink);
						if(!(((Next==MLTKINT) || (Next==MLTKREAL)) && MLGetDouble(stdlink, tExtraSpec++))) { rad.Send.ErrorMessage("Radia::Error065"); return;}
					}
				}
				else { rad.Send.ErrorMessage("Radia::Error065"); return;}
			}
			else 
			{ 
				//MLDisownSymbol(stdlink, FunName);
				MLWrapDeleteSymbol(stdlink, FunName); //OC091015
				rad.Send.ErrorMessage("Radia::Error065"); return;
			}
		}
	}

// General Subdivision Options
	const int AmOfOptions = 3;
	const char* OptionNames[] = {0,0,0};
	const char* OptionValues[] = {0,0,0};
	int OptionCount = 0;
	for(int ii=0; ii<AmOfOptions; ii++)
	{
		if(!StartReadingGeneralOptions)
		{
			Next = MLGetNext(stdlink);
			if(Next==MLTKFUNC)
			{
				ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
				if((!ReadOK) || strcmp(FunName, "Rule") || (ArgNum!=2)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
				//MLDisownSymbol(stdlink, FunName);
				MLWrapDeleteSymbol(stdlink, FunName); //OC091015

				if(!MLGetSymbol(stdlink, OptionNames + OptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
				if(!MLGetSymbol(stdlink, OptionValues + (OptionCount++))) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			}
			else { MLNewPacket(stdlink); break;}
		}
		else
		{
			if(!MLGetSymbol(stdlink, OptionNames + OptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			if(!MLGetSymbol(stdlink, OptionValues + (OptionCount++))) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			StartReadingGeneralOptions = 0;
		}
	}

	rad.SubdivideElement_g3d(ElemKey, SubdivArray, 6, TypeExtraSpec, ExtraSpec, LenExtraSpec, OptionNames, OptionValues, OptionCount);

	for(int k=0; k<AmOfOptions; k++) 
	{
		//if(OptionNames[k] != 0) MLDisownSymbol(stdlink, OptionNames[k]);
		if(OptionNames[k] != 0) MLWrapDeleteSymbol(stdlink, OptionNames[k]); //OC091015
		//if(OptionValues[k] != 0) MLDisownSymbol(stdlink, OptionValues[k]);
		if(OptionValues[k] != 0) MLWrapDeleteSymbol(stdlink, OptionValues[k]); //OC091015
	}
#endif
}

//-------------------------------------------------------------------------

void SubdivideElementG3DOpt(int ElemKey, double* SubdivArray, char TypeExtraSpec, double* ExtraSpec, int LenExtraSpec, const char* Opt1, const char* Opt2, const char* Opt3)
{
	const char* OptionNames[] = {0,0,0};
	const char* OptionValues[] = {0,0,0};
	int OptionCount = 0;

	char CharBuf1[200], CharBuf2[200], CharBuf3[200];
	if(Opt1 != 0)
	{
		if(*Opt1 != '\0') 
		{
			strcpy(CharBuf1, Opt1);
			char *pEndOptName = strrchr(CharBuf1, '-');
			if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			*pEndOptName = '\0';
			OptionNames[0] = CharBuf1;
			OptionValues[0] = strrchr(Opt1, '>') + 1;
			OptionCount++;
		}
	}
	if(Opt2 != 0)
	{
		if(*Opt2 != '\0') 
		{
			strcpy(CharBuf2, Opt2);
			char *pEndOptName = strrchr(CharBuf2, '-');
			if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			*pEndOptName = '\0';
			OptionNames[OptionCount] = CharBuf2;
			OptionValues[OptionCount] = strrchr(Opt2, '>') + 1;
			OptionCount++;
		}
	}
	if(Opt3 != 0)
	{
		if(*Opt3 != '\0') 
		{
			strcpy(CharBuf3, Opt3);
			char *pEndOptName = strrchr(CharBuf3, '-');
			if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			*pEndOptName = '\0';
			OptionNames[OptionCount] = CharBuf3;
			OptionValues[OptionCount] = strrchr(Opt3, '>') + 1;
			OptionCount++;
		}
	}

	rad.SubdivideElement_g3d(ElemKey, SubdivArray, 6, TypeExtraSpec, ExtraSpec, LenExtraSpec, OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

void CutElementG3D()
{
#ifdef __MATHEMATICA__
	int ElemKey;
	MLGetInteger(stdlink, &ElemKey);

	double PointOnPlane[3];
	double PlaneNormal[3];

	const char *FunName;
	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;
	int ReadOK = 0;

	int Next = MLGetNext(stdlink);
	if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	else { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if((!ReadOK) || strcmp(FunName, "List") || ArgNum!=3) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015
	for(int kk=0; kk<3; kk++)
	{
		Next = MLGetNext(stdlink);
		if((Next==MLTKINT) || (Next==MLTKREAL))
		{
			if(!MLGetDouble(stdlink, PointOnPlane + kk)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		}
		else { rad.Send.ErrorMessage("Radia::Error000"); return;}
	}

	Next = MLGetNext(stdlink);
	if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	else { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if((!ReadOK) || strcmp(FunName, "List") || ArgNum!=3) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015
	for(int jj=0; jj<3; jj++)
	{
		Next = MLGetNext(stdlink);
		if((Next==MLTKINT) || (Next==MLTKREAL))
		{
			if(!MLGetDouble(stdlink, PlaneNormal + jj)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		}
		else { rad.Send.ErrorMessage("Radia::Error000"); return;}
	}

// Optional parameters
	const int AmOfOptions = 2;
	const char* OptionNames[] = {0, 0};
	const char* OptionValues[] = {0, 0};
	int OptionCount = 0;
	for(int ii=0; ii<AmOfOptions; ii++)
	{
		Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC)
		{
			ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
			if((!ReadOK) || strcmp(FunName, "Rule") || (ArgNum!=2)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			if(!MLGetSymbol(stdlink, OptionNames + OptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			if(!MLGetSymbol(stdlink, OptionValues + (OptionCount++))) { rad.Send.ErrorMessage("Radia::Error062"); return;}

			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC091015
		}
		else { MLNewPacket(stdlink); break;}
	}

	rad.CutElement_g3d(ElemKey, PointOnPlane, 3, PlaneNormal, 3, OptionNames, OptionValues, OptionCount);

	for(int k=0; k<AmOfOptions; k++) 
	{
		//if(OptionNames[k] != 0) MLDisownSymbol(stdlink, OptionNames[k]);
		if(OptionNames[k] != 0) MLWrapDeleteSymbol(stdlink, OptionNames[k]); //OC091015
		//if(OptionValues[k] != 0) MLDisownSymbol(stdlink, OptionValues[k]);
		if(OptionValues[k] != 0) MLWrapDeleteSymbol(stdlink, OptionValues[k]); //OC091015
	}
#endif
}

//-------------------------------------------------------------------------

void CutElementG3DOpt1(int ElemKey, double x, double y, double z, double nx, double ny, double nz, const char* Opt)
{
	double PointOnPlane[] = {x,y,z};
	double PlaneNormal[] = {nx,ny,nz};

	char CharBuf[200];
	strcpy(CharBuf, Opt);
	char *pEndOptName = strrchr(CharBuf, '-');
	if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
	*pEndOptName = '\0';

	const char* OptionNames[] = {CharBuf};
	const char* OptionValues[] = {0};
	OptionValues[0] = strrchr(Opt, '>') + 1;

	rad.CutElement_g3d(ElemKey, PointOnPlane, 3, PlaneNormal, 3, OptionNames, OptionValues, 1);
}

void CutElementG3DOpt0(int ElemKey, double x, double y, double z, double nx, double ny, double nz)
{
	double PointOnPlane[] = {x,y,z};
	double PlaneNormal[] = {nx,ny,nz};

	const char* OptionNames[] = {0};
	const char* OptionValues[] = {0};

	rad.CutElement_g3d(ElemKey, PointOnPlane, 3, PlaneNormal, 3, OptionNames, OptionValues, 0);
}

void CutElementG3DOpt(int ElemKey, double x, double y, double z, double nx, double ny, double nz, const char* Opt1)
{
	double PointOnPlane[] = {x,y,z};
	double PlaneNormal[] = {nx,ny,nz};

	const char* OptionNames[] = {0};
	const char* OptionValues[] = {0};
	int OptionCount = 0;

	char CharBuf1[200];
	if(*Opt1 != '\0') 
	{
		strcpy(CharBuf1, Opt1);
		char *pEndOptName = strrchr(CharBuf1, '-');
		if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
		*pEndOptName = '\0';
		OptionNames[0] = CharBuf1;
		OptionValues[0] = strrchr(Opt1, '>') + 1;
		OptionCount++;
	}
	rad.CutElement_g3d(ElemKey, PointOnPlane, 3, PlaneNormal, 3, OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

void SubdivideElementG3DByParPlanes()
{
#ifdef __MATHEMATICA__
	int ElemKey;
	if(!MLGetInteger(stdlink, &ElemKey)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	double SubdivArray[15];

#if(MLVERSION >= 3)
	const char *FunName;
	const char *LabOrLocFrame;
#else
	char *FunName;
	char *LabOrLocFrame;
#endif

	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;
	double DoubleValue;

	int Next = MLGetNext(stdlink);
	int ReadOK = 0;
	if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	else { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if((!ReadOK) || strcmp(FunName, "List") || (ArgNum>3) || (ArgNum==0)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	//MLDisownSymbol(stdlink, FunName);
	MLWrapDeleteSymbol(stdlink, FunName); //OC091015

	int AmOfSubdivDirections = ArgNum;

	for(int kk=0; kk<AmOfSubdivDirections; kk++)
	{
		Next = MLGetNext(stdlink);
		if(Next!=MLTKFUNC) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		if((!ReadOK) || strcmp(FunName, "List") || (ArgNum!=3)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015

		Next = MLGetNext(stdlink);
		if(Next!=MLTKFUNC) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		if((!ReadOK) || strcmp(FunName, "List") || (ArgNum!=3)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		//MLDisownSymbol(stdlink, FunName);
		MLWrapDeleteSymbol(stdlink, FunName); //OC091015

		int kk_mu_5 = kk*5;
		for(int jj=0; jj<5; jj++)
		{
			Next = MLGetNext(stdlink);
			if((Next==MLTKINT) || (Next==MLTKREAL)) { if(MLGetDouble(stdlink, &DoubleValue)) SubdivArray[kk_mu_5 + jj] = DoubleValue;}
			else { rad.Send.ErrorMessage("Radia::Error000"); return;}
		}
	}
	if(!MLGetString(stdlink, &LabOrLocFrame)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	rad.SubdivideElement_g3dByParPlanes(ElemKey, SubdivArray, AmOfSubdivDirections, LabOrLocFrame);
	//MLDisownString(stdlink, LabOrLocFrame);
	MLWrapDeleteString(stdlink, LabOrLocFrame); //OC091015
#endif
}

//-------------------------------------------------------------------------

void GeometricalVolume(int ElemKey)
{
	rad.ComputeGeometricalVolume(ElemKey);
}

//-------------------------------------------------------------------------

void GeometricalLimits(int ElemKey)
{
	rad.ComputeGeometricalLimits(ElemKey);
}

//-------------------------------------------------------------------------

void FldCmpMetForSubdRecMag(int ElemKey, int Switch, int SubLevel)
{
	rad.FieldCompMethForSubdividedRecMag(ElemKey, Switch, SubLevel);
}

//-------------------------------------------------------------------------

void SetLocMgnInSbdRecMag() // May be removed
{
#ifdef __MATHEMATICA__
	int ElemKey;
	MLGetInteger(stdlink, &ElemKey);

	TVector3d* ArrayOfVectIndx = NULL;
	TVector3d* ArrayOfMagn = NULL;
	int Len;

	long* Dims;
	char** Heads;
	long Depth;
	double* DataPtr;

	int ArrayReadOK = MLGetDoubleArray(stdlink, &DataPtr, &Dims, &Heads, &Depth);
	if(!ArrayReadOK || Dims[2]!=3 || Dims[1]!=2 || Depth!=3) 
	{ 
		rad.Send.ErrorMessage("Radia::Error000"); return;
	}

	Len = Dims[0];
	//try
	//{
		ArrayOfVectIndx = new TVector3d[Len];
		ArrayOfMagn = new TVector3d[Len];
	//}
	//catch (radTException* radExceptionPtr)
	//{
	//	rad.Send.ErrorMessage(radExceptionPtr->what()); return;
	//}
	//catch (...)
	//{
	//	rad.Send.ErrorMessage("Radia::Error999"); return;
	//}

	for(int i=0; i<Len; i++)
	{
		int i_mu_6 = i*6;
		ArrayOfVectIndx[i] = TVector3d(DataPtr[i_mu_6], DataPtr[i_mu_6+1], DataPtr[i_mu_6+2]);
		ArrayOfMagn[i] = TVector3d(DataPtr[i_mu_6+3], DataPtr[i_mu_6+4], DataPtr[i_mu_6+5]);
	}
	//MLDisownDoubleArray(stdlink, DataPtr, Dims, Heads, Depth);
	MLWrapDeleteDoubleArray(stdlink, DataPtr, Dims, Heads, Depth); //OC091015

	rad.SetLocMgnInSbdRecMag(ElemKey, ArrayOfVectIndx, ArrayOfMagn, Len);

	delete[] ArrayOfVectIndx;
	delete[] ArrayOfMagn;
#endif
}

//-------------------------------------------------------------------------

void Translation(double vx, double vy, double vz)
{
	//double TranslArray[] = {radCR.Double(vx), radCR.Double(vy), radCR.Double(vz)};
	double TranslArray[] = {vx, vy, vz};
	rad.SetTranslation(TranslArray, 3);
}

//-------------------------------------------------------------------------

void Rotation(double xc, double yc, double zc, 
			  double vx, double vy, double vz, double Angle)
{
	//double PoiOnAxis[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	//double AxisVect[] = {radCR.Double(vx), radCR.Double(vy), radCR.Double(vz)};

	double PoiOnAxis[] = {xc, yc, zc};
	double AxisVect[] = {vx, vy, vz};

	rad.SetRotation(PoiOnAxis, 3, AxisVect, 3, Angle);
}

//-------------------------------------------------------------------------

void PlaneSym(double xc, double yc, double zc, 
			  double nx, double ny, double nz)
{
	//double PoiOnPlane[] = {radCR.Double(xc), radCR.Double(yc), radCR.Double(zc)};
	//double PlaneNormal[] = {radCR.Double(nx), radCR.Double(ny), radCR.Double(nz)};

	double PoiOnPlane[] = {xc, yc, zc};
	double PlaneNormal[] = {nx, ny, nz};

	rad.SetPlaneSym(PoiOnPlane, 3, PlaneNormal, 3, 1);
}

//-------------------------------------------------------------------------

void FieldInversion()
{
	rad.SetFieldInversion();
}

//-------------------------------------------------------------------------

void CombineTransformLeft(int ThisElemKey, int AnotherElemKey)
{
	rad.CombineTransformations(ThisElemKey, AnotherElemKey, 'L');
}

//-------------------------------------------------------------------------

void CombineTransformRight(int ThisElemKey, int AnotherElemKey)
{
	rad.CombineTransformations(ThisElemKey, AnotherElemKey, 'R');
}

//-------------------------------------------------------------------------

void ApplySymmetry(int g3dElemKey, int TransElemKey, int Multiplicity)
{
	rad.ApplySymmetry(g3dElemKey, TransElemKey, Multiplicity);
}

//-------------------------------------------------------------------------

void TransformObject(int g3dElemKey, int TransElemKey)
{
	rad.ApplySymmetry(g3dElemKey, TransElemKey, 1);	
}

//-------------------------------------------------------------------------

void LinearMaterial(double KsiPar, double KsiPer, double Mrx, double Mry, double Mrz)
{
	double KsiArray[] = {KsiPar, KsiPer};
	double RemMagnArray[] = {Mrx, Mry, Mrz};
	rad.SetLinearMaterial(KsiArray, 2, RemMagnArray, 3);
}

//-------------------------------------------------------------------------

void LinearMaterial2(double KsiPar, double KsiPer, double Mr)
{
	double KsiArray[] = {KsiPar, KsiPer};
	rad.SetLinearMaterial(KsiArray, 2, &Mr, 1);
}

//-------------------------------------------------------------------------

void MaterialStd(char* MatName, double Mr)
{
	rad.SetMaterialStd(MatName, Mr);
}

//-------------------------------------------------------------------------

void NonlinearIsotropMaterial(double Ms1, double Ms2, double Ms3, 
							  double ks1, double ks2, double ks3)
{
	double Ms[] = {Ms1, Ms2, Ms3};
	double ks[] = {ks1, ks2, ks3};
	rad.SetNonlinearIsotropMaterial(Ms, 3, ks, 3);
}

//-------------------------------------------------------------------------

void NonlinearIsotropMaterial2(double ks1, double Ms1, double ks2, double Ms2, double ks3, double Ms3)
{
	double Ms[] = {Ms1, Ms2, Ms3};
	double ks[] = {ks1, ks2, ks3};
	rad.SetNonlinearIsotropMaterial(Ms, 3, ks, 3);
}

//-------------------------------------------------------------------------

void NonlinearIsotropMaterial3()
{
	int lenArrayOfPoints2d;
	TVector2d* ArrayOfPoints2d = NULL;
	if(!rad.Send.GetArrayOfVector2d(ArrayOfPoints2d, lenArrayOfPoints2d)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	rad.SetNonlinearIsotropMaterial(ArrayOfPoints2d, lenArrayOfPoints2d);

	if(ArrayOfPoints2d != NULL) delete[] ArrayOfPoints2d;
}

//-------------------------------------------------------------------------

void NonlinearIsotropMaterial3Opt(double** HandM_Array, long LenHandM_Array)
{
	TVector2d* ArrayOfPoints2d = new TVector2d[LenHandM_Array];
	if(ArrayOfPoints2d == 0) { rad.Send.ErrorMessage("Radia::Error900"); return;}
	TVector2d* tArrayOfPoints2d = ArrayOfPoints2d;
	double** tHandM_Array = HandM_Array;

	for(long i=0; i<LenHandM_Array; i++)
	{
		tArrayOfPoints2d->x = **tHandM_Array;
		(tArrayOfPoints2d++)->y = (*tHandM_Array)[1];
		tHandM_Array++;
	}

	rad.SetNonlinearIsotropMaterial(ArrayOfPoints2d, (int)LenHandM_Array);

	if(ArrayOfPoints2d != 0) delete[] ArrayOfPoints2d;
}

//-------------------------------------------------------------------------

void NonlinearLaminatedMaterialML()
{
#ifdef __MATHEMATICA__

	//const char *FunName;
	//long ArgNum = 0;
	//int ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
    //if((!ReadOK) || strcmp(FunName, "List")) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	int lenArrayOfPoints2d;
	TVector2d* ArrayOfPoints2d = NULL;
	if(!rad.Send.GetArrayOfVector2d(ArrayOfPoints2d, lenArrayOfPoints2d)) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	double PackFactor = 1;
	if(!rad.Send.GetDouble(PackFactor)) { rad.Send.ErrorMessage("Radia::Error000"); return;}

    double *dN=0;
    long Len_dN=0;
	int Next = MLGetNext(stdlink);
	if(Next == MLTKFUNC) 
	{
        if(!rad.Send.GetArrayOfDouble(dN, Len_dN)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		if(Len_dN != 3) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	}

    rad.SetNonlinearLaminatedMaterial(ArrayOfPoints2d, lenArrayOfPoints2d, PackFactor, dN);

	if(ArrayOfPoints2d != NULL) delete[] ArrayOfPoints2d;
	if(dN != 0) delete[] dN;
#endif
}

//-------------------------------------------------------------------------

void NonlinearLaminatedMaterialFrm(double* pKsiMs1, double* pKsiMs2, double* pKsiMs3, double PackFactor, double* dN)
{
	int lenArrayOfPoints2d = 0;
	if(pKsiMs1 != 0) lenArrayOfPoints2d++;
	if(pKsiMs2 != 0) lenArrayOfPoints2d++;
	if(pKsiMs3 != 0) lenArrayOfPoints2d++;

	if(lenArrayOfPoints2d == 0) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	TVector2d* ArrayOfPoints2d = new TVector2d[lenArrayOfPoints2d];

	if(pKsiMs1 != 0)
	{
        ArrayOfPoints2d->x = pKsiMs1[0];
        ArrayOfPoints2d->y = pKsiMs1[1];
	}
	if(pKsiMs2 != 0)
	{
        ArrayOfPoints2d[1].x = pKsiMs2[0];
        ArrayOfPoints2d[1].y = pKsiMs2[1];
	}
	if(pKsiMs3 != 0)
	{
        ArrayOfPoints2d[2].x = pKsiMs3[0];
        ArrayOfPoints2d[2].y = pKsiMs3[1];
	}
    rad.SetNonlinearLaminatedMaterial(ArrayOfPoints2d, lenArrayOfPoints2d, PackFactor, dN);

	if(ArrayOfPoints2d != NULL) delete[] ArrayOfPoints2d;
}

//-------------------------------------------------------------------------

void NonlinearLaminatedMaterialTab(double* pFlatMatDef, int AmOfMatPts, double PackFactor, double* dN)
{
	if((pFlatMatDef == 0) || (AmOfMatPts <= 3))
	{
        rad.Send.ErrorMessage("Radia::Error088");
	}
	TVector2d* ArrayOfPoints2d = new TVector2d[AmOfMatPts];

	TVector2d *tArrayOfPoints2d = ArrayOfPoints2d;
	double *tFlatMatDef = pFlatMatDef;
	for(int i=0; i<AmOfMatPts; i++)
	{
		tArrayOfPoints2d->x = *(tFlatMatDef++);
		(tArrayOfPoints2d++)->y = *(tFlatMatDef++);
	}
    rad.SetNonlinearLaminatedMaterial(ArrayOfPoints2d, AmOfMatPts, PackFactor, dN);

	if(ArrayOfPoints2d != NULL) delete[] ArrayOfPoints2d;
}

//-------------------------------------------------------------------------

void NonlinearAnisotropMaterial()
{
#ifdef __MATHEMATICA__
	double KsiPar[4], KsiPer[4], MsPar[3], MsPer[3], Hci[4]; //Hc[2];
	double* Ksi[] = {KsiPar, KsiPer};
	double* Ms[] = {MsPar, MsPer};

	char DependenceIsNonlinear[2];

	for(int i=0; i<4; i++) Hci[i] = 0;

	const char *FunName;
	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;
	int Next, ReadOK;

	for(int k=0; k<2; k++)
	{
		double *CurKsi = Ksi[k], *CurMs = Ms[k], *CurHci = Hci;
		Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC) //{{ksi1,ms1,hc1},{ksi2,ms2,hc2},{ksi3,ms3,hc3},{ksi0,hc0}} or {{ksi1,ms1},{ksi2,ms2},{ksi3,ms3},ksi0,hc:0}
		{
			DependenceIsNonlinear[k] = 1;
			//long PartArgNum; //OC240907 (port to math6)
			int PartArgNum;
			ReadOK = MLGetFunction(stdlink, &FunName, &PartArgNum);
			if((!ReadOK) || strcmp(FunName, "List") || (!(PartArgNum==4 || PartArgNum==5))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC091015

			for(int j=0; j<3; j++)
			{
				Next = MLGetNext(stdlink);
				if(Next!=MLTKFUNC) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
				if((!ReadOK) || strcmp(FunName, "List") || ((ArgNum!=2) && (ArgNum!=3))) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				//MLDisownSymbol(stdlink, FunName);
				MLWrapDeleteSymbol(stdlink, FunName); //OC091015

				Next = MLGetNext(stdlink);
				if(!(Next==MLTKINT || Next==MLTKREAL)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				if(!MLGetDouble(stdlink, CurKsi++)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				Next = MLGetNext(stdlink);
				if(!(Next==MLTKINT || Next==MLTKREAL)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				if(!MLGetDouble(stdlink, CurMs++)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				if(ArgNum==3)
				{
					if(k==0)
					{
						Next = MLGetNext(stdlink);
						if(!(Next==MLTKINT || Next==MLTKREAL)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
						if(!MLGetDouble(stdlink, CurHci++)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
					}
					else { rad.Send.ErrorMessage("Radia::Error000"); return;}
				}
			}

			Next = MLGetNext(stdlink);
			if(Next==MLTKFUNC)
			{
				ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
				if((!ReadOK) || strcmp(FunName, "List") || (ArgNum!=2)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				//MLDisownSymbol(stdlink, FunName);
				MLWrapDeleteSymbol(stdlink, FunName); //OC091015

				Next = MLGetNext(stdlink);
				if(!(Next==MLTKINT || Next==MLTKREAL)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				if(!MLGetDouble(stdlink, CurKsi)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				Next = MLGetNext(stdlink);
				if(!(Next==MLTKINT || Next==MLTKREAL)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				if(!MLGetDouble(stdlink, CurHci)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
			}
			else
			{
				if(!(Next==MLTKINT || Next==MLTKREAL)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				if(!MLGetDouble(stdlink, CurKsi)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
			}

			if(PartArgNum==5)
			{
				Next = MLGetNext(stdlink);
				if(!(Next==MLTKINT || Next==MLTKREAL)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				//if(!MLGetDouble(stdlink, &Hc[k])) { rad.Send.ErrorMessage("Radia::Error000"); return;}
				if(!MLGetDouble(stdlink, CurHci)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
			}
			//else Hc[k] = 0.;
			else 
			{
				if(k == 0) *CurHci = 0.;
			}
		}
		else if((Next==MLTKINT) || (Next==MLTKREAL)) //ksi0
		{
			DependenceIsNonlinear[k] = 0;
			if(!MLGetDouble(stdlink, CurKsi)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		}
		else { rad.Send.ErrorMessage("Radia::Error000"); return;}
	}
	//rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, DependenceIsNonlinear);
	rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hci, 4, DependenceIsNonlinear);
#endif
}

//-------------------------------------------------------------------------

void NonlinearAnisotropMaterialOpt0(double* pDataPar, int lenDataPar, double* pDataPer, int lenDataPer)
{
	rad.SetNonlinearAnisotropMaterial0(pDataPar, lenDataPar, pDataPer, lenDataPer);
}

//-------------------------------------------------------------------------

void NonlinearAnisotropMaterialOpt1(double** Par, double** Per)
{
	double KsiPar[4], KsiPer[4], MsPar[3], MsPer[3], Hc[2];
	double* Ksi[] = {KsiPar, KsiPer};
	double* Ms[] = {MsPar, MsPer};

	char DependenceIsNonlinear[] = {1,1};

	double **tPar = Par, **tPer = Per;
	for(int i=0; i<3; i++)
	{
		KsiPar[i] = (*tPar)[0]; MsPar[i] = (*tPar)[1];
		KsiPer[i] = (*tPer)[0]; MsPer[i] = (*tPer)[1];
		tPar++; tPer++;
	}
	KsiPar[3] = (*tPar)[0]; Hc[0] = (*tPar)[1];
	KsiPer[3] = (*tPer)[0]; Hc[1] = 0; //Hc[1] = (*tPer)[1];

	//rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, DependenceIsNonlinear);
	rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, 2, DependenceIsNonlinear);
}

//-------------------------------------------------------------------------

void NonlinearAnisotropMaterialOpt2(double** Par, double Per)
{
	double KsiPar[4], KsiPer[4], MsPar[3], MsPer[3], Hc[2];
	double* Ksi[] = {KsiPar, KsiPer};
	double* Ms[] = {MsPar, MsPer};

	char DependenceIsNonlinear[] = {1,0};

	double **tPar = Par;
	for(int i=0; i<3; i++)
	{
		KsiPar[i] = (*tPar)[0]; MsPar[i] = (*tPar)[1];
		tPar++;
	}
	KsiPar[3] = (*tPar)[0]; Hc[0] = (*tPar)[1];
	KsiPer[0] = Per; Hc[1] = 0;

	//rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, DependenceIsNonlinear);
	rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, 2, DependenceIsNonlinear);
}

//-------------------------------------------------------------------------

void NonlinearAnisotropMaterialOpt3(double Par, double** Per)
{
	double KsiPar[4], KsiPer[4], MsPar[3], MsPer[3], Hc[2];
	double* Ksi[] = {KsiPar, KsiPer};
	double* Ms[] = {MsPar, MsPer};

	char DependenceIsNonlinear[] = {0,1};

	double **tPer = Per;
	for(int i=0; i<3; i++)
	{
		KsiPer[i] = (*tPer)[0]; MsPer[i] = (*tPer)[1];
		tPer++;
	}
	KsiPar[0] = Par; Hc[0] = 0;
	KsiPer[3] = (*tPer)[0]; Hc[1] = 0; //Hc[1] = (*tPer)[1];

	//rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, DependenceIsNonlinear);
	rad.SetNonlinearAnisotropMaterial(Ksi, Ms, Hc, 2, DependenceIsNonlinear);
}

//-------------------------------------------------------------------------

void ApplyMaterial(int g3dRelaxElemKey, int MaterElemKey)
{
	rad.ApplyMaterial(g3dRelaxElemKey, MaterElemKey);
}

//-------------------------------------------------------------------------

void MvsH(int g3dRelaxOrMaterElemKey, char* MagnChar, double hx, double hy, double hz)
{
	double H[] = {hx, hy, hz};
	rad.ComputeMvsH(g3dRelaxOrMaterElemKey, MagnChar, H, 3);
}

//-------------------------------------------------------------------------

void PreRelax(int ElemKey, int SrcElemKey)
{
	rad.PreRelax(ElemKey, SrcElemKey);
}

//-------------------------------------------------------------------------

void ShowInteractMatrix(int InteractElemKey)
{
	rad.ShowInteractMatrix(InteractElemKey);
}

//-------------------------------------------------------------------------

void ShowInteractVector(int InteractElemKey, char* FieldVectID)
{
	rad.ShowInteractVector(InteractElemKey, FieldVectID);
}

//-------------------------------------------------------------------------

void ManualRelax(int InteractElemKey, int MethNo, int IterNumber, double RelaxParam)
{
	rad.MakeManualRelax(InteractElemKey, MethNo, IterNumber, RelaxParam);
}

//-------------------------------------------------------------------------

void AutoRelax()
{
#ifdef __MATHEMATICA__

	int InteractElemKey, MaxIterNumber, MethNo = 4;
	if(!MLGetInteger(stdlink, &InteractElemKey)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	double PrecOnMagnetiz, dMethNo = 0;
	if(!MLGetDouble(stdlink, &PrecOnMagnetiz)) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	if(!MLGetInteger(stdlink, &MaxIterNumber)) { rad.Send.ErrorMessage("Radia::Error000"); return;}

	int numOptions = 0;
	const char* arOptionNames[] = {0};
	const char* arOptionValues[] = {0};
	int parCount = 0, auxInt = -1, Next = 0;
	for(int i=0; i<2; i++)
	{//eventual absence of MethNo, Opt
		if(Next <= 0) Next = MLGetNext(stdlink);		
		if((Next==MLTKINT) || (Next==MLTKREAL))
		{
			if(!MLGetInteger(stdlink, &auxInt)) 
			{
				if(MLGetReal(stdlink, &dMethNo)) auxInt = (int)(dMethNo + 0.00001);
				else { rad.Send.ErrorMessage("Radia::Error000"); return;}
			}
			if(auxInt >= 0) MethNo = auxInt;
			parCount += 1;
			Next = 0;
		}
		else
		{
			numOptions = 1; //it will be eventually modified in AuxGetMathematicaOptions
			//if(AuxGetMathematicaOptions(Next, arOptionNames, arOptionValues, numOptions)) parCount += numOptions;
			if(!AuxGetMathematicaOptions(Next, arOptionNames, arOptionValues, numOptions)) return; //OC240612
			if(numOptions > 0) parCount += numOptions; //OC240612
		}

		if(parCount == 0) { rad.Send.ErrorMessage("Radia::Error000"); return;}
		else if(parCount >= 2) break;
	}

	rad.MakeAutoRelax(InteractElemKey, PrecOnMagnetiz, MaxIterNumber, MethNo, arOptionNames, arOptionValues, numOptions);

	for(int k=0; k<numOptions; k++) 
	{
		//if(arOptionNames[k] != 0) MLDisownSymbol(stdlink, arOptionNames[k]);
		if(arOptionNames[k] != 0) MLWrapDeleteSymbol(stdlink, arOptionNames[k]); //OC091015
		//if(arOptionValues[k] != 0) MLDisownSymbol(stdlink, arOptionValues[k]);
		if(arOptionValues[k] != 0) MLWrapDeleteSymbol(stdlink, arOptionValues[k]); //OC091015
	}

#endif
}

//-------------------------------------------------------------------------

void AutoRelaxOpt(int InteractElemKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo, const char* Opt1)
{
	char CharBuf1[200];
	const char* OptionNames[] = {CharBuf1};
	const char* OptionValues[] = {0};
	const char* NonParsedOpts[] = {Opt1};
	int OptionCount = 1;
	AuxParseOptionNamesAndValues(NonParsedOpts, OptionNames, OptionValues, OptionCount);

	rad.MakeAutoRelax(InteractElemKey, PrecOnMagnetiz, MaxIterNumber, MethNo, OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

void UpdateSourcesForRelax(int InteractElemKey)
{
	rad.UpdateSourcesForRelax(InteractElemKey);
}

//-------------------------------------------------------------------------

void SolveGen(int ObjKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo)
{
	if(MethNo == 0) MethNo = 4; //Default method
	rad.SolveGen(ObjKey, PrecOnMagnetiz, MaxIterNumber, MethNo);
}

//-------------------------------------------------------------------------

void CompCriterium(double InAbsPrecB, double InAbsPrecA, double InAbsPrecB_int, 
				   double InAbsPrecFrc, double InAbsPrecTrjCoord, double InAbsPrecTrjAngle)
{
	rad.SetCompCriterium(InAbsPrecB, InAbsPrecA, InAbsPrecB_int, InAbsPrecFrc, InAbsPrecTrjCoord, InAbsPrecTrjAngle);
}

//-------------------------------------------------------------------------

void CompPrecision()
{
#ifdef __MATHEMATICA__
	const char* ValNames[] = {0,0,0,0,0,0,0,0,0};
	double Values[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

	const char *FunName;
	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;
	int Next, ReadOK;

	int ValCount = 0;
	for(int i=0; i<9; i++)
	{
		Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC)
		{
			ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
			if((!ReadOK) || strcmp(FunName, "Rule") || (ArgNum!=2)) { rad.Send.ErrorMessage("Radia::Error057"); return;}

			if(!MLGetSymbol(stdlink, &(ValNames[ValCount]))) { rad.Send.ErrorMessage("Radia::Error057"); return;}
			if(!MLGetDouble(stdlink, &Values[ValCount++])) { rad.Send.ErrorMessage("Radia::Error057"); return;}
			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC091015
		}
		else { MLNewPacket(stdlink); break;}
	}
	rad.SetCompPrecisions(ValNames, Values, ValCount);

	//for(int k=0; k<9; k++) if(ValNames[k] != 0) MLDisownSymbol(stdlink, ValNames[k]);
	for(int k=0; k<9; k++) if(ValNames[k] != 0) MLWrapDeleteSymbol(stdlink, ValNames[k]); //OC091015
#endif
}

//-------------------------------------------------------------------------

void CompPrecisionOpt(const char* Opt1, const char* Opt2, const char* Opt3, const char* Opt4, const char* Opt5, const char* Opt6, const char* Opt7, const char* Opt8)
//void CompPrecisionOpt(const char* Opt1, const char* Opt2, const char* Opt3, const char* Opt4, const char* Opt5, const char* Opt6, const char* Opt7, const char* Opt8, const char* Opt9)
{
	const char* OptionNames[] = {0,0,0,0,0,0,0,0,0};
	double OptionValues[] = {0,0,0,0,0,0,0,0,0};
	int OptionCount = 0;

	int MaxAmOfOptions = 8;
	char CharBuf1[200], CharBuf2[200], CharBuf3[200], CharBuf4[200], CharBuf5[200], CharBuf6[200], CharBuf7[200], CharBuf8[200];
	char *TotCharBuf[] = {CharBuf1, CharBuf2, CharBuf3, CharBuf4, CharBuf5, CharBuf6, CharBuf7, CharBuf8};
	const char *InOpt[] = {Opt1, Opt2,  Opt3, Opt4, Opt5, Opt6, Opt7, Opt8};
	for(int i=0; i<MaxAmOfOptions; i++)
	{
		if(InOpt[i] != 0) //OC240108
		{
			if(*(InOpt[i]) != '\0') 
			{
				strcpy(TotCharBuf[i], InOpt[i]);
				char *pEndOptName = strrchr(TotCharBuf[i], '-');
				if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
				*pEndOptName = '\0';
				OptionNames[i] = TotCharBuf[i];
				OptionValues[i] = atof(strrchr(InOpt[i], '>') + 1);
				if(OptionValues[i] == 0) { rad.Send.ErrorMessage("Radia::Error057"); return;}
				OptionCount++;
			}
		}
	}
	rad.SetCompPrecisions(OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

// Maybe to be removed later
void MultipoleThresholds(double a0, double a1, double a2, double a3)
{
	double MltplThresh[] = {a0, a1, a2, a3};
	rad.SetMltplThresh(MltplThresh); 
}

//-------------------------------------------------------------------------

void Field(int ElemKey, char* FieldChar, double x1, double y1, double z1, 
		   double x2, double y2, double z2, int Np, char* ShowArgFlag, double StrtArg)
{
	double StObsPoi[] = {radCR.Double(x1), radCR.Double(y1), radCR.Double(z1)};
	double FiObsPoi[] = {radCR.Double(x2), radCR.Double(y2), radCR.Double(z2)};
	rad.ComputeField(ElemKey, FieldChar, StObsPoi, 3, FiObsPoi, 3, Np, ShowArgFlag, radCR.Double(StrtArg));
}

//-------------------------------------------------------------------------
/**
void FieldArbitraryPointsStruct()
{
	int ElemKey;
	if(!rad.Send.GetInteger(ElemKey)) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	char* FieldChar;
	if(!rad.Send.GetString((const char*&)FieldChar)) { rad.Send.ErrorMessage("Radia::Error000"); return;};
	//if(!rad.Send.GetString(FieldChar)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	radTVectorOfVector3d VectorOfVector3d;
	radTVectInputCell VectInputCell;
	//if(!rad.Send.GetArbitraryListOfVector3d(VectorOfVector3d, VectInputCell)) { rad.Send.ErrorMessage("Radia::Error000"); return;};

	int resListRead = rad.Send.GetArbitraryListOfVector3d(VectorOfVector3d, VectInputCell);
	if(!resListRead) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	else if(resListRead == 2) return; //OC29092010 calculation should not be performed because symbol was sent (?)
	//else if(resListRead == 2) //OCTEST
	//{
	//	MLPutSymbol(stdlink, "x");
	//	return; //OC29092010 calculation should not be performed because symbol was sent (?)
	//}

	for(int i=0; i<(int)(VectorOfVector3d.size()); i++) //OC050504
	{
	TVector3d vSmallRand(radCR.Double(0), radCR.Double(0), radCR.Double(0));
	VectorOfVector3d[i] += vSmallRand;
	}

	rad.ComputeField(ElemKey, FieldChar, VectorOfVector3d, VectInputCell);
	rad.Send.DisownString(FieldChar);
}
**/

//This version is derived from the one of D. Hidas (18/04/2016).
//It works with Math. v.10 Plot[] (that eventually sends a symbolic variable to radFld).
//Attempts to place "getting" of a List of Points to a separate function in radTSend 
//lead to re-appearing the 
void FieldArbitraryPointsStruct(int ElemKey, char* FieldChar) 
{
#ifdef __MATHEMATICA__
	radTVectorOfVector3d VectorOfVector3d;
	radTVectInputCell VectInputCell;
	TVector3d vec3d;
	const char *FunName, *locFunName;
	int ArgNum, locArgNum;

	if(!MLGetFunction(stdlink, &FunName, &ArgNum)) { rad.Send.ErrorMessage("Radia::Error000"); return; }
	if(strcmp(FunName, "List")) { rad.Send.ErrorMessage("Radia::Error000"); return; }
	MLWrapDeleteSymbol(stdlink, FunName);

	int Next = MLGetNext(stdlink);
	if((ArgNum == 3) && ((Next == MLTKINT) || (Next == MLTKREAL) || (Next == MLTKSYM)))
	{
		if(!MLGetReal(stdlink, &(vec3d.x))) { rad.Send.ErrorMessage("Radia::Error000"); return; }
		if(!MLGetReal(stdlink, &(vec3d.y))) { rad.Send.ErrorMessage("Radia::Error000"); return; }
		if(!MLGetReal(stdlink, &(vec3d.z))) { rad.Send.ErrorMessage("Radia::Error000"); return; }

		VectorOfVector3d.push_back(vec3d);
		VectInputCell.push_back(radTInputCell('P', 3));
	}
	else 
	{
		VectInputCell.push_back(radTInputCell('L', ArgNum));

		for(int i = 0; i < ArgNum; ++i) 
		{
			if(!MLGetFunction(stdlink, &locFunName, &locArgNum)) { rad.Send.ErrorMessage("Radia::Error000"); return; }
			if(strcmp(locFunName, "List")) { rad.Send.ErrorMessage("Radia::Error000"); return; }
			MLWrapDeleteSymbol(stdlink, locFunName);

			int locNext = MLGetNext(stdlink);
			if((locArgNum == 3) && ((locNext == MLTKINT) || (locNext == MLTKREAL) || (locNext == MLTKSYM)))
			{
				if(!MLGetReal(stdlink, &(vec3d.x))) { rad.Send.ErrorMessage("Radia::Error000"); return; }
				if(!MLGetReal(stdlink, &(vec3d.y))) { rad.Send.ErrorMessage("Radia::Error000"); return; }
				if(!MLGetReal(stdlink, &(vec3d.z))) { rad.Send.ErrorMessage("Radia::Error000"); return; }

				VectorOfVector3d.push_back(vec3d);
				VectInputCell.push_back(radTInputCell('P', 3));
			}
		}
	}

	for(int i = 0; i<(int)(VectorOfVector3d.size()); i++)
	{
		TVector3d vSmallRand(radCR.Double(0), radCR.Double(0), radCR.Double(0));
		VectorOfVector3d[i] += vSmallRand;
	}

	rad.ComputeField(ElemKey, FieldChar, VectorOfVector3d, VectInputCell);
#endif
}

//-------------------------------------------------------------------------

void FieldArbitraryPointsArray(long ElemKey, const char* FieldChar, double** Points, long LenPoints)
{
	char LocStr[50];
	strcpy(LocStr, FieldChar);
	rad.ComputeField((int)ElemKey, LocStr, Points, LenPoints);
}

//-------------------------------------------------------------------------

void FieldInt(int ElemKey, char* IntID, char* FieldIntChar, double x1, double y1, double z1, double x2, double y2, double z2)
{
	double StPoi[] = {radCR.Double(x1), radCR.Double(y1), radCR.Double(z1)};
	double FiPoi[] = {radCR.Double(x2), radCR.Double(y2), radCR.Double(z2)};
	rad.ComputeFieldInt(ElemKey, IntID, FieldIntChar, StPoi, 3, FiPoi, 3);
}

//-------------------------------------------------------------------------

void FieldForce(int ElemKey, int ShapeElemKey)
{
	rad.ComputeFieldForce(ElemKey, ShapeElemKey);
}

//-------------------------------------------------------------------------

void FieldEnergy(int DestElemKey, int SourceElemKey, int kx, int ky, int kz)
{
	int SubdArray[] = {kx, ky, kz};
	rad.ComputeFieldEnergy(DestElemKey, SourceElemKey, SubdArray, 3);
}

//-------------------------------------------------------------------------

void FieldForceThroughEnergy(int DestElemKey, int SourceElemKey, char* ForceComponID, int kx, int ky, int kz)
{
	int SubdArray[] = {kx, ky, kz};
	rad.ComputeFieldForceThroughEnergy(DestElemKey, SourceElemKey, ForceComponID, SubdArray, 3);
}

//-------------------------------------------------------------------------

void FieldTorqueThroughEnergy(int DestElemKey, int SourceElemKey, char* TorqueComponID, double x0, double y0, double z0, int kx, int ky, int kz)
{
	int SubdArray[] = {kx, ky, kz};
	double TorqueCenPo[] = {x0, y0, z0};
	rad.ComputeFieldTorqueThroughEnergy(DestElemKey, SourceElemKey, TorqueComponID, SubdArray, 3, TorqueCenPo, 3);
}

//-------------------------------------------------------------------------

void ParticleTrajectory(int ElemKey, double E, double x0, double dxdy0, double z0, double dzdy0, double y0, double y1, int Np)
{
	rad.ComputeParticleTrajectory(ElemKey, E, radCR.Double(x0), dxdy0, radCR.Double(z0), dzdy0, radCR.Double(y0), radCR.Double(y1), Np);
}

//-------------------------------------------------------------------------

void FocusingPotential(int ElemKey, double x1, double y1, double z1, double x2, double y2, double z2, int Np)
{
	double StPoi[] = {radCR.Double(x1), radCR.Double(y1), radCR.Double(z1)};
	double FiPoi[] = {radCR.Double(x2), radCR.Double(y2), radCR.Double(z2)};
	rad.ComputeFocusPotent(ElemKey, StPoi, 3, FiPoi, 3, Np);
}

//-------------------------------------------------------------------------

void FocusingKickPer(int ElemKey, double x1, double y1, double z1, double nsx, double nsy, double nsz, double per, double nper, double n1x, double n1y, double n1z, double r1, int np1, double r2, int np2, const char* Comment, int nharm, int ns, double d1, double d2, const char* strKickUnit, double energyGeV, const char* strOutFormat)
{
	double P1[] = {radCR.Double(x1), radCR.Double(y1), radCR.Double(z1)};
	double Nlong[] = {radCR.Double(nsx), radCR.Double(nsy), radCR.Double(nsz)};
	double N1[] = {radCR.Double(n1x), radCR.Double(n1y), radCR.Double(n1z)};
	rad.ComputeFocusKickPer(ElemKey, P1, Nlong, per, nper, N1, r1, np1, r2, np2, Comment, nharm, ns, d1, d2, strKickUnit, energyGeV, strOutFormat);
}

//-------------------------------------------------------------------------

void FocusingKickPerFormStrRep(double* pKickData1, double* pKickData2, double* pBtE2Int, double* pCoordDir1, double* pCoordDir2, int np1, int np2, double per, int nper, const char* Comment)
{
	rad.ComposeFocusKickPerFormStrRep(pKickData1, pKickData2, pBtE2Int, pCoordDir1, pCoordDir2, np1, np2, per, nper, Comment);
}

//-------------------------------------------------------------------------

//void FocusingKick(int ElemKey, double x1,double y1,double z1, double nsx,double nsy,double nsz, double* dsArr,long Len_dsArr,int ns, double n1x,double n1y,double n1z, double r1,int np1,double r2,int np2, const char* Comment, double d1,double d2)
void FocusingKickML()
{
#ifdef __MATHEMATICA__

	int ElemKey=0;
	if(!rad.Send.GetInteger(ElemKey)) return;

    double *pAuxData=0;
	long lenAuxData=0;

	if(!rad.Send.GetArrayOfDouble(pAuxData, lenAuxData)) return;
	if(lenAuxData < 3) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	double P1[] = {radCR.Double(pAuxData[0]), radCR.Double(pAuxData[1]), radCR.Double(pAuxData[2])};
	delete[] pAuxData; pAuxData = 0; lenAuxData = 0;

	if(!rad.Send.GetArrayOfDouble(pAuxData, lenAuxData)) return;
	if(lenAuxData < 3) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	double Nlong[] = {radCR.Double(pAuxData[0]), radCR.Double(pAuxData[1]), radCR.Double(pAuxData[2])};
	delete[] pAuxData; pAuxData = 0; lenAuxData = 0;

    double *ArrLongDist=0;
	long lenArrLongDist=0;
	if(!rad.Send.GetArrayOfDouble(ArrLongDist, lenArrLongDist)) return;

	int ns;
	if(!rad.Send.GetInteger(ns)) return;

	if(!rad.Send.GetArrayOfDouble(pAuxData, lenAuxData)) return;
	if(lenAuxData < 3) { rad.Send.ErrorMessage("Radia::Error000"); return;}
	double Ntr1[] = {radCR.Double(pAuxData[0]), radCR.Double(pAuxData[1]), radCR.Double(pAuxData[2])};
	delete[] pAuxData; pAuxData = 0; lenAuxData = 0;

	double r1=0, r2=0;
	int np1=0, np2=0;
	if(!rad.Send.GetDouble(r1)) return;
	if(!rad.Send.GetInteger(np1)) return;
	if(!rad.Send.GetDouble(r2)) return;
	if(!rad.Send.GetInteger(np2)) return;

	char *pCommentStr=0;
	//if(!rad.Send.GetString(pCommentStr)) return;
	if(!rad.Send.GetString((const char*&)pCommentStr)) return;

	double d1=0, d2=0;
	if(!rad.Send.GetDouble(d1)) return;
	if(!rad.Send.GetDouble(d2)) return;

	rad.ComputeFocusKick(ElemKey, P1, Nlong, ArrLongDist, (int)lenArrLongDist, ns, Ntr1, r1, np1, r2, np2, pCommentStr, d1, d2);

	delete[] ArrLongDist; ArrLongDist = 0; lenArrLongDist = 0;
	rad.Send.DisownString(pCommentStr);

#endif
}

//-------------------------------------------------------------------------

void ShimSignature(int ElemKey, char* FldID, double vx, double vy, double vz, double x1, double y1, double z1, double x2, double y2, double z2, int Np, double vix, double viy, double viz)
{
	double V[] = {vx, vy, vz};
	double StPoi[] = {radCR.Double(x1), radCR.Double(y1), radCR.Double(z1)};
	double FiPoi[] = {radCR.Double(x2), radCR.Double(y2), radCR.Double(z2)};
	double Vi[] = {radCR.Double(vix), radCR.Double(viy), radCR.Double(viz)};
	rad.ComputeShimSignature(ElemKey, FldID, V, StPoi, FiPoi, Np, Vi);
}

//-------------------------------------------------------------------------

void TolForConvergence(double AbsRandMagnitude, double RelRandMagnitude, double ZeroRandMagnitude)
{
	rad.SetTolForConvergence(AbsRandMagnitude, RelRandMagnitude, ZeroRandMagnitude);
}

//-------------------------------------------------------------------------

void RandomizationOnOrOff(char* OnOrOff)
{
	rad.RandomizationOnOrOff(OnOrOff);
}

//-------------------------------------------------------------------------

void PhysicalUnits()
{
	rad.SetAndShowPhysUnits();
}

//-------------------------------------------------------------------------

//void DumpElem(int ElemKey, char* strFormat)
//{
//	rad.DumpElem(ElemKey, strFormat);
//}

void DumpElem()
{
#ifdef __MATHEMATICA__ //OC070613

	int *arKeys = 0;
	bool arKeysAllocInMathLink = false;
	long nElem = 0;
	int next = MLGetNext(stdlink);
	if(next == MLTKFUNC)
	{
		//if(!MLGetIntegerList(stdlink, &arKeys, &nElem)) { rad.Send.ErrorMessage("Radia::Error000"); MLDisownIntegerList(stdlink, arKeys, nElem); return;}
		if(!MLGetIntegerList(stdlink, &arKeys, &nElem)) { rad.Send.ErrorMessage("Radia::Error000"); MLWrapDeleteIntegerList(stdlink, arKeys, nElem); return;} //OC091015
		arKeysAllocInMathLink = true;
	}
	else if((next == MLTKREAL) || (next == MLTKINT))
	{
		arKeys = new int[1];
		if(!MLGetInteger(stdlink, arKeys)) { rad.Send.ErrorMessage("Radia::Error000"); if(arKeys != 0) delete[] arKeys; return;}
		nElem = 1;
	}

	const char *strFormat = 0;
	next = MLGetNext(stdlink);
	if(next == MLTKSTR)
	{
		if(!MLGetString(stdlink, &strFormat)) 
		{ 
			rad.Send.ErrorMessage("Radia::Error000"); 
			if(arKeys != 0) 
			{
				//if(nElem > 1) MLDisownIntegerList(stdlink, arKeys, nElem); 
				//if(arKeysAllocInMathLink) MLDisownIntegerList(stdlink, arKeys, nElem); 
				if(arKeysAllocInMathLink) MLWrapDeleteIntegerList(stdlink, arKeys, nElem); //OC091015
				else delete[] arKeys;
			}
			//if(strFormat != 0) MLDisownString(stdlink, strFormat); return;
			if(strFormat != 0) MLWrapDeleteString(stdlink, strFormat); return; //OC091015
		}
	}
	else 
	{
		if(arKeys != 0) 
		{
			//if(nElem > 1) MLDisownIntegerList(stdlink, arKeys, nElem); 
			//if(arKeysAllocInMathLink) MLDisownIntegerList(stdlink, arKeys, nElem); 
			if(arKeysAllocInMathLink) MLWrapDeleteIntegerList(stdlink, arKeys, nElem); //OC091015
			else delete[] arKeys;
		}
	}

	rad.DumpElem(arKeys, (int)nElem, strFormat, arKeysAllocInMathLink);

	if(arKeys != 0) 
	{
		//if(nElem > 1) MLDisownIntegerList(stdlink, arKeys, nElem); 
		//if(arKeysAllocInMathLink) MLDisownIntegerList(stdlink, arKeys, nElem); 
		if(arKeysAllocInMathLink) MLWrapDeleteIntegerList(stdlink, arKeys, nElem); //OC091015
		else delete[] arKeys;
	}
	//if(strFormat != 0) MLDisownString(stdlink, strFormat);
	if(strFormat != 0) MLWrapDeleteString(stdlink, strFormat); //OC091015
#endif
}

//-------------------------------------------------------------------------

void DumpElemOpt(int* arKeys, int nKeys, const char* AscOrBin)
{
	rad.DumpElem(arKeys, nKeys, AscOrBin, true); //?
}

//-------------------------------------------------------------------------

void DumpElemParse()
{
#ifdef __MATHEMATICA__ //OC070613

	const unsigned char *bstr=0;
	int bstrLen=0;
	if((!MLGetByteString(stdlink, &bstr, &bstrLen, 0)) || (bstr == 0) || (bstrLen == 0))
	{
		rad.Send.ErrorMessage("Radia::Error500"); return;        
	}

	rad.DumpElemParse(bstr, bstrLen);
	
	if((bstr != 0) && (bstrLen > 0)) MLReleaseByteString(stdlink, bstr, bstrLen);
#endif
}

//-------------------------------------------------------------------------

void DumpElemParseOpt(const unsigned char* sBytes, int nBytes)
{
	rad.DumpElemParse(sBytes, nBytes); //?
}

//-------------------------------------------------------------------------

void GenDump()
{
	rad.GenDump();
}

//-------------------------------------------------------------------------

void GraphicsForElemWithoutSymChilds(int ElemKey)
{
	rad.GraphicsForElem_g3d(ElemKey, 0);
}

//-------------------------------------------------------------------------

//void GraphicsForElemWithSymChilds(int ElemKey)
//{
//	rad.GraphicsForElem_g3d(ElemKey, 1);
//}

void GraphicsForElemWithSymChildsExt()
{
#ifdef __MATHEMATICA__
	int ElemKey;
	MLGetInteger(stdlink, &ElemKey);

// Optional parameters
	const char* arOptionNames[] = {0,0,0};
	const char* arOptionValues[] = {0,0,0};
	int Next = 0, AmOfOptions = 1;
	//AuxGetMathematicaOptions(Next, arOptionNames, arOptionValues, AmOfOptions);
	if(!AuxGetMathematicaOptions(Next, arOptionNames, arOptionValues, AmOfOptions)) return; //OC240612

	rad.GraphicsForElem_g3d(ElemKey, 1, arOptionNames, arOptionValues, AmOfOptions);
	//rad.GraphicsForElem_g3d(ElemKey, 1);

	for(int k=0; k<AmOfOptions; k++) 
	{
		//if(arOptionNames[k] != 0) MLDisownSymbol(stdlink, arOptionNames[k]);
		if(arOptionNames[k] != 0) MLWrapDeleteSymbol(stdlink, arOptionNames[k]); //OC240612
		//if(arOptionValues[k] != 0) MLDisownSymbol(stdlink, arOptionValues[k]);
		if(arOptionValues[k] != 0) MLWrapDeleteSymbol(stdlink, arOptionValues[k]); //OC240612
	}
#endif
}

//-------------------------------------------------------------------------

void GraphicsForAllWithoutSymChilds()
{
	rad.GraphicsForAll_g3d(0);
}

//-------------------------------------------------------------------------

void GraphicsForAllWithSymChilds()
{
	rad.GraphicsForAll_g3d(1);
}

//-------------------------------------------------------------------------

void ApplyDrawAttrToElem(int ElemKey, double R, double G, double B, double Thickness)
{
	double RGB_col[] = {R, G, B};
	rad.ApplyDrawAttrToElem_g3d(ElemKey, RGB_col, 3, Thickness);
}

//-------------------------------------------------------------------------

void ApplyColorToElem(int ElemKey, double R, double G, double B)
{
	double RGB_col[] = {R, G, B};
	rad.ApplyDrawAttrToElem_g3d(ElemKey, RGB_col, 3);
}

//-------------------------------------------------------------------------

void RemoveDrawAttrFromElem(int ElemKey)
{
	rad.RemoveDrawAttrFromElem_g3d(ElemKey);
}

//-------------------------------------------------------------------------

void QuickDraw3D_Viewer()
{
#ifdef __MATHEMATICA__
	int ElemKey;
	MLGetInteger(stdlink, &ElemKey);

	const char *FunName;
	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;
	int ReadOK = 0;

// Optional parameters
	const int AmOfOptions = 3;
	const char* OptionNames[] = {0,0,0};
	const char* OptionValues[] = {0,0,0};
	int OptionCount = 0;
	for(int ii=0; ii<AmOfOptions; ii++)
	{
		int Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC)
		{
			ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
			if((!ReadOK) || strcmp(FunName, "Rule") || (ArgNum!=2)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			if(!MLGetSymbol(stdlink, OptionNames + OptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			if(!MLGetSymbol(stdlink, OptionValues + (OptionCount++))) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC240612
		}
		else { MLNewPacket(stdlink); break;}
	}

	rad.GoQuickDraw3D_Viewer(ElemKey, OptionNames, OptionValues, OptionCount);

	for(int k=0; k<AmOfOptions; k++) 
	{
		//if(OptionNames[k] != 0) MLDisownSymbol(stdlink, OptionNames[k]);
		if(OptionNames[k] != 0) MLWrapDeleteSymbol(stdlink, OptionNames[k]); //OC240612
		//if(OptionValues[k] != 0) MLDisownSymbol(stdlink, OptionValues[k]);
		if(OptionValues[k] != 0) MLWrapDeleteSymbol(stdlink, OptionValues[k]); //OC240612
	}
#endif
}

//-------------------------------------------------------------------------

void OpenGL_3D_Viewer()
{
#ifdef __MATHEMATICA__
	int ElemKey;
	MLGetInteger(stdlink, &ElemKey);

	const char *FunName;
	//long ArgNum; //OC240907 (port to math6)
	int ArgNum;
	int ReadOK = 0;

// Optional parameters
	const int AmOfOptions = 3;
	const char* OptionNames[] = {0,0,0};
	const char* OptionValues[] = {0,0,0};
	int OptionCount = 0;
	for(int ii=0; ii<AmOfOptions; ii++)
	{
		int Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC)
		{
			ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
			if((!ReadOK) || strcmp(FunName, "Rule") || (ArgNum!=2)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			if(!MLGetSymbol(stdlink, OptionNames + OptionCount)) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			if(!MLGetSymbol(stdlink, OptionValues + (OptionCount++))) { rad.Send.ErrorMessage("Radia::Error062"); return;}
			//MLDisownSymbol(stdlink, FunName);
			MLWrapDeleteSymbol(stdlink, FunName); //OC240612
		}
		else { MLNewPacket(stdlink); break;}
	}

	rad.GoOpenGL_3D_Viewer(ElemKey, OptionNames, OptionValues, OptionCount);

	for(int k=0; k<AmOfOptions; k++) 
	{
		//if(OptionNames[k] != 0) MLDisownSymbol(stdlink, OptionNames[k]);
		if(OptionNames[k] != 0) MLWrapDeleteSymbol(stdlink, OptionNames[k]); //OC240612
		//if(OptionValues[k] != 0) MLDisownSymbol(stdlink, OptionValues[k]);
		if(OptionValues[k] != 0) MLWrapDeleteSymbol(stdlink, OptionValues[k]); //OC240612
	}
#endif
}

//-------------------------------------------------------------------------

void QuickDraw3D_ViewerOpt(int ElemKey, const char* Opt1, const char* Opt2, const char* Opt3)
{
	char CharBuf1[200], CharBuf2[200], CharBuf3[200];
	const char* OptionNames[] = {CharBuf1, CharBuf2, CharBuf3};
	const char* OptionValues[] = {0,0,0};
	const char* NonParsedOpts[] = {Opt1, Opt2, Opt3};
	int OptionCount = 3;
	//AuxParseOptionNamesAndValues(3, NonParsedOpts, OptionNames, OptionValues, OptionCount);
	AuxParseOptionNamesAndValues(NonParsedOpts, OptionNames, OptionValues, OptionCount);

/**
	if(Opt1 != 0)
	{
		if(*Opt1 != '\0') 
		{
			if(!AuxSetOptionNameAndValue(Opt1, (char*)(OptionNames[OptionCount]), (char**)(OptionValues + OptionCount))) return;
			OptionCount++;
		}
	}
	if(Opt2 != 0)
	{
		if(*Opt2 != '\0') 
		{
			if(!AuxSetOptionNameAndValue(Opt2, (char*)(OptionNames[OptionCount]), (char**)(OptionValues + OptionCount))) return;
			OptionCount++;
		}
	}
	if(Opt2 != 0)
	{
		if(*Opt3 != '\0') 
		{
			if(!AuxSetOptionNameAndValue(Opt3, (char*)(OptionNames[OptionCount]), (char**)(OptionValues + OptionCount))) return;
			OptionCount++;
		}
	}
**/
	rad.GoQuickDraw3D_Viewer(ElemKey, OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

void OpenGL_3D_ViewerOpt(int ElemKey, const char* Opt1, const char* Opt2, const char* Opt3)
{
	char CharBuf1[200], CharBuf2[200], CharBuf3[200];
	const char* OptionNames[] = {CharBuf1, CharBuf2, CharBuf3};
	const char* OptionValues[] = {0,0,0};
	const char* NonParsedOpts[] = {Opt1, Opt2, Opt3};
	int OptionCount = 3;
	//AuxParseOptionNamesAndValues(3, NonParsedOpts, OptionNames, OptionValues, OptionCount);
	AuxParseOptionNamesAndValues(NonParsedOpts, OptionNames, OptionValues, OptionCount);

	rad.GoOpenGL_3D_Viewer(ElemKey, OptionNames, OptionValues, OptionCount);
}

//-------------------------------------------------------------------------

void QuickDraw3D_ViewerOpt3(int ElemKey, const char* Opt1, const char* Opt2, const char* Opt3)
{
	char CharBuf1[200], CharBuf2[200], CharBuf3[200];
	strcpy(CharBuf1, Opt1);
	strcpy(CharBuf2, Opt2);
	strcpy(CharBuf3, Opt3);

	char *pEndOptName = strrchr(CharBuf1, '-');
	if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
	*pEndOptName = '\0';

	pEndOptName = strrchr(CharBuf2, '-');
	if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
	*pEndOptName = '\0';

	pEndOptName = strrchr(CharBuf3, '-');
	if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
	*pEndOptName = '\0';

	const char* OptionNames[] = {CharBuf1, CharBuf2, CharBuf3};
	const char* OptionValues[] = {0,0,0};
	OptionValues[0] = strrchr(Opt1, '>') + 1;
	OptionValues[1] = strrchr(Opt2, '>') + 1;
	OptionValues[2] = strrchr(Opt3, '>') + 1;
	rad.GoQuickDraw3D_Viewer(ElemKey, OptionNames, OptionValues, 3);
}

void QuickDraw3D_ViewerOpt2(int ElemKey, const char* Opt1, const char* Opt2)
{
	char CharBuf1[200], CharBuf2[200];
	strcpy(CharBuf1, Opt1);
	strcpy(CharBuf2, Opt2);

	char *pEndOptName = strrchr(CharBuf1, '-');
	if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
	*pEndOptName = '\0';

	pEndOptName = strrchr(CharBuf2, '-');
	if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
	*pEndOptName = '\0';

	const char* OptionNames[] = {CharBuf1, CharBuf2, 0};
	const char* OptionValues[] = {0,0,0};
	OptionValues[0] = strrchr(Opt1, '>') + 1;
	OptionValues[1] = strrchr(Opt2, '>') + 1;
	rad.GoQuickDraw3D_Viewer(ElemKey, OptionNames, OptionValues, 2);
}

void QuickDraw3D_ViewerOpt1(int ElemKey, const char* Opt1)
{
	char CharBuf[200];
	strcpy(CharBuf, Opt1);
	char *pEndOptName = strrchr(CharBuf, '-');
	if(pEndOptName == NULL) { rad.Send.ErrorMessage("Radia::Error062"); return;}
	*pEndOptName = '\0';

	const char* OptionNames[] = {CharBuf, 0, 0};
	const char* OptionValues[] = {0,0,0};
	OptionValues[0] = strrchr(Opt1, '>') + 1;
	rad.GoQuickDraw3D_Viewer(ElemKey, OptionNames, OptionValues, 1);
}

void QuickDraw3D_ViewerOpt0(int ElemKey)
{
	const char* OptionNames[] = {0, 0, 0};
	const char* OptionValues[] = {0,0,0};
	rad.GoQuickDraw3D_Viewer(ElemKey, OptionNames, OptionValues, 0);
}

//-------------------------------------------------------------------------

void DeleteElement(int ElemKey)
{
	rad.DeleteElement(ElemKey);
}

//-------------------------------------------------------------------------

void DeleteAllElements1()
{
	rad.DeleteAllElements(1);
}

//-------------------------------------------------------------------------

void DeleteAllElements2()
{
	rad.DeleteAllElements(2);
}

//-------------------------------------------------------------------------

void InterruptTime(double t)
{
	radYield.YieldInit(t);
	rad.Send.Double(t);
}

//-------------------------------------------------------------------------

void RadiaVersion()
{
	rad.ReturnVersionID();
}

//-------------------------------------------------------------------------

void ReturnInput(double Input, int NumTimes)
{
	rad.ReturnInput(Input, NumTimes);
}

//-------------------------------------------------------------------------

void MemAllocMethForIntrctMatr(char* TotOrParts)
{
	rad.SetMemAllocMethForIntrctMatr(TotOrParts);
}

//-------------------------------------------------------------------------

void OutCommandsInfo()
{
	rad.Send.OrdinaryMessage("Radia::usage");
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

#ifdef __GCC__

//-------------------------------------------------------------------------

void StartProf(int, int, int) {}
void StopProf() {}

//-------------------------------------------------------------------------

#ifdef __MATHEMATICA__

//int main()
//{
//	int argc;
//	char* argv[] = {"-linkmode", "parentconnect", 0};
//	return MLMain(argc, argv);
//}
//To enable passing default link parameters by Mathematica 
//(modified after changing default protocol in Mathematica 7):
int main(int argc, char* argv[])
{
     return MLMain(argc, argv);
}

#endif
#endif

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

#ifdef WIN32

//-------------------------------------------------------------------------

void StartProf(int, int, int) {}
void StopProf() {}

//-------------------------------------------------------------------------

HINSTANCE hinstCurrentRadia;
HINSTANCE hinstPreviousRadia;
LPSTR lpszCmdLineRadia;
int nCmdShowRadia;

#ifdef _WITH_QD3D
extern "C" {
extern void FreeQD3DViewerDLL();
extern void FreeQD3D_DLL();
}
#endif

//-------------------------------------------------------------------------

#ifdef __MATHEMATICA__

int WINAPI WinMain(HINSTANCE hinstCurrent, HINSTANCE hinstPrevious, LPSTR lpszCmdLine, int nCmdShow)
{
	char buff[512];
	char FAR * buff_start = buff;

#ifdef NDEBUG_EXTRA
	char FAR * argv[32];
	char FAR * FAR * argv_end = argv + 32;
#else
	char FAR * argv[] = {"-linkmode", "listen", "-linkname", "RadiaLink"};
	char FAR * FAR * argv_end = argv + 4;
#endif

	if( !MLInitializeIcon(hinstCurrent, nCmdShow)) return 1;

#ifdef NDEBUG_EXTRA
	MLScanString( argv, &argv_end, &lpszCmdLine, &buff_start);
#endif

	hinstCurrentRadia = hinstCurrent;
	hinstPreviousRadia = hinstPrevious;
	lpszCmdLineRadia = lpszCmdLine;
	nCmdShowRadia = nCmdShow;

	int MLMainFinishedOK = MLMain((int)(argv_end - argv), argv);

#ifdef _WITH_QD3D
	FreeQD3DViewerDLL();
	FreeQD3D_DLL();
#endif

	return MLMainFinishedOK;
}

#endif
#endif

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

#ifdef __MWERKS__
#ifdef _MAC_OS

//-------------------------------------------------------------------------

int profile=0;

void StartProf(int flag,int NumFunc, int Depth) 
{
/**
	if (profile==0) {
	if (flag==0)
	ProfilerInit(collectSummary,bestTimeBase,(long)NumFunc,(long)Depth);
	if (flag==1)
	ProfilerInit(collectDetailed,bestTimeBase,(long)NumFunc,(long)Depth);
	profile=1;
	MLPutString(stdlink,"Profiler Started");
	}
	else {
	MLPutString(stdlink,"Profiler Already Started");
	}
**/
}

//-------------------------------------------------------------------------

void StopProf() 
{
/**
	if (profile==1) {
	long fun,sta;
	ProfilerDump("\pProf");
	ProfilerGetDataSizes(&fun,&sta);
	ProfilerTerm();
	MLPutFunction(stdlink, "List", 2);
	MLPutLongInteger(stdlink,fun);
	MLPutLongInteger(stdlink,sta);
	profile=0;}
	else {
	MLPutString(stdlink,"Profiler not Started");
	}
**/
}

//-------------------------------------------------------------------------

#ifdef _WITH_QD3D
extern "C" {
extern void QuitQD3D();
}
#endif

//-------------------------------------------------------------------------

#ifdef __MATHEMATICA__

int main(int argc, char* argv[])
{
	set_new_handler(&radMemAllocHandler);
	int MLMainFinishedOK = MLMain(argc, argv);

#ifdef _WITH_QD3D
	QuitQD3D();
#endif

	return MLMainFinishedOK;
}

#endif

//-------------------------------------------------------------------------

#endif
#endif